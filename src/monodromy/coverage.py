"""monodromy/coverage.py.

Routines for converting a family of native gates to a minimal set of
minimum- cost circuits whose union covers the entire monodromy polytope.
"""

import heapq
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Optional, Tuple

import numpy as np
from qiskit.circuit import Instruction
from qiskit.quantum_info import Operator
from qiskit.transpiler.exceptions import TranspilerError

from monodromy.approximate import polytope_approx_decomp_fidelity
from monodromy.coordinates import unitary_to_monodromy_coordinate
from monodromy.static.examples import exactly

from .coordinates import monodromy_alcove, monodromy_alcove_c2, rho_reflect
from .elimination import cylinderize, project
from .io.base import CircuitPolytopeData
from .polytopes import Polytope, trim_polytope_set
from .static.examples import everything_polytope, identity_polytope
from .static.qlr_table import qlr_polytope

# duration cost ratio of single qubit gates to two qubit gates
# e.g. 50ns to 100ns, -> 1/2
Q1_GATE_COST = 0


@dataclass
class CircuitPolytope(Polytope, CircuitPolytopeData):
    """A polytope describing the alcove coverage of a particular circuit
    type."""

    def __gt__(self, other):
        return (self.cost > other.cost) or (
            self.cost == other.cost and self.volume > other.volume
        )

    def __ge__(self, other):
        return (self.cost > other.cost) or (
            self.cost == other.cost and self.volume >= other.volume
        )

    def __lt__(self, other):
        return (self.cost < other.cost) or (
            self.cost == other.cost and self.volume < other.volume
        )

    def __le__(self, other):
        return (self.cost < other.cost) or (
            self.cost == other.cost and self.volume <= other.volume
        )


def _operation_to_circuit_polytope(
    operation: Instruction,
    cost=1,
    name=None,
    instruction=None,
    single_qubit_cost=Q1_GATE_COST,
) -> CircuitPolytope:
    """The operation_to_circuit_polytope() function takes a qiskit.Instruction
    object and returns a CircuitPolytope object that represents the unitary of
    the operation.

    :param operation: A qiskit.Instruction object.
    :return: A CircuitPolytope object
    """

    gd = operation.to_matrix()
    b_polytope = exactly(
        *(
            Fraction(x).limit_denominator(10_000)
            for x in unitary_to_monodromy_coordinate(gd)[:-1]
        )
    )
    convex_polytope = deduce_qlr_consequences(
        target="c",
        a_polytope=identity_polytope,
        b_polytope=b_polytope,
        c_polytope=everything_polytope,
    )

    if name is None:  # try to come up with something
        op_name = (
            f"{operation.name}({operation.params[0]:.5f})"
            if operation.params
            and len(operation.params) > 0
            and not isinstance(operation.params[0], np.ndarray)
            else f"{operation.name}"
        )
    else:
        op_name = name  # Use provided name if available

    return CircuitPolytope(
        operations=[op_name],
        instructions=[instruction],
        # cost=cost + 2 * single_qubit_cost, # * 2 for first and last layers
        cost=cost
        + single_qubit_cost,  # NOTE new convention only counts 1 of the exterior
        convex_subpolytopes=convex_polytope.convex_subpolytopes,
    )


def gates_to_coverage(
    *gates: Instruction,
    costs=None,
    names=None,
    instructions=None,
    sort=True,
    single_qubit_cost=Q1_GATE_COST,
) -> List[CircuitPolytope]:
    """Calculates coverage given a basis gate set."""
    for gate in gates:
        assert gate.num_qubits == 2, "Basis gate must be a 2Q gate."

    # default costs for all gates are 1, except SWAP which is 0 (virtual-SWAP)
    if costs is None:
        costs = [1 if gate.name != "swap" else 0 for gate in gates]

    if names is None:
        names = [None for _ in range(len(costs))]

    if instructions is None:
        instructions = [None for _ in range(len(costs))]

    operations = [
        _operation_to_circuit_polytope(
            gate, cost=c, name=n, instruction=i, single_qubit_cost=single_qubit_cost
        )
        for gate, c, n, i in zip(gates, costs, names, instructions)
    ]
    coverage_set = build_coverage_set(operations, single_qubit_cost=single_qubit_cost)

    if sort:
        return sorted(coverage_set, key=lambda k: k.cost)

    return coverage_set


# want to separate these functions from existing code
# tricky because the hacky way to perform approx cost uses decomposition
# cost? -> trial and error converge -> return k when success
# but I would prefer if we would instead do cost -> decomposition at k.
# this works for exact cost, but not approx cost. which means
# (a) need to not make a mess of the existing code,
# (b) refactor it later once have a better nearest coordinate function.
# see https://github.com/evmckinney9/monodromy/issues/1
# be careful, we need to use approx_degree to represent total fidelity
# idea is that threshold of accepting changes as function of instruction set
# more expensive ansatz as a higher threshold of accepting approximations


# TODO?
def coverage_lookup_decomposition():
    raise NotImplementedError


def convert_gate_to_monodromy_coordinate(
    gate: Instruction, use_fast_settings: bool = True
) -> List[Fraction]:
    """Converts a gate to a monodromy coordinate."""
    # convert gate to monodromy coordinate
    if (
        hasattr(gate, "_monodromy_coord")
        and gate._monodromy_coord is not None
        and use_fast_settings
    ):
        return gate._monodromy_coord
    else:
        try:
            return unitary_to_monodromy_coordinate(gate.to_matrix())
        except AttributeError:
            return unitary_to_monodromy_coordinate(Operator(gate).data)


def coverage_lookup_cost(
    coverage_set: List[CircuitPolytope],
    target: Instruction,
    error_model=None,
    allow_approx=False,
    use_fast_settings: bool = True,
) -> Tuple[float, float]:
    """Calculates the cost of an operation.

    Find the cost of an operation by iterating through the coverage set, sorted by cost.
    Args:
        coverage_set (List[CircuitPolytope]): The coverage set to search
        target (Instruction): The operation to find the cost of
        error_model (ErrorModel): Model to use when calculating fidelity from time costs
        allow_approx (bool): Whether to allow approximation of the operation
        use_fast_settings (bool): Whether to use fast settings for polytope
    Returns:
        float: The cost of the operation
        float: Expected fidelity of the operation,
            total fidelity = fidelity(circuit cost) * fidelity(decomposition)
    """

    target_coords = convert_gate_to_monodromy_coordinate(target, use_fast_settings)

    # iterate through coverage set to find exact decomposition polytope
    # XXX assume has been sorted by cost
    exact_polytope, exact_decomp_fidelity = None, None
    for polytope_index, circuit_polytope in enumerate(coverage_set):
        # assume if last polytope, then must be contained
        if polytope_index == len(coverage_set) - 1:
            exact_polytope, exact_decomp_fidelity = circuit_polytope, 1.0
        if circuit_polytope.has_element(
            target_coords, use_fast_settings=use_fast_settings
        ):
            exact_polytope, exact_decomp_fidelity = circuit_polytope, 1.0
            break

    # if not allowing approximations, return exact decomposition
    if not allow_approx or error_model is None:
        polytope_fid = error_model.fidelity(exact_polytope.cost) if error_model else 1.0
        exact_fid = polytope_fid * exact_decomp_fidelity
        return exact_polytope.cost, exact_fid

    # else set approx_degree to be set from polytope cost
    approx_degree = error_model.fidelity(exact_polytope.cost) * exact_decomp_fidelity

    # check special case, no approximation improvement possible
    if approx_degree == 0.0:
        return exact_polytope.cost, 1.0

    # iterate through coverage set again
    # look for approximate decomp that outperforms approx_degree infidelity
    for polytope_index, circuit_polytope in enumerate(coverage_set):
        if polytope_index == len(coverage_set) - 1:
            return circuit_polytope.cost, error_model.fidelity(exact_polytope.cost)

        # case when approx can no longer do better
        if circuit_polytope.cost >= exact_polytope.cost:
            return exact_polytope.cost, approx_degree

        polytope_fid = error_model.fidelity(circuit_polytope.cost)
        approx_decomp_fid = polytope_approx_decomp_fidelity(circuit_polytope, target)
        if approx_decomp_fid * polytope_fid > approx_degree:
            return circuit_polytope.cost, approx_decomp_fid * polytope_fid

    raise TranspilerError("Operation not found in coverage set.")


def deduce_qlr_consequences(
    target: str,
    a_polytope: Polytope,
    b_polytope: Polytope,
    c_polytope: Polytope,
    extra_polytope: Optional[Polytope] = None,
) -> Polytope:
    """Produces the consequences for `target` for a family of a-, b-, and
    c-inequalities.

    `target` can take on the values 'a', 'b', or 'c'.
    """

    coordinates = {
        "a": [0, 1, 2, 3],
        "b": [0, 4, 5, 6],
        "c": [0, 7, 8, 9],
    }
    assert target in coordinates.keys()

    if extra_polytope is None:
        extra_polytope = everything_polytope

    p = extra_polytope.intersect(qlr_polytope)
    p = p.union(rho_reflect(p, coordinates[target]))
    # impose the "large" alcove constraints
    for value in coordinates.values():
        p = p.intersect(cylinderize(monodromy_alcove, value))

    # also impose whatever constraints we were given besides
    p = p.intersect(cylinderize(a_polytope, coordinates["a"]))
    p = p.intersect(cylinderize(b_polytope, coordinates["b"]))
    p = p.intersect(cylinderize(c_polytope, coordinates["c"]))

    # restrict to the A_{C_2} part of the target coordinate
    p = p.intersect(cylinderize(monodromy_alcove_c2, coordinates[target]))

    # lastly, project away the non-target parts
    p = p.reduce()
    for index in range(9, 0, -1):
        if index in coordinates[target]:
            continue
        p = project(p, index)
        p = p.reduce()

    return p


def prereduce_operation_polytopes(
    operations: List[CircuitPolytope],
    target_coordinate: str = "c",
    background_polytope: Optional[Polytope] = None,
    chatty: bool = False,
) -> Dict[str, CircuitPolytope]:
    """Specializes the "b"-coordinates of the monodromy polytope to a
    particular operation, then projects them away."""

    coordinates = {
        "a": [0, 1, 2, 3],
        "b": [0, 4, 5, 6],
        "c": [0, 7, 8, 9],
    }
    prereduced_operation_polytopes = {}

    for operation in operations:
        if chatty:
            print(f"Prereducing QLR relations for {'.'.join(operation.operations)}")
        p = (
            background_polytope
            if background_polytope is not None
            else everything_polytope
        )
        p = p.intersect(qlr_polytope)
        p = p.union(rho_reflect(p, coordinates[target_coordinate]))
        for value in coordinates.values():
            p = p.intersect(cylinderize(monodromy_alcove, value))
        p = p.intersect(cylinderize(operation, coordinates["b"]))
        p = p.intersect(
            cylinderize(monodromy_alcove_c2, coordinates[target_coordinate])
        )

        # project away the operation part
        p = p.reduce()
        for index in [6, 5, 4]:
            p = project(p, index)
            p = p.reduce()
        prereduced_operation_polytopes[operation.operations[-1]] = p

    return prereduced_operation_polytopes


def build_coverage_set(
    operations: List[CircuitPolytope],
    chatty: bool = False,
    single_qubit_cost=Q1_GATE_COST,
) -> List[CircuitPolytope]:
    """Given a set of `operations`, thought of as members of a native gate set,
    this emits a list of circuit shapes built as sequences of those operations
    which is:

    + Exhaustive: Every two-qubit unitary is covered by one of the
    circuit designs in the list.
    + Irredundant: No circuit design is completely contained within other designs
    in the list which are of equal or lower cost.

    If `chatty` is toggled, emits progress messages.
    """

    # assert that operations.operation are unique
    if len(set([operation.operations[-1] for operation in operations])) != len(
        operations
    ):
        raise ValueError("Operations must be unique.")

    # start by generating precalculated operation polytopes
    prereduced_operation_polytopes = prereduce_operation_polytopes(operations)

    # a collection of polytopes explored so far, and their union
    total_polytope = CircuitPolytope(
        convex_subpolytopes=identity_polytope.convex_subpolytopes,
        operations=[],
        instructions=[],
        cost=0.0,
    )
    necessary_polytopes = [total_polytope]

    # a priority queue of sequences to be explored
    to_be_explored = []
    for operation in operations:
        heapq.heappush(to_be_explored, operation)

    # a set of polytopes waiting to be reduced, all of equal cost
    to_be_reduced = []
    waiting_cost = 0

    # main loop: dequeue the next cheapest gate combination to explore
    while 0 < len(to_be_explored):
        next_polytope = heapq.heappop(to_be_explored)

        # if this cost is bigger than the old cost, flush
        if next_polytope.cost > waiting_cost:
            to_be_reduced = trim_polytope_set(
                to_be_reduced, fixed_polytopes=[total_polytope]
            )
            necessary_polytopes += to_be_reduced
            for new_polytope in to_be_reduced:
                total_polytope = total_polytope.union(new_polytope).reduce()
            to_be_reduced = []
            waiting_cost = next_polytope.cost

        if chatty:
            print(f"Considering {'·'.join(next_polytope.operations)};\t", end="")

        # find the ancestral polytope
        tail_polytope = next(
            (
                p
                for p in necessary_polytopes
                if p.operations == next_polytope.operations[:-1]
            ),
            None,
        )
        # if there's no ancestor, skip
        if tail_polytope is None:
            if chatty:
                print("no ancestor, skipping.")
            continue

        # take the head's polytopes, adjoin the new gate (& its rotation),
        # calculate new polytopes, and add those polytopes to the working set
        # TODO: This used to be part of a call to `deduce_qlr_consequences`, which
        #       we split up for efficiency.  See GH #13.
        new_polytope = prereduced_operation_polytopes[next_polytope.operations[-1]]
        new_polytope = new_polytope.intersect(
            cylinderize(tail_polytope, [0, 1, 2, 3], parent_dimension=7)
        )
        new_polytope = new_polytope.reduce()
        for index in [3, 2, 1]:
            new_polytope = project(new_polytope, index).reduce()
        # specialize it from a Polytope to a CircuitPolytope
        new_polytope = CircuitPolytope(
            operations=next_polytope.operations,
            instructions=next_polytope.instructions,
            cost=next_polytope.cost,
            convex_subpolytopes=new_polytope.convex_subpolytopes,
        )

        to_be_reduced.append(new_polytope)

        if chatty:
            print(f"Cost {next_polytope.cost} ", end="")
            volume = new_polytope.volume
            if volume.dimension == 3:
                volume = volume.volume / monodromy_alcove_c2.volume.volume
                print(f"and Euclidean volume {float(100 * volume):6.2f}%")
            else:
                print(f"and Euclidean volume {0:6.2f}%")

        # if this polytope is NOT of maximum volume,
        if monodromy_alcove_c2.volume > new_polytope.volume:
            # add this polytope + the continuations to the priority queue
            # NOTE subtract 1Q cost, because each successive layer +1, but each circuitpolytope is +2
            for operation in operations:
                heapq.heappush(
                    to_be_explored,
                    CircuitPolytope(
                        operations=next_polytope.operations + operation.operations,
                        instructions=next_polytope.instructions
                        + operation.instructions,
                        # cost=next_polytope.cost + operation.cost - single_qubit_cost,
                        cost=next_polytope.cost
                        + operation.cost,  # NOTE updated cost convention
                        convex_subpolytopes=operation.convex_subpolytopes,
                    ),
                )
        else:
            # the cheapest option that gets us to 100% is good enough.
            break

    # one last flush
    necessary_polytopes += trim_polytope_set(
        to_be_reduced, fixed_polytopes=[total_polytope]
    )
    for new_polytope in to_be_reduced:
        total_polytope = total_polytope.union(new_polytope).reduce()

    return necessary_polytopes


def print_coverage_set(necessary_polytopes):
    print("Percent volume of A_C2\t | Cost\t | Sequence name")
    for gate in necessary_polytopes:
        vol = gate.volume
        if vol.dimension == 3:
            vol = vol.volume / monodromy_alcove_c2.volume.volume
        else:
            vol = Fraction(0)
        print(
            f"{float(100 * vol):6.2f}% = "
            f"{str(vol.numerator): >4}/{str(vol.denominator): <4} "
            f"\t | {float(gate.cost):4.2f}"
            f"\t | {'.'.join(gate.operations)}"
        )
