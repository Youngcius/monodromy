import logging
from functools import lru_cache

import numpy as np
import retworkx
from qiskit.circuit import Instruction
from qiskit.dagcircuit import DAGCircuit, DAGOpNode
from qiskit.transpiler import AnalysisPass, TransformationPass
from qiskit.transpiler.exceptions import TranspilerError
from qiskit.transpiler.passes import Collect2qBlocks, ConsolidateBlocks

from monodromy.coverage import (
    convert_gate_to_monodromy_coordinate,
    coverage_lookup_cost,
    gates_to_coverage,
    print_coverage_set,
)


class MonodromyPass(AnalysisPass):
    """MonodromyDepth class extends the AnalysisPass to perform cost analysis
    on a given CircuitDAG with respect to a specified 2-qubit basis gate. This
    basis gate is crucial in calculating the minimum execution cost of 2-qubit
    blocks within the CircuitDAG.

    This class is particularly useful for quantum circuit optimization
    where the cost associated with the execution of certain gates is a
    crucial factor in the overall performance of the quantum computer.

    This class requires the Collect2qBlocks and ConsolidateBlocks passes
    to decompose the CircuitDAG into 2-qubit blocks and consolidate them
    respectively.
    """

    _coverage_cache = {}  # Class level cache

    def __init__(
        self,
        basis_gate: Instruction,
        gate_cost: float = 1.0,
        consolidate: bool = False,
        consolidator: TransformationPass = None,
        use_fast_settings: bool = True,
    ):
        """Initializes the MonodromyDepth pass.

        Args:
            basis_gate (Instruction): The 2-qubit basis gate to be used in
                the cost analysis.
            gate_cost (float, optional): The cost of the basis gate. Defaults
                to 1.0.
            consolidate (bool, optional): Whether to consolidate the 2-qubit
                blocks in the CircuitDAG. Defaults to False.
            consolidator (TransformationPass, optional): The pass to use for
                consolidation. Defaults to None.

            XXX For correctness: consolidate should be True. It defaults to False under
            the assumption that the pass manager has already included consolidation.
        """
        super().__init__()
        assert basis_gate.num_qubits == 2, "Basis gate must be a 2Q gate."

        if (
            consolidate or consolidator
        ):  # default False: assume pass manager has already done this
            if consolidator is not None:
                self.requires = [consolidator]
            else:
                self.requires = [
                    Collect2qBlocks(),
                    ConsolidateBlocks(force_consolidate=True),
                ]

        self.basis_gate = basis_gate
        self.gate_cost = gate_cost
        self.chatty = True
        self.use_fast_settings = use_fast_settings

        # Use the basis_gate as a key for the cache
        basis_gate_key = str(self.basis_gate)
        if basis_gate_key in MonodromyDepth._coverage_cache:
            self.coverage_set = MonodromyDepth._coverage_cache[basis_gate_key]
        else:
            self.coverage_set = self._gate_set_to_coverage()
            MonodromyDepth._coverage_cache[basis_gate_key] = self.coverage_set

    @lru_cache(maxsize=None)
    def _gate_set_to_coverage(self):
        """The gate_set_to_coverage() function takes the basis gate and creates
        a CircuitPolytope object that represents all the possible 2Q unitaries
        that can be formed by piecing together different instances of the basis
        gate.

        The coverage set considers the relative durations of the gates
        in the basis set so costs are already accounted in calls to
        coverage_lookup_cost.

        :return: A CircuitPolytope object
        """
        if self.chatty:
            logging.info("==== Working to build a set of covering polytopes ====")

        # TODO, here could add functionality for multiple basis gates
        # just need to fix the cost function to account for relative durations
        coverage_set = gates_to_coverage(
            self.basis_gate, costs=[self.gate_cost], sort=True
        )

        # TODO: add some warning or fail condition if the coverage set fails to coverage
        # one way, (but there may be a more direct way) is to check if expected haar == 0
        if self.chatty:
            logging.info("==== Done. Here's what we found: ====")
            logging.info(print_coverage_set(coverage_set))

        return coverage_set


class MonodromyCountSwaps(MonodromyPass):
    """MonodromyCountSwaps counts the number of SWAP gates in a given
    CircuitDAG."""

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        # check if gate has coordinate [0.25, 0.25, 0.25, -0.75], if so, it is a SWAP
        swap_count = 0
        for gate_node in dag.two_qubit_ops():
            t_c = convert_gate_to_monodromy_coordinate(gate_node.op)
            # need to check if close, ie [0.25, 0.25, 0.249999, -0.7499999]
            if np.allclose(t_c, [0.25, 0.25, 0.25, -0.75]):
                swap_count += 1
        self.property_set["total_swaps"] = swap_count
        return dag


class MonodromyTotal(MonodromyPass):
    """MonodromyTotal calculates the total number of gates in a given
    CircuitDAG, with respect to a specified 2-qubit basis gate."""

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        total_cost = 0
        for gate_node in dag.two_qubit_ops():
            total_cost += coverage_lookup_cost(self.coverage_set, gate_node.op)[0]
        self.property_set["monodromy_total"] = total_cost


class MonodromyDepth(MonodromyPass):
    """MonodromyDepth calculates the depth of a given CircuitDAG with respect
    to a specified 2-qubit basis gate."""

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        """Find longest_path by iterating over edges and using a weight
        function."""
        # from qiskit.converters import dag_to_circuit
        # print(dag_to_circuit(dag).draw(fold=-1))

        SCALE_FACTOR = 1000  # scale factor to convert floats to integers

        # Keeping track of the last iteration saves nearly 50% of the time
        last_lookup = (None, None)

        def weight_fn(_1, node, _2) -> int:
            """Weight function for longest path algorithm.

            Since this iterates over edges, every node is counted twice.
            (Since every 2Q gate will have 2 incoming edges).
            """
            if self.use_fast_settings:
                nonlocal last_lookup
                if node == last_lookup[0]:
                    return last_lookup[1]
            target_node = dag._multi_graph[node]
            if not isinstance(target_node, DAGOpNode):
                return 0
            elif target_node.op.name in ["barrier", "measure"]:
                return 0
            elif len(target_node.qargs) == 1:
                return 0
            elif len(target_node.qargs) > 2:
                raise TranspilerError("Operation not supported.")
            else:
                float_cost = coverage_lookup_cost(
                    self.coverage_set,
                    target_node.op,
                    use_fast_settings=self.use_fast_settings,
                )[0]
                int_cost = int(float_cost * SCALE_FACTOR)
                last_lookup = (node, int_cost)
                return int_cost

        longest_path_length = retworkx.dag_longest_path_length(
            dag._multi_graph, weight_fn=weight_fn
        )

        # convert back to float
        longest_path_length /= SCALE_FACTOR

        if self.chatty:
            logging.info(f"Longest path length: {longest_path_length}")

        self.property_set["monodromy_depth"] = longest_path_length
        return dag
