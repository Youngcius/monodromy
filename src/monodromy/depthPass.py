import logging
from functools import lru_cache

import retworkx
from qiskit.circuit import Instruction
from qiskit.dagcircuit import DAGCircuit, DAGOpNode
from qiskit.transpiler import AnalysisPass, TransformationPass
from qiskit.transpiler.exceptions import TranspilerError
from qiskit.transpiler.passes import Collect2qBlocks, ConsolidateBlocks

from monodromy.coverage import (
    coverage_lookup_operation,
    gates_to_coverage,
    print_coverage_set,
)


class MonodromyDepth(AnalysisPass):
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
            if consolidator is None:
                self.requires = [consolidator]
            else:
                self.requires = [
                    Collect2qBlocks(),
                    ConsolidateBlocks(force_consolidate=True),
                ]

        self.basis_gate = basis_gate
        self.gate_cost = gate_cost
        self.chatty = True

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

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        """The run() method is the main entry point for the AnalysisPass. It
        takes a CircuitDAG as input and returns an updated CircuitDAG. This
        method applies the basis gate to the CircuitDAG, computes the cost of
        the applied gate, and updates the CircuitDAG accordingly.

        :param dag: The CircuitDAG to be analyzed.
        :return: An updated CircuitDAG.
        """
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
                float_cost = coverage_lookup_operation(
                    self.coverage_set, target_node.op
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
