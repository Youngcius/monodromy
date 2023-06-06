# `monodromy`

Computations in the monodromy polytope for quantum gate sets

## Fork

Qiskit AnalysisPass which determines the depth of a circuit. I want to use this in various transpiler passes as a subroutine, without actually having to perform decomposition. This tool lets me count the iterate the cost of a circuit without actually having to decompose it. Simple procedure, given a basis gate(s):

1. force consolidate unitary blocks
2. look up the cost using the monodromy polytope
3. return the longest path, weighted against the decomposition cost

```python
from monodromy.depthPass import MonodromyDepthPass
from qiskit.transpiler.passmanager import PassManager
from qiskit.circuit.library import CXGate
from qiskit import QuantumCircuit

pm = PassManager()
pm.append(MonodromyDepthPass(basis_gate=CXGate()))

qc= QuantumCircuit(4)
qc.swap(0,1)
qc.cx(0,1)
qc.cx(1,2)
qc.swap(0,1)
qc.rz(0.5, 0)
qc.cx(0,1)
qc.ry(0.7, 1)
qc.cx(1,2)
qc.swap(0,1)
qc.cx(2,3)

pm.run(qc)
expected_value = 9
assert pm.property_set["monodromy_depth"] == expected_value, "Monodromy depth not calculated correctly!"

```

#### Change log:

1. Small fix in setup.py for installation
2. WIP Write consolidate + depth as a Qiskit AnalysisPass

The original work was done in this repository (https://github.com/Pitt-JonesLab/slam_decomposition), which became a [mess](https://github.com/Pitt-JonesLab/slam_decomposition/blob/main/src/slam/utils/polytopes/polytope_wrap.py).

---

## Overview

[Fixed-Depth Two-Qubit Circuits and the Monodromy Polytope](https://arxiv.org/abs/1904.10541) described a technique for determining whether a given two-qubit unitary can be written as a circuit with a prescribed sequence of two-qubit interactions, interleaved with arbitrary single-qubit unitaries.
This python package is a computationally effective implementation of that technique.

## Installation

1. This package comes with a `requirements.txt` file.
   Begin by installing those requirements, using `pip -r requirements.txt`.
2. Install [`lrs`](http://cgm.cs.mcgill.ca/~avis/C/lrs.html).
   Typically, this means downloading the source, building it, and placing the generated executable somewhere in the search path for your python process.
   **NOTE:** We require either version ≥0.7.2 _or_ version ≥0.7.1 of `lrsgmp` (renamed to `lrs`).
3. _Optionally,_ install [`lrcalc`](https://sites.math.rutgers.edu/~asbuch/lrcalc/). Typically, this means downloading the source, building it, then building the Cython bindings, and installing the resulting package. _This is not necessary: `lrcalc` is needed only for the curious user who wants to regenerate the contents of `qlr_table.py`._

## Configuration

If you're unable to put `lrs` in your search path, you can specify its location via the environment variable `LRS_PATH`.

## Usage

Give it a whirl with `scripts/demo.py`.

## Notes

We've designed the package around polytopes with rational coordinates.
In practice this suits the use cases of quantum computer scientists fine, but it is easy to imagine use cases outside of this.
`lrs` makes strong assumptions about the kind of arithmetic used, so if one were to want to process irrational polytopes, one would be obligated to move away from `lrs`.

`lrs` is licensed under GPLv2, which makes it insuitable for inclusion in some projects.
We may someday want to remove (or make optional) this dependency from this project.
To this end, we have set up a `backend` which describes a contained and minimal set of calls we make of any computational library.
