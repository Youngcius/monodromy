"""monodromy/render.py.

Utilities for rendering polytopes.
"""

from typing import List

import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial as ss
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from qiskit.circuit import Instruction
from weylchamber import WeylChamber

from monodromy.coordinates import monodromy_to_positive_canonical_polytope
from monodromy.coverage import CircuitPolytope
from monodromy.haar import gates_to_coverage


def _plot_polytope(circuit_polytope, w, color="red"):
    polytope_vertices = (
        monodromy_to_positive_canonical_polytope(circuit_polytope).reduce().vertices
    )
    polytope_vertices = np.array([[float(b) for b in a] for a in polytope_vertices[0]])

    left_vertices = polytope_vertices
    right_vertices = np.array([[1 - a[0], a[1], a[2]] for a in polytope_vertices])

    for vertices in [left_vertices, right_vertices]:
        # delete duplicates that might exist
        vertices = np.unique(vertices, axis=0)

        if len(vertices) < 3:
            w.ax.scatter3D(*zip(*vertices), color=color)
        elif len(vertices) == 3:
            triangle = Poly3DCollection([vertices])
            triangle.set_facecolor(color)
            triangle.set_edgecolor("k")
            triangle.set_alpha(0.5)
            w.ax.add_collection3d(triangle)
        else:
            # TODO use Qbk:0Bk:0 - drop dimension k from the input points
            hull = ss.ConvexHull(vertices, qhull_options="QJ")
            faces = Poly3DCollection([vertices[simplex] for simplex in hull.simplices])
            faces.set_facecolor(color)
            faces.set_alpha(0.2)
            faces.set_edgecolor("k")
            w.ax.add_collection3d(faces)


def _plot_coverage_set(coverage_set, overlap=False):
    """Plot a set of 3D polytopes.

    Args:
        coverage_set (list): a list of CircuitPolytope objects.
        overlap (bool): If True, all polytopes are drawn on the same plot. If False, each polytope is drawn in a separate subplot.
    """
    colors = [
        "red",
        "green",
        "blue",
        "orange",
        "purple",
        "cyan",
        "black",
        "pink",
        "brown",
        "grey",
    ]

    # Preprocess coverage_set to organize CircuitPolytope objects based on their cost
    organized_set = {}
    for circuit_polytope in coverage_set:
        cost = circuit_polytope.cost
        if cost == 0:
            continue
        if cost not in organized_set:
            organized_set[cost] = []
        organized_set[cost].append(circuit_polytope)

    if overlap:
        fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
        w = WeylChamber()
        w.labels = {}
        w.render(ax)
        for i, (cost, polytopes) in enumerate(organized_set.items()):
            color = colors[i % len(colors)]
            for circuit_polytope in polytopes:
                _plot_polytope(circuit_polytope, color=color, w=w)
    else:
        n = len(organized_set)
        fig, axs = plt.subplots(
            1, n, subplot_kw={"projection": "3d"}, figsize=(n * 5, 5)
        )  # Adjust size to avoid crowding
        for i, (cost, polytopes) in enumerate(organized_set.items()):
            ax = axs[i] if n > 1 else axs
            w = WeylChamber()
            w.labels = {}
            w.render(ax)
            color = colors[i % len(colors)]
            for circuit_polytope in polytopes:
                _plot_polytope(circuit_polytope, color=color, w=w)
            w.ax.set_title(f"Cost: {cost}")

    plt.show()
    # save fig
    fig.savefig("coverage_set.svg")


def gates_to_coverage_plot(*gates: Instruction, costs=None, overlap=False):
    """Plot the coverage set of a gate.

    Args:
        gate (Instruction): a gate.
        overlap (bool): If True, all polytopes are drawn on the same plot. If False, each polytope is drawn in a separate subplot.
    """
    coverage_set = gates_to_coverage(*gates, costs=costs, sort=True)
    _plot_coverage_set(coverage_set, overlap=overlap)
    return coverage_set


def polytopes_to_mathematica(necessary_polytopes: List[CircuitPolytope]):
    output = ""
    output += "polytopeData = {"
    for n, gate_polytope in enumerate(necessary_polytopes):
        output += "{"
        output += f"{float(gate_polytope.cost)}, "
        for i, polytope in enumerate(gate_polytope.convex_subpolytopes):
            output += "{"
            vertices = polytope.vertices
            for j, vertex in enumerate(vertices):
                output += "{"
                output += f"{vertex[0]}, {vertex[1]}, {vertex[2]}"
                if 1 + j != len(vertices):
                    output += "}, "
                else:
                    output += "}"
            if 1 + i != len(gate_polytope.convex_subpolytopes):
                output += "}, "
            else:
                output += "}"
        if 1 + n != len(necessary_polytopes):
            output += "}, "
        else:
            output += "}"
    output += "};"

    output += """
corners = {{0, 0, 0}, {1/4, 1/4, 1/4}, {1/4, 1/4, -(1/4)}, {3/8, 3/
    8, -(1/8)}, {3/8, -(1/8), -(1/8)}, {1/2, 0, 0}};
names = {{0, 0, 0} -> "I", {1/4, 1/4, -1/4} -> "CZ", {1/2, 0, 0} ->
    "ISWAP", {1/4, 1/4, 1/4} -> "SWAP", {3/8, 3/8, -1/8} -> Sqrt[
    SWAP], {3/8, -1/8, -1/8} -> Sqrt[SWAP]'};

(* tune this loop bound to skip degenerate solids *)
skipTo := 6
(* tune these scalars to get a picture w/o Z-fighting *)
OffsetCoords[coord_, n_] := ((1 + 0.5^n) (coord - {0.25, 0.1, -0.1}))

Module[{directives = {}, n, depth, vertices, maxdepth},
  maxdepth = Max[First /@ polytopeData];
  For[n = skipTo, n <= Length[polytopeData], n++,
   (* inject new color settings *)
   depth = polytopeData[[n, 1]];
   vertices = Rest[polytopeData[[n]]];
   directives = Join[directives, {
      Hue[(depth - skipTo)/maxdepth],
      Opacity[1.0 - (depth - skipTo)/maxdepth]
      }];
   (* add new polyhedra *)
   directives = Join[directives,
     With[{mesh = ConvexHullMesh[OffsetCoords[#, n] & /@ #]},
        GraphicsComplex[
         MeshCoordinates[mesh], {EdgeForm[], MeshCells[mesh, 2]}]] & /@
       vertices];
   ];
  directives];
Show[Graphics3D[Lighting -> "Neutral", Boxed -> False],
 Graphics3D@(Join @@
    Table[{Sphere[OffsetCoords[corner, skipTo], 0.02],
      Text[corner /. names,
       OffsetCoords[corner, 5 skipTo] +
        0.05*If[Norm[corner] == 0, 0, corner/Norm[corner]]]}, {corner,
       corners}]),
 Graphics3D[%, Boxed -> False, ViewPoint -> {0, 1, 1},
  Lighting -> {{"Ambient", White}}]]
"""

    return output
