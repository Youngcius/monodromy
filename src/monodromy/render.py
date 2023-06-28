"""
monodromy/render.py

Utilities for rendering polytopes.
"""

from typing import List

from monodromy.coverage import CircuitPolytope

import matplotlib.pyplot as plt
from monodromy.coordinates import monodromy_to_positive_canonical_polytope
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from weylchamber import WeylChamber
import scipy.spatial as ss

def _plot_polytope(circuit_polytope, w, color='red'):
    polytope_vertices = monodromy_to_positive_canonical_polytope(circuit_polytope).reduce().vertices
    polytope_vertices = np.array([[float(b) for b in a] for a in polytope_vertices[0]])

    # seems that if (0,0,0) is included we need to manually add (1,0,0) [both are identity]
    if np.any(np.all(polytope_vertices == [0,0,0], axis=1)):
        polytope_vertices = np.append(polytope_vertices, [(1,0,0)], axis=0)

    if len(polytope_vertices) == 1:
        w.ax.scatter3D(*zip(*polytope_vertices), color=color)
    elif len(polytope_vertices) == 3:
        triangle = Poly3DCollection([polytope_vertices])
        triangle.set_facecolor(color)
        triangle.set_edgecolor('k')
        triangle.set_alpha(0.5)
        w.ax.add_collection3d(triangle)
    elif len(polytope_vertices) > 3:
        hull = ss.ConvexHull(polytope_vertices)
        faces = Poly3DCollection([polytope_vertices[simplex] for simplex in hull.simplices])
        faces.set_facecolor(color)
        faces.set_alpha(0.2)
        faces.set_edgecolor('k')
        w.ax.add_collection3d(faces)

def plot_coverage_set(coverage_set, overlap=True):
    """Plot a set of 3D polytopes.
    
    Args:
        coverage_set (list): a list of CircuitPolytope objects.
        overlap (bool): If True, all polytopes are drawn on the same plot. If False, each polytope is drawn in a separate subplot.
    """
    colors = ['red', 'green', 'blue', 'orange', 'purple', 'cyan', 'black', 'pink', 'brown', 'grey']
    coverage_set = coverage_set[1:] # remove the identity

    if overlap:
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection':'3d'})
        w = WeylChamber()
        w.labels = {}
        w.render(ax)
        for i, circuit_polytope in enumerate(coverage_set):
            _plot_polytope(circuit_polytope, color=colors[i % len(colors)], w=w)
    else:
        n = len(coverage_set)
        fig, axs = plt.subplots(1, n, subplot_kw={'projection':'3d'}, figsize=(n*5, 5))  # Adjust size to avoid crowding
        for i, ax in enumerate(axs):
            w = WeylChamber()
            w.labels = {}
            w.render(ax)
            _plot_polytope(coverage_set[i], color=colors[i % len(colors)], w=w)

    plt.show()

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
