from __future__ import annotations

from itertools import combinations
from typing import TYPE_CHECKING, List

import numpy as np

from monodromy.coordinates import average_infidelity

if TYPE_CHECKING:
    from monodromy.coverage import CircuitPolytope

from cvxopt import matrix, solvers


def _quadratic_programming_nearest_point(polytope, target):
    """Compute the nearest point to a given target in a specified polytope,
    using quadratic programming.

    Args:
        polytope: The polytope to search.
        target: The point to which we want to find the nearest point in the polytope.

    Returns:
        The point in the polytope that is nearest to the target, according to quadratic programming.
    """

    # Create the P matrix
    P = matrix(np.eye(len(target)))

    # Create the q vector
    q = matrix(-2 * target)

    # Create the G matrix and h vector
    G = matrix(polytope.inequalities)
    h = matrix([eq[-1] for eq in polytope.inequalities])

    # Create the A matrix and b vector
    A = matrix(polytope.equalities)
    b = matrix([eq[-1] for eq in polytope.equalities])

    # Solve the QP problem
    sol = solvers.qp(P, q, G, h, A, b)

    return np.array(sol["x"]).reshape((-1,))


def _affine_subspace_projection(point, inequalities, equalities=[]):
    """Project a point onto an affine subspace defined by an intersection of
    hyperplanes, which are defined by inequalities and equalities.

    Args:
        point: A numpy array representing the point to project.
        inequalities: A list of lists representing the inequalities defining the polytope. Each list's
                      last element is the constant term, and the others are the coefficients of the variables.
        equalities: A list of lists representing the equalities defining the polytope. Each list's
                    last element is the constant term, and the others are the coefficients of the variables.

    Returns:
        The projection of the point onto the affine subspace.
    """
    # Convert inequalities into equalities (ignoring the "greater than or equal to" part)
    # and concatenate them with the given equalities
    equations = [ineq[1:] + [-ineq[0]] for ineq in inequalities] + equalities

    # Form the matrix A from the equations
    A = np.array([eq[:-1] for eq in equations])

    # Here, we compute the 'b' vector as A @ point rather than from the equations
    b = A @ point

    # Compute the projection
    # This involves solving the normal equation: (A^T A) x = A^T b
    x = np.linalg.solve(A.T @ A, A.T @ b)

    return x


def _nearest(circuit_polytope, target, method="affine"):
    """Finds the nearest point to the target within a CircuitPolytope, which is
    defined as the union of multiple ConvexPolytopes.

    Args:
        circuit_polytope: The CircuitPolytope object.
        target: A list representing the target point.
        method: A string representing the projection method ('affine' or 'quadratic').

    Returns:
        The point within the CircuitPolytope that is nearest to the target point.
    """
    nearest_points = []

    # Iterate over all convex polytopes in the circuit polytope
    for polytope in circuit_polytope.convex_subpolytopes:
        if method == "affine":
            valid_projections = []

            # Iterate over all subsets of the inequalities defining the convex polytope
            for inequalities in combinations(
                polytope.inequalities, len(polytope.inequalities)
            ):
                # Calculate the projection of the target point onto the subspace defined by this subset
                projection = _affine_subspace_projection(
                    target, inequalities, polytope.equalities
                )

                # If the projection is in the convex polytope, it's a valid candidate for the nearest point
                if polytope.has_element(projection):
                    valid_projections.append(projection)

            # Find the valid projection that's closest to the target point in this convex polytope
            if valid_projections:
                distances = [
                    np.linalg.norm(projection - target)
                    for projection in valid_projections
                ]
                nearest_points.append(valid_projections[np.argmin(distances)])

        elif method == "quadratic":
            projection = _quadratic_programming_nearest_point(polytope, target)
            # Add a check for equalities as well, to make sure the solution satisfies all the constraints
            if polytope.has_element(projection) and all(
                np.isclose(np.dot(eq[:-1], projection), eq[-1])
                for eq in polytope.equalities
            ):
                nearest_points.append(projection)

    # Among the nearest points in each convex polytope, find the one that's closest to the target
    if nearest_points:
        overall_distances = [np.linalg.norm(point - target) for point in nearest_points]
        return nearest_points[np.argmin(overall_distances)]
    else:
        return None


def polytope_approx_contains(
    polytope: CircuitPolytope,
    point: List[float],
    approximation_degree=0.0,
    method="affine",
):
    """Checks if a polytope contains a point, with an optional degree of
    approximation.

    Args:
        polytope: The CircuitPolytope object.
        point: A numpy array representing the point to check.
        approximation_degree: The allowed degree of approximation.
        method: Method for computing nearest point. Options are "affine" (default) and "quadratic".

    Returns:
        True if the polytope contains the point or the infidelity between the point and the polytope
        is less than or equal to the approximation_degree, False otherwise.
    """
    # If approximation_degree is 0.0, then we do an exact check
    if approximation_degree == 0.0:
        return polytope.has_element(point)

    # find the point in the polytope that minimizes the distance -> maximizes fidelity
    nearest = _nearest(polytope, point, method)
    # Compute infidelity between the nearest point in the polytope and the point
    # If the infidelity is less than or equal to approximation_degree, then we say that the point is in the polytope
    return average_infidelity(nearest, point) <= approximation_degree
