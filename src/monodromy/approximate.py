"""TODO: implement approximate decomposition, e.g. instead of satisfying has_element,
the target gate just needs to be sufficiently close to the polytope.
find the point in the polytope that minimizes the distance -> maximizes fidelity.
"""
from itertools import combinations

import numpy as np

from monodromy.coordinates import average_infidelity


def _affine_subspace_projection(point, equations):
    """Project a point onto an affine subspace defined by an intersection of
    hyperplanes.

    Args:
        point: A numpy array representing the point to project.
        equations: A list of lists representing the equations of the hyperplanes. Each list's last
                   element is the constant term, and the others are the coefficients of the variables.

    Returns:
        The projection of the point onto the affine subspace.
    """
    # Form the matrix A and the vector b from the equations
    A = np.array([eq[:-1] for eq in equations])
    b = np.array([eq[-1] for eq in equations])

    # Compute the projection
    # This involves solving the normal equation: (A^T A) x = A^T b
    x = np.linalg.solve(A.T @ A, A.T @ b)

    return x


def _nearest(polytope, target):
    """Finds the nearest point to the target within a polytope.

    Args:
        polytope: A list of lists representing the inequalities defining the polytope. Each list
                  represents an equation of a facet of the polytope, with the last element being
                  the constant term and the others being the coefficients of the variables.
        target: A numpy array representing the target point.

    Returns:
        The point within the polytope that is nearest to the target point.
    """
    valid_projections = []

    # Iterate over all subsets of the equations defining the polytope
    for i in range(1, len(polytope) + 1):
        for subset in combinations(polytope, i):
            # Calculate the projection of the target point onto the subspace defined by this subset
            projection = _affine_subspace_projection(target, subset)

            # If the projection is in the polytope, it's a valid candidate for the nearest point
            if polytope.has_element(projection):
                valid_projections.append(projection)

    # Find the valid projection that's closest to the target point
    distances = [
        np.linalg.norm(projection - target) for projection in valid_projections
    ]
    nearest_point = valid_projections[np.argmin(distances)]

    return nearest_point


def polytope_approx_contains(polytope, point, approximation_degree=0.0):
    """Checks if a polytope contains a point.

    Args:
        polytope: A list of lists representing the inequalities defining the polytope. Each list
                  represents an equation of a facet of the polytope, with the last element being
                  the constant term and the others being the coefficients of the variables.
        point: A numpy array representing the point to check.

    Returns:
        True if the polytope contains the point, False otherwise.
    """

    if approximation_degree == 0.0:
        return polytope.has_element(point)

    # find the point in the polytope that minimizes the distance -> maximizes fidelity
    nearest = _nearest(polytope, point)
    if average_infidelity(nearest, point) <= approximation_degree:
        return True
    else:
        return False
