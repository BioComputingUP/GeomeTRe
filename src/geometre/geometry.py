import numpy as np
import logging
from numpy.linalg import norm
from scipy.spatial.transform import Rotation
from sklearn.decomposition import PCA
from skimage.measure import CircleModel
from scipy.optimize import minimize

import tmtools

logger = logging.getLogger(__name__)


def orthogonalize(v, n):
    """Orthogonalize vector v to vector n."""
    n_orth = n / norm(n)
    v_orth = v - np.dot(v, n_orth) * n_orth
    v_orth /= norm(v_orth)
    return v_orth


def dihedral_angle(v1, v2, n):
    """Calculate the dihedral angle between two vectors in a plane defined by a normal vector."""
    n = n / norm(n)
    v1_to_orth = v1 / norm(v1)
    v2_to_orth = v2 / norm(v2)
    v1_orth = orthogonalize(v1_to_orth, n)
    v2_orth = orthogonalize(v2_to_orth, n)
    dirn = np.cross(v1_orth, v2_orth)
    dirn /= norm(dirn)
    return get_angle(v1_orth, v2_orth), np.sign(np.dot(n, dirn))


def get_angle(v1, v2):
    """Calculate the angle between two vectors."""
    return np.arccos(np.dot(v1, v2) / (norm(v1) * norm(v2)))


def widest_circle(c, data):
    """Find the widest circular crown within units."""
    nearest = -np.inf
    farthest = np.inf
    for unit in data:
        distances = [norm(ca - c) for ca in unit]
        is_farthest = max(distances)
        if is_farthest < farthest:
            farthest = is_farthest
        is_nearest = min(distances)
        if is_nearest > nearest:
            nearest = is_nearest
    return nearest - farthest  # Minimize the opposite of the width


def widest_circle_fit(units, centers, window=6):
    """Calculate curvature of repeat units"""
    num_units = len(units)
    index_list = []
    centers_list = []
    score_list = []

    for i in range(max(num_units - window + 1, 1)):
        min_index = i
        max_index = min(i + window, num_units)
        data_to_fit = units[min_index:max_index]
        pca_centers = centers[min_index:max_index]

        # Find plane of rotation of units, and project them onto it
        pca = PCA(n_components=2)
        pca.fit(pca_centers)

        # Project centroids along the plane defined by the first 2 PCA components
        pca_centers = pca.transform(pca_centers)

        # Project all point of the units
        data_transformed = []
        for unit in data_to_fit:
            data_transformed.append(pca.transform(unit))

        # Find the widest crown in the 2D plane starting from the centers of the centroids but
        # considering all the projected coordinates of the C-alpha
        circle = CircleModel()
        circle.estimate(pca_centers)
        res = minimize(widest_circle, circle.params[0:2], args=(data_transformed))

        center = res.x
        centers_list.append(pca.inverse_transform(center))
        index_list.append([*range(min_index, max_index)])

        # Fun is the value of the objective function at x
        # Here we are changing the sign for keeping the smallest crown
        score_list.append(-res.fun)

    # Select the widest crown
    def_centers = np.empty((num_units - 1, 3))
    best_score = np.full(num_units - 1, np.inf)
    for center, indexes, score in zip(centers_list, index_list, score_list):  # For each unit pair, select center corresponding to the widest crown
        act_indexes = indexes[:-1]
        score_to_confront = best_score[act_indexes]
        centers_to_confront = def_centers[act_indexes]
        mask = score_to_confront > score
        centers_to_confront[mask] = np.vstack([center for _ in range(len(centers_to_confront[mask]))])
        score_to_confront[mask] = score
        best_score[act_indexes] = score_to_confront
        def_centers[act_indexes] = centers_to_confront

    logger.debug("Widest circle fit completed.")
    return def_centers


def build_ref_axes(geometric_centers, rot_centers):
    """Build reference axes for pitch and twist."""
    num_centers = len(geometric_centers)
    pitch_axis = []
    twist_axis = []
    rots = []

    # Pitch vectors
    for i in range(num_centers - 1):
        vec_1 = rot_centers[i] - geometric_centers[i]
        vec_1 /= norm(vec_1)  # Obtain the versor instead of the vector
        vec_2 = rot_centers[i] - geometric_centers[i + 1]
        vec_2 /= norm(vec_2)
        pitch_axis.append((vec_1, vec_2))

    # Twist vectors
    # Find the tangent of the rotational circle passing through the two centroids
    for i in range(num_centers - 1):
        twist_vect = geometric_centers[i + 1] - geometric_centers[i]
        vec_1 = orthogonalize(twist_vect, pitch_axis[i][0])
        vec_2 = orthogonalize(twist_vect, pitch_axis[i][1])
        twist_axis.append((vec_1, vec_2))

    # Obtain the curvature (yaw) axis with cross product between the twist and pitch axis
    for i in range(num_centers - 1):
        try:
            rots.append(
                (
                    Rotation.from_matrix(
                        np.array(
                            [
                                twist_axis[i][0],
                                pitch_axis[i][0],
                                np.cross(twist_axis[i][0], pitch_axis[i][0]),
                            ]
                        )
                    ),
                    Rotation.from_matrix(
                        np.array(
                            [
                                twist_axis[i][1],
                                pitch_axis[i][1],
                                np.cross(twist_axis[i][1], pitch_axis[i][1]),
                            ]
                        )
                    ),
                )
            )
        except Exception as e:
            logger.warning(f"Error constructing rotation matrix at index {i}: {e}")
            rots.append(Rotation.from_matrix(np.eye(3)))
    return pitch_axis, twist_axis, rots


def get_unit_rotation(coords, seqs, rotations):
    """Align 2 units using TM-align and return rotation."""
    coords_1 = coords[0]
    coords_2 = coords[1]
    coords_1 = rotations[0].apply(coords_1)
    coords_2 = rotations[1].apply(coords_2)
    alignment = tmtools.tm_align(coords_1, coords_2, seqs[0], seqs[1])
    return alignment
