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
    logger.debug(f"Orthogonalizing vector {v} to {n}.")
    n_orth = n / norm(n)
    v_orth = v - np.dot(v, n_orth) * n_orth
    v_orth /= norm(v_orth)
    return v_orth


def dihedral_angle(v1, v2, n):
    """Calculate the dihedral angle between two vectors in a plane defined by a normal vector."""
    logger.debug(f"Calculating dihedral angle between {v1} and {v2} with normal {n}.")
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
    logger.debug(f"Calculating angle between vectors {v1} and {v2}.")
    return np.arccos(np.dot(v1, v2) / (norm(v1) * norm(v2)))


def create_list(list_indexes):
    """Used to process the unit_def argument."""
    logger.debug(f"Creating list from indexes: {list_indexes}.")
    list_indexes = list_indexes.split(',')
    list_indexes = list(dict.fromkeys(list_indexes))
    len_unit = []
    for item in list_indexes:
        item = item.strip().split('_')
        item = (int(item[0]), int(item[1]))
        if len(len_unit) > 0:
            if item[0] > len_unit[-1][1]:
                len_unit.append(item)
        else:
            len_unit.append(item)
    logger.debug(f"Processed list indexes: {len_unit}.")
    return len_unit


def widest_circle(c, data):
    """Find the widest circular crown within units."""
    logger.debug(f"Finding widest circle for center {c} and data set.")
    nearest = np.NINF
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
    """Alternative method for curvature."""
    logger.debug("Starting widest circle fit.")
    num_units = len(units)
    index_list = []
    centers_list = []
    score_list = []

    for i in range(max(num_units - window + 1, 1)):
        min_index = i
        max_index = min(i + window, num_units)
        data_to_fit = units[min_index:max_index]
        pca_centers = centers[min_index:max_index]

        pca = PCA(n_components=2)  # Find plane of rotation of units, and project them onto it
        pca.fit(pca_centers)
        pca_centers = pca.transform(pca_centers)

        data_transformed = []
        for unit in data_to_fit:
            data_transformed.append(pca.transform(unit))

        circle = CircleModel()
        circle.estimate(pca_centers)
        res = minimize(widest_circle, circle.params[0:2], args=(data_transformed))  # Find widest crown in the 2D plane

        center = res.x
        centers_list.append(pca.inverse_transform(center))
        index_list.append([*range(min_index, max_index)])
        score_list.append(np.std([norm(center - geo_center) for geo_center in pca_centers]))

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
    logger.debug("Building reference axes for pitch and twist.")
    num_centers = len(geometric_centers)
    pitch_axis = []
    twist_axis = []
    rots = []

    for i in range(num_centers - 1):
        vec_1 = rot_centers[i] - geometric_centers[i]
        vec_1 /= norm(vec_1)
        vec_2 = rot_centers[i] - geometric_centers[i + 1]
        vec_2 /= norm(vec_2)
        pitch_axis.append((vec_1, vec_2))

    for i in range(num_centers - 1):
        twist_vect = geometric_centers[i + 1] - geometric_centers[i]
        vec_1 = orthogonalize(twist_vect, pitch_axis[i][0])
        vec_2 = orthogonalize(twist_vect, pitch_axis[i][1])
        twist_axis.append((vec_1, vec_2))

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
    logger.debug("Reference axes built.")
    return pitch_axis, twist_axis, rots


def get_unit_rotation(coords, seqs, rotations):
    """Align 2 units using CEAlign and return rotation."""
    logger.debug("Aligning units using CEAlign.")
    coords_1 = coords[0]
    coords_2 = coords[1]
    coords_1 = rotations[0].apply(coords_1)
    coords_2 = rotations[1].apply(coords_2)
    alignment = tmtools.tm_align(coords_1, coords_2, seqs[0], seqs[1])
    logger.debug("Alignment completed.")
    return alignment

def pymol_drawing(filepath, geometric_centers, rot_centers, twist_axis, pitch_axis, rots, units_rots, units_coords):

    import pymol
    from pymol import cmd

    from sklearn.decomposition import PCA

    """Draw geometrical properties using PyMOL."""
    logger.info(f"Drawing geometrical properties for file: {filepath}.")

    # Perform PCA on the first unit to find a principal axis for visualization
    draw_pca = PCA()
    draw_pca.fit(units_coords[0])
    unit_vector = draw_pca.components_[0]

    num_centers = len(geometric_centers)
    pymol.finish_launching()
    cmd.load(filepath)
    cmd.hide('all')

    # Place pseudoatoms to draw distances and angles
    for i in range(num_centers):
        cmd.pseudoatom('geo_centers', pos=tuple(geometric_centers[i]))
        # Uncommented lines preserved from original function:
        # cmd.pseudoatom('ref_1', pos=tuple(geometric_centers[i] + 6 * unit_vector))
        # cmd.pseudoatom('ref_2', pos=tuple(geometric_centers[i] - 6 * unit_vector))
        # cmd.select('unit_1', selection='model ref_1 and name PS{}'.format(str(i + 1)))
        # cmd.select('unit_2', selection='model ref_2 and name PS{}'.format(str(i + 1)))
        # cmd.distance('unit_vector', selection1='unit_1', selection2='unit_2')

        if i < num_centers - 1:
            unit_vector = units_rots[i] @ (rots[i][0].apply(unit_vector))
            unit_vector = rots[i][1].apply(unit_vector, inverse=True)

    # Uncommented lines preserved from original function:
    # for i in range(len(rot_centers)):
    #     cmd.pseudoatom('rot_centers', pos=tuple(rot_centers[i]))

    # Draw rotation angles and protein geometry
    for i in range(num_centers - 1):
        cmd.pseudoatom('twist_ref', pos=tuple(geometric_centers[i] + 9 * twist_axis[i][0]))
        cmd.pseudoatom('pitch_ref', pos=tuple(geometric_centers[i] + 9 * pitch_axis[i][0]))
        cmd.pseudoatom('yaw_ref', pos=tuple(geometric_centers[i] + 9 * np.cross(pitch_axis[i][0], twist_axis[i][0])))
        cmd.select('point1', selection='model geo_centers and name PS{}'.format(str(i + 1)))
        # cmd.select('point2', selection='model geo_centers and name PS{}'.format(str(i + 2)))
        # cmd.select('rot_center', selection='model rot_centers and name PS{}'.format(str(i + 1)))
        cmd.select('twist_point', selection='model twist_ref and name PS{}'.format(str(i + 1)))
        cmd.select('pitch_point', selection='model pitch_ref and name PS{}'.format(str(i + 1)))
        cmd.select('yaw_point', selection='model yaw_ref and name PS{}'.format(str(i + 1)))
        # cmd.angle('rot_angle', selection1='point1', selection2='rot_center', selection3='point2')
        # cmd.distance('superaxis', selection1='point1', selection2='point2')
        cmd.distance('twist_axis', selection1='point1', selection2='twist_point')
        cmd.distance('pitch_axis', selection1='point1', selection2='pitch_point')
        cmd.distance('yaw_axis', selection1='point1', selection2='yaw_point')

    # Place pseudoatoms for the last geometric center
    cmd.pseudoatom('twist_ref', pos=tuple(geometric_centers[-1] + 6 * twist_axis[-1][1]))
    cmd.pseudoatom('pitch_ref', pos=tuple(geometric_centers[-1] + 6 * pitch_axis[-1][1]))
    cmd.pseudoatom('yaw_ref', pos=tuple(geometric_centers[-1] + 6 * np.cross(pitch_axis[-1][1], twist_axis[-1][1])))
    cmd.select('point1', selection='model geo_centers and name PS{}'.format(str(num_centers)))
    # cmd.select('point2', selection='model geo_centers and name PS{}'.format(str(i + 2)))
    # cmd.select('rot_center', selection='model rot_centers and name PS{}'.format(str(i + 1)))
    cmd.select('twist_point', selection='model twist_ref and name PS{}'.format(str(num_centers)))
    cmd.select('pitch_point', selection='model pitch_ref and name PS{}'.format(str(num_centers)))
    cmd.select('yaw_point', selection='model yaw_ref and name PS{}'.format(str(num_centers)))
    # cmd.angle('rot_angle', selection1='point1', selection2='rot_center', selection3='point2')
    # cmd.distance('superaxis', selection1='point1', selection2='point2')
    cmd.distance('twist_axis', selection1='point1', selection2='twist_point')
    cmd.distance('pitch_axis', selection1='point1', selection2='pitch_point')
    cmd.distance('yaw_axis', selection1='point1', selection2='yaw_point')

    # Uncommented lines preserved from original function:
    # cmd.select('point1_1', selection='model ref_1 and name PS{}'.format(str(i + 1)))
    # cmd.select('point1_2', selection='model ref_2 and name PS{}'.format(str(i + 1)))
    # cmd.select('point2_1', selection='model ref_1 and name PS{}'.format(str(i + 2)))
    # cmd.select('point2_2', selection='model ref_2 and name PS{}'.format(str(i + 2)))
    # cmd.distance('dist_1', selection1='point1_1', selection2='point2_1')
    # cmd.distance('dist_2', selection1='point1_2', selection2='point2_2')

    # cmd.color('orange', 'unit_vector')
    # cmd.color('red', 'dist_1')
    # cmd.color('red', 'dist_2')
    cmd.color('red', 'twist_axis')
    cmd.color('orange', 'pitch_axis')
    cmd.color('white', 'yaw_axis')
    cmd.hide('labels')
    cmd.deselect()
