import numpy as np
import logging
import pymol
from pymol import cmd
from sklearn.decomposition import PCA

logger = logging.getLogger(__name__)


def pymol_drawing(filepath, geometric_centers, rot_centers, twist_axis, pitch_axis, rots, units_rots, units_coords, units_ids, chain):
    """Draw geometrical properties using PyMOL."""
    logger.info(f"Drawing geometrical properties for file: {filepath}.")

    # Perform PCA on the first unit to find a principal axis for visualization
    draw_pca = PCA()
    draw_pca.fit(units_coords[0])
    unit_vector = draw_pca.components_[0]

    num_centers = len(geometric_centers)
    pymol.finish_launching(['pymol', '-q'])
    cmd.load(filepath)
    cmd.hide('all')

    # Color the units
    for i, pos in enumerate(units_ids):
        cmd.color("red" if i%2 == 0 else "blue", f"resi {pos[0]}-{pos[1]} and chain {chain}")

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
    for i in range(len(rot_centers)):
        cmd.pseudoatom('rot_centers', pos=tuple(rot_centers[i]))

    # Orient the repeated region
    region_sele = f'resi {units_ids[0][0]}-{units_ids[-1][1]} and chain {chain}'
    cmd.orient(region_sele)
    cmd.show('cartoon', region_sele)

    # Draw rotation angles and protein geometry
    for i in range(num_centers - 1):

        # Length of the vectors
        l = np.sqrt(np.sum((geometric_centers[i+1] - geometric_centers[i])**2, axis=0)) * 1.5


        cmd.pseudoatom('twist_ref', pos=tuple(geometric_centers[i] + l * twist_axis[i][0]))
        cmd.pseudoatom('pitch_ref', pos=tuple(geometric_centers[i] + l * pitch_axis[i][0]))
        cmd.pseudoatom('curvature_ref', pos=tuple(geometric_centers[i] + l * np.cross(pitch_axis[i][0], twist_axis[i][0])))
        if i > 0:
            cmd.select('point1', selection='model geo_centers and name PS{}'.format(str(i + 1)))
            # cmd.select('point2', selection='model geo_centers and name PS{}'.format(str(i + 2)))
            # cmd.select('rot_center', selection='model rot_centers and name PS{}'.format(str(i + 1)))
            cmd.select('twist_point', selection='model twist_ref and name PS{}'.format(str(i + 1)))
            cmd.select('pitch_point', selection='model pitch_ref and name PS{}'.format(str(i + 1)))
            cmd.select('curvature_point', selection='model curvature_ref and name PS{}'.format(str(i + 1)))
            # cmd.angle('rot_angle', selection1='point1', selection2='rot_center', selection3='point2')
            # cmd.distance('superaxis', selection1='point1', selection2='point2')
            cmd.distance('twist_axis', selection1='point1', selection2='twist_point')
            cmd.distance('pitch_axis', selection1='point1', selection2='pitch_point')
            cmd.distance('curvature_axis', selection1='point1', selection2='curvature_point')

    # Place pseudoatoms for the last geometric center
    l = np.sqrt(np.sum((geometric_centers[-1] - geometric_centers[-2]) ** 2, axis=0)) * 1.5
    cmd.pseudoatom('twist_ref', pos=tuple(geometric_centers[-1] + l * twist_axis[-1][1]))
    cmd.pseudoatom('pitch_ref', pos=tuple(geometric_centers[-1] + l * pitch_axis[-1][1]))
    cmd.pseudoatom('curvature_ref', pos=tuple(geometric_centers[-1] + l * np.cross(pitch_axis[-1][1], twist_axis[-1][1])))
    cmd.select('point1', selection='model geo_centers and name PS{}'.format(str(num_centers)))
    # cmd.select('point2', selection='model geo_centers and name PS{}'.format(str(i + 2)))
    # cmd.select('rot_center', selection='model rot_centers and name PS{}'.format(str(i + 1)))
    cmd.select('twist_point', selection='model twist_ref and name PS{}'.format(str(num_centers)))
    cmd.select('pitch_point', selection='model pitch_ref and name PS{}'.format(str(num_centers)))
    cmd.select('curvature_point', selection='model curvature_ref and name PS{}'.format(str(num_centers)))
    # cmd.angle('rot_angle', selection1='point1', selection2='rot_center', selection3='point2')
    # cmd.distance('superaxis', selection1='point1', selection2='point2')
    cmd.distance('twist_axis', selection1='point1', selection2='twist_point')
    cmd.distance('pitch_axis', selection1='point1', selection2='pitch_point')
    cmd.distance('curvature_axis', selection1='point1', selection2='curvature_point')

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
    cmd.color('green', 'pitch_axis')
    cmd.color('blue', 'curvature_axis')
    cmd.hide('labels')



    cmd.deselect()