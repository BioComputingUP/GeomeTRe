import numpy as np
import logging
import argparse
import pymol
from pymol import cmd
from sklearn.decomposition import PCA

logger = logging.getLogger(__name__)

def pymol_drawing(filepath, geometric_centers, rot_centers, twist_axis, pitch_axis, rots, units_rots, units_coords):

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


#Command-line Support
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate PyMOL visualization from saved .npy data.")
    parser.add_argument("--input", required=True, help="Path to the input PDB file.")
    parser.add_argument("--data", required=True, help="Path to the .npy file containing geometry data.")

    args = parser.parse_args()
    
    # Load saved PyMOL data
    try:
        pymol_data = np.load(args.data, allow_pickle=True).item()
    except Exception as e:
        logger.error(f"Failed to load PyMOL data: {e}")
        exit(1)

    pymol_drawing(args.input, **pymol_data)

#run command
#python pymol_drawing.py --input my_structure.pdb --data saved_data.npy
