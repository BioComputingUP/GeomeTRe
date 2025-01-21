from pathlib import Path
import pandas as pd
import numpy as np
import logging
from Bio.PDB import PDBParser, FastMMCIFParser, Polypeptide
from Bio.SeqUtils import seq1
from scipy.spatial.transform import Rotation

import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Suppress PDBConstructionWarnings
from .geometry import create_list, widest_circle_fit, get_unit_rotation,get_angle,build_ref_axes, pymol_drawing

# Use the shared logger
logger = logging.getLogger(__name__)

def geometre(filepath, chain, units_ids, o_path, ins_ids=None, draw=False):
    """Calculate geometrical parameters for repeats."""
    logging.info(f"Processing file: {filepath}, chain: {chain}")

    units_ids = create_list(units_ids)
    file_type = Path(filepath).suffix.lower()

    if file_type == '.cif':
        parser = FastMMCIFParser(QUIET=True)
        structure = parser.get_structure('structure', Path(filepath))
    elif file_type == '.pdb':
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', Path(filepath))
    else:
        raise ValueError(f"Unsupported file type: {file_type}. Provide a '.pdb' or '.cif' file.")


    # Ensure structure is loaded
    if structure is None:
        logging.error("Structure could not be initialized. Check file type and content.")
        return None, None

    chain_s = structure[0][chain]

    # Handle insertions
    ins_ids = []
    if ins_ids:
        ins_ids = create_list(ins_ids)

    # If insertions are present, we make sure to remove them from the structure
    if len(ins_ids) > 0:
        units = []
        to_remove = []
        for limits in units_ids:
            a, b = limits
            unit_ins = [ins for ins in ins_ids if a <= ins[0] <= b or a <= ins[1] <= b]
            unit = []
            for residue in chain_s:
                res_id = residue.get_id()[1]
                res_in_ins = np.any([ins[0] <= res_id <= ins[1] for ins in unit_ins])
                if a <= res_id <= b and not res_in_ins and Polypeptide.is_aa(residue):
                    unit.append(residue)
            if len(unit) > 0:
                units.append(unit)
            else:
                to_remove.append(limits)
        for ids in to_remove:
            units_ids.remove(ids)
    else:
        units = []
        to_remove = []
        for limits in units_ids:
            a, b = limits
            unit = []
            for residue in chain_s:
                if a <= residue.get_id()[1] <= b:
                    if Polypeptide.is_aa(residue):
                        unit.append(residue)
            if len(unit) > 0:
                units.append(unit)
            else:
                to_remove.append(limits)
        for ids in to_remove:
            units_ids.remove(ids)

    # Extract sequences and coordinates
    units_seqs = [seq1(''.join(res.get_resname() for res in unit)) for unit in units]
    units_coords = [[res['CA'].get_coord() for res in unit] for unit in units]

    # Calculate geometrical centers and rotations
    geometric_centers = [np.mean(coords, axis=0) for coords in units_coords]
    num_centers = len(geometric_centers)
    rot_centers = widest_circle_fit(units_coords, geometric_centers)
    logging.info("Geometric centers and rotation centers calculated.")
    rot_angles = [
        get_angle(geometric_centers[i] - rot_centers[i], geometric_centers[i + 1] - rot_centers[i])
        for i in range(num_centers - 1)
    ]
    pitch_axis, twist_axis, rots = build_ref_axes(geometric_centers, rot_centers)
    logging.info("Reference axes built for geometry calculations.")

    # Calculate rotations and scores
    units_rots, tmscores = [], []
    for i in range(num_centers - 1):
        alignment = get_unit_rotation(units_coords[i:i + 2], units_seqs[i:i + 2], rots[i])
        units_rots.append(alignment.u)
        tmscores.append(alignment.tm_norm_chain1)

    # Decompose rotation into pitch, twist, and yaw
    pitchlist, twistlist, twist_handednesslist, pitch_handednesslist, yawlist = [], [], [], [], []
    for i in range(num_centers - 1):
        rotation = units_rots[i]
        twist, pitch, yaw = Rotation.from_matrix(rotation).as_euler('xyz')
        twistlist.append(abs(twist))
        twist_handednesslist.append(np.sign(twist))
        pitchlist.append(abs(pitch))
        pitch_handednesslist.append(np.sign(pitch))
        yawlist.append(abs(yaw))

    # Compute statistics
    stats = [
        np.nanmean(rot_angles), np.nanstd(rot_angles),
        np.nanmean(twistlist), np.nanstd(twistlist),
        np.nanmean(twist_handednesslist), np.nanstd(twist_handednesslist),
        np.nanmean(pitchlist), np.nanstd(pitchlist),
        np.nanmean(pitch_handednesslist), np.nanstd(pitch_handednesslist),
        np.nanmean(tmscores), np.nanstd(tmscores),
        np.nanmean(yawlist), np.nanstd(yawlist)
    ]
    # DataFrame output
    rot_angles.extend(stats[0:2])
    twistlist.extend(stats[2:4])
    twist_handednesslist.extend(stats[4:6])
    pitchlist.extend(stats[6:8])
    pitch_handednesslist.extend(stats[8:10])
    tmscores.extend(stats[10:12])
    yawlist.extend(stats[12:14])

    rot_angles.insert(0,0)
    twistlist.insert(0,0)
    twist_handednesslist.insert(0,0)
    pitch_handednesslist.insert(0,0)
    pitchlist.insert(0,0)
    tmscores.insert(0,0)
    yawlist.insert(0,0)

    # Prepare DataFrame
    num_centers = len(rot_angles) - 2
    starts = [unit[0] for unit in units_ids]
    starts.append('mean')
    starts.append('std deviation')
    ends = [unit[1] for unit in units_ids]
    ends.append('-')
    ends.append('-')

    pdb = filepath.split('/')[-1][:4]
    pdbs = [pdb for _ in range(num_centers + 2)]
    chains = [chain for _ in range(num_centers + 2)]

    d = {
        'pdb_id': pdbs,
        'chain': chains,
        'unit_start': starts,
        'unit_end': ends,
        'curvature': rot_angles,
        'twist': twistlist,
        'twist_hand': twist_handednesslist,
        'pitch': pitchlist,
        'pitch_hand': pitch_handednesslist,
        'TM-score': tmscores,
        'yaw': yawlist
    }
    df = pd.DataFrame(data=d)

    # Save or display output
    output_path = Path(o_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)  # Ensure the parent directory exists
    pd.DataFrame([stats]).to_csv(output_path, index=False)
    logging.info(f"Results saved to {output_path}")

    if not df.empty:
        df.to_csv(output_path, index=False, float_format='%.4f')
        logging.info(f"Output successfully saved to {output_path}")
    else:
        logging.warning(f" The DataFrame is empty and no data to save. No file was created at {output_path}.")

    # Drawing with PyMOL
    if draw:
        if not df.empty:
            # Call the PyMOL drawing function
            pymol_drawing(filepath, geometric_centers, rot_centers, twist_axis, pitch_axis, rots, units_rots,
                          units_coords)

            logger.info(f"PyMOL visualization saved.")
        else:
            logging.error(f"Missing required inputs for PYMOL visualization for file {filepath}, chain: {chain}: {e}")

    return df, stats

