from pathlib import Path
import pandas as pd
import numpy as np
import logging
from Bio.PDB import PDBParser, FastMMCIFParser, Polypeptide
from Bio.SeqUtils import seq1
from scipy.spatial.transform import Rotation

from geometry import widest_circle_fit, get_unit_rotation, get_angle, build_ref_axes

# Suppress PDBConstructionWarnings
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Use the shared logger
logger = logging.getLogger(__name__)

def compute(filepath, chain, units_ids, ins_ids=None):
    """Calculate geometrical parameters for repeats."""
    logger.info(f"Processing file: {filepath}, chain: {chain}")

    units_ids = [int(el) for ele in units_ids.split(',') for el in ele.split('_')]
    units_ids = [units_ids[i:i+2] for i in range(0,len(units_ids),2)]

    file_type = Path(filepath).suffix.lower()

    if file_type == '.cif':
        parser = FastMMCIFParser(QUIET=True)
    elif file_type == '.pdb':
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported file type: {file_type}. Provide a '.pdb' or '.cif' file.")

    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            structure = parser.get_structure('structure', Path(filepath))
    except:
        logger.error(f"Structure parse error {filepath}")
        return None, None

    # Ensure structure is loaded
    if structure is None:
        logger.error(f"Empty structure error {filepath}")
        return None, None

    chain_s = structure[0][chain]

    # If insertions are present, we make sure to remove them from the structure
    if ins_ids is not None:
        ins_ids = [int(el) for ele in ins_ids.split(',') for el in ele.split('_')]
        ins_ids = [ins_ids[i:i + 2] for i in range(0, len(ins_ids), 2)]

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

    # Find the center of the circle
    rot_centers = widest_circle_fit(units_coords, geometric_centers)
    logger.debug("Geometric centers and rotation centers calculated.")

    # Calculate rotation angle (yaw angle) for each pair of units
    rot_angles = [
        get_angle(geometric_centers[i] - rot_centers[i], geometric_centers[i + 1] - rot_centers[i])
        for i in range(num_centers - 1)
    ]

    # Calculate the axes
    pitch_axis, twist_axis, rots = build_ref_axes(geometric_centers, rot_centers)
    logger.debug("Reference axes built for geometry calculations.")

    # Calculate rotations and scores
    units_rots, tmscores = [], []
    for i in range(num_centers - 1):
        # TM-align
        # We provide rots to superimpose the reference systems of the two consecutive units
        alignment = get_unit_rotation(units_coords[i:i + 2], units_seqs[i:i + 2], rots[i])
        units_rots.append(alignment.u)  # Rotation matrix
        tmscores.append(alignment.tm_norm_chain1)  # TM-score

    # Decompose rotation into pitch, twist, and yaw
    pitchlist, twistlist, twist_handednesslist, pitch_handednesslist, yawlist = [], [], [], [], []
    for i in range(num_centers - 1):
        rotation = units_rots[i]
        # The perfect units should have a yaw of zero here, because they start from the same reference system.
        # It is used here a sanity check
        twist, pitch, yaw = Rotation.from_matrix(rotation).as_euler('xyz')
        twistlist.append(abs(twist))
        twist_handednesslist.append(np.sign(twist))
        pitchlist.append(abs(pitch))
        pitch_handednesslist.append(np.sign(pitch))
        yawlist.append(abs(yaw))

    # Compute statistics
    # Curvature, twist, twist_h, pitch, puitch_h, yaw (curvature sanity check)
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

    obj = {"geometric_centers": np.array(geometric_centers),
        "rot_centers": np.array(rot_centers),
        "twist_axis": np.array(twist_axis),
        "pitch_axis": np.array(pitch_axis),
        "rots": np.array(rots),
        "units_rots": np.array(units_rots),
        "units_coords": np.array(units_coords, dtype="object"),
    }

    df = pd.DataFrame(data={'pdb_id': pdbs,
        'chain': chains,
        'unit_start': starts,
        'unit_end': ends,
        'curvature': rot_angles,
        'twist': twistlist,
        'twist_hand': twist_handednesslist,
        'pitch': pitchlist,
        'pitch_hand': pitch_handednesslist,
        'tmscore': tmscores,
        'yaw': yawlist
    })

    return df, obj
