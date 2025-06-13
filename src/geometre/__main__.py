#!/usr/bin/python

from argparse import ArgumentParser
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBList
from concurrent.futures import ThreadPoolExecutor, as_completed

from geometre.process import compute


def command_line():
    arg_parser = ArgumentParser(description="Calculate repeat protein geometrical properties.")
    subparsers = arg_parser.add_subparsers(dest='mode')

    # Single file mode
    single_parser = subparsers.add_parser('single', help='Process a single PDB/CIF file')
    single_parser.add_argument('filepath', help='Input structure file')
    single_parser.add_argument('chain', help='Chain ID')
    single_parser.add_argument('out_file', help='Output CSV file path')
    single_parser.add_argument('unit_def', help='Units (10_50,51_100,...)')
    single_parser.add_argument('-ins_def', help='Insertions (same as units)')

    # Batch mode
    batch_parser = subparsers.add_parser('batch', help='Batch process files from a TSV input')
    batch_parser.add_argument('filepath', help='TSV file with list of structures and units')
    batch_parser.add_argument('out_file', help='Output CSV file for results')
    batch_parser.add_argument('-pdb_dir', help='PDB files directory for download. '
                                                         'If not provided assume file path in the input TSV')
    batch_parser.add_argument('-threads', type=int, default=1,
                              help='Number of threads for parallel processing')

    drawing_parser = subparsers.add_parser('pymol', help='Generate PyMOL visualization from saved .npy data.')
    drawing_parser.add_argument("pdb_filepath", help="Path to the input PDB file.")
    drawing_parser.add_argument("npy_filepath", help="Path to the .npy file containing geometry data.")

    # Logging setup
    arg_parser.add_argument('-l', help='Log file path. Default standard error')
    return arg_parser.parse_args()


def batch_compute(tsv_path, pdb_dir, output_path, threads=4):
    df_input = pd.read_csv(tsv_path, sep="\t",
                           dtype={"pdb_file": str, "chain": str, "units": str, "insertion": str})

    if pdb_dir is not None:
        # Need to download structures, pdb_file column contains only identifiers
        pdbl = PDBList()
        pdbl.download_pdb_files([ele.split("/")[-1][:4] for ele in df_input['pdb_file'].unique()], pdir=pdb_dir)  # Download mmCIF by default

    # Iterate over all structures
    df_list = []
    columns = ['pdb_id', 'chain', 'region_start', 'region_end', 'curvature', 'twist', 'twist_hand', 'pitch', 'pitch_hand', 'tmscore', 'yaw']

    with ThreadPoolExecutor(max_workers=threads) as executor:
        fs = {}
        for i, row in df_input.iterrows():
            structure_file = row["pdb_file"] if pdb_dir is None else "{}/{}.cif".format(pdb_dir, row["pdb_file"].split("/")[-1][:4])
            ins_ids = None if pd.isnull(row['insertion']) else row["insertion"]

            fs[executor.submit(compute, structure_file, row["chain"], row["units"], ins_ids)] = structure_file

        for future in as_completed(fs):
            try:
                df_, obj_ = future.result()
            except Exception as exc:
                logger.error(f'Execution error {fs[future]}')
            else:
                if df_ is not None and not df_.empty:
                    df_[['region_start', 'region_end']] = [df_.iloc[0]['unit_start'], df_.iloc[-3]['unit_end']]
                    df_ = pd.merge(df_.loc[df_['unit_start'] == 'mean', columns],
                                   df_.loc[df_['unit_start'] == 'std', columns],
                                   on=['pdb_id', 'chain'], suffixes=('_mean', '_std')).round(3)
                    df_list.append(df_)
                else:
                    logger.info(f'Empty output {fs[future]}')

    pd.concat(df_list).to_csv(output_path, index=False, sep=",")
    return


def main():
    """CLI entry point for geometre"""
    args = command_line()

    # Set up logging before processing
    global logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Add stdout handler
    h = logging.StreamHandler()
    h.setFormatter(fmt)
    logger.addHandler(h)

    # Add file handler
    if args.l:
        f = logging.FileHandler(args.l, mode='w')
        f.setFormatter(fmt)
        logger.addHandler(f)

    logger.debug("Process started.")

    # Single mode execution
    if args.mode == 'single':
        df, obj = compute(filepath=args.filepath,
                     chain=args.chain,
                     units_ids=args.unit_def,
                     ins_ids=args.ins_def)

        # Save or display output
        if df is not None:
            output_path = Path(args.out_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_path, index=False, float_format='%.4f')
            logger.debug(f"Results saved to {output_path}")

            # Save PyMOL data separately for future visualization
            pymol_output_path = output_path.with_suffix('.npy')
            np.save(pymol_output_path, obj)
            logger.debug(f"PyMOL-compatible data saved to {pymol_output_path}")

        # logger.info(f"Single mode results saved to: {result_path}")

    elif args.mode == 'batch':
        logger.info(f"Running batch mode with arguments: {args.mode}, {args.filepath}, {args.out_file}, {args.threads}, {args.pdb_dir}")
        batch_compute(args.filepath,
                      args.pdb_dir,
                      args.out_file,
                      threads=args.threads)
        # logger.info(f"Batch mode results saved to: {result_path}")

    # Load saved PyMOL data and launch pymol
    elif args.mode == 'pymol':
        from geometre.draw import pymol_drawing
        pymol_data = np.load(args.npy_filepath, allow_pickle=True).item()
        pymol_drawing(args.pdb_filepath, **pymol_data)

if __name__ == "__main__":
    main()
