from argparse import ArgumentParser, FileType
import logging
import numpy as np
import pandas as pd
from Bio.PDB import PDBList

from geometre.process import compute


def batch_compute(tsv_path, pdb_dir, output_path, num_threads=4, file_format="cif", ):
    df_input = pd.read_csv(tsv_path, sep="\t",
                           dtype={"pdb_file": str, "chain": str, "units": str, "insertion": str}).loc[0:10]

    if pdb_dir is not None:
        # Need to download structures, pdb_file column contains only identifiers
        pdbl = PDBList()
        pdbl.download_pdb_files([ele.split("/")[-1][:4] for ele in df_input['pdb_file'].unique()], pdir=pdb_dir)  # Download mmCIF by default

    # Iterate over all structures
    df_list = []
    for i, row in df_input.iterrows():
        structure_file = row["pdb_file"] if pdb_dir is None else "{}/{}.cif".format(pdb_dir, row["pdb_file"].split("/")[-1][:4])
        df, stats = compute(filepath=structure_file, chain= row["chain"], units_ids=row["units"], ins_ids=row["insertion"])
        df_list.append(df)
    df_out = pd.concat(df_list)

    # Combine mean and std on columns
    columns = ['pdb_id', 'chain', 'curvature', 'twist', 'pitch', 'TM-score', 'yaw']
    df_out = pd.merge(df_out.loc[df_out['unit_start'] == 'mean', columns], df_out.loc[df_out['unit_start'] == 'mean', columns],
             on=['pdb_id', 'chain'], suffixes=('_mean', '_std')).round(3)

    df_out.to_csv(output_path, index=False, sep=",")
    return


def main():
    arg_parser = ArgumentParser(description="Calculate repeat protein geometrical properties.")
    subparsers = arg_parser.add_subparsers(dest='mode')

    # Single file mode
    single_parser = subparsers.add_parser('single', help='Process a single PDB/CIF file')
    single_parser.add_argument('filepath', type=str, help='Path to input PDB or CIF file')
    single_parser.add_argument('chain', help='Chain ID')
    single_parser.add_argument('unit_def', help='Unit limits (e.g., 10_50,51_100)')
    single_parser.add_argument('-ins', default='', help='Insertions (optional)')
    single_parser.add_argument('-o', default='output.csv', help='Output CSV file path')

    # Batch mode
    batch_parser = subparsers.add_parser('batch', help='Batch process files from a TSV input')
    batch_parser.add_argument('filepath', help='Input TSV file for batch processing')
    batch_parser.add_argument('-pdb_dir', type=str, help='Directory containing PDB files, existing or to download (supports .gz files)')
    batch_parser.add_argument('-o', help='Output CSV file for results')
    batch_parser.add_argument('-format', choices=['cif', 'pdb'], default='cif',
                               help='Choose file format for downloading (default: cif)')
    batch_parser.add_argument('-threads', type=int, default=4, help='Number of threads for parallel processing')

    drawing_parser = subparsers.add_parser('pymol', help='Generate PyMOL visualization from saved .npy data.')
    drawing_parser.add_argument("pdb_filepath", help="Path to the input PDB file.")
    drawing_parser.add_argument("filepath", help="Path to the .npy file containing geometry data.")

    # Logging setup
    # TODO log should not default to a file but to stderr only
    arg_parser.add_argument('-l', help='Log file path (default: ./process.log)', default='./process.log')
    args = arg_parser.parse_args()

    # Set up logging before processing
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(args.l, mode="w"),
            logging.StreamHandler()
        ]
    )
    logger = logging.getLogger(__name__)
    logger.info("Process started.")

    # Single mode execution
    if args.mode == 'single':
        df, stats = compute(
            filepath=args.filepath,
            chain=args.chain,
            units_ids=args.unit_def,
            ins_ids=args.ins,
            o_path=args.o
        )
        # logger.info(f"Single mode results saved to: {result_path}")

    elif args.mode == 'batch':
        logger.info(f"Running batch mode with arguments: {args.mode}, {args.o}, {args.format}, {args.threads}, {args.pdb_dir}")
        batch_compute(
            args.filepath,
            args.pdb_dir,
            args.o,
            num_threads=args.threads,
            file_format=args.format
        )
        # logger.info(f"Batch mode results saved to: {result_path}")

    # Load saved PyMOL data
    elif args.mode == 'pymol':
        from geometre.draw import pymol_drawing
        pymol_data = np.load(args.filepath, allow_pickle=True).item()
        pymol_drawing(args.pdb_filepath, **pymol_data)

    return


if __name__ == '__main__':
    main()
