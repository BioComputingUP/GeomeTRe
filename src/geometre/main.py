import argparse
import logging
from .repeats_geometry import repeats_geometry
from .batch_processing import batch_repeats_geometry
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Suppress PDBConstructionWarnings
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

# Set up logging to save logs in a file
logging.basicConfig(
    level=logging.WARNING,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("process.log",mode="a"),
    ]
)
logger = logging.getLogger(__name__)
logger.info("Batch process is started.")

def main():
    arg_parser = argparse.ArgumentParser(description="Calculate repeat protein geometrical properties.")
    subparsers = arg_parser.add_subparsers(dest='mode', required=True)

    # Single file mode
    single_parser = subparsers.add_parser('single', help='Process a single PDB/CIF file')
    single_parser.add_argument('filepath', help='Path to input PDB or CIF file')
    single_parser.add_argument('chain', help='Chain ID')
    single_parser.add_argument('unit_def', help='Unit limits (e.g., 10_50,51_100)')
    single_parser.add_argument('-ins', default='', help='Insertions (optional)')
    single_parser.add_argument('-o', default='', help='Output CSV file path')
    single_parser.add_argument('--draw', action='store_true', help='Generate PyMOL visualization')

    # Batch mode
    batch_parser = subparsers.add_parser('batch', help='Batch process files from a TSV input')
    batch_parser.add_argument('--batch', required=True, help='Input TSV file for batch processing')
    batch_parser.add_argument('--output', required=True, help='Output CSV file for results')
    batch_parser.add_argument('--format', choices=['cif', 'pdb'], default='cif',
                               help='Choose file format for downloading (default: cif)')
    batch_parser.add_argument('--threads', type=int, default=4, help='Number of threads for parallel processing')
    batch_parser.add_argument('--pdb_dir', type=str, help='Directory containing local PDB files (supports .gz files).')

    args = arg_parser.parse_args()

    logging.info(f"Arguments parsed: {args}")

    if args.mode == 'single':
        logging.info(f"Running in single mode with arguments: {args.filepath}, {args.chain}, {args.unit_def}, {args.ins}, {args.o}, {args.draw}")
        repeats_geometry(
            filepath=args.filepath,
            chain=args.chain,
            units_ids=args.unit_def,
            ins_ids=args.ins,
            o_path=args.o,
            draw=args.draw,
            batch=False
        )
        logging.info(f"Single mode results saved to: {args.o}")
    elif args.mode == 'batch':
        logging.info(f"Running in batch mode with arguments: {args.batch}, {args.output}, {args.format}, {args.threads}, {args.pdb_dir}")
        batch_repeats_geometry(
            tsv_path=args.batch,
            output_path=args.output,
            file_format=args.format,
            num_threads=args.threads,
            pdb_dir=args.pdb_dir
        )
        logging.info(f"Batch mode results saved to: {args.output}")

if __name__ == '__main__':
    main()
