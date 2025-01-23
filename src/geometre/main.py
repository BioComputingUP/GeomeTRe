from argparse import ArgumentParser, FileType
import logging
from geometre.wrapper import GeomeTRe


def main():
    arg_parser = ArgumentParser(description="Calculate repeat protein geometrical properties.")
    subparsers = arg_parser.add_subparsers(dest='mode')

    # Single file mode
    single_parser = subparsers.add_parser('single', help='Process a single PDB/CIF file')
    single_parser.add_argument('filepath', type=FileType("rt", encoding="UTF-8"), help='Path to input PDB or CIF file')
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
    arg_parser.add_argument('-l', help='Log file path (default: ./process.log)', default='./process.log')
    args = arg_parser.parse_args()

    # Initialize GeomeTRe
    geometre = GeomeTRe(draw_enabled=args.draw if 'draw' in args else False)

    if args.mode == 'single':
        geometre.single(
            pdb_filepath=args.filepath.name,
            chain=args.chain,
            units_def=args.unit_def,
            output_path=args.o,
            insertion=args.ins
        )
        logging.info(f"Single mode results saved to: {args.o}")

    elif args.mode == 'batch':
        # Set up logging to save logs in a file
        log_file = args.l
        logging.basicConfig(
            level=logging.WARNING,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file, mode="w"),
                logging.StreamHandler()
            ]
        )
        logger = logging.getLogger(__name__)
        logger.info("Batch process is started.")
        logging.info(f"Arguments parsed: {args}")
        logging.info(f"Running in batch mode with arguments: {args.batch}, {args.output}, {args.format}, {args.threads}, {args.pdb_dir}")
        geometre.batch(
            tsv_filepath=args.batch,
            output_path=args.output,
            file_format=args.format,
            num_threads=args.threads,
            pdb_dir=args.pdb_dir
        )
        logging.info(f"Batch mode results saved to: {args.output}")

if __name__ == '__main__':
    main()
