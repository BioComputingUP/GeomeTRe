import os
from single_processing import geometre
from batch_processing import batch_repeats_geometry


class GeomeTRe:
    def __init__(self, draw_enabled=False):
        """
        Initialize the GeomeTRe API.

        Args:
            draw_enabled (bool): Whether to enable PyMOL drawing during processing.
        """
        self.draw_enabled = draw_enabled

    def single(self, pdb_filepath, chain, units_def, insertion="", output_path=""):
        """
        Process a single PDB or CIF file.

        Args:
            pdb_filepath (str): Path to the input PDB or CIF file.
            chain (str): Chain ID to analyze.
            units_def (str): Repeat unit definitions (e.g., "10_50,51_100").
            insertion (str): Optional insertions (default: "").
            output_path (str): Path to save the output CSV file.

        Returns:
            result: DataFrame of computed properties.
            stats: Summary statistics.
        """
        return geometre(
            filepath=pdb_filepath,
            chain=chain,
            units_ids=units_def,
            o_path=output_path,
            ins_ids=insertion,
            draw=self.draw_enabled
        )

    def batch(self, tsv_filepath, output_path, num_threads=4, file_format="cif", pdb_dir=None):
        """
        Process multiple files in batch mode.

        Args:
            tsv_filepath (str): Path to the input TSV file.
            output_path (str): Path to save the output CSV file.
            num_threads (int): Number of threads for parallel processing (default: 4).
            file_format (str): Format of input files ("cif" or "pdb").
            pdb_dir (str): Directory containing local PDB files (required).

        Raises:
            ValueError: If pdb_dir is not provided or does not exist.
        """
        if not pdb_dir:
            raise ValueError("The 'pdb_dir' argument is required for batch processing.")

        if not os.path.exists(pdb_dir):
            raise ValueError(f"The specified pdb_dir '{pdb_dir}' does not exist.")

        return batch_repeats_geometry(
            tsv_path=tsv_filepath,
            output_path=output_path,
            num_threads=num_threads,
            file_format=file_format,
            pdb_dir=pdb_dir
        )

