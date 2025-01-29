import os
import logging
import pandas as pd
from .single_processing import geometre
from .batch_processing import batch_repeats_geometry

logger = logging.getLogger(__name__)

class GeomeTRe:
    def __init__(self, draw_enabled=False, output_dir="./results"):
        """
        Initialize the GeomeTRe API.
        Args:
            draw_enabled (bool): Enable PyMOL visualization.
            output_dir (str): Directory to save output files (default: "./results").
        """
        self.draw_enabled = draw_enabled
        self.output_dir = os.path.abspath(output_dir)

        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)
        logger.info(f"Initialized GeomeTRe with output directory: {self.output_dir}")

    def compute(self, mode, **kwargs):
        """
        Method for both single and batch processing.
        Args:
            mode (str): "single" or "batch"
            **kwargs: Parameters for `geometre` (single) or `batch_repeats_geometry` (batch).
        Returns:
            str: Path to the output file.
        """
        if mode not in ["single", "batch"]:
            raise ValueError("Invalid mode. Choose 'single' or 'batch'.")

        # Default output filenames based on mode
        if mode == "single":
            kwargs["filepath"] = kwargs.pop("filepath", None)  # Rename argument
            if "filepath" not in kwargs or kwargs["filepath"] is None:
                raise ValueError("Missing required argument: 'filepath' for single mode.")

            # Automatically generate filename based on PDB ID
            output_path = os.path.join(self.output_dir, f"{os.path.basename(kwargs['filepath'])}_result.csv")
            kwargs["o_path"] = output_path  # Fix argument name

        else:  # Batch mode
            output_path = os.path.join(self.output_dir, "batch_results.csv")
            kwargs["output_path"] = output_path  # Assign default batch output path

        # Ensure the directory exists
        os.makedirs(self.output_dir, exist_ok=True)

        # Select the correct function
        processing_function = geometre if mode == "single" else batch_repeats_geometry
        logger.info(f"Starting {mode} processing: Output -> {output_path}")

        # Call function with updated kwargs
        processing_function(**kwargs)

        return output_path  # Return final saved file path

    def display_results(self, file_path):
        """
        Display processed results from a CSV file.
        Args:
            file_path (str): Path to the output CSV file.
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File '{file_path}' not found.")

        df = pd.read_csv(file_path)
        print(df.head())  # Display only the first few rows
        return df