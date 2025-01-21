from concurrent.futures import ThreadPoolExecutor
import tempfile
import os
import pandas as pd
from .single_processing import geometre
import requests
import gzip
import shutil
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

def get_structure_file(pdb_id, temp_dir, file_format="cif", pdb_dir=None):
    """
    Retrieve a structure file from a local directory if available, otherwise download it.
    Handle both uncompressed and .gz files.
    """
    try:
        if pdb_dir:
            # Check for both uncompressed and .gz files in the local directory
            possible_extensions = [".pdb", ".pdb.gz", ".ent", ".ent.gz"]
            for ext in possible_extensions: #this is stupid
                local_file = os.path.join(pdb_dir, f"{pdb_id}{ext}")
                if os.path.exists(local_file):
                    logger.info(f"Using local structure file: {local_file}")
                    print(local_file)
                    return local_file

        # If not found locally, download the file
        file_path = os.path.join(temp_dir, f"{pdb_id}.{file_format}")
        url = f"https://files.rcsb.org/download/{pdb_id}.{file_format}"
        logger.debug(f"Attempting to download file from: {url}")
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        with open(file_path, "wb") as f:
            f.write(response.content)
        logger.info(f"Downloaded structure file to: {file_path}")

        # Return the downloaded file path
        return file_path

    except Exception as e:
        logger.error(f"Failed to retrieve structure file for {pdb_id}: {e}")
        raise


def process_entry(row, temp_dir, output_file, file_format, pdb_dir=None, processed_files=None):
    try:
        pdb_file = row["pdb_file"]
        chain = row["chain"]
        units = row["units"]
        insertions = row["insertion"] if pd.notna(row["insertion"]) else ""

        # Track processed files
        if processed_files is None:
            processed_files = {}

        # Check if the structure file is already processed
        if pdb_file in processed_files:
            structure_file = processed_files[pdb_file]
        else:
            # Retrieve structure file (download or locate locally)
            structure_file = get_structure_file(
                pdb_id=Path(pdb_file).stem,  # Extract base filename without extensions
                temp_dir=temp_dir,
                file_format=file_format,
                pdb_dir=pdb_dir
            )
            if not structure_file:
                raise FileNotFoundError(f"No structure file found for {pdb_file} with valid extensions.")

            # Handle decompression if necessary
            if structure_file.endswith(".gz"):
                decompressed_file = os.path.join(temp_dir, Path(structure_file).stem)  # Remove .gz
                if not os.path.exists(decompressed_file):
                    with gzip.open(structure_file, "rb") as gz_file:
                        with open(decompressed_file, "wb") as out_file:
                            out_file.write(gz_file.read())
                    logger.info(f"Decompressed file: {decompressed_file}")
                structure_file = decompressed_file

            # Mark the file as processed
            processed_files[pdb_file] = structure_file

        # Process the geometry for the current chain
        geometre(
            filepath=structure_file,
            chain=chain,
            units_ids=units,
            o_path=output_file,
            ins_ids=insertions,
            draw=False
        )
        # Verify output file
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            logger.info(f"Processed and saved: {output_file}")
        else:
            logger.warning(f"Output file not created: {output_file}")
    except Exception as e:
        logger.error(f"Error processing {row['pdb_file']}: {e}")


def merge_results(temp_files, output_path):
    """
    Merge temporary output files into a single output file.
    """
    try:
        dfs = []
        for temp_file in temp_files:
            try:
                temp_df = pd.read_csv(temp_file, sep=",")
                dfs.append(temp_df)
            except Exception as e:
                logger.warning(f"Failed to load {temp_file}: {e}")

        if not dfs:
            logger.warning("No data files to merge.")
            return None

        final_df = pd.concat(dfs, ignore_index=True)
        final_df.to_csv(output_path, index=False, sep=",")
        logger.info(f"Merged results saved to: {output_path}")
        return output_path
    except Exception as e:
        logger.error(f"Error merging results: {e}")
        return None


def modify_batch_columns(input_file, output_file):
    """
    Modify the batch mode output by merging 'mean' and 'std deviation' rows into one line per pdb_id and chain??.
    """
    try:
        # Read the input file
        data = pd.read_csv(input_file, sep=",")
        # Filter out rows with 'mean' or 'std deviation'
        mean_std_rows = data[data['unit_start'].isin(['mean', 'std deviation'])]
        # Remove these rows from the main data to focus on actual unit data
        data_clean = data[~data['unit_start'].isin(['mean', 'std deviation'])]

        final_rows = []

        # Process each pdb_id and chain group
        for (pdb_id, chain), group in data_clean.groupby(['pdb_id', 'chain']):
            mean_row = mean_std_rows[(mean_std_rows['pdb_id'] == pdb_id) & (mean_std_rows['chain'] == chain) & (
                        mean_std_rows['unit_start'] == 'mean')]
            std_row = mean_std_rows[(mean_std_rows['pdb_id'] == pdb_id) & (mean_std_rows['chain'] == chain) & (
                        mean_std_rows['unit_start'] == 'std deviation')]

            if not mean_row.empty and not std_row.empty:
                # Combine the 'mean' and 'std deviation' into one row
                result_row = {
                    'pdb_id': pdb_id,
                    'chain': chain,
                    'curv_mean': round(mean_row['curvature'].values[0], 6),
                    'curv_std': round(std_row['curvature'].values[0], 6),
                    'twist_mean': round(mean_row['twist'].values[0], 6),
                    'twist_std': round(std_row['twist'].values[0], 6),
                    'twist_sign_mean': round(mean_row['twist_hand'].values[0], 6),
                    'twist_sign_std': round(std_row['twist_hand'].values[0], 6),
                    'pitch_mean': round(mean_row['pitch'].values[0], 6),
                    'pitch_std': round(std_row['pitch'].values[0], 6),
                    'pitch_sign_mean': round(mean_row['pitch_hand'].values[0], 6),
                    'pitch_sign_std': round(std_row['pitch_hand'].values[0], 6),
                    'tmscores_mean': round(mean_row['TM-score'].values[0], 6),
                    'tmscores_std': round(std_row['TM-score'].values[0], 6),
                    'yaw_mean': round(mean_row['yaw'].values[0], 6),
                    'yaw_std': round(std_row['yaw'].values[0], 6),
                }
                final_rows.append(result_row)

        final_df = pd.DataFrame(final_rows)
        final_df.to_csv(output_file, index=False, sep=",")
        logger.info(f"Batch mode aggregated results saved to: {output_file}")

    except Exception as e:
        print(f"Error modifying batch output columns: {e}")


def batch_repeats_geometry(tsv_path, output_path, num_threads=4, file_format="cif", pdb_dir=None):
    """
    Batch processes single_processing in parallel and modifies the output columns.
    """
    try:
        data = pd.read_csv(tsv_path, sep="\t", dtype={"pdb_file": str, "chain": str, "units": str, "insertion": str})
        temp_files = []
        processed_files = {}

        with tempfile.TemporaryDirectory() as temp_dir:
            with ThreadPoolExecutor(max_workers=num_threads) as executor:
                futures = []
                for idx, row in data.iterrows():
                    logger.debug(f"Submitting row {idx}: {row.to_dict()}")
                    temp_file = os.path.join(temp_dir, f"output_{idx}.csv")
                    temp_files.append(temp_file)
                    futures.append(
                        executor.submit(
                            process_entry, row, temp_dir, temp_file, file_format, pdb_dir, processed_files
                        )
                    )

                for future in futures:
                    try:
                        future.result()
                    except Exception as e:
                        logger.info(f"Row {idx} processed successfully.")

            logger.info(f"Temporary files created: {temp_files}")
            merged_file = merge_results(temp_files, output_path)

            if merged_file:
                logger.info(f"Batch processing completed. Results saved to: {merged_file}")
                modify_batch_columns(merged_file, output_path)
            else:
                logger.error("No merged file generated.")
    except Exception as e:
        logger.error(f"Error in batch processing: {e}")