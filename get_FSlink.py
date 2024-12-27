#!/public/home/tianyi/tools/miniforge3/envs/biopython/bin/python
"""
FASTA Sequence ID Extractor (get_FSlink.py)

Overview:
    A high-performance utility for extracting and cataloging sequence identifiers 
    from multiple FASTA files in a directory. Designed for bioinformatics 
    preprocessing and large-scale genomic data management.

Key Features:
    - Concurrent processing of multiple FASTA files
    - Supports .fasta, .fna, and .fa file extensions
    - Memory-efficient ID extraction
    - Parallel processing with progress tracking
    - Comprehensive logging

Output:
    Generates a CSV file with two columns:
    1. file_id: Base name of the source FASTA file
    2. seq_id: Sequence identifier from FASTA header

Usage:
    python get_FSlink.py <input_directory> <output_csv> <log_directory>

Example:
    python get_FSlink.py /data/genomes output_ids.csv /logs

Performance:
    - Uses ThreadPoolExecutor for concurrent file processing
    - Limits concurrent workers to min(CPU cores, 8)
    - Provides real-time progress tracking with tqdm

Dependencies:
    - Python 3.7+
    - concurrent.futures
    - tqdm
    - log_utils module

Typical Use Cases:
    - Creating sequence ID catalogs
    - Preprocessing genomic data
    - Extracting identifiers for downstream analysis

Error Handling:
    - Logs individual file processing errors
    - Continues processing other files if one fails
    - Provides summary of total IDs extracted

Author: Tianyi Chang
Version: 1.0.0
Last Updated: 2024-12-27
"""

import argparse
from pathlib import Path
import csv
import os
import logging
from typing import Generator, Tuple
import concurrent.futures
import log_utils
from tqdm import tqdm

def extract_ids_from_fasta(file_path: Path, logger: logging.Logger) -> Generator[Tuple[str, str], None, None]:
    """
    Memory-efficient generator to extract file and sequence IDs from a FASTA file.
    
    Args:
        file_path (Path): Path to the FASTA file
        logger (logging.Logger): Logger for error reporting
    
    Yields:
        Tuple[str, str]: (file_id, seq_id) pairs
    """
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Remove '>' and strip whitespace
                    header = line.strip()[1:].split()[0]
                    yield (file_path.stem, header)
    except Exception as e:
        logger.error(f"Error processing {file_path}: {str(e)}")

def process_fasta_files(indir: Path, logger: logging.Logger) -> Generator[Tuple[str, str], None, None]:
    """
    Process all FASTA files in the input directory using concurrent processing.
    
    Args:
        indir (Path): Directory containing FASTA files
        logger (logging.Logger): Logger for error reporting
    
    Yields:
        Tuple[str, str]: (file_id, seq_id) pairs
    """
    # Find all FASTA files
    fasta_files = list(indir.glob('*.fasta')) + list(indir.glob('*.fna')) + list(indir.glob('*.fa'))
    logger.info(f"Found {len(fasta_files)} FASTA files to process")
    
    # Use ThreadPoolExecutor for concurrent processing
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(os.cpu_count(), 8)) as executor:
        # Submit all files for processing
        future_to_file = {
            executor.submit(extract_ids_from_fasta, file_path, logger): file_path 
            for file_path in fasta_files
        }
        
        # Use tqdm for progress tracking
        with tqdm(total=len(future_to_file), desc="Processing FASTA files", unit="file") as pbar:
            # Yield results as they complete
            for future in concurrent.futures.as_completed(future_to_file):
                file_path = future_to_file[future]
                try:
                    # Yield all IDs from this file
                    yield from future.result()
                    pbar.update(1)
                except Exception as e:
                    logger.error(f"Error processing {file_path}: {str(e)}")
                    pbar.update(1)

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(description="Extract file and sequence IDs from FASTA files")
    parser.add_argument("indir", type=str, help="Directory containing FASTA files")
    parser.add_argument("outfile", type=str, help="Path to output CSV file")
    parser.add_argument("logdir", type=str, help="Directory for log files")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Convert paths
    indir = Path(args.indir)
    outfile = Path(args.outfile)
    logdir = Path(args.logdir)
    
    # Setup logging
    log_utils.setup_logging(logdir)
    log_utils.log_message(f"Starting ID extraction from {indir}")
    
    try:
        # Ensure output directory exists
        outfile.parent.mkdir(parents=True, exist_ok=True)
        
        # Process files and write to CSV
        with open(outfile, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            # Write header
            csv_writer.writerow(['file_id', 'seq_id'])
            
            # Use tqdm to track overall progress
            total_ids = 0
            with tqdm(desc="Extracting IDs", unit="id") as pbar:
                # Process files and write IDs
                for file_id, seq_id in process_fasta_files(indir, logging.getLogger()):
                    csv_writer.writerow([file_id, seq_id])
                    total_ids += 1
                    pbar.update(1)
        
        log_utils.log_message(f"ID extraction completed. Output written to {outfile}")
        print(f"\nTotal IDs extracted: {total_ids}")
    
    except Exception as e:
        log_utils.log_message(f"Error during ID extraction: {str(e)}", logging.ERROR)
        raise

if __name__ == "__main__":
    main()
