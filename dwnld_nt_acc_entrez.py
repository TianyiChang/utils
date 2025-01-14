#!/public/home/tianyi/tools/miniforge3/envs/biopython/bin/python

"""
NCBI Nucleotide Sequence Downloader

This script provides a flexible utility to download nucleotide sequences from NCBI 
using the Entrez E-utilities. It supports multiple input methods for specifying 
accession numbers:

1. Direct input of space-separated accession numbers
2. TSV file with a single column of accession numbers
3. CSV file with a single column of accession numbers

Features:
- Validates NCBI nucleotide accession numbers
- Parallel downloading of sequences
- Robust error handling and logging
- Configurable output and log directories

Usage Examples:
    # Download sequences using direct accession input
    python dwnld_nt_acc_entrez.py "NC_012589.1 NC_012589.2"

    # Download sequences from a TSV/CSV file
    python dwnld_nt_acc_entrez.py accessions.tsv
    python dwnld_nt_acc_entrez.py accessions.csv

    # Specify custom output and log directories
    python dwnld_nt_acc_entrez.py accessions.tsv -o output_seqs -l logs

Dependencies:
- pandas
- tqdm
- NCBI E-utilities (efetch)

Note: Requires Entrez command-line tools to be installed and configured.
"""

import logging
import os
import re
import pandas as pd
import subprocess
import concurrent.futures
from tqdm import tqdm
import argparse
from typing import List, Set, Union

# Import log_utils
from log_utils import setup_logging, log_message


def is_valid_ncbi_accession(accession: str) -> bool:
    """
    Validate if a string is a valid NCBI nucleotide accession.
    
    Checks for common NCBI nucleotide accession formats:
    - Starts with letters (e.g., NC_, NZ_, XM_, XR_, etc.)
    - Followed by numbers and letters
    - Followed by version number after a dot
    
    Args:
        accession (str): Accession to validate
    
    Returns:
        bool: True if valid NCBI nucleotide accession, False otherwise
    """
    # Regex pattern for NCBI nucleotide accessions
    # Matches formats like: 
    # NC_012345.1, NZ_ABCD01000001.1, XM_123456.7, 
    # NZ_ABCD01000001Z01, etc.
    pattern = r'^[A-Z]{2}_?[A-Z0-9]+\.\d{1,2}$'
    
    # Check if the accession matches the pattern and is not empty
    return bool(re.match(pattern, str(accession).strip()))

def read_accessions(input_source: Union[str, List[str]]) -> Set[str]:
    """
    Read unique accessions from various input sources.
    
    Args:
        input_source (Union[str, List[str]]): 
            - Path to a TSV/CSV file 
            - Space-separated string of accessions
            - List of accessions
    
    Returns:
        Set of unique valid NCBI nucleotide accession numbers
    """
    # If input is a space-separated string, split it
    if isinstance(input_source, str):
        # Check if input is a file path
        if os.path.isfile(input_source):
            # Determine file type
            file_ext = os.path.splitext(input_source)[1].lower()
            
            try:
                # Read file based on extension, skip headers
                if file_ext in ['.tsv', '.csv']:
                    # Read first column, skip headers
                    df = pd.read_csv(input_source, sep=',' if file_ext == '.csv' else '\t', 
                                     header=None, names=['accession'])
                    
                    # Filter for valid NCBI nucleotide accessions
                    accessions = set(
                        acc for acc in df['accession'].astype(str) 
                        if is_valid_ncbi_accession(acc)
                    )
                    
                    return accessions
                else:
                    raise ValueError(f"Unsupported file type: {file_ext}")
            except Exception as e:
                log_message(f"Error reading file {input_source}: {e}", level=logging.ERROR)
                return set()
        else:
            # Treat as space-separated string of accessions
            return set(
                acc for acc in input_source.split() 
                if is_valid_ncbi_accession(acc)
            )
    
    # If input is already a list, filter for valid accessions
    return set(
        acc for acc in input_source 
        if is_valid_ncbi_accession(acc)
    )

def download_sequence(accession: str, output_folder: str) -> bool:
    """
    Download nucleotide sequence for a given accession.
    
    Args:
        accession (str): Nucleotide accession number
        output_folder (str): Folder to save downloaded sequences
    
    Returns:
        bool: True if download successful, False otherwise
    """
    # Sanitize accession to remove any whitespace
    accession = accession.strip()
    
    output_file = os.path.join(output_folder, f"{accession}.fasta")
    
    try:
        # Use subprocess to run efetch command
        cmd = [
            'efetch', 
            '-db', 'nucleotide', 
            '-id', accession, 
            '-format', 'fasta'
        ]
        
        with open(output_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
        
        # Check if the command was successful
        if result.returncode == 0 and os.path.getsize(output_file) > 0:
            return True
        else:
            os.remove(output_file)  # Remove empty or failed download
            return False
    
    except Exception as e:
        # Remove the file if an exception occurs
        if os.path.exists(output_file):
            os.remove(output_file)
        log_message(f"Error downloading {accession}: {e}", level=logging.ERROR)
        return False

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description='Download nucleotide sequences from NCBI.')
    parser.add_argument('input', help='Input source of accessions (file path or space-separated accessions)')
    parser.add_argument('-o', '--output', default='whole_seqs', 
                        help='Output folder for downloaded sequences (default: whole_seqs)')
    parser.add_argument('-l', '--log', default='../logs', 
                        help='Log folder (default: ../logs)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Create necessary directories
    os.makedirs(args.output, exist_ok=True)
    
    # Setup logging
    log_file = setup_logging(args.log)
    
    # Get unique accessions
    accessions = read_accessions(args.input)
    log_message(f"Total unique accessions found: {len(accessions)}")
    
    # Prepare for parallel download
    failed_accessions = []
    successful_downloads = 0
    
    # Use concurrent futures for parallel download
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        # Create a list of futures
        futures = {
            executor.submit(download_sequence, acc, args.output): acc 
            for acc in accessions
        }
        
        # Use tqdm for progress tracking
        for future in tqdm(concurrent.futures.as_completed(futures), 
                           total=len(futures), 
                           desc="Downloading Sequences",
                           unit = "accession"):
            accession = futures[future]
            try:
                success = future.result()
                if success:
                    successful_downloads += 1
                else:
                    failed_accessions.append(accession)
            except Exception as e:
                failed_accessions.append(accession)
    
    # Log results
    log_message(f"Successful downloads: {successful_downloads}")
    log_message(f"Failed downloads: {len(failed_accessions)}")
    log_message("Failed Accessions: " + ", ".join(failed_accessions))
    
    print(f"Download complete. Check {log_file} for details.")

if __name__ == "__main__":
    main()
