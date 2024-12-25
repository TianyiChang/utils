#!/public/home/tianyi/tools/miniforge3/envs/biopython/bin/python

import os
import sys
import csv
import time
import argparse
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from log_utils import setup_logging, log_message

def read_accessions(input_source):
    """
    Read nuccore accessions from a space-separated list or CSV file.
    
    Args:
        input_source (str): Path to CSV or space-separated list of accessions
    
    Returns:
        list: List of nuccore accessions
    """
    if os.path.isfile(input_source):
        try:
            with open(input_source, 'r') as f:
                # Check if it's a CSV file
                if input_source.lower().endswith('.csv'):
                    reader = csv.DictReader(f)
                    return [row['nt_acc'].strip() for row in reader]
                else:
                    # Assume space-separated text file
                    return [acc.strip() for acc in f.read().split()]
        except Exception as e:
            log_message(f"Error reading input file: {e}", logging.ERROR)
            sys.exit(1)
    else:
        # Assume space-separated input
        return input_source.split()

def get_assembly_accession(nuccore_acc, max_retries=3):
    """
    Convert a nuccore accession to its corresponding assembly accession.
    
    Args:
        nuccore_acc (str): Nuccore accession number
        max_retries (int): Maximum number of retry attempts
    
    Returns:
        tuple: (nuccore_acc, assembly_acc) or (nuccore_acc, None)
    """
    for attempt in range(max_retries):
        try:
            # Use a more robust command with error handling
            cmd = f'esearch -db nuccore -query "{nuccore_acc}" | elink -target assembly | esummary | xtract -pattern DocumentSummary -element AssemblyAccession'
            
            # Run the command with timeout and error handling
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
            
            if result.returncode == 0:
                assembly_acc = result.stdout.strip()
                if assembly_acc:
                    return (nuccore_acc, assembly_acc)
                else:
                    log_message(f"No assembly accession found for {nuccore_acc}", logging.WARNING)
            else:
                log_message(f"Command failed for {nuccore_acc}: {result.stderr}", logging.ERROR)
            
            # Wait before retrying
            time.sleep(2 ** attempt)
        
        except subprocess.TimeoutExpired:
            log_message(f"Timeout occurred for {nuccore_acc}, attempt {attempt + 1}", logging.WARNING)
        except Exception as e:
            log_message(f"Error processing {nuccore_acc}: {e}", logging.ERROR)
        
        # Brief pause between retries to avoid overwhelming NCBI servers
        time.sleep(1)
    
    log_message(f"Failed to get assembly accession for {nuccore_acc} after {max_retries} attempts", logging.ERROR)
    return (nuccore_acc, None)

def main():
    parser = argparse.ArgumentParser(description='Convert Nuccore Accessions to Assembly Accessions')
    parser.add_argument('input', help='Space-separated list of accessions or path to CSV with "nt_acc" column')
    parser.add_argument('output', help='Path to output CSV file')
    parser.add_argument('log_dir', help='Directory for log files')
    
    args = parser.parse_args()
    
    # Setup logging
    log_file = setup_logging(args.log_dir)
    log_message(f"Log file created at: {log_file}")
    
    # Read input accessions
    nuccore_accs = read_accessions(args.input)
    log_message(f"Found {len(nuccore_accs)} nuccore accessions to process")
    
    # Process accessions with concurrent processing and progress bar
    results = []
    with ThreadPoolExecutor(max_workers=4) as executor:
        # Submit all jobs
        futures = {executor.submit(get_assembly_accession, acc): acc for acc in nuccore_accs}
        
        # Process results as they complete
        for future in tqdm(as_completed(futures), 
                           total=len(futures), 
                           desc="Converting Accessions", 
                           unit="accession"):
            nuccore_acc, asm_acc = future.result()
            if asm_acc:
                results.append({'nt_acc': nuccore_acc, 'asm_acc': asm_acc})
            else:
                log_message(f"Skipping {nuccore_acc} due to no assembly accession", logging.WARNING)
    
    # Write results to CSV
    try:
        with open(args.output, 'w', newline='') as csvfile:
            fieldnames = ['nt_acc', 'asm_acc']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        
        log_message(f"Converted {len(results)} accessions. Output saved to {args.output}")
    
    except Exception as e:
        log_message(f"Error writing output file: {e}", logging.ERROR)
        sys.exit(1)

if __name__ == '__main__':
    main()
