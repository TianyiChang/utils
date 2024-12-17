#!/public/home/tianyi/tools/miniforge3/envs/biopython/bin/python

import os
import sys
import argparse
import subprocess
import csv
import logging

# Import local logging utilities
from log_utils import setup_logging, log_message

def run_hmmsearch(
    path2hmm, 
    path2query, 
    path2out, 
    log_dir='../logs', 
    cpu=16, 
    threshold=24.5,
    no_convert=False
):
    """
    Run HMMER hmmsearch with specified parameters
    
    Args:
        path2hmm (str): Path to HMM file
        path2query (str): Path to query sequence file
        path2out (str): Path to output tblout file
        log_dir (str): Directory to store log files
        cpu (int): Number of CPUs to use
        threshold (float): Bit score threshold
        no_convert (bool): Skip automatic conversion of tblout
    
    Returns:
        bool: True if successful, False otherwise
    """
    # Set up logging
    setup_logging(log_dir)
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(path2out), exist_ok=True)
    
    # HMMER command
    cmd = [
        "hmmsearch",
        "--cpu", str(cpu),
        "-T", str(threshold),
        "--tblout", path2out,
        path2hmm,
        path2query
    ]
    
    try:
        # Log command details
        log_message(f"Running HMMER hmmsearch with parameters:")
        log_message(f"  HMM file: {path2hmm}")
        log_message(f"  Query file: {path2query}")
        log_message(f"  Output file: {path2out}")
        log_message(f"  CPUs: {cpu}")
        log_message(f"  Bit score threshold: {threshold}")
        
        # Run the command
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Log output
        if result.stdout:
            log_message(result.stdout)
        if result.stderr:
            log_message(result.stderr, level=logging.WARNING)
        
        log_message(f"HMMER hmmsearch completed. Results saved to {path2out}")
        
        # Automatically convert unless specified not to
        if not no_convert:
            convert_path = f"{os.path.splitext(path2out)[0]}.tsv"
            parse_hmmsearch_tblout(path2out, convert_path)
            log_message(f"Converted tblout to TSV: {convert_path}")
        
        return True
    
    except subprocess.CalledProcessError as e:
        log_message(f"HMMER hmmsearch failed: {e}", level=logging.ERROR)
        log_message(f"Command output: {e.stdout}", level=logging.ERROR)
        log_message(f"Command error: {e.stderr}", level=logging.ERROR)
        return False

def parse_hmmsearch_tblout(input_file, output_file=None, output_format='tsv'):
    """
    Convert HMMER hmmsearch tblout to TSV or CSV
    
    Args:
        input_file (str): Path to input hmmsearch tblout file
        output_file (str, optional): Path to output file. If None, print to stdout
        output_format (str): Output format, either 'tsv' or 'csv'
    """
    # Define the columns to extract from the tblout file
    columns = [
        'target_name',     # 1
        'target_accession', # 2
        'query_name',       # 3
        'query_accession',  # 4
        'full_evalue',      # 6
        'full_score',       # 5
        'full_bias',        # 7
        'domain_evalue',    # 8
        'domain_score',     # 9
        'domain_bias',      # 10
    ]

    # Prepare output
    if output_file:
        output_handle = open(output_file, 'w')
        if output_format == 'tsv':
            writer = csv.writer(output_handle, delimiter='\t')
        else:
            writer = csv.writer(output_handle)
    else:
        output_handle = sys.stdout
        writer = csv.writer(output_handle, delimiter='\t')

    # Write header
    writer.writerow(columns)

    # Parse tblout file
    with open(input_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            # Split the line
            parts = line.split()
            
            # Extract relevant columns
            row_data = [
                parts[0],   # target_name
                parts[1],   # target_accession
                parts[2],   # query_name
                parts[3],   # query_accession
                parts[5],   # full_evalue
                parts[4],   # full_score
                parts[6],   # full_bias
                parts[7],   # domain_evalue
                parts[8],   # domain_score
                parts[9],   # domain_bias
            ]
            
            # Write row
            writer.writerow(row_data)

    # Close output file if it was opened
    if output_file:
        output_handle.close()

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Run HMMER hmmsearch and convert tblout")
    
    # Input and output file arguments
    parser.add_argument('path2hmm', help='Path to HMM file')
    parser.add_argument('path2query', help='Path to query sequence file')
    parser.add_argument('path2out', help='Path to output tblout file')
    
    # Optional arguments
    parser.add_argument('--cpu', type=int, default=16, help='Number of CPUs to use')
    parser.add_argument('-T', '--threshold', type=float, default=24.5, help='Bit score threshold')
    parser.add_argument('--log-dir', default='../logs', help='Directory to store log files')
    parser.add_argument('--no-convert', action='store_true', help='Skip automatic conversion of tblout')
    parser.add_argument('-f', '--format', choices=['tsv', 'csv'], default='tsv', 
                        help='Output format for conversion (default: tsv)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Run HMMER search
    success = run_hmmsearch(
        args.path2hmm, 
        args.path2query, 
        args.path2out, 
        log_dir=args.log_dir,
        cpu=args.cpu, 
        threshold=args.threshold,
        no_convert=args.no_convert
    )
    
    # If conversion is not skipped and a specific format is requested
    if not args.no_convert and args.format != 'tsv':
        convert_path = f"{os.path.splitext(args.path2out)[0]}.{args.format}"
        parse_hmmsearch_tblout(args.path2out, convert_path, output_format=args.format)
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
