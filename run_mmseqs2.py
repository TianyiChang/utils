#!/public/home/tianyi/tools/miniforge3/envs/biopython/bin/python

import os
import sys
import argparse
import subprocess
import logging

# Import local logging utilities
from log_utils import setup_logging, log_message

def run_mmseqs_linclust(
    input_file, 
    output_file, 
    log_dir='../logs', 
    threads=8, 
    min_seq_id=0.4,  # Default sequence identity threshold
    coverage=0.8,    # Default coverage threshold
    cov_mode=2,      # Default coverage mode
    cluster_mode=1   # Default clustering mode
):
    """
    Run MMseqs2 easy-linclust to cluster protein sequences
    
    Args:
        input_file (str): Path to input protein sequence file
        output_file (str): Path to output (name prefix)
        log_dir (str): Directory to store log files
        threads (int): Number of threads to use
        min_seq_id (float): Minimum sequence identity threshold (0.0-1.0)
        coverage (float): Coverage threshold (0.0-1.0)
        cov_mode (int): Coverage mode (0-3)
        cluster_mode (int): Clustering mode (0-3)
    """
    # Set up logging
    setup_logging(log_dir)
    
    # Prepare output and temporary directories
    output_base = os.path.splitext(output_file)[0]
    tmp_dir = f"{output_base}_tmp"
    
    # Ensure directories exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # MMseqs2 command with configurable parameters
    cmd = [
        "mmseqs", "easy-linclust", 
        input_file,          # Input file
        output_base,         # Output base name
        tmp_dir,             # Temporary directory
        "--threads", str(threads),
        "--min-seq-id", str(min_seq_id),  # Sequence identity threshold
        "-c", str(coverage),              # Coverage threshold
        "--cov-mode", str(cov_mode),      # Coverage mode
        "--cluster-mode", str(cluster_mode)  # Clustering mode
    ]
    
    try:
        # Run the command
        log_message(f"Running MMseqs2 linclust with parameters:")
        log_message(f"  Input: {input_file}")
        log_message(f"  Output: {output_file}")
        log_message(f"  Min Seq ID: {min_seq_id}")
        log_message(f"  Coverage: {coverage}")
        log_message(f"  Coverage Mode: {cov_mode}")
        log_message(f"  Cluster Mode: {cluster_mode}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Log output
        if result.stdout:
            log_message(result.stdout)
        if result.stderr:
            log_message(result.stderr, level=logging.WARNING)
        
        log_message(f"MMseqs2 linclust completed. Representative sequences saved to {output_file}")
    
    except subprocess.CalledProcessError as e:
        log_message(f"MMseqs2 linclust failed: {e}", level=logging.ERROR)
        log_message(f"Command output: {e.stdout}", level=logging.ERROR)
        log_message(f"Command error: {e.stderr}", level=logging.ERROR)
        raise

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Run MMseqs2 linclust for protein sequence clustering")
    
    # Input and output file arguments
    parser.add_argument("input", help="Input FASTA file path")
    parser.add_argument("output", help="Output representative sequences file path")
    
    # MMseqs2 parameter arguments
    parser.add_argument("--min-seq-id", type=float, default=0.4, 
                        help="Minimum sequence identity threshold (0.0-1.0)")
    parser.add_argument("-c", "--coverage", type=float, default=0.8, 
                        help="Coverage threshold (0.0-1.0)")
    parser.add_argument("--cov-mode", type=int, default=2, 
                        help="Coverage mode (0-3)")
    parser.add_argument("--cluster-mode", type=int, default=1, 
                        help="Clustering mode (0-3)")
    
    # Additional optional arguments
    parser.add_argument("--threads", type=int, default=8, 
                        help="Number of threads to use")
    parser.add_argument("--log-dir", default="../logs", 
                        help="Directory to store log files")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Run MMseqs2 linclust
    run_mmseqs_linclust(
        input_file=args.input,
        output_file=args.output,
        log_dir=args.log_dir,
        threads=args.threads,
        min_seq_id=args.min_seq_id,
        coverage=args.coverage,
        cov_mode=args.cov_mode,
        cluster_mode=args.cluster_mode
    )

if __name__ == "__main__":
    sys.exit(main())
