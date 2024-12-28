#!/public/home/tianyi/tools/miniforge3/envs/biopython/bin/python

"""
NUCmer Alignment Pipeline

This script provides a flexible and parallel interface for running NUCmer 
sequence alignments with advanced filtering options.

Features:
- Support for single query/reference or batch processing via CSV
- Configurable delta-filter modes (single2single, many2many)
- Translocation filtering options
- Parallel processing with configurable workers
- Comprehensive logging and progress tracking

Usage:
    python run_nucmer.py -q QUERY -r REF [options]
    python run_nucmer.py -c COMPARISON_CSV [options]

Options:
    -q, --query         Query sequence file/directory
    -r, --ref           Reference sequence file/directory
    -c, --csv           CSV file with 'query' and 'ref' columns
    --single2single     Use single-to-single alignment mode
    --many2many         Use many-to-many alignment mode
    --translo           Enable translocation filtering
    --no_translo        Disable translocation filtering
    --min_pident        Minimum percent identity (default: 90)
    --min_aln_len       Minimum alignment length (default: 100)
    -p, --processes     Number of parallel workers (default: 8)
    -o, --outdir        Output directory
    --logdir            Log directory
"""

import sys
import csv
import argparse
import logging
import subprocess
from pathlib import Path
from typing import List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

# Import logging and progress utilities
from log_utils import setup_logging, log_message
from tqdm import tqdm


def parse_comparison_csv(csv_path: Path) -> List[Tuple[str, str]]:
    """
    Parse CSV file with query and reference paths.
    
    Args:
        csv_path (Path): Path to the CSV file
    
    Returns:
        List of (query, reference) path tuples
    """
    comparisons = []
    try:
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                comparisons.append((row['query'], row['ref']))
    except Exception as e:
        log_message(f"Error parsing CSV: {e}", level=logging.ERROR)
        raise
    return comparisons


def run_nucmer_alignment(
    query: str, 
    ref: str, 
    outdir: Path, 
    min_pident: int, 
    min_aln_len: int, 
    delta_mode: str, 
    translo_mode: Optional[str] = None
) -> Tuple[bool, str]:
    """
    Run NUCmer alignment with delta-filter and show-coords.
    
    Args:
        query (str): Path to query sequence
        ref (str): Path to reference sequence
        outdir (Path): Output directory for results
        min_pident (int): Minimum percent identity
        min_aln_len (int): Minimum alignment length
        delta_mode (str): Alignment mode (single2single/many2many)
        translo_mode (Optional[str]): Translocation filtering mode
    
    Returns:
        Tuple of (success_status, error_message)
    """
    try:
        # Ensure inputs are strings
        query = str(query)
        ref = str(ref)
        
        # Create output subdirectories
        delta_orig_dir = outdir / 'delta' / 'orig'
        delta_filtered_dir = outdir / 'delta' / 'filtered'
        coord_dir = outdir / 'coord'
        delta_orig_dir.mkdir(parents=True, exist_ok=True)
        delta_filtered_dir.mkdir(parents=True, exist_ok=True)
        coord_dir.mkdir(parents=True, exist_ok=True)

        # Generate unique prefix based on input files
        prefix = f"{Path(ref).stem}_vs_{Path(query).stem}"
        
        # Paths for output files
        out_delta_orig = delta_orig_dir / prefix
        out_delta_filtered = delta_filtered_dir / f"{prefix}_filtered.delta"
        out_coord = coord_dir / f"{prefix}_coords.txt"

        # NUCmer alignment
        nucmer_cmd = [
            'nucmer', 
            ref, 
            query, 
            '-p', str(out_delta_orig), 
            '-t', '1'
        ]
        
        # Delta-filter command
        delta_filter_cmd = [
            'delta-filter', 
            '-i', str(min_pident),
            '-l', str(min_aln_len)
        ]
        
        # Add mode-specific flags
        if delta_mode == 'many2many':
            delta_filter_cmd.append('-m')
        
        if translo_mode == 'translo':
            delta_filter_cmd.append('-1')
        elif translo_mode == 'no_translo':
            delta_filter_cmd.append('-g')
        
        delta_filter_cmd.extend([
            str(f"{out_delta_orig}.delta"),
            '>', 
            str(out_delta_filtered)
        ])
        
        # Show-coords command
        show_coords_cmd = [
            'show-coords', 
            '-r', 
            str(out_delta_filtered), 
            '>', 
            str(out_coord)
        ]

        # Execute commands
        subprocess.run(nucmer_cmd, check=True)
        subprocess.run(' '.join(delta_filter_cmd), shell=True, check=True)
        subprocess.run(' '.join(show_coords_cmd), shell=True, check=True)

        return True, ""
    
    except subprocess.CalledProcessError as e:
        return False, f"Subprocess error: {e}"
    except Exception as e:
        return False, f"Unexpected error: {e}"
        

def main():
    parser = argparse.ArgumentParser(description='Run NUCmer Alignments')
    
    # Input arguments
    parser.add_argument('-q', '--query', help='Query sequence file')
    parser.add_argument('-r', '--ref', help='Reference sequence file')
    parser.add_argument('-c', '--csv', type=Path, help='CSV file with query and ref paths')
    
    # Alignment mode selection
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--single2single', action='store_true', help='Single-to-single alignment')
    mode_group.add_argument('--many2many', action='store_true', help='Many-to-many alignment')
    
    # Translocation mode group (required only for single2single)
    translo_group = parser.add_mutually_exclusive_group()
    translo_group.add_argument('--translo', 
                                action='store_true',
                                help='Enable translocation filtering (required for single2single)')
    translo_group.add_argument('--no_translo', 
                                action='store_true',
                                help='Disable translocation filtering (required for single2single)')
    
    # Other parameters
    parser.add_argument('--min_pident', type=float, default=90.0, help='Minimum percent identity (default: 90)')
    parser.add_argument('--min_aln_len', type=float, default=100.0, help='Minimum alignment length (default: 100)')
    parser.add_argument('-p', '--processes', type=int, default=8, help='Number of parallel workers (default: 8)')
    parser.add_argument('-o', '--outdir', type=Path, required=True, help='Output directory')
    parser.add_argument('--logdir', type=Path, required=True, help='Log directory')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate inputs
    if args.csv:
        # If CSV is used, no need for separate query and ref
        comparisons = parse_comparison_csv(args.csv)
        if args.query or args.ref:
            parser.error("When using -c/--csv, -q/--query and -r/--ref are not required")
    else:
        # If not using CSV, both query and ref are required
        if not (args.query and args.ref):
            parser.error("When not using -c/--csv, both -q/--query and -r/--ref are required")
        comparisons = [(args.query, args.ref)]

    # Determine translocation mode
    if args.single2single:
        # For single2single, translo mode must be specified
        if not (args.translo or args.no_translo):
            parser.error("For single2single mode, you must specify either --translo or --no_translo")
        translo_mode = 'translo' if args.translo else 'no_translo'
    else:
        # For many2many, translocation mode is banned
        if args.translo or args.no_translo:
            parser.error("For many2many mode (translocation is enabled by default) and you cannot specify --translo or --no_translo")
        translo_mode = None

    # Setup logging
    setup_logging(args.logdir)
    log_message(f"NUCmer Alignment Process Started: {' '.join(sys.argv)}")
    
    # Determine alignment mode
    delta_mode = 'many2many' if args.many2many else 'single2single'
    
    # Create output directory
    args.outdir.mkdir(parents=True, exist_ok=True)
    
    # Parallel processing
    successful_jobs = 0
    failed_jobs = 0
    
    with ProcessPoolExecutor(max_workers=args.processes) as executor:
        futures = {
            executor.submit(
                run_nucmer_alignment, 
                query, 
                ref, 
                args.outdir, 
                args.min_pident, 
                args.min_aln_len, 
                delta_mode, 
                translo_mode
            ): (query, ref) for query, ref in comparisons
        }
        
        with tqdm(total=len(futures), desc='NUCmer Alignments', unit="alignments") as pbar:
            for future in as_completed(futures):
                query, ref = futures[future]
                try:
                    success, error_msg = future.result()
                    if success:
                        successful_jobs += 1
                        # log_message(f"Successfully aligned {query} vs {ref}")
                    else:
                        failed_jobs += 1
                        log_message(f"Failed to align {query} vs {ref}: {error_msg}", level=logging.ERROR)
                    pbar.update(1)
                except Exception as e:
                    failed_jobs += 1
                    log_message(f"Unexpected error aligning {query} vs {ref}: {e}", level=logging.ERROR)
                    pbar.update(1)
    
    # Final summary
    log_message("\nNUCmer Alignment Summary:")
    log_message(f"Total Jobs: {len(futures)}")
    log_message(f"Successful Alignments: {successful_jobs}")
    log_message(f"Failed Alignments: {failed_jobs}")

if __name__ == '__main__':
    main()
