#!/public/home/tianyi/tools/miniforge3/envs/biopython/bin/python

"""
Parallel FTP Genome Downloader

This script provides a robust, parallel downloading mechanism for genomic files 
from FTP sources. It supports:
- Parallel downloads with configurable worker count
- Checkpoint and resume functionality
- Optional file decompression
- Detailed logging

Key Features:
- Multithreaded download using ThreadPoolExecutor
- Progress tracking with tqdm
- Checkpoint mechanism to track download status
- Error handling and logging

Usage:
    python wget_gb_ftp_asm.py <ftp_list> <outdir> <logdir> <max_workers> [--decompress]

Dependencies:
    - Python 3.7+
    - tqdm
    - log_utils module
"""

import os
import sys
import re
import json
import argparse
import subprocess
import logging
import shutil
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Any

# Import logging utilities from log_utils
from log_utils import setup_logging, log_message

# Import tqdm for progress bar
from tqdm import tqdm


def load_previous_checkpoint(checkpoint_file: str) -> Dict[str, Any]:
    """
    Load previous checkpoint if exists, otherwise return empty checkpoint data.
    
    Handles potential file reading and JSON decoding errors gracefully.
    
    Args:
        checkpoint_file (str): Path to the checkpoint file
    
    Returns:
        Dict[str, Any]: Checkpoint data with 'total_jobs' and 'jobs' keys
    """
    try:
        if os.path.exists(checkpoint_file):
            with open(checkpoint_file, 'r') as f:
                return json.load(f)
    except (IOError, json.JSONDecodeError):
        pass
    
    return {'total_jobs': 0, 'jobs': {}}


def download_file(ftp_path: str, checkpoint_dir: str, decompress: bool = False) -> Dict[str, str]:
    """
    Download a file using wget and optionally decompress it using pigz.
    
    Supports downloading files from FTP sources with optional decompression.
    Provides detailed error tracking and status reporting.
    
    Args:
        ftp_path (str): Full FTP path to download
        checkpoint_dir (str): Directory to save downloaded files
        decompress (bool, optional): Whether to decompress .gz files. Defaults to False.
    
    Returns:
        Dict[str, str]: Download status with details including success/failure information
    """
    filename = os.path.basename(ftp_path)
    local_path = os.path.join(checkpoint_dir, filename)
    
    try:
        # Download file using wget
        wget_cmd = ['wget', '-O', local_path, ftp_path]
        result = subprocess.run(wget_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            return {
                'status': 'failed',
                'ftp_path': ftp_path,
                'error': result.stderr
            }
        
        # Decompress if requested and file is .gz
        if decompress and filename.endswith('.gz'):
            pigz_cmd = ['pigz', '-d', local_path]
            decomp_result = subprocess.run(pigz_cmd, capture_output=True, text=True)
            
            if decomp_result.returncode != 0:
                return {
                    'status': 'decompress_failed',
                    'ftp_path': ftp_path,
                    'error': decomp_result.stderr
                }
            
            # Update local_path to decompressed file
            local_path = local_path.rstrip('.gz')
        
        return {
            'status': 'success',
            'ftp_path': ftp_path,
            'local_path': local_path
        }
    
    except Exception as e:
        return {
            'status': 'failed',
            'ftp_path': ftp_path,
            'error': str(e)
        }


def main():
    """
    Main function to orchestrate parallel FTP genome downloads.
    
    Handles argument parsing, setup of download environment, 
    parallel download processing, and result tracking.
    
    Workflow:
    1. Parse command-line arguments
    2. Setup logging and checkpoint directories
    3. Generate download paths
    4. Execute parallel downloads
    5. Track and log download progress
    6. Generate download summary
    """
    parser = argparse.ArgumentParser(description='Download files via FTP with parallel processing')
    parser.add_argument('-i', dest='ftp_list', required=True, help='Path to file with FTP base paths')
    parser.add_argument('-o', dest='outdir', required=True, help='Output directory for final files')
    parser.add_argument('--logdir', dest='logdir', required=True, help='Directory for logs and checkpoints')
    parser.add_argument('-t', dest='threads', type=int, default=4, help='Number of parallel workers (default: 4)')
    parser.add_argument('-d', dest='decompress', action='store_true', help='Decompress .gz files')
    parser.add_argument('--full_path', dest='full_path', action='store_true', help='Provided paths are: /path/to/sequence.fna.gz')
    
    args = parser.parse_args()

    ftp_list = args.ftp_list
    outdir = args.outdir
    logdir = args.logdir
    max_workers = args.threads
    decompress = args.decompress
    full_path = args.full_path

    # Create directories if they don't exist
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(logdir, exist_ok=True)

    # Clean up previous checkpoint directories
    for folder in os.listdir(logdir):
        if folder.startswith('checkpoint_'):
            shutil.rmtree(os.path.join(logdir, folder))
            
    # Create timestamped checkpoint directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    checkpoint_dir = os.path.join(logdir, f'checkpoint_{timestamp}')
    os.makedirs(checkpoint_dir, exist_ok=True)
    
    # Setup logging
    setup_logging(logdir)
    log_message(f"Process started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    # Log the command-line arguments
    log_message(f"Command-line arguments: {' '.join(sys.argv)}\n")
    log_message(f"Parsed arguments:")
    log_message(f"  Input FTP list: {args.ftp_list}")
    log_message(f"  Output directory: {args.outdir}")
    log_message(f"  Log directory: {args.logdir}")
    log_message(f"  Parallel workers: {args.threads}")
    log_message(f"  Decompress files: {args.decompress}")
    log_message(f"  Use full path: {args.full_path}\n")
    
    # Read FTP base paths and generate full paths
    with open(ftp_list, 'r') as f:
        ftp_paths = [line.strip() for line in f if line.strip()]
    
    # Generate full paths based on the full_path flag
    if full_path:
        full_paths = ftp_paths
    else:
        full_paths = [f"{ftp_path}/{re.sub(r'.*/', '', ftp_path)}_genomic.fna.gz" for ftp_path in ftp_paths]
    
    # Load previous checkpoint
    checkpoint_file = os.path.join(logdir, 'latest_checkpoint.json')
    checkpoint_data = load_previous_checkpoint(checkpoint_file)
    
    # Filter out already downloaded paths
    pending_paths = [
        path for path in full_paths 
        if path not in checkpoint_data['jobs'] or 
           checkpoint_data['jobs'][path].get('status') != 'success'
    ]

    log_message(f"Total paths: {len(full_paths)}, Pending downloads: {len(pending_paths)}")
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        with tqdm(total=len(pending_paths), desc='Downloading Genomes', unit='file') as pbar:
            futures = {
                executor.submit(download_file, path, checkpoint_dir, decompress): path 
                for path in pending_paths
            }
            
            for future in as_completed(futures):
                ftp_path = futures[future]
                try:
                    result = future.result()
                    checkpoint_data['jobs'][ftp_path] = result
                    
                    # Update checkpoint file in real-time
                    with open(checkpoint_file, 'w') as f:
                        json.dump(checkpoint_data, f, indent=2)
                    
                    if result['status'] == 'success':
                        # Move file immediately after successful download
                        local_file = result['local_path']
                        dest_file = os.path.join(outdir, os.path.basename(local_file))
                        os.rename(local_file, dest_file)
                        
                        log_message(f"Successfully downloaded and moved: {ftp_path}")
                        pbar.update(1)
                    else:
                        log_message(f"Failed to download: {ftp_path}. Error: {result.get('error', 'Unknown error')}", level=logging.ERROR)
                        pbar.update(1)
                
                except Exception as e:
                    log_message(f"Unexpected error downloading {ftp_path}: {e}", level=logging.ERROR)
                    pbar.update(1)

    # Clean up checkpoint directories
    for folder in os.listdir(logdir):
        if folder.startswith('checkpoint_'):
            shutil.rmtree(os.path.join(logdir, folder))
            
    # Prepare summary of jobs
    successful_jobs = [
        job for job, details in checkpoint_data['jobs'].items() 
        if details['status'] == 'success'
    ]
    failed_jobs = [
        job for job, details in checkpoint_data['jobs'].items() 
        if details['status'] != 'success'
    ]
    
    # Log summary
    summary_message = [
        "Download Process Summary:",
        f"Total Jobs: {len(full_paths)}",
        f"Successful Downloads: {len(successful_jobs)}",
        f"Failed Downloads: {len(failed_jobs)}"
    ]
    log_message('\n'.join(summary_message))


if __name__ == '__main__':
    main()
    