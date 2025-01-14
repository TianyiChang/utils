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
import re
import json
import argparse
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import logging utilities from log_utils
from log_utils import setup_logging, log_message

# Import tqdm for progress bar
from tqdm import tqdm


def load_checkpoint(checkpoint_file):
    """
    Load checkpoint data from a JSON file.
    
    Args:
        checkpoint_file (str): Path to the checkpoint JSON file
    
    Returns:
        dict: Checkpoint data with 'total_jobs' and 'jobs' keys
    """
    try:
        if os.path.exists(checkpoint_file):
            with open(checkpoint_file, 'r') as f:
                return json.load(f)
        else:
            # Return default checkpoint structure if file doesn't exist
            return {'total_jobs': 0, 'jobs': {}}
    except (json.JSONDecodeError, IOError):
        # Handle potential file reading or parsing errors
        log_message(f"Error reading checkpoint file: {checkpoint_file}", level=logging.ERROR)
        return {'total_jobs': 0, 'jobs': {}}


def save_checkpoint(checkpoint_file, checkpoint_data):
    """
    Save checkpoint data to a JSON file.
    
    Args:
        checkpoint_file (str): Path to save the checkpoint JSON file
        checkpoint_data (dict): Checkpoint data to save
    """
    try:
        # Ensure the directory exists
        os.makedirs(os.path.dirname(checkpoint_file), exist_ok=True)
        
        # Write checkpoint data
        with open(checkpoint_file, 'w') as f:
            json.dump(checkpoint_data, f, indent=2)
        
        log_message(f"Checkpoint saved to {checkpoint_file}")
    except Exception as e:
        log_message(f"Error saving checkpoint: {e}", logging.ERROR)


def determine_pending_jobs(full_paths, checkpoint_data, decompress):
    """
    Determine which jobs need to be processed based on checkpoint and decompress flag
    
    Args:
        full_paths (list): List of full FTP paths
        checkpoint_data (dict): Existing checkpoint data
        decompress (bool): Whether decompression is requested
    
    Returns:
        list: Paths that need to be processed
    """
    pending_paths = []
    
    for path in full_paths:
        # If path not in checkpoint, needs processing
        if path not in checkpoint_data['jobs']:
            pending_paths.append(path)
            continue
        
        job_status = checkpoint_data['jobs'][path].get('status')
        
        if decompress:
            # For decompress mode, check decompress status
            if job_status != 'decompress_success':
                pending_paths.append(path)
        else:
            # For non-decompress mode, only check download status
            if job_status != 'download_success':
                pending_paths.append(path)
    
    return pending_paths


def download_and_process_file(ftp_path, local_path, checkpoint_data, decompress=False, timeout=120):
    """
    Download and optionally decompress a file with timeout
    
    Args:
        ftp_path (str): Full FTP path to download
        local_path (str): Local directory to save the file
        checkpoint_data (dict): Existing checkpoint data
        decompress (bool): Whether to decompress .gz files
        timeout (int): Timeout in seconds for wget download (default: 120 seconds)
    
    Returns:
        dict: Download and processing status
    """
    # Check existing job status in checkpoint
    existing_job = checkpoint_data['jobs'].get(ftp_path, {})
    
    try:
        # Determine if download is needed
        if not existing_job or existing_job.get('status') in ['download_failed', 'decompress_failed', 'unexpected_error']:
            
            # Ensure output directory exists
            os.makedirs(os.path.dirname(local_path), exist_ok=True)

            try:
                # Download file using wget with timeout
                wget_cmd = ['timeout', str(timeout), 'wget', '-O', local_path, ftp_path]
                subprocess.run(
                    wget_cmd, 
                    capture_output=True, 
                    text=True, 
                    check=True
                )
            except subprocess.CalledProcessError as wget_error:
                log_message(f"wget failed for {ftp_path}:\n{wget_error.stderr}", level=logging.ERROR)
                return {
                    'status': 'download_failed',
                    'ftp_path': ftp_path,
                    'error': f"wget failed with return code {wget_error.returncode}: {wget_error.stderr}"
                }
            except subprocess.TimeoutExpired:
                log_message(f"Download timed out after {timeout} seconds for {ftp_path}", level=logging.ERROR)
                return {
                    'status': 'download_failed',
                    'ftp_path': ftp_path,
                    'error': f'Download timed out after {timeout} seconds'
                }
            
            # Verify download
            if not os.path.exists(local_path) or os.path.getsize(local_path) == 0:
                log_message(f"Download failed or empty file for {ftp_path}", level=logging.ERROR)
                return {
                    'status': 'download_failed',
                    'ftp_path': ftp_path,
                    'error': f'Download failed or empty file for {ftp_path}'
                }
            
            # Successful download
            job_result = {
                'status': 'download_success',
                'ftp_path': ftp_path,
                'local_path': local_path
            }
        else:
            # Use existing local path from previous checkpoint
            job_result = existing_job
            local_path = job_result.get('local_path', os.path.join(local_path, os.path.basename(ftp_path)))
        
        # Handle decompression if requested
        if decompress and local_path.endswith('.gz'):
            decompressed_path = local_path.rstrip('.gz')
            
            # Only decompress if not already successfully decompressed
            if existing_job.get('status') != 'decompress_success':
                try:
                    # Decompress, overwriting if partial
                    pigz_cmd = ['pigz', '-d', '-f', local_path]
                    subprocess.run(pigz_cmd, check=True, capture_output=True, text=True)
                    
                    # Update result for successful decompression
                    job_result.update({
                        'status': 'decompress_success',
                        'local_path': decompressed_path
                    })
                
                except subprocess.CalledProcessError as decomp_error:
                    log_message(f"Decompression failed for {ftp_path}:\n{decomp_error}", level=logging.ERROR)
                    job_result = {
                        'status': 'decompress_failed',
                        'ftp_path': ftp_path,
                        'error': f'Decompression failed for {ftp_path}'
                    }
            else:
                # Already successfully decompressed
                job_result['local_path'] = decompressed_path
        
        return job_result
    
    except subprocess.CalledProcessError as download_error:
        log_message(f"Download failed for {ftp_path}: {download_error}", level=logging.ERROR)
        return {
            'status': 'download_failed',
            'ftp_path': ftp_path,
            'error': str(download_error)
        }
    except Exception as e:
        log_message(f"Unexpected error for {ftp_path}: {e}", level=logging.ERROR)
        return {
            'status': 'unexpected_error',
            'ftp_path': ftp_path,
            'error': str(e)
        }


def main():
    parser = argparse.ArgumentParser(description='Download genomic sequences from FTP paths')
    parser.add_argument('ftp_list', help='File containing FTP paths')
    parser.add_argument('outdir', help='Output directory for downloaded files')
    parser.add_argument('--logdir', dest='logdir', required=True, help='Directory for logs and checkpoints')
    parser.add_argument('-t', dest='threads', type=int, default=4, help='Number of parallel workers (default: 4)')
    parser.add_argument('-d', dest='decompress', action='store_true', help='Decompress .gz files')
    
    args = parser.parse_args()
    
    # Setup directories and logging
    outdir = args.outdir
    logdir = args.logdir
    max_workers = args.threads
    decompress = args.decompress

    # Create directories if they don't exist
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(logdir, exist_ok=True)

    # Setup logging
    log_file = setup_logging(logdir)
    log_message(f"Log file created at: {log_file}")
    
    # Log input parameters
    log_message(f"  Input FTP list: {args.ftp_list}")
    log_message(f"  Output directory: {args.outdir}")
    log_message(f"  Log directory: {args.logdir}")
    log_message(f"  Parallel workers: {args.threads}")
    log_message(f"  Decompress files: {args.decompress}\n")
        
    # Checkpoint file management
    checkpoint_file = os.path.join(logdir, 'latest_checkpoint.json')
    checkpoint_data = load_checkpoint(checkpoint_file)
    
    # Read FTP base paths and generate full paths
    with open(args.ftp_list, 'r') as f:
        ftp_paths = [line.strip() for line in f if line.strip()]
    
    # Generate full paths with the specified pattern
    full_paths = [f"{ftp_path}/{re.sub(r'.*/', '', ftp_path)}_genomic.fna.gz" for ftp_path in ftp_paths]
    
    # Determine pending paths based on checkpoint and decompress flag
    pending_paths = determine_pending_jobs(full_paths, checkpoint_data, decompress)
    
    # Track successful downloads and decompressions
    successful_downloads = 0
    successful_decompressions = 0
    
    # Create progress bar
    progress_bar = tqdm(total=len(pending_paths), 
                        desc="Processing Files", 
                        unit="file")
    
    # Concurrent download and processing
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Prepare futures
        futures = {
            executor.submit(download_and_process_file, 
                            path, 
                            os.path.join(outdir, os.path.basename(path)), 
                            checkpoint_data, 
                            decompress): path 
            for path in pending_paths
        }
        
        # Process results
        for future in as_completed(futures):
            result = future.result()
            
            # Update checkpoint
            checkpoint_data['jobs'][result['ftp_path']] = result
            
            # Save checkpoint in real-time
            try:
                with open(checkpoint_file, 'w') as f:
                    json.dump(checkpoint_data, f, indent=2)
            except Exception as e:
                log_message(f"Error saving real-time checkpoint: {e}", level=logging.ERROR)
            
            # Count successes
            if result['status'] == 'download_success':
                successful_downloads += 1
            elif result['status'] == 'decompress_success':
                successful_downloads += 1
                successful_decompressions += 1
            
            # Update progress bar
            progress_bar.update(1)
        
        # Close progress bar
        progress_bar.close()
    
    # Log summary
    log_message(f"Total paths: {len(full_paths)}")
    log_message(f"Pending paths: {len(pending_paths)}")
    log_message(f"Successful downloads: {successful_downloads}")
    log_message(f"Successful decompressions: {successful_decompressions}")
    
    # Final checkpoint save (optional, but ensures final state is saved)
    save_checkpoint(checkpoint_file, checkpoint_data)

if __name__ == "__main__":
    main()
