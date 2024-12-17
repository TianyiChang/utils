#!/public/home/tianyi/tools/miniforge3/envs/biopython/bin/python

import sys
import argparse
import concurrent.futures
from pathlib import Path
import logging
from tqdm import tqdm
import io
import multiprocessing
import json
import os
import hashlib

# Import local logging utilities
from log_utils import setup_logging, log_message

def find_fasta_files(input_paths):
    """
    Find all FASTA files in the given input paths.
    
    Args:
        input_paths (list): List of paths to search for FASTA files
    
    Returns:
        list: List of paths to FASTA files
    """
    fasta_files = []
    for path in input_paths:
        path_obj = Path(path)
        if path_obj.is_dir():
            fasta_files.extend(list(path_obj.rglob('*.fa*')) + 
                               list(path_obj.rglob('*.fna')))
        elif path_obj.is_file() and path_obj.suffix in ['.fa*', '.fna']:
            fasta_files.append(path_obj)
    
    return fasta_files

def read_fasta_sequences(file_path):
    """
    Generator to read FASTA sequences one at a time.
    
    Args:
        file_path (Path): Path to input FASTA file
    
    Yields:
        tuple: (header, sequence) for each FASTA entry
    """
    with io.BufferedReader(open(file_path, 'rb')) as infile:
        header = None
        sequence = []
        
        for line in map(lambda x: x.decode('utf-8').strip(), infile):
            if line.startswith('>'):
                if header and sequence:
                    yield (header, ''.join(sequence))
                    sequence = []
                header = line
            else:
                sequence.append(line)
        
        # Last sequence
        if header and sequence:
            yield (header, ''.join(sequence))

def concat_fasta(input_file, output_file, checkpoint_file):
    """
    Concatenate sequences from a FASTA file with minimal memory usage and checkpointing.
    
    Args:
        input_file (Path): Path to input FASTA file
        output_file (Path): Path to output concatenated file
        checkpoint_file (Path): Path to checkpoint file
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Calculate file hash for unique identification
        file_hash = hashlib.md5(str(input_file).encode()).hexdigest()
        
        # Check if this file was already processed
        if os.path.exists(checkpoint_file):
            with open(checkpoint_file, 'r') as cp_file:
                checkpoints = json.load(cp_file)
                if file_hash in checkpoints.get('processed_files', {}):
                    log_message(f"Skipping already processed file: {input_file}")
                    return True
        
        # Process the file
        with open(output_file, 'a') as outfile:
            for header, sequence in read_fasta_sequences(input_file):
                outfile.write(f"{header}\n{sequence}\n")
        
        # Update checkpoint
        if not os.path.exists(checkpoint_file):
            checkpoints = {'processed_files': {}}
        else:
            with open(checkpoint_file, 'r') as cp_file:
                checkpoints = json.load(cp_file)
        
        checkpoints['processed_files'][file_hash] = {
            'path': str(input_file),
            'timestamp': str(Path(input_file).stat().st_mtime)
        }
        
        with open(checkpoint_file, 'w') as cp_file:
            json.dump(checkpoints, cp_file, indent=2)
        
        return True
    except Exception as e:
        log_message(f"Error processing {input_file}: {e}", logging.ERROR)
        return False

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description='Concatenate FASTA files')
    parser.add_argument('input_paths', nargs='+', help='Paths to input FASTA files or directories')
    parser.add_argument('output_file', help='Path to output concatenated FASTA file')
    parser.add_argument('--log_dir', default='logs', help='Directory for log files')
    parser.add_argument('--max_workers', type=int, 
                        default=max(1, multiprocessing.cpu_count() - 1), 
                        help='Maximum number of concurrent workers (default: number of CPUs - 1)')
    parser.add_argument('--checkpoint_dir', default='checkpoints', 
                        help='Directory to store checkpoint files')
    
    args = parser.parse_args()
    
    # Setup logging
    log_file = setup_logging(args.log_dir)
    log_message(f"Logging to {log_file}")
    
    # Ensure checkpoint directory exists
    Path(args.checkpoint_dir).mkdir(parents=True, exist_ok=True)
    
    # Create checkpoint file path
    checkpoint_file = Path(args.checkpoint_dir) / f"concat_fasta_checkpoint_{Path(args.output_file).name}.json"
    
    # Find FASTA files
    fasta_files = find_fasta_files(args.input_paths)
    
    # Log summary of files to process
    log_message(f"Found {len(fasta_files)} FASTA files to process:")
    for path in args.input_paths:
        dir_files = [f for f in fasta_files if Path(f).is_relative_to(Path(path))]
        log_message(f"  {path}: {len(dir_files)} files")
    
    # Ensure output file is empty before concatenation if no checkpoint exists
    if not checkpoint_file.exists():
        Path(args.output_file).unlink(missing_ok=True)
    
    # Concurrent processing with progress bar
    successful_files = 0
    failed_files = 0
    failed_file_list = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {executor.submit(concat_fasta, file, args.output_file, checkpoint_file): file for file in fasta_files}
        
        with tqdm(total=len(futures), desc="Processing FASTA Files", unit="file") as pbar:
            for future in concurrent.futures.as_completed(futures):
                file = futures[future]
                try:
                    result = future.result()
                    if result:
                        successful_files += 1
                    else:
                        failed_files += 1
                        failed_file_list.append(str(file))
                    pbar.update(1)
                except Exception as e:
                    log_message(f"Unexpected error processing {file}: {e}", logging.ERROR)
                    failed_files += 1
                    failed_file_list.append(str(file))
                    pbar.update(1)
    
    # Final logging and summary
    log_message(f"Concatenation complete.")
    log_message(f"Total files processed: {len(fasta_files)}")
    log_message(f"Successful files: {successful_files}")
    log_message(f"Failed files: {failed_files}")
    
    if failed_files > 0:
        log_message("Failed files:", logging.WARNING)
        for failed_file in failed_file_list:
            log_message(f"  {failed_file}", logging.WARNING)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
