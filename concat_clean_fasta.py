#!/public/home/tianyi/tools/miniforge3/envs/biopython/bin/python

import os
import re
import uuid
import csv
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from tqdm import tqdm

# Import local logging utilities
from log_utils import setup_logging, log_message
import logging

def generate_file_unicode_mapping(input_files):
    """
    Generate unique identifiers for input files.
    
    Args:
        input_files (list): List of input file paths
    
    Returns:
        list: List of dictionaries with 'unicode' and 'path' keys
    """
    file_unicode_mappings = []
    for file_path in input_files:
        file_unicode_mappings.append({
            'unicode': str(uuid.uuid4()),
            'path': str(file_path)
        })
    return file_unicode_mappings

def add_unicode_to_fasta_headers(content, unicode_mapping):
    """
    Add unicode to fasta headers at the beginning.
    
    Args:
        content (str): Fasta file content
        unicode_mapping (str): Unicode to add to headers
    
    Returns:
        str: Modified fasta content with unicodes at the beginning of headers
    """
    lines = content.split('\n')
    modified_lines = []
    for line in lines:
        if line.startswith('>'):
            # Remove the '>' and add unicode at the beginning
            header = line[1:]
            line = f">{unicode_mapping}_unicode_{header}"
        modified_lines.append(line)
    return '\n'.join(modified_lines)
    
def find_fasta_files(input_paths):
    """
    Find all FASTA files in the given input paths using concurrent processing.
    
    Args:
        input_paths (list): List of paths to search for FASTA files
    
    Returns:
        list: List of paths to FASTA files
    """
    
    def find_files_in_path(path):
        path_str = str(path)
        # Handle wildcard paths
        if '*' in path_str:
            matching_files = list(Path(path_str.rsplit('*', 1)[0]).glob(path_str.split('*')[-1]))
            return [f for f in matching_files if f.suffix in ['.fa', '.fasta', '.faa', '.fas', '.fna']]
        
        path_obj = Path(path)
        if path_obj.is_dir():
            return list(path_obj.rglob('*.fa*')) + list(path_obj.rglob('*.fna'))
        elif path_obj.is_file() and path_obj.suffix in ['.fa', '.fasta', '.faa', '.fas', '.fna']:
            return [path_obj]
        return []

    with ThreadPoolExecutor() as executor:
        # Submit search tasks for each input path
        futures = {executor.submit(find_files_in_path, path): path for path in input_paths}
        
        # Collect results
        fasta_files = []
        for future in as_completed(futures):
            fasta_files.extend(future.result())
    
    return list(set(fasta_files))  # Remove duplicates

def read_fasta_file(file_path):
    """Read the contents of a fasta file."""
    try:
        with open(file_path, 'r') as f:
            return file_path, f.read()
    except Exception as e:
        log_message(f'Failed to read {file_path}: {e}', logging.ERROR)
        return file_path, ''

def ensure_newline_between_sequences(content):
    """Ensure there is a newline between sequences."""
    return re.sub(r'([autcgAUTCG]{3,})\s{0,}>', '\\1\n>', content)

def clean_fasta_chunk(content):
    """Clean a chunk of FASTA content."""
    content = ensure_newline_between_sequences(content)
    cleaned_content = re.sub(r'\n{2,}', '\n', content).strip()
    return cleaned_content

def chunk_fasta_content(content, chunk_size):
    """Yield chunks of the content, ensuring no split at '>'."""
    start = 0
    while start < len(content):
        end = start + chunk_size
        if end >= len(content):
            yield content[start:]
            break
        
        end = content.rfind('>', start, end)
        if end == -1:
            end = len(content)
        
        yield content[start:end]
        start = end

def clean_fasta(input_file, output_file):
    """Clean the combined FASTA file."""
    with open(input_file, 'r') as infile:
        content = infile.read()

    chunk_size = 100_000
    chunks = chunk_fasta_content(content, chunk_size)

    max_workers = os.cpu_count() or 1  # Default to 1 if None
    with ProcessPoolExecutor(max_workers=max_workers * 2) as executor:
        cleaned_chunks = list(executor.map(clean_fasta_chunk, chunks))

    with open(output_file, 'w') as outfile:
        outfile.write('\n'.join(cleaned_chunks) + '\n')

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Combine and clean FASTA files from specified directories or files.')
    parser.add_argument('input_paths', nargs='+', help='One or more directories containing FASTA files or specific FASTA files.')
    parser.add_argument('output_file', help='The path to the combined fasta file.')
    parser.add_argument('--log_dir', default='logs', help='Directory for log files')
    parser.add_argument('--clean', action='store_true', help='Clean the concatenated fasta file (optional)')
    parser.add_argument('--no-unicode', action='store_true', help='Disable adding unicodes to fasta headers')
    parser.add_argument('--max_workers', type=int, 
                        default=max(1, (os.cpu_count() or 1) - 1), 
                        help='Maximum number of concurrent workers (default: number of CPUs - 1)')

    args = parser.parse_args()

    # Setup logging using log_utils
    log_file = setup_logging(args.log_dir)
    log_message(f"Logging to {log_file}")

    # Find all FASTA files
    all_files = find_fasta_files(args.input_paths)
    
    # Resolve full path for output file
    output_file_path = Path(args.output_file).resolve()
    log_message(f"Output file will be saved to: {output_file_path}", logging.INFO)
    
    # Generate unicode mappings if not disabled
    file_unicode_mappings = []
    if not args.no_unicode:
        file_unicode_mappings = generate_file_unicode_mapping(all_files)
    
    # Log summary of files to process
    log_message(f"Found {len(all_files)} FASTA files to process:", logging.INFO)
    for path in args.input_paths:
        dir_files = [f for f in all_files if Path(f).is_relative_to(Path(path))]
        log_message(f"  {path}: {len(dir_files)} files", logging.INFO)

    # Track processing with progress bar
    combined_content = []
    successful_files = 0
    failed_files = 0
    failed_file_list = []

    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        # Submit all tasks
        futures = {executor.submit(read_fasta_file, file): file for file in all_files}
        
        # Use tqdm for progress tracking
        with tqdm(total=len(futures), desc="Processing FASTA Files", unit="file") as pbar:
            for future in as_completed(futures):
                file = futures[future]
                try:
                    _, content = future.result()
                    if content:  # Only add non-empty content
                        # Add unicode to headers if not disabled
                        if not args.no_unicode:
                            unicode_mapping = next((m['unicode'] for m in file_unicode_mappings if m['path'] == str(file)), '')
                            content = add_unicode_to_fasta_headers(content, unicode_mapping)
                        
                        combined_content.append(content.strip())
                        successful_files += 1
                    else:
                        failed_files += 1
                        failed_file_list.append(file)
                    pbar.update(1)
                except Exception as e:
                    log_message(f"Failed to read file {file}: {e}", logging.ERROR)
                    failed_files += 1
                    failed_file_list.append(file)
                    pbar.update(1)

    # Create a folder if output file path contains relative path
    output_file_path.parent.mkdir(parents=True, exist_ok=True)

    # Write the combined content to a file
    with open(output_file_path, 'w') as f:
        f.write('\n'.join(combined_content))
        log_message(f'Successfully wrote combined sequences to {output_file_path}', logging.INFO)

    # Export unicode-file link CSV if unicodes were generated
    if file_unicode_mappings and not args.no_unicode:
        csv_path = output_file_path.parent / 'unicode_file_link.csv'
        with open(csv_path, 'w', newline='') as csvfile:
            fieldnames = ['unicode', 'path']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for mapping in file_unicode_mappings:
                writer.writerow(mapping)
        log_message(f'Exported unicode-file link to {csv_path}', logging.INFO)

    # Final logging and summary
    log_message(f"Total files processed: {len(all_files)}", logging.INFO)
    log_message(f"Successful files: {successful_files}", logging.INFO)
    log_message(f"Failed files: {failed_files}", logging.INFO)

    if failed_files > 0:
        log_message("Failed files:", logging.WARNING)
        for failed_file in failed_file_list:
            log_message(f"  {failed_file}", logging.WARNING)

    # Clean the file if --clean is specified
    if args.clean:
        log_message('Cleaning the combined fasta file...', logging.INFO)
        clean_fasta(output_file_path, output_file_path)
        log_message('Finished cleaning the combined fasta file', logging.INFO)
    else:
        log_message('Skipping cleaning step as --clean was not specified', logging.INFO)

    log_message("Processing complete!", logging.INFO)
    log_message(f"Output file location: {output_file_path}", logging.INFO)

if __name__ == '__main__':
    main()
