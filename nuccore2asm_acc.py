#!/public/home/tianyi/tools/miniforge3/envs/biopython/bin/python

import os
import sys
import csv
import time
import argparse
import subprocess
import logging
import re
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
    Convert a nuccore accession to its assembly information.
    
    Args:
        nuccore_acc (str): Nuccore accession number
        max_retries (int): Maximum number of retry attempts
    
    Returns:
        tuple: (asm_acc, asm_name, asm_species_taxid, asm_ftp_path) or (None, None, None, None)
    """
    for attempt in range(max_retries):
        try:
            # Get all available information
            cmd = f'esearch -db nuccore -query "{nuccore_acc}" | elink -target assembly | esummary | xtract -pattern DocumentSummary -element "*"'
            
            # Run the command with timeout and error handling
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
            
            if result.returncode == 0:
                # Extract all possible elements
                output = result.stdout.strip()
                
                # Attempt to find AssemblyAccession
                asm_acc_match = re.findall(r'<AssemblyAccession>(.*?)</AssemblyAccession>', output)
                asm_acc = asm_acc_match[0] if asm_acc_match else None
                
                # Attempt to find AssemblyName
                asm_name_match = re.findall(r'<AssemblyName>(.*?)</AssemblyName>', output)
                asm_name = asm_name_match[0] if asm_name_match else None
                
                # Attempt to find SpeciesTaxid
                species_taxid_match = re.findall(r'<SpeciesTaxid>(.*?)</SpeciesTaxid>', output)
                asm_species_taxid = species_taxid_match[0] if species_taxid_match else None
                
                # Attempt to find FtpPath_RefSeq
                ftp_path_match = re.findall(r'<FtpPath_RefSeq>(.*?)</FtpPath_RefSeq>', output)
                asm_ftp_path = ftp_path_match[0] if ftp_path_match else None
                
                # Log warnings for missing information
                if not asm_acc:
                    log_message(f"Warning: No valid AssemblyAccession found for {nuccore_acc}", logging.WARNING)
                
                if not asm_name:
                    log_message(f"Warning: No valid AssemblyName found for {nuccore_acc}", logging.WARNING)

                if not asm_species_taxid:
                    log_message(f"Warning: No valid SpeciesTaxid found for {nuccore_acc}", logging.WARNING)
                
                if not asm_ftp_path:
                    log_message(f"Warning: No valid FtpPath_RefSeq found for {nuccore_acc}", logging.WARNING)
                
                # Generate FTP path if not already found and both asm_acc and asm_name are present
                if not asm_ftp_path and asm_acc and asm_name:
                    try:
                        # Extract components from assembly accession
                        prefix = re.match(r'([^_]+)', asm_acc).group(1)
                        level1 = re.search(r'_(\d{3})', asm_acc).group(1)
                        level2 = re.search(r'_\d{3}(\d{3})', asm_acc).group(1)
                        level3 = re.search(r'_\d{6}(\d{3})', asm_acc).group(1)
                        
                        # Construct FTP path
                        asm_ftp_path = '/'.join([
                            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all',
                            prefix,
                            level1,
                            level2,
                            level3,
                            f"{asm_acc}_{asm_name}"
                        ])
                    
                    except (AttributeError, IndexError) as e:
                        log_message(f"Failed to generate FTP path for {nuccore_acc}: {e}", logging.WARNING)
                
                return (asm_acc, asm_name, asm_species_taxid, asm_ftp_path)
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
    
    log_message(f"Failed to get assembly information for {nuccore_acc} after {max_retries} attempts", logging.ERROR)
    return (None, None, None, None)


def main():
    parser = argparse.ArgumentParser(description='Convert Nuccore Accessions to Assembly Information')
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
        # Create futures for all jobs
        futures = {executor.submit(get_assembly_accession, acc): acc for acc in nuccore_accs}
        
        # Create tqdm progress bar
        progress_bar = tqdm(total=len(futures), 
                            desc="Converting Accessions", 
                            unit="accession")
        
        # Process results as they complete
        for future in as_completed(futures):
            try:
                asm_acc, asm_name, asm_species_taxid, asm_ftp_path = future.result()
                
                # Only add results with ftp_path
                if asm_ftp_path:
                    result_entry = {
                        'nt_acc': futures[future],
                        'asm_acc': asm_acc if asm_acc is not None else '',
                        'asm_name': asm_name if asm_name is not None else '',
                        'asm_ftp_path': asm_ftp_path,
                        'asm_species_taxid': asm_species_taxid if asm_species_taxid is not None else ''
                    }
                    
                    results.append(result_entry)
                
                # Update progress bar
                progress_bar.update(1)
            
            except Exception as e:
                log_message(f"Error processing future: {e}", logging.ERROR)
                progress_bar.update(1)
        
        # Close the progress bar
        progress_bar.close()
        
    # Write results to TSV
    try:
        with open(args.output, 'w', newline='') as tsvfile:
            # Dynamically create fieldnames
            fieldnames = ['nt_acc', 'asm_acc', 'asm_name', 'asm_ftp_path']
            if any('asm_species_taxid' in result for result in results):
                fieldnames.append('asm_species_taxid')
            
            # Use csv writer with tab delimiter
            writer = csv.writer(tsvfile, delimiter='\t')
            
            # Write header
            writer.writerow(fieldnames)
            
            # Write data rows
            for result in results:
                # Ensure consistent order of fields
                row = [result.get(field, '') for field in fieldnames]
                writer.writerow(row)
        
        log_message(f"Converted {len(results)} accessions. Output saved to {args.output}")
    
    except Exception as e:
        log_message(f"Error writing output file: {e}", logging.ERROR)
        sys.exit(1)

if __name__ == '__main__':
    main()
