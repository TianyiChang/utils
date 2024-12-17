import os
import shutil
import glob
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Split Fasta into chunks.')

# Add command-line arguments
parser.add_argument('-i', dest='path2infile', help='Path to input Fasta')
parser.add_argument('-o', dest='output_path', help='Output directory')
parser.add_argument('--scripts', dest='path2scripts', nargs='*', help='Path to a list of scripts to move')
parser.add_argument('--filesize', dest='char_count_per_file', type=int, help='Total number of characters for each file')
parser.add_argument('--foldersize', dest='file_count_per_folder', type=int, help='Total number of files for each subfolder')

# Parse the command-line arguments
args = parser.parse_args()

# Extract values from the args namespace
path2infile = args.path2infile
output_path = args.output_path
path2scripts = args.path2scripts
char_count_per_file = args.char_count_per_file
file_count_per_folder = args.file_count_per_folder

current_output_file = None
current_size = 0
output_count = 1

# 1. Split the file into batches
with open(path2infile, 'r') as infile:
    current_output_file = open(
        f'{output_path}/split_batch_{output_count}.fasta', 'w')

    for line in infile:
        if line.startswith('>'):
            if current_size > char_count_per_file:
                current_output_file.close()
                output_count += 1
                current_size = 0
                current_output_file = open(
                    f'{output_path}/split_batch_{output_count}.fasta', 'w')
                current_output_file.write(line)
            else:
                current_output_file.write(line)
        else:
            current_output_file.write(line)
            current_size += len(line)

if current_output_file:
    current_output_file.close()

# 2. Create subfolders each contain N batch files
subfolder_count = round(output_count / file_count_per_folder)
batch_files = os.listdir(output_path)

if subfolder_count == 0:
    subfolder_count = 1

def move_files_to_subfolder(start_index, end_index, i):
    subfolder_qseq = f'chunk_{i + 1}/query_seqs'
    subfolder_scripts = f'chunk_{i + 1}/scripts'
    
    # Create subfolder if it doesn't exist
    subfolder_qseq_path = os.path.join(output_path, subfolder_qseq)
    subfolder_scripts_path = os.path.join(output_path, subfolder_scripts)

    os.makedirs(subfolder_qseq_path, exist_ok=True)
    os.makedirs(subfolder_scripts_path, exist_ok=True)

    # Move batch sequences to each subfolder
    files_to_move = batch_files[start_index:end_index]
    for file_name in files_to_move:
        file_path = os.path.join(output_path, file_name)
        new_file_path = os.path.join(subfolder_qseq_path, file_name)
        shutil.move(file_path, new_file_path)

    # Copy scripts to the subfolder
    for script in path2scripts:
        shutil.copy(script, subfolder_scripts_path)

# Using ThreadPoolExecutor to handle subfolder file moves concurrently
with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
    futures = []
    for i in range(subfolder_count):
        start_index = i * file_count_per_folder
        end_index = min((i + 1) * file_count_per_folder, output_count)
        futures.append(executor.submit(move_files_to_subfolder, start_index, end_index, i))

    # Wait for all threads to complete
    for future in futures:
        future.result()

# Move remaining batch files to chunk_1
if glob.glob(f'{output_path}/*.fasta'):
    subprocess.run(f'mv {output_path}/*.fasta {output_path}/chunk_1/query_seqs/',
                   shell=True, capture_output=True, text=True)
