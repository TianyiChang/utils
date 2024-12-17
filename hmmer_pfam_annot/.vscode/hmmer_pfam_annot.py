import os
import subprocess
import re
import shutil
import glob
import logging
import argparse
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor

from my_logger import logger
from my_logger import formatter

# Parse command-line options using argparse
parser = argparse.ArgumentParser(description='hmmer_pfam_annot')
parser.add_argument('-i', type=str, dest='INPUT', help='A query file or a path to all the query files (query should only be protein sequences)', required=True)
parser.add_argument('-o', type=str, dest='OUTPUT_DIR', help='Path to the output directory', required=True)
parser.add_argument('-r', type=str, dest='DATABASE_DIR', help='Path to the directory containing Pfam-A database', required=True)
parser.add_argument('-t', type=int, dest='THREADS', help='The number of threads', default=4)
parser.add_argument('--no_infile_name', action='store_true', help='Skip the creation of FSlink records')
# parser.add_argument('--annot', action='store_true', help='Annotate MGE proteins using the Pfam-A database (default is off)')

# Store values into variables
args = parser.parse_args()

#===================#
# Passing arguments #
#===================#

INPUT = args.INPUT
OUTPUT_DIR = args.OUTPUT_DIR
DATABASE_DIR = args.DATABASE_DIR
THREADS = args.THREADS
NO_INFILE_NAME = args.no_infile_name

#=============#
# Set up logs #
#=============#

os.makedirs(f'{OUTPUT_DIR}/log', exist_ok=True)

# get the date
current_time = datetime.now().strftime("%Y%m%d%H%M%S")

# create a filehandler
fh = logging.FileHandler(f'{OUTPUT_DIR}/log/hmmer_pfam_annot_{current_time}.log')

# set the imported formatter
fh.setFormatter(formatter)

# add fh to the imported logger
logger.addHandler(fh)

# Log the command-line arguments
logger.info(f'Command-line arguments: {vars(args)}')

#=======================#
# Install denpendencies #
#=======================#

# locate the scripts
SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))

# customize messages
def custom_main_message(header, text, color_code = '\033[1;33m'):
    return f'\n{color_code}{header} {text} {header} {color_code}\033[0m\n'

print(custom_main_message('===', 'Setting up environments'))
logger.info('Step1: Setting up environment')

# check if conda is in PATH
try:
    subprocess.run(['conda'],
                   check=True, capture_output=True, text=True)
    logger.info('Conda detected in PATH')
except FileNotFoundError:
    logger.error('Conda is not in PATH. Please ensure that Conda has been installed and added to PATH')
    exit(1)

# check if mamba is in PATH
try:
    subprocess.run(['mamba'],
                   check=True, capture_output=True, text=True)
    logger.info('Mamba detected in PATH')
except FileNotFoundError:
    logger.info('Installing Mamba')
    try:
        subprocess.run(
            ['conda', 'install', '-y', '-c', 'conda-forge', 'mamba=1.5.2'],
            check=True, capture_output=True, text=True
            )
        logger.info('Mamba has been installed')
    except:
        logger.error('Mamba cannot be installed. To avoid this error, you can manually install Mamba and add it to PATH')
        exit(1)


# install dependencies
def install_conda_env(env_name, yml_name, message):
    env_dir = os.path.join(SCRIPT_DIR, f'../conda/{env_name}')
    
    if os.path.exists(env_dir):
        logger.info(f'{message} detected')
    else:
        logger.info(f'Installing {message}')
        try:
            subprocess.run(
                ['mamba', 'env', 'create', '-f',
                 f'{SCRIPT_DIR}/../../envs/{yml_name}.yml',
                 '--prefix', env_dir],
                check=True, capture_output=True, text=True
            )
            logger.info(f'{message} has been installed')
        except:
            logger.error(f'{message} cannot be installed')
            exit(1)

install_conda_env(env_name='hmmer_pfam_annot_R_pcks', yml_name='R_pcks', message='R packages')
install_conda_env(env_name='hmmer_pfam_annot_snakemake', yml_name='snakemake', message='Snakemake')
install_conda_env(env_name='hmmer_pfam_annot_python3', yml_name='python3_12', message='Python v3.12')
install_conda_env(env_name='hmmer_pfam_annot_hmmer', yml_name='hmmer', message='HMMER')

#===============#
# All functions #
#===============#

def run_smk(path2smk, rerun_triggers=None):
    try:
        cmd = f'''
        source {CONDA_PATH}
        conda activate {SCRIPT_DIR}/../conda/hmmer_pfam_annot_snakemake

        snakemake -s {path2smk} \
        --configfile={CONFIG_DIR}/config.yaml \
        --cores {THREADS} \
        --latency-wait 120'''

        if rerun_triggers is not None:
            cmd += f' --rerun-triggers {rerun_triggers}'
        
        subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    
    # continue the running of incomplete intermediate files
    except subprocess.CalledProcessError:
        try:
            cmd_unlock = cmd + f' --unlock'
            subprocess.run(cmd_unlock, shell=True, capture_output=True, text=True, check=True)
            cmd_unlock = cmd + f' --rerun-incomplete'
            subprocess.run(cmd_unlock, shell=True, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(e.stderr)
            exit(1)


def run_smk_in_chunks(working_dir, smk_file):

    qfile_chunks = glob.glob(working_dir, recursive=True)

    qfile_chunk_dirs = {os.path.dirname(qfile_chunk) for qfile_chunk in qfile_chunks}

    count_qfile_chunk_dirs = len(qfile_chunk_dirs)
    job_number = 1

    for qfile_chunk_dir in qfile_chunk_dirs:
        os.chdir(qfile_chunk_dir)
        logger.info(f'Running chunk {job_number} of {count_qfile_chunk_dirs}')
        #! use mtime as the only rerun-trigger
        run_smk(smk_file, rerun_triggers='mtime')
        job_number += 1


def read_fasta_file(file_path):
    """Read the contents of a fasta file"""
    try:
        with open(file_path, 'r') as f:
            logger.info(f'Successfully read {file_path}')
            content = f.read().strip().splitlines()

            return '\n'.join(content)
    except Exception as e:
        logger.warning(f'Failed to read {file_path}: {e}')
        return ''


def get_FSlink_record(input, outdir):
    """Create a record contains file name and sequence header link for each file"""
    filename = os.path.basename(input)
    filename_wo_extension = os.path.splitext(filename)[0]
    path2out = os.path.join(outdir, f'{filename_wo_extension}_fslink.tsv')

    try:
        with open(input, 'r') as in_fh, open(path2out, 'w') as out_fh:
            out_fh.write('cds_id' + '\t' + 'filename' + '\n')
            for line in in_fh:
                if line.startswith('>'):
                    cds_id = re.sub('>', '', line).strip().split()[0]
                    out_fh.write(cds_id + '\t' + filename + '\n')
    except Exception as e:
        logger.warning(f'Failed to create FSlink for {input}: {e}')


#========================#
# Create smk CONFIG yaml #
#========================#

default_rerun_triggers = ['mtime', 'params', 'input', 'software-env', 'code']

# get the conda info
conda_info = subprocess.run(['conda', 'info'], check=True, capture_output=True, text=True)
con_env_line = [line for line in conda_info.stdout.split('\n') if 'base environment' in line][0]
CONDA_PATH = con_env_line.split()[3]

# handling unusal cases
if 'anaconda' in CONDA_PATH or 'miniconda' in CONDA_PATH:
    CONDA_PATH = re.sub(r'(^.*/(miniconda|anaconda)[^/]+).*$', r'\1/etc/profile.d/conda.sh', CONDA_PATH)
else:
    CONDA_PATH = CONDA_PATH + '/etc/profile.d/conda.sh'

# create YAML content
yaml_content = f'''\
path2query_seqs: {INPUT}
path2outputs: {OUTPUT_DIR}
path2refdb: {DATABASE_DIR}
path2scripts: {SCRIPT_DIR}
path2conda: {CONDA_PATH}
'''

# write YAML content to config.yaml file
CONFIG_DIR = os.path.join(OUTPUT_DIR, 'config')
os.makedirs(CONFIG_DIR, exist_ok=True)

with open(f'{CONFIG_DIR}/config.yaml', 'w') as file:
    file.write(yaml_content)

#====================#
# Protein annotation #
#====================#

logger.info('Processing started...')

checkpoint_dir = os.path.join(OUTPUT_DIR, 'checkpoints')
fslink_dir = os.path.join(checkpoint_dir, 'FSlink')
os.makedirs(checkpoint_dir, exist_ok=True)
os.makedirs(fslink_dir, exist_ok=True)

# For multiple query files in directory
if os.path.isdir(INPUT):
    files = [os.path.join(INPUT, f) for f in os.listdir(INPUT)
                if os.path.isfile(os.path.join(INPUT, f)) and f.endswith(('.fasta', '.faa', '.fa', '.fna'))]

    if not NO_INFILE_NAME:
        print(custom_main_message('===', 'Creating FSlink records'))
        logger.info(f'Create FSlink for {len(files)} FASTA files')
        for file in files:
            get_FSlink_record(file, fslink_dir)

    logger.info(f"Concatenate FASTA files in directory: {INPUT}")
    combined_content = []

    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = [executor.submit(read_fasta_file, file) for file in files]

        for future in futures:
            try:
                content = future.result()
                if content:  # Only add non-empty content
                    combined_content.append(content.strip())
            except Exception as e:
                logging.error(f"read_fasta_file executor failed: {e}")
    
    # Write the combined content to a file
    with open(f'{checkpoint_dir}/query.faa', 'w') as combined_file:
        combined_file.write('\n'.join(combined_content) + '\n')

# For single query file
elif os.path.isfile(INPUT) and INPUT.endswith(('.fasta', '.faa', '.fa', '.fna')):
    shutil.copy(INPUT, f'{checkpoint_dir}/query.faa')
    logger.info(f"Using provided FASTA file: {INPUT}")

    if not NO_INFILE_NAME:
        print(custom_main_message('===', 'Creating FSlink records'))
        logger.info(f'Create FSlink for single FASTA file')
        get_FSlink_record(INPUT, os.path.join(checkpoint_dir, 'FSlink'))
else:
    logger.error(f"{INPUT} is not a valid directory or FASTA file")
    exit(1)

print(custom_main_message('===', 'Splitting coding sequences into chunks'))
logger.info('Splitting coding sequences into chunks')
os.makedirs(f'{OUTPUT_DIR}/checkpoints/split_combined_mge_aa_into_chunks', exist_ok=True)
os.chdir(f'{OUTPUT_DIR}/checkpoints/split_combined_mge_aa_into_chunks')
run_smk(
    f'{SCRIPT_DIR}/split_combined_mge_aa_into_chunks.smk',
    rerun_triggers='mtime')

print(custom_main_message('===', 'Annotating proteins using the Pfam database'))
logger.info('Annotating proteins using the Pfam database')
run_smk_in_chunks(
    f'{OUTPUT_DIR}/checkpoints/candid_mge_aa_chunks/**/*.smk',
    'mge_protein_annot.smk')

print(custom_main_message('===', 'Joining annotation tables'))
logger.info('Joining annotation tables')
os.makedirs(f'{OUTPUT_DIR}/checkpoints/join_annot', exist_ok=True)
os.chdir(f'{OUTPUT_DIR}/checkpoints/join_annot')

run_smk(
    f'{SCRIPT_DIR}/join_annot.smk',
    rerun_triggers='mtime')

logger.info('Pipeline finished')
