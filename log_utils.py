#!/public/home/tianyi/tools/miniforge3/bin/python

import os
import sys
import logging
import colorlog
from datetime import datetime

def setup_logging(log_dir):
    """Set up logging configuration with color and timestamp in filename
    
    Args:
        log_dir (str): Directory where log files will be stored
        
    Returns:
        str: Path to the created log file
    """
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    log_file = os.path.join(log_dir, f'{script_name}_{timestamp}.log')
    
    # Create logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # Remove any existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Console handler with color
    console_handler = colorlog.StreamHandler()
    console_handler.setFormatter(colorlog.ColoredFormatter(
        '%(log_color)s%(message)s%(reset)s',
        log_colors={
            'DEBUG': 'cyan',
            'INFO': 'green',
            'WARNING': 'yellow',
            'ERROR': 'red',
            'CRITICAL': 'red,bg_white',
        }
    ))
    
    # File handler without color
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))
    
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    
    return log_file

def log_message(message, level=logging.INFO):
    """Helper function to log messages with proper color
    
    Args:
        message (str): Message to log
        level (int): Logging level (default: logging.INFO)
    """
    logger = logging.getLogger()
    if level == logging.INFO:
        logger.info(message)
    elif level == logging.ERROR:
        logger.error(message)
    elif level == logging.WARNING:
        logger.warning(message)
