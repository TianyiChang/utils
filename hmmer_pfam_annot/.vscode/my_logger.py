import logging
import os
import subprocess

# locate this script
SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))

# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# create a console handler
ch = logging.StreamHandler()
# ch.setLevel(logging.WARNING)

# create formatter
formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",
                              "%Y-%m-%d %H:%M:%S")

# add formatter to handler
ch.setFormatter(formatter)
# fh.setFormatter(formatter)

# add handler to logger
logger.addHandler(ch)
# logger.addHandler(fh)

# check if colorlog in installed
try:
    from colorlog import ColoredFormatter
except ModuleNotFoundError:
    logger.info('Installing colorlog via pip install')
    try:
        subprocess.run(['pip', 'install', 'colorlog'], capture_output=True)
    except FileExistsError:
        logger.warning('pip is not found in PATH. No worries, the pipeline will use the default logging module')
    except subprocess.CalledProcessError:
        logger.warning('colorlog cannot be installed via pip. No worries, the pipeline will use the default logging module')


# update console Formatter
try:
    from colorlog import ColoredFormatter
    colored_formatter = ColoredFormatter(
        "%(log_color)s%(asctime)s [%(levelname)s] %(message)s",
                                "%Y-%m-%d %H:%M:%S",
        reset=True,
        log_colors={
            'DEBUG': 'cyan',
            'INFO': 'green',
            'WARNING': 'yellow',
            'ERROR': 'red',
            'CRITICAL': 'red,bg_white',
        })

    ch.setFormatter(colored_formatter)
    logger.addHandler(ch)

    logger.info('colorlog has been installed')
except:
    pass
