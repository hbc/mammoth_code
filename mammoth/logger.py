import logging
import requests
requests.packages.urllib3.disable_warnings()

from mammoth import utils
try:
    from colorlog import ColoredFormatter
except:
    pass
import os


def getLogger(name="mammoth"):
    return logging.getLogger(__name__)

def set_format(frmt, frmt_col=None, datefmt=None):
    if frmt_col:
        try:
            formatter = ColoredFormatter(frmt_col, datefmt=datefmt)
            return formatter
        except:
            pass
    formatter = logging.Formatter(frmt)
    return formatter

def initialize_logger(output_dir, debug, level=False):
    NOTE = 15
    COLOR_FORMAT = "%(log_color)s%(asctime)s%(levelname)s-%(name)s(%(lineno)d)%(reset)s: %(message)s"
    COLOR_FORMAT_INFO = "%(log_color)s%(asctime)s %(levelname)s%(reset)s: %(message)s"
    DATE_FRT = '%m/%d/%Y %I:%M:%S %p'
    FORMAT = "%(levelname)s-%(name)s(%(lineno)d): %(message)s"
    FORMAT_INFO = "%(levelname)s %(message)s"
    if output_dir:
        if not os.path.exists(output_dir):
            utils.safe_makedir(output_dir)

    numeric_level = getattr(logging, "DEBUG", None)
    if not debug:
        numeric_level = getattr(logging, "INFO", None)
    logger = logging.getLogger("mammoth")
    logger.setLevel(numeric_level)
    # create console handler and set level to info
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    if level:
        handler.setLevel(logging.DEBUG)

    formatter = set_format(FORMAT_INFO, COLOR_FORMAT_INFO, datefmt=DATE_FRT)
    handler.setFormatter(formatter)
    logger.addHandler(handler)


    if output_dir:
        # create error file handler and set level to error
        handler = logging.FileHandler(os.path.join(output_dir, "error.log"), "w", encoding=None)
        handler.setLevel(logging.ERROR)
        formatter = set_format(FORMAT, COLOR_FORMAT, datefmt=DATE_FRT)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        # create debug file handler and set level to debug
        handler = logging.FileHandler(os.path.join(output_dir, "run.log"), "w")
        handler.setLevel(numeric_level)
        formatter = logging.Formatter(FORMAT)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        requests_log = logging.getLogger("requests")
        requests_log.setLevel(logging.DEBUG)
        requests_log.addHandler(handler)
        requests_log.propagate = True
