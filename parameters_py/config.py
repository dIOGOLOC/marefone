"""
--------------------------------------------------------------------------------
         Module that parses global parameters from a configuration file
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 07/2022


Description:
Module that parses global parameters from a configuration file at first import,
to make them available to the other parts of the program.

More information in:
https://wiki.python.org/moin/ConfigParserExamples

Input:
Configuration file, wherein global paths and parameters are defined.

Outputs:
The module provides a parser for simple configuration files consisting of groups
of named values.

"""

import configparser
import os
import glob
import json


def select_and_parse_config_file(basedir='.', ext='cnf', verbose=True):
    """
    Reads a configuration file and returns an instance of ConfigParser:
    First, looks for files in *basedir* with extension *ext*.
    Asks user to select a file if several files are found,
    and parses it using ConfigParser module.
    @rtype: L{ConfigParser.ConfigParser}
    """
    config_files = glob.glob(os.path.join(basedir, u'*.{}'.format(ext)))


    if not config_files:
        raise Exception("No configuration file found!")

    if len(config_files) == 1:
        # only one configuration file
        config_file = config_files[0]
    else:
        print("Select a configuration file:")
        for i, f in enumerate(config_files, start=1):
            print("{} - {}".format(i, f))
        res = int(input(''))
        config_file = config_files[res - 1]

    if verbose:
        print("Reading configuration file: {}".format(config_file))

    conf = configparser.ConfigParser(allow_no_value=True)
    conf.read(config_file)

    return conf

# ==========================
# parsing configuration file
# ==========================

config = select_and_parse_config_file(basedir='.', ext='cnf', verbose=True)

# -----
# lang
# -----

#choose between portuguese (br) or english (en):
LABEL_LANG = config.get('lang', 'LABEL_LANG')

# -----
# paths
# -----

##
# inputs
##

# directory of raw files
DIR_DATA_FILES = config.get('paths', 'DIR_DATA_FILES')

# directory of raw status files
DIR_STATUS_FILES = config.get('paths', 'DIR_STATUS_FILES')

#Shapefile  boundary area
BOUNDARY_STATES_SHP = config.get('paths', 'BOUNDARY_STATES_SHP')

##
# outputs
##

#Directory to save Figures
OUTPUT_FIGURE_DIR = config.get('paths', 'OUTPUT_FIGURE_DIR')

#Directory to save processed data
OUTPUT_OUT_DIR = config.get('paths', 'OUTPUT_OUT_DIR')

# -------
# process
# -------

#Number worker processes (see https://docs.python.org/3/library/multiprocessing.html)
NUM_PROCESS = config.getint('process', 'NUM_PROCESS')