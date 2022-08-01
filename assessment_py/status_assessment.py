'''
--------------------------------------------------------------------------------
      Function to get and plot information from status files (Seaglider format)
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 02/2022


Description:
This code will retrieve and plot information from status files recorded by the
Segaglider.


More information in:
https://hugepdf.com/download/1ka-seaglider-users-guide-download_pdf

'''

import matplotlib.pyplot as plt
import pandas as pd
import obspy
from obspy import read,read_inventory, UTCDateTime, Stream
import os
import glob
import json
import numpy as np
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
import matplotlib.dates as mdates
import matplotlib as mpl
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from obspy.signal import PPSD
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# ============================
# Function to find status data
# ============================

def filelist(basedir):
    """
    Returns the list of .nc files in *basedir*
    """
    
    files_nc = []
    for root, dirs, files in os.walk(basedir):
        for f in files:
            if f.endswith('.nc'):
                files_nc.append(os.path.join(root, f))

    files_ncs = sorted(files_nc)


    return files_ncs


