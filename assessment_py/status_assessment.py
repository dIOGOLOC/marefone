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
import time
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool
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

from parameters_py.config import (
					OUTPUT_FIGURE_DIR,DIR_STATUS_FILES,NUM_PROCESS,LABEL_LANG,OUTPUT_OUT_DIR
					)

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


# ==========================
# Function to plot wav files
# ==========================

def date_extract(ielement):
    '''
    Function to extract the UTCDateTime from the *.wav file names.
    '''
    # splitting subdir/basename
    subdir, filename = os.path.split(str(ielement))
    name_file = filename.split('.wav')[0]
    if 'pa' and 'au' in name_file:
        wav_file_name = name_file.split('_')
        if 'pa' in wav_file_name[0] and len(wav_file_name) > 3:
            return UTCDateTime(year=int('20'+wav_file_name[2][:2]),month=int(wav_file_name[2][2:4]),day=int(wav_file_name[2][4:6]),hour=int(wav_file_name[3][:2]),minute=int(wav_file_name[3][2:4]),second=int(wav_file_name[3][4:6]))
        if len(wav_file_name) > 3 and 'pa' in wav_file_name[2]:
            return UTCDateTime(year=int('20'+wav_file_name[0][:2]),month=int(wav_file_name[0][2:4]),day=int(wav_file_name[0][4:6]),hour=int(wav_file_name[1][:2]),minute=int(wav_file_name[1][2:4]),second=int(wav_file_name[1][4:6]))
         
# -----------------------------------------------------------------------
         
def check_24_hour(input):
    '''
    Function to check how many hour have in a single day.
    '''
    if len(input) > 1:
        temp_24 = [1 if h in input else 0 for h in np.arange(24)]

    if len(input) == 1:
        temp_24 = [1 if input[0] == h else 0 for h in np.arange(24)]
    
    if len(input) == 0:
        temp_24 = list(np.zeros_like(np.arange(24)))

    return temp_24

# -----------------------------------------------------------------------

def create_date_day_array(input):
    '''
    Function to create a list of daily datetime arrays.
    '''
    lst_temp = [dattime_l.hour for dattime_l in input[1] if input[0].date() == dattime_l.datetime.date()]
    lst_temp_24 = check_24_hour(list(set(lst_temp)))

    return lst_temp_24
    
# -----------------------------------------------------------------------

def READ_PLOT_DATE_file(FIG_FOLDER_OUTPUT,TXT_FILE):

    if LABEL_LANG == 'br':
        print('Procurando por dados no arquivo TXT: '+TXT_FILE)

    else:
        print('Looking for data in the TXT file: '+TXT_FILE)

    data_lis = np.genfromtxt(TXT_FILE,delimiter=',',dtype='str')

    #----------------------------
    pool = Pool(processes=NUM_PROCESS)
    date_lst = []
    for result in tqdm(pool.imap(func=date_extract, iterable=data_lis), total=len(data_lis)):
        if isinstance(result, type(None)) == False: 
            date_lst.append(result)
    #----------------------------
    datatime_initial = min(date_lst).datetime
    datatime_final = max(date_lst).datetime
    datetime_lista = np.arange(datatime_initial, datatime_final, datetime.timedelta(days=1)).astype(datetime.datetime)

    #Contador de horas:
    input_create_date_day = []
    for k,l in enumerate(datetime_lista):
        input_create_date_day.append([l,date_lst])

    #--

    pool = Pool(processes=NUM_PROCESS)
    DAY_lst_dates = []
    for result in tqdm(pool.imap(func=create_date_day_array, iterable=input_create_date_day), total=len(input_create_date_day)):
        if isinstance(result, type(None)) == False: 
            DAY_lst_dates.append(result)
    #----------------------------
   
    data_x_axis = np.array(DAY_lst_dates).T

    xlim_initial = mdates.date2num(datatime_initial)
    xlim_final = mdates.date2num(datatime_final)

    # ====================================
    # Function to plot DATA availability
    # ====================================

    month_1 = DayLocator(interval=1)   # every month
    months = MonthLocator(interval=3)  # every 3-month
    yearsFmt = DateFormatter('%Y-%m-%d')

    #Matplotlib parameters
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(40,10))

    im = ax.imshow(data_x_axis,origin='lower',extent = [xlim_initial,xlim_final,0,24],cmap=plt.cm.Greens,interpolation=None, vmin=0, vmax=1)
    ax.set_xlim(datatime_initial,datatime_final)
    ax.set_aspect(10)
    ax.yaxis.set_major_locator(MultipleLocator(4))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(months)
    ax.xaxis.set_major_formatter(yearsFmt)
    ax.xaxis.set_minor_locator(month_1)
    ax.tick_params(which='minor', length=4)
    ax.tick_params(which='major', length=10)
    ax.set_ylim(0,24)
    ax.set_ylabel('Hora do dia')
    ax.grid(b=True, which='major', color='k', linestyle='-')
    ax.grid(b=True, which='minor', color='k', linestyle='-')
    plt.setp(ax.xaxis.get_majorticklabels(), fontsize=10, rotation=30)

    os.makedirs(FIG_FOLDER_OUTPUT,exist_ok=True)
    fig.savefig(FIG_FOLDER_OUTPUT+'COMPLETENESS_'+str(datatime_initial.year)+'_'+str(datatime_final.year)+'.png',dpi=300)