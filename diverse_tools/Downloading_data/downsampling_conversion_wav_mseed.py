#!/usr/bin/env python
# coding: utf-8

import os
import glob
from datetime import datetime,timedelta,date
from obspy import read,UTCDateTime,Trace
import numpy as np
from scipy.io import wavfile

import pandas as pd
from multiprocessing import Pool, RLock, freeze_support
from tqdm.auto import tqdm
import time

import matplotlib.dates as mdates
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ======
# Config
# ======

wav_file_dir = '/home/diogoloc/dados_posdoc/gliders_project/gliders_data/WAV_DATA/'

FOLDER_OUTPUT = '/home/diogoloc/dados_posdoc/gliders_project/OUTPUT/'

# ========
# Function
# ========

def dataframe_extraction_from_wavfile(i):
    '''
    i: .wav file file.
    '''

    try:
          
        subdir, filename = os.path.split(i)
        mergulho = filename.split('.wav')[0].split('_')[0].split('a')[1]
        stream_number = filename.split('.wav')[0].split('_')[1]

        year_month_day = filename.split('.wav')[0].split('_')[2]
        hour_minute_second = filename.split('.wav')[0].split('_')[3]

        year = int('20'+year_month_day[:2])
        month = int(year_month_day[2:4])
        day = int(year_month_day[4:])

        hour = int(hour_minute_second[:2])
        minute = int(hour_minute_second[2:4])
        second = int(hour_minute_second[4:])

        d = UTCDateTime(datetime(year,month,day,hour,minute,second).isoformat())
        #----------------------------
        #Starting Dataframe

        st = read(i,format='WAV',headonly=True)
        starttime = d.datetime
        endtime = (d+(st[0].stats.npts/st[0].stats.sampling_rate)).datetime
        n_minutes = round((endtime - starttime).total_seconds() / 60.0)
        hour_day = starttime.hour
        sampling_rate = st[0].stats.sampling_rate
        delta = st[0].stats.delta
        npts = st[0].stats.npts

        df = pd.DataFrame([[filename],[mergulho],[stream_number],[sampling_rate],[delta],[npts],[starttime],[endtime],[hour_day],[n_minutes]], index=['filename', 'mergulho', 'stream_number','sampling_rate','delta','npts','starttime','endtime','hour_day','number_of_minutes']).T
        #Ending Dataframe
        #----------------------------

        return df
    
    except:
        pass
#----------------------------

def downsampling_function(i, sampling_rate=100):
    """
    Decimate a trace to achieve the desired sampling rate, sr.

    NOTE: data will be detrended and a cosine taper applied before
    decimation, in order to avoid edge effects when applying the lowpass
    filter before decimating.

    source: https://quakemigrate.readthedocs.io/en/latest/_modules/quakemigrate/util.html#decimate

    Parameters:
    -----------
    i : .WAV file path
        Stream to be decimated.
    sampling_rate : int
        Output sampling rate.

    Returns:
    --------
    trace : `obspy.Trace.stats` object
        Decimated trace.

    """
    try:

        subdir, filename = os.path.split(i)
        mergulho = filename.split('.wav')[0].split('_')[0].split('a')[1]
        stream_number = filename.split('.wav')[0].split('_')[1]

        year_month_day = filename.split('.wav')[0].split('_')[2]
        hour_minute_second = filename.split('.wav')[0].split('_')[3]

        year = int('20'+year_month_day[:2])
        month = int(year_month_day[2:4])
        day = int(year_month_day[4:])

        hour = int(hour_minute_second[:2])
        minute = int(hour_minute_second[2:4])
        second = int(hour_minute_second[4:])

        d = UTCDateTime(datetime(year,month,day,hour,minute,second).isoformat())

        #----------------------------
        #Collecting wav data

        sampleratetr, datatr = wavfile.read(i)

        tr = Trace(data=datatr)
        tr.stats.sampling_rate = sampleratetr
        tr.stats.starttime = d 

        starttime = d.datetime
        endtime = (d+(tr.stats.npts/tr.stats.sampling_rate)).datetime
        n_minutes = round((endtime - starttime).total_seconds() / 60.0)
        hour_day = starttime.hour

        # Work on a copy of the trace
        trace = tr.copy()

        # Detrend and apply cosine taper
        trace.detrend('linear')
        trace.detrend('demean')
        trace.taper(type='cosine', max_percentage=0.05)

        # Zero-phase Butterworth-lowpass filter at Nyquist frequency
        trace.filter("lowpass", freq=float(sampling_rate) / 2.000001, corners=2,zerophase=True)
        trace.decimate(factor=int(trace.stats.sampling_rate / sampling_rate), strict_length=False, no_filter=True)

        OUTPUT_TRACE = FOLDER_OUTPUT+'/MSEED/'+d.strftime("%Y")+'/'+d.strftime("%Y-%m-%d")+'/'
        os.makedirs(OUTPUT_TRACE,exist_ok=True)
        trace.write(OUTPUT_TRACE+filename.split('.')[0]+'.mseed', format='MSEED')

        return trace.stats.sampling_rate

    except:
        pass
#----------------------------

def check_datetime_in_period(datetime_lst,dataf):
    '''
    Function to check if the dates in data set are inside the chosen time period
    
    '''
    array_to_plot_by_xlim = []
    for x,c in enumerate(datetime_lst):
        lista_temp = []
        for t,y in enumerate(dataf['DATETIME'].values):
            if y == c.date():
                lista_temp.append(np.array(dataf[dataf['DATETIME'] == y]['NUMBER_HOUR'].tolist()[0]))
        array_to_plot_by_xlim.append(lista_temp)
   
    data_x_axis = []
    for x,c in enumerate(array_to_plot_by_xlim):
        if c != []:
            data_x_axis.append(c[0])
        else:
            data_x_axis.append(np.zeros_like(np.arange(24)))

    data_x_axis = np.array(data_x_axis).T

    return data_x_axis
#----------------------------

# =======
# Program
# =======

start_time = time.time()

wav_files_lst = sorted(glob.glob(wav_file_dir+'*/*/*.wav'))

files_datetime_input = sorted(wav_files_lst)

# ================================
# Calculating waveforms parameters
# ================================

df_lst = []

# create and configure the process pool
with Pool() as pool:
    # execute tasks
    for result in tqdm(pool.imap_unordered(dataframe_extraction_from_wavfile, files_datetime_input),total=len(files_datetime_input), desc='WAV files processing'):
        df_lst.append(result)
# process pool is closed automatically

dataframe_final = pd.concat(df_lst, ignore_index=True)

dataframe_final['DayMonthYear'] = dataframe_final['starttime'].dt.strftime("%Y-%m-%d")

day_date_lst = sorted(list(set(dataframe_final['DayMonthYear'].values)))

# creating the array to plot
dataframe_lista = []

for day_date in  tqdm(day_date_lst,total=len(day_date_lst), desc='Calculatig minutes per hour'):
    NUMBER_HOUR = []
    df_day = dataframe_final[dataframe_final['DayMonthYear'] == day_date]
    for g,h in enumerate(np.arange(24)):
        df_hour = df_day[df_day['hour_day'] == h]
        n_minutes_day = df_hour['number_of_minutes'].sum()
        NUMBER_HOUR.append(n_minutes_day)
    
    dataframe_lista.append(pd.DataFrame([date.fromisoformat(day_date),NUMBER_HOUR], index=['DATETIME','NUMBER_HOUR']).T)

df_to_plot = pd.concat(dataframe_lista, ignore_index=True)

# ==========================================================
# Calculating datetime between INITIAL_DATE and  FINAL_DATE
# ==========================================================

datatime_initial = df_to_plot['DATETIME'].values[0]

datatime_final = df_to_plot['DATETIME'].values[-1]

datetime_lista = np.arange(datatime_initial, datatime_final, timedelta(days=1)).astype(datetime)

xlim_initial = mdates.date2num(datatime_initial)
xlim_final = mdates.date2num(datatime_final)

# ==========================
# Plotting DATA availability
# ==========================
#x axis parameters

days1 = DayLocator(interval=5)   # every day
months = MonthLocator()  # every month
monthsFmt = DateFormatter('%b-%y')
daysFmt = DateFormatter('%-d')

#Matplotlib parameters
fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(20,5))

data_x_axis = check_datetime_in_period(datetime_lista,df_to_plot)
datetime_pcolormesh = np.arange(datatime_initial, datatime_final+timedelta(days=1), timedelta(days=1)).astype(datetime)

im = ax.pcolormesh(datetime_pcolormesh,np.arange(25),data_x_axis,cmap='gist_heat_r', vmin=0, vmax=60,shading='flat',ec='k')
ax.set_xlim(datatime_initial,datatime_final)
ax.yaxis.set_major_locator(MultipleLocator(4))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(monthsFmt)
ax.xaxis.set_minor_locator(days1)
ax.xaxis.set_minor_formatter(daysFmt)
ax.tick_params(which='minor', length=2)
ax.tick_params(which='major', length=20)
ax.set_ylim(0,24)
ax.set_ylabel('Hora do Dia',fontsize=15)
ax.grid(visible=True, which='major', color='k', linestyle='-')
ax.grid(visible=True, which='minor', color='k', linestyle='-')

plt.setp(ax.xaxis.get_majorticklabels(), fontsize=20)
plt.setp(ax.xaxis.get_minorticklabels(), fontsize=15)

#criando a localização da barra de cores:
axins = inset_axes(ax,
                    width="10%",  # width = 15% of parent_bbox width
                    height="2.5%",  # height : 2.5%
                    loc='upper left',
                    bbox_to_anchor=(0.85, 0.05, 1, 1),
                    bbox_transform=ax.transAxes,
                    borderpad=0,
                    )
cbar = fig.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top',ticks=[0,20,40,60],label='Minutos/hora')

os.makedirs(FOLDER_OUTPUT+'/FIGURAS/',exist_ok=True)
fig.savefig(FOLDER_OUTPUT+'/FIGURAS/'+'COMPLETENESS_'+datatime_initial.strftime("%Y.%m.%d")+datatime_final.strftime("%Y.%m.%d")+'.png',dpi=300)

# =========================
# Converting wav into mseed
# =========================

sr_lst = []
# create and configure the process pool
with Pool(processes=3) as pool:
    # execute tasks
    for result in tqdm(pool.imap_unordered(downsampling_function, files_datetime_input),total=len(files_datetime_input), desc='Converting WAV --> MSEED'):
        sr_lst.append(result)
# process pool is closed automatically

print(list(set(sr_lst)))

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print(str(len(files_datetime_input))+' waveforms processed!')
print('\n')