#!/usr/bin/env python
# coding: utf-8

import os
import glob
from datetime import datetime,timedelta,date
from obspy import read,UTCDateTime,Trace
import numpy as np

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

mseed_files = '/run/media/dIOGOLOC/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/OUTPUT/MSEED/'

FOLDER_OUTPUT = '/home/dIOGOLOC/dados_posdoc/gliders_project/OUTPUT/'

# ==========================================================
# Calculating datetime between INITIAL_DATE and  FINAL_DATE
# ==========================================================

datatime_initial = datetime.strptime('2015-01-01', "%Y-%m-%d").date() 

datatime_final = datetime.strptime('2022-12-31', "%Y-%m-%d").date() 

datetime_lista = np.arange(datatime_initial, datatime_final, timedelta(days=1)).astype(datetime)

datetime_lista_months = sorted(list(set([i.strftime("%Y-%m") for i in datetime_lista])))

datetime_lista_years = sorted(list(set([i.strftime("%Y") for i in datetime_lista])))

xlim_initial = mdates.date2num(datatime_initial)
xlim_final = mdates.date2num(datatime_final)

# ========
# Function
# ========

def dataframe_extraction_from_mseedfile(i):
        '''
        i: .mseed file.
        '''
                  
        subdir, filename_mseed = os.path.split(i)
        filename = filename_mseed.split('.mseed')[0]
        if 'pa' in filename.split('_')[0]:
            mergulho = filename.split('_')[0].split('a')[1]
            stream_number = filename.split('_')[1]

            year_month_day = filename.split('_')[2]
            hour_minute_second = filename.split('_')[3]

            year = int('20'+year_month_day[:2])
            month = int(year_month_day[2:4])
            day = int(year_month_day[4:])

            hour = int(hour_minute_second[:2])
            minute = int(hour_minute_second[2:4])
            second = int(hour_minute_second[4:])

            d = UTCDateTime(datetime(year,month,day,hour,minute,second).isoformat())


        if 'pa' in filename.split('_')[2]:

            mergulho = filename.split('_')[2].split('a')[1]
            stream_number = filename.split('_')[3]

            year_month_day = filename.split('_')[0]
            hour_minute_second = filename.split('_')[1]

            year = int('20'+year_month_day[:2])
            month = int(year_month_day[2:4])
            day = int(year_month_day[4:])

            hour = int(hour_minute_second[:2])
            minute = int(hour_minute_second[2:4])
            second = int(hour_minute_second[4:])

            d = UTCDateTime(datetime(year,month,day,hour,minute,second).isoformat())
        
        #----------------------------
        #Starting Dataframe

        starttime = d.datetime

        df = pd.DataFrame([[filename],[mergulho],[stream_number],[starttime],[str(year)],[year_month_day[2:4]],[starttime.strftime("%b")]], index=['filename', 'mergulho', 'stream_number','starttime','year','month','name_month']).T
        #Ending Dataframe
        #----------------------------
        
        return df

#----------------------------

# =======
# Program
# =======

start_time = time.time()

mseed_files_lst = sorted(glob.glob(mseed_files+'*/*/*.mseed'))
mseed_files_lst = mseed_files_lst
# ================================
# Calculating waveforms parameters
# ================================

df_lst = []

# create and configure the process pool
with Pool(processes=6) as pool:
    # execute tasks
    for result in tqdm(pool.imap_unordered(dataframe_extraction_from_mseedfile, mseed_files_lst),total=len(mseed_files_lst), desc='MSEED files processing'):
        df_lst.append(result)
# process pool is closed automatically

dataframe_final = pd.concat(df_lst, ignore_index=True)
month_date_lst = sorted(list(set(dataframe_final['month'].values)))

# creating the array to plot
dataframe_lista = []
for h in tqdm(np.arange(1,13),total=len(np.arange(1,13)), desc='Creating the dataframe:'):
    try:

        df_month = dataframe_final[dataframe_final['month'] == str(h).zfill(2)]
            
        MONTH_LST = df_month['name_month'].values
        df_temp = [str(h).zfill(2),list(set(MONTH_LST))[0]]
        df_temp_index = ['number_month','name_month']
        for i in datetime_lista_years:
            if df_month['year'].str.contains(i).any():

                df_year = df_month[df_month['year'] == i]
                NUMBER_MINUTES_LST = df_year['year'].value_counts().values

                df_temp.append(NUMBER_MINUTES_LST[0])
                df_temp_index.append(i)
            else:
                df_temp.append(0)
                df_temp_index.append(i) 

        dataframe_lista.append(pd.DataFrame(df_temp, index=df_temp_index).T)
    except:
        pass

df_to_plot = pd.concat(dataframe_lista, ignore_index=True)

name_months = df_to_plot['name_month'].values
data_x_axis = df_to_plot[datetime_lista_years].values.astype(float).T
#-------------------------

# ==========================
# Plotting DATA availability
# ==========================
#x axis parameters

days = DayLocator(interval=1)  # every 1 day
months = MonthLocator(interval=1)  # every 1 month
monthsFmt = DateFormatter('%b-%y')

#Matplotlib parameters
fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(16,9))

im = ax.pcolormesh(data_x_axis,cmap='viridis', shading ='flat',ec='k')

# Get the dimensions of the array
rows, cols = data_x_axis.shape

# Loop over each cell and add the cell value as text
for i in range(rows):
    for j in range(cols):
        cell_value = int(data_x_axis[i, j])
        text_color = 'white' if cell_value < np.max(data_x_axis)/2 else 'black'
        ax.text(j + 0.5, i + 0.5, str(cell_value), color=text_color,ha='center', va='center')

# Set the tick locations and labels
# Set the x and y axis tick locations and labels
ax.set_xticks(np.arange(data_x_axis.shape[1]) + 0.5, name_months,fontsize=15)
ax.set_yticks(np.arange(data_x_axis.shape[0]) + 0.5, datetime_lista_years,fontsize=15)

ax.tick_params(which='minor', length=2)
ax.tick_params(which='major', length=10)
ax.set_aspect(1)

plt.setp(ax.xaxis.get_majorticklabels(), fontsize=12,rotation=30)

#criando a localização da barra de cores:
axins = inset_axes(ax,
                    width="15%",  # width = 15% of parent_bbox width
                    height="2.5%",  # height : 2.5%
                    loc='upper left',
                    bbox_to_anchor=(0.0, 0.05, 1, 1),
                    bbox_transform=ax.transAxes,
                    borderpad=0,
                    )
cbar = fig.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top',label='Arquivos/Mês')
os.makedirs(FOLDER_OUTPUT+'/FIGURAS/',exist_ok=True)
fig.savefig(FOLDER_OUTPUT+'/FIGURAS/'+'COMPLETENESS_'+datatime_initial.strftime("%Y")+'_'+datatime_final.strftime("%Y")+'_compact_mseed.png',dpi=300)
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')
