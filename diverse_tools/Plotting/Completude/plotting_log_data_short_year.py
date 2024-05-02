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

log_file_folder = '/run/media/dIOGOLOC/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/gliders_data/LOG_data/'

FOLDER_OUTPUT = '/home/dIOGOLOC/dados_posdoc/gliders_project/OUTPUT/'

# ==========================================================
# Calculating datetime between INITIAL_DATE and  FINAL_DATE
# ==========================================================

datatime_initial = datetime.strptime('2015-01-01', "%Y-%m-%d").date() 

datatime_final = datetime.strptime('2019-12-31', "%Y-%m-%d").date() 

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

def dm(x):
    south = False
    
    if x<0:
        south = True
        x = abs(x)
        
    degrees = int(x) // 100
    minutes = x - 100*degrees
    
    x = degrees + minutes/60
    if south:
        x = -x
    return x

#----------------------------

# =======
# Program
# =======

start_time = time.time()

log_files = sorted(glob.glob(log_file_folder+'*/*.log'))

log_files_lst = sorted(log_files)

# ================================
# Calculating waveforms parameters
# ================================

pandas_lst = []
for file in log_files_lst:
    glider_dic = {'file_location':[],'id':[],"campanha":[],'mergulho':[],"time_d_UTC":[],"time_d":[],"time_s":[],"time_s_UTC":[],"lat_s":[],"lon_s":[],"lat_d":[],"lon_d":[],"surface_drift_direction_deg_s":[],"surface_drift_speed_knots_s":[],"surface_drift_direction_deg_d":[],"surface_drift_speed_knots_d":[]}

    with open(file,'r') as f:
        for x in f:
            
            if '$MISSION,' in x:
                glider_dic['file_location'].append(file.split('/')[0])
                glider_dic['campanha'].append(float(x.split('$MISSION,')[1].split('\n')[0]))

            if '$DIVE,' in x:
                glider_dic['mergulho'].append(int(x.split('$DIVE,')[1].split('\n')[0]))
                
            if '$ID,' in x:
                glider_dic['id'].append(int(x.split('$ID,')[1].split('\n')[0]))
    
            if '$GPS2,' in x:
                glider_dic['time_d_UTC'].append(UTCDateTime(year=int('20'+x.split('$GPS2,')[1].split(',')[0][4:6]), month=int(x.split('$GPS2,')[1].split(',')[0][2:4]), day=int(x.split('$GPS2,')[1].split(',')[0][0:2]),hour=int(x.split('$GPS2,')[1].split(',')[1][0:2]), minute=int(x.split('$GPS2,')[1].split(',')[1][2:4]), second=int(x.split('$GPS2,')[1].split(',')[1][4:6])).datetime)
                glider_dic['time_d'].append(mdates.date2num(UTCDateTime(year=int('20'+x.split('$GPS2,')[1].split(',')[0][4:6]), month=int(x.split('$GPS2,')[1].split(',')[0][2:4]), day=int(x.split('$GPS2,')[1].split(',')[0][0:2]),hour=int(x.split('$GPS2,')[1].split(',')[1][0:2]), minute=int(x.split('$GPS2,')[1].split(',')[1][2:4]), second=int(x.split('$GPS2,')[1].split(',')[1][4:6])).datetime))
                glider_dic['lat_d'].append(dm(float(x.split('$GPS2,')[1].split(',')[2])))
                glider_dic['lon_d'].append(dm(float(x.split('$GPS2,')[1].split(',')[3])))
                glider_dic['surface_drift_direction_deg_d'].append(float(x.split('$GPS2,')[1].split(',')[-3]))
                glider_dic['surface_drift_speed_knots_d'].append(float(x.split('$GPS2,')[1].split(',')[-4]))         

            if '$GPS,' in x:
                glider_dic['time_s_UTC'].append(UTCDateTime(year=int('20'+x.split('$GPS,')[1].split(',')[0][4:6]), month=int(x.split('$GPS,')[1].split(',')[0][2:4]), day=int(x.split('$GPS,')[1].split(',')[0][0:2]),hour=int(x.split('$GPS,')[1].split(',')[1][0:2]), minute=int(x.split('$GPS,')[1].split(',')[1][2:4]), second=int(x.split('$GPS,')[1].split(',')[1][4:6])).datetime)
                glider_dic['time_s'].append(mdates.date2num(UTCDateTime(year=int('20'+x.split('$GPS,')[1].split(',')[0][4:6]), month=int(x.split('$GPS,')[1].split(',')[0][2:4]), day=int(x.split('$GPS,')[1].split(',')[0][0:2]),hour=int(x.split('$GPS,')[1].split(',')[1][0:2]), minute=int(x.split('$GPS,')[1].split(',')[1][2:4]), second=int(x.split('$GPS,')[1].split(',')[1][4:6])).datetime))
                glider_dic['lat_s'].append(dm(float(x.split('$GPS,')[1].split(',')[2])))
                glider_dic['lon_s'].append(dm(float(x.split('$GPS,')[1].split(',')[3])))
                glider_dic['surface_drift_direction_deg_s'].append(float(x.split('$GPS,')[1].split(',')[-3]))
                glider_dic['surface_drift_speed_knots_s'].append(float(x.split('$GPS,')[1].split(',')[-4]))
    if glider_dic['time_s'] != []:
        df = pd.DataFrame.from_dict(glider_dic)
        pandas_lst.append(df)                

dataframe_final = pd.concat(pandas_lst, ignore_index=True)

dataframe_final['year'] = dataframe_final['time_d_UTC'].dt.strftime("%Y")
dataframe_final['month'] = dataframe_final['time_d_UTC'].dt.strftime("%m")
dataframe_final['name_month'] = dataframe_final['time_d_UTC'].dt.strftime("%b")

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
fig.savefig(FOLDER_OUTPUT+'/FIGURAS/'+'COMPLETENESS_'+datatime_initial.strftime("%Y")+'-'+datatime_final.strftime("%Y")+'compact_log.png',dpi=300)
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')
