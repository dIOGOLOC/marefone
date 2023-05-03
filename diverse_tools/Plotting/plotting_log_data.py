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

log_file_folder = '/home/diogoloc/dados_posdoc/gliders_project/gliders_data/LOG_data/'

FOLDER_OUTPUT = '/home/diogoloc/dados_posdoc/gliders_project/OUTPUT/'

# ========
# Function
# ========

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

dataframe_final['DayMonthYear'] = dataframe_final['time_d_UTC'].dt.strftime("%Y-%m-%d")
dataframe_final['hour_day'] = [i.hour for i in dataframe_final['time_d_UTC']]
dataframe_final['minute_day'] = [i.minute for i in dataframe_final['time_d_UTC']]

day_date_lst = sorted(list(set(dataframe_final['DayMonthYear'].values)))

# creating the array to plot
dataframe_lista = []

for day_date in  tqdm(day_date_lst,total=len(day_date_lst), desc='Calculatig minutes per hour'):
    NUMBER_HOUR = []
    df_day = dataframe_final[dataframe_final['DayMonthYear'] == day_date]
    for g,h in enumerate(np.arange(24)):
        df_hour = df_day[df_day['hour_day'] == h]
        n_minutes_day = df_hour['minute_day'].tolist()
        NUMBER_HOUR.append(len(n_minutes_day))
    
    dataframe_lista.append(pd.DataFrame([date.fromisoformat(day_date),NUMBER_HOUR], index=['DATETIME','NUMBER_HOUR']).T)

df_to_plot = pd.concat(dataframe_lista, ignore_index=True)

#--------------------------

campanha_dic_str = {
    'C01':['12/11/15','19/12/15'],'C02':['08/01/16','02/02/16'],'C03':['02/02/16','18/04/16'],'C05':['08/07/16','16/08/16'],
    'C06':['17/08/16','15/09/16'],'C07':['16/09/16','15/10/16'],'C08':['21/10/16','20/11/16'],'C09':['20/11/16','13/01/17'],
    'C10':['14/01/17','16/02/17'],'C11':['17/02/17','24/03/17'],'C12':['25/03/17','29/04/17'],'C13':['30/04/17','02/05/17'],
    'C15':['03/06/17','10/07/17'],'C16':['12/07/17','16/08/17'],'SG612a':['17/08/17','04/10/17'],'SG612b':['04/10/17','28/10/17'],
    'SG671':['29/10/17','09/12/17'],'20':['10/12/17','20/01/18'],'21':['21/01/18','22/02/18'],'22':['22/02/18','05/04/18'],
    '23':['10/04/18','17/05/18'],'24':['22/05/18','06/07/18'],'25':['06/07/18','22/07/18'],
    '01':['23/07/18','28/08/18'],'02':['27/08/18','03/10/18'],'03':['02/10/18','14/11/18'],'04':['13/11/18','23/12/18'],'06':['10/02/19','14/03/19'],
    '07':['23/03/19','01/05/19'],'08':['02/05/19','16/06/19'],'09':['29/06/19','31/07/19'],'10':['02/08/19','28/08/19'],'11':['04/09/19','10/10/19'],
                    }

campanha_dic_dates = []
for i in campanha_dic_str.values():
    campanha_dic_dates.append(UTCDateTime(year=int('20'+i[0].split('/')[2]), month=int(i[0].split('/')[1]), day=int(i[0].split('/')[0])).datetime.date())


campanha_dic_labels = []
for i in campanha_dic_str.keys():
    campanha_dic_labels.append(i)


# ==========================================================
# Calculating datetime between INITIAL_DATE and  FINAL_DATE
# ==========================================================

datatime_initial = datetime.strptime('2015-06-01', "%Y-%m-%d").date() 

datatime_final = datetime.strptime('2019-12-31', "%Y-%m-%d").date() 

datetime_lista = np.arange(datatime_initial, datatime_final, timedelta(days=1)).astype(datetime)

xlim_initial = mdates.date2num(datatime_initial)
xlim_final = mdates.date2num(datatime_final)

# ==========================
# Plotting DATA availability
# ==========================
#x axis parameters

months1 = MonthLocator(interval=1)  # every 1 month
months = MonthLocator(interval=6)  # every 6 month
monthsFmt = DateFormatter('%b-%y')
months1Fmt = DateFormatter('%b')

#Matplotlib parameters
fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(20,5))

data_x_axis = check_datetime_in_period(datetime_lista,df_to_plot)
datetime_pcolormesh = np.arange(datatime_initial, datatime_final+timedelta(days=1), timedelta(days=1)).astype(datetime)

im = ax.pcolormesh(datetime_pcolormesh,np.arange(25),data_x_axis,cmap='nipy_spectral_r', vmin=0, vmax=60,shading='flat',ec='none')

for dc,date_camp in enumerate(campanha_dic_dates):
    ax.axvline(x = date_camp, color='k', ls=':',lw=2)
    ax.text(date_camp, 24.5, campanha_dic_labels[dc], horizontalalignment='center',rotation=45) 

ax.set_xlim(datatime_initial,datatime_final)
ax.yaxis.set_major_locator(MultipleLocator(4))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(monthsFmt)
ax.xaxis.set_minor_locator(months1)
ax.xaxis.set_minor_formatter(months1Fmt)
ax.tick_params(which='minor', length=2)
ax.tick_params(which='major', length=20)
ax.set_ylim(0,24)
ax.set_aspect(20.0)
ax.set_ylabel('Hora do Dia',fontsize=15)
ax.grid(visible=True, which='major', color='k', linestyle='-')
ax.grid(visible=True, which='minor', color='k', linestyle='-')

plt.setp(ax.xaxis.get_majorticklabels(), fontsize=20)
plt.setp(ax.xaxis.get_minorticklabels(), fontsize=10,rotation=30)

#criando a localização da barra de cores:
axins = inset_axes(ax,
                    width="7%",  # width = 15% of parent_bbox width
                    height="2.5%",  # height : 2.5%
                    loc='upper left',
                    bbox_to_anchor=(0.0, 0.05, 1, 1),
                    bbox_transform=ax.transAxes,
                    borderpad=0,
                    )
cbar = fig.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top',ticks=[0,30,60],label='Minutos/hora')

os.makedirs(FOLDER_OUTPUT+'/FIGURAS/',exist_ok=True)
fig.savefig(FOLDER_OUTPUT+'/FIGURAS/'+'COMPLETENESS_'+datatime_initial.strftime("%Y_%m_%d")+datatime_final.strftime("%Y_%m_%d")+'_log.png',dpi=300)

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')