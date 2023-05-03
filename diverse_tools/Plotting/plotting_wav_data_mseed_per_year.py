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

YEAR = 2019

# ==========================================================
# Calculating datetime between INITIAL_DATE and  FINAL_DATE
# ==========================================================

datatime_initial = datetime.strptime(str(YEAR)+'-01-01', "%Y-%m-%d").date() 

datatime_final = datetime.strptime(str(YEAR)+'-12-31', "%Y-%m-%d").date() 

datetime_lista = np.arange(datatime_initial, datatime_final, timedelta(days=1)).astype(datetime)

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
        hour_day = starttime.hour
        minute_day = starttime.minute
   
        df = pd.DataFrame([[filename],[mergulho],[stream_number],[starttime],[hour_day],[minute_day]], index=['filename', 'mergulho', 'stream_number','starttime','hour_day','minute_day']).T
        #Ending Dataframe
        #----------------------------
        
        return df

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

mseed_files_lst = sorted(glob.glob(mseed_files+'*/*/*.mseed'))

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

dataframe_final['DayMonthYear'] = dataframe_final['starttime'].dt.strftime("%Y-%m-%d")

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
campanha_dic_dates_id = []
for id_camp,i in enumerate(campanha_dic_str.values()):
    date_campanha = UTCDateTime(year=int('20'+i[0].split('/')[2]), month=int(i[0].split('/')[1]), day=int(i[0].split('/')[0])).datetime.date()
    datai = [date_campanha for c in datetime_lista if i == c]
    if len(datai) == 1:
        campanha_dic_dates.append(datai[0])
        campanha_dic_dates_id.append(id_camp)
    else:
        pass


campanha_dic_labels = []
for id_label,i in enumerate(campanha_dic_str.keys()):
    if id_label in campanha_dic_dates_id:
        campanha_dic_labels.append(i)

# ==========================
# Plotting DATA availability
# ==========================
#x axis parameters

days = DayLocator(interval=1)  # every 1 day
months = MonthLocator(interval=1)  # every 1 month
monthsFmt = DateFormatter('%b-%y')

#Matplotlib parameters
fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(16,4))

data_x_axis = check_datetime_in_period(datetime_lista,df_to_plot)
datetime_pcolormesh = np.arange(datatime_initial, datatime_final+timedelta(days=1), timedelta(days=1)).astype(datetime)

im = ax.pcolormesh(datetime_pcolormesh,np.arange(25),data_x_axis,cmap='nipy_spectral_r', vmin=0, vmax=60,shading='flat',ec='none')

for dc,date_camp in enumerate(campanha_dic_dates):
    ax.axvline(x = date_camp, color='k', ls='--',lw=2)
    ax.text(date_camp, 24.5, campanha_dic_labels[dc], horizontalalignment='center',rotation=45)
    
ax.set_xlim(datatime_initial,datatime_final)
ax.yaxis.set_major_locator(MultipleLocator(4))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(monthsFmt)
ax.xaxis.set_minor_locator(days)
#ax.xaxis.set_minor_formatter(dayFmt)
ax.tick_params(which='minor', length=2)
ax.tick_params(which='major', length=10)
ax.set_ylim(0,24)
ax.set_aspect(2)
ax.set_ylabel('Hora do Dia',fontsize=15)
ax.grid(visible=True, which='major', color='k', linestyle='-')
ax.grid(visible=True, which='minor', color='k', linestyle='-')

plt.setp(ax.xaxis.get_majorticklabels(), fontsize=12,rotation=30)

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
fig.savefig(FOLDER_OUTPUT+'/FIGURAS/'+'COMPLETENESS_'+datatime_initial.strftime("%Y")+'_mseed.png',dpi=300)

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')