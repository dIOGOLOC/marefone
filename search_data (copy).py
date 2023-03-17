
'''

--------------------------------------------------------------------------------
                 Plotting the COMPLETENESS of the raw data
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 04/2022


Project: ---


Description:
This code will retrieve information about the filename of the raw data for each
file and plot a mosaic the the Data availability.


Inputs:
dive .wav file

Outputs:
Images of the COMPLETENESS of files per hour (format: PDF)

'''

import sys
import os
import glob
from datetime import datetime
import pandas as pd

# ======
# Inputs
# ======

HD_NUMER = 'HD_11'

INPUT_HD_FOLDER = '/mnt/medata02/HDs_01.04.22/'+HD_NUMER+'/'

# ========
# Function
# ========

def find_files(folder_hd):
    files_wav = sorted(glob.glob(folder_hd+'**/**/**/**/**/*.wav'))

    return files_wav

# =======
# Program
# =======

#command = str(find_files(INPUT_HD_FOLDER))

wav_files = find_files(INPUT_HD_FOLDER)

files_datetime = []

for i in wav_files:
    # splitting subdir/basename
    subdir, filename = os.path.split(i)
    year_month_day = filename.split('.wav')[0].split('_')[2]
    hour_minute_second = filename.split('.wav')[0].split('_')[3]

    year = int('20'+year_month_day[:2])
    month = int(year_month_day[2:4])
    day = int(year_month_day[4:])

    hour = int(hour_minute_second[:2])
    minute = int(hour_minute_second[2:4])
    second = int(hour_minute_second[4:])

    d = datetime(year,month,day,hour,minute,second).isoformat()
    files_datetime.append(d)

files_datetime = sorted(files_datetime)
for dt in files_datetime:
    print(dt)


# ajeitar para amanhã

 ==============================
# Function to store data in dic
# ==============================

def plot_date_file(FIG_FOLDER_OUTPUT,directory_data,XML_FILE):

    data_lista = []
    if LABEL_LANG == 'br':
        print('Procurando por dados no diretório: '+directory_data)

    else:
        print('Looking for data in the directory: '+directory_data)

    for root, dirs, files in os.walk(directory_data):
        for name in files:
            data_lista.append(os.path.join(root, name))

    data_lista = sorted(data_lista)

    dataframe_lista = []
    #create a empty dataframe with pandas
    for i,j in enumerate(data_lista):
        if LABEL_LANG == 'br':
            print("Extraindo os data do cabeçalho: "+str(i+1)+" of "+str(len(data_lista)))

        else:
            print("Extracting data from header: "+str(i+1)+" of "+str(len(data_lista)))

        #Reading header from data
        st = obspy.read(j)

        #----------------------------
        #Dataframe starting

        network = st[0].stats.network
        station = st[0].stats.station
        channel = st[0].stats.channel

        time_lst = []

        for t,trace in enumerate(st):
            starttime = trace.stats.starttime
            endtime = trace.stats.endtime
            time_lst.append(np.arange(starttime,endtime,60))

        flat_time_lst = [item for sublist in time_lst for item in sublist]

        DATETIME = str(st[0].stats.starttime.year)+','+str(st[0].stats.starttime.month)+','+str(st[0].stats.starttime.day)

        #Contador da lista de horas
        time_flat_time_lst = [[]]*24
        for g,h in enumerate(np.arange(24)):
            lst_time = []
            for x,c in enumerate(flat_time_lst):
                if c.hour == h:
                    lst_time.append(c.hour)
            time_flat_time_lst[g] = lst_time

        NUMBER_HOUR = [[]]*24
        for q,w in enumerate(time_flat_time_lst):
            NUMBER_HOUR[q] = len(w)


        dataframe_lista.append(pd.DataFrame([[network],[station],[channel],[DATETIME],[NUMBER_HOUR]], index=['NETWORK', 'STATION', 'CHANNEL', 'DATETIME','NUMBER_HOUR']).T)
        #Dataframe ending
        #----------------------------

    df = pd.concat(dataframe_lista, ignore_index=True)

    #Sorting according to station

    station_lista = list(set(df['STATION']))

    for i,j in enumerate(station_lista):
        df_sta = df[df['STATION'] == j]

        channel_lista = list(set(df_sta['CHANNEL']))
        channel_lista = sorted(channel_lista)

        # ==========================================================
        # Calculating datetime between INITIAL_DATE and  FINAL_DATE
        # ==========================================================

        datatime_initial = datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day)

        datatime_final = datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day)

        datetime_lista = np.arange(datatime_initial, datatime_final, datetime.timedelta(days=1)).astype(datetime.datetime)


        xlim_initial = mdates.date2num(datatime_initial)
        xlim_final = mdates.date2num(datatime_final)

        #----------------------------
        #Function to check if the dates in data set are inside the period chosen (INITIAL_DATE to FINAL_DATE)

        def check_datetime_in_period(datetime_lst,df_DATETIME,df_NUMBER_HOUR):

            array_to_plot_by_xlim = []
            for x,c in enumerate(datetime_lst):
                lista_temp = []
                for t,y in enumerate(df_DATETIME):
                        if datetime.datetime(obspy.UTCDateTime(y).year,obspy.UTCDateTime(y).month,obspy.UTCDateTime(y).day) == c:
                                lista_temp.append(df_NUMBER_HOUR[df_DATETIME[df_DATETIME == y].index[0]])
                array_to_plot_by_xlim.append(lista_temp)

            data_x_axis = []
            for x,c in enumerate(array_to_plot_by_xlim):
                if c != []:
                    data_x_axis.append(c[0][::-1])
                else:
                    data_x_axis.append(np.zeros_like(np.arange(24)))

            data_x_axis = np.array(data_x_axis).T

            return data_x_axis

        # ====================================
        # Function to plot DATA availability
        # ====================================

        #x axis parameters

        days1 = DayLocator(interval=1)   # every day
        days5 = DayLocator(interval=int(len(datetime_lista)*5/100))   # every day
        months = MonthLocator()  # every month
        yearsFmt = DateFormatter('%Y-%m-%d')

        days1.MAXTICKS = 10000


        #Matplotlib parameters
        fig, ax = plt.subplots(nrows=len(channel_lista), ncols=1,sharex=True,sharey=True,figsize=(40,15))
        fig.suptitle(j,fontsize=25,y=0.9)
        for k,l in enumerate(channel_lista):

            df_ch = df_sta[df_sta['CHANNEL'] == l]

            data_x_axis = check_datetime_in_period(datetime_lista,df_ch['DATETIME'],df_ch['NUMBER_HOUR'])

            im = ax[k].imshow(data_x_axis,extent = [xlim_initial,xlim_final,0,24],cmap=plt.cm.Greens,interpolation=None, vmin=0, vmax=60)
            ax[k].set_xlim(datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day),datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day))
            ax[k].yaxis.set_major_locator(MultipleLocator(4))
            ax[k].yaxis.set_minor_locator(MultipleLocator(1))
            ax[k].xaxis.set_major_locator(days5)
            ax[k].xaxis.set_major_formatter(yearsFmt)
            ax[k].xaxis.set_minor_locator(days1)
            ax[k].tick_params(which='minor', length=4)
            ax[k].tick_params(which='major', length=10)
            ax[k].set_ylim(0,24)
            ax[k].set_ylabel(l,fontsize=15)
            ax[k].grid(b=True, which='major', color='k', linestyle='-')
            ax[k].grid(b=True, which='minor', color='k', linestyle='-')


        plt.setp(ax[k].xaxis.get_majorticklabels(), fontsize=10, rotation=30)
        if LABEL_LANG == 'br':
            ax[-1].set_xlabel('Data', fontsize=20)

        else:
            ax[-1].set_xlabel('Time', fontsize=20)

        #criando a localização da barra de cores:
        axins = inset_axes(ax[0],
                           width="10%",  # width = 10% of parent_bbox width
                           height="5%",  # height : 50%
                           loc='upper left',
                           bbox_to_anchor=(0.85, 0.1, 1, 1),
                           bbox_transform=ax[0].transAxes,
                           borderpad=0,
                           )
        cbar = fig.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top',ticks=[0,30,60],label='Files per hour')
        cbar.ax.set_xticklabels(['0%','50%','100%'])

        os.makedirs(FIG_FOLDER_OUTPUT,exist_ok=True)
        fig.savefig(FIG_FOLDER_OUTPUT+j+'_'+'COMPLETENESS_'+str(obspy.UTCDateTime(INITIAL_DATE).year)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).month)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).day)+'_'+str(obspy.UTCDateTime(FINAL_DATE).year)+'_'+str(obspy.UTCDateTime(FINAL_DATE).month)+'_'+str(obspy.UTCDateTime(FINAL_DATE).day)+'.pdf',dpi=500)
        #plt.show()
