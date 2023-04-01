#!/usr/bin/env python
# coding: utf-8

import os
import glob
from datetime import datetime,timedelta,date
from obspy import read,UTCDateTime,Trace
import numpy as np
from scipy.io import wavfile
import functools

import pandas as pd
from multiprocessing import Pool, RLock, freeze_support
from tqdm.auto import tqdm
import time

from paramiko import SSHClient,AutoAddPolicy,Transport,SFTPClient

# ======
# Config
# ======

host = '....'

port = 22

username = '....'

password= '....'

wav_file_txt = '/home/diogoloc/dados_posdoc/Gliders_DATA/file_wav.txt'

main_directory_hds = '/mnt/medata02/HDs_01.04.22/'

FOLDER_OUTPUT = '/home/diogoloc/dados_posdoc/Gliders_DATA/OUTPUT/'

# ========
# Function
# ========

def download_downsampling_function(file_wav, sampling_rate=100,host=host,port=port,username=username,password=password):
        """
        Download wav file via SSH and decimate the trace to achieve the desired sampling rate, sr.

        NOTE: data will be detrended and a cosine taper applied before
        decimation, in order to avoid edge effects when applying the lowpass
        filter before decimating.

        source: https://quakemigrate.readthedocs.io/en/latest/_modules/quakemigrate/util.html#decimate

        Parameters:
        -----------
        file_wav : .WAV file path
            Stream to be decimated.
        sampling_rate : int
            Output sampling rate.

        Returns:
        --------
        trace : `obspy.Trace.stats` object
            Decimated trace.

        """
    #try:

        #----------------------------
        #Collecting wav data
        subdir = '/'.join(file_wav.split('/')[4:-1])+'/'
        filename = file_wav.split('/')[-1].split("'")[0]

        # Connecting via SSH
        client = SSHClient()
        client.set_missing_host_key_policy(AutoAddPolicy())
        client.connect(hostname=host,port=port,username=username,password=password)
    
        # read the file using SFTP
        sftp = client.open_sftp()
        sftp.chdir(path=main_directory_hds)
        sftp.chdir(path=subdir)
        sftp.get(remotepath=filename,localpath='/home/diogoloc/dados_posdoc/Gliders_DATA/DATA_DOWNLOAD_GLIDER/tmp/'+filename)
           
        # close the connections
        sftp.close()
        client.close()

        # Retrieving header informations
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

        sampleratetr, datatr = wavfile.read('/home/diogoloc/dados_posdoc/Gliders_DATA/DATA_DOWNLOAD_GLIDER/tmp/'+filename)

        tr = Trace(data=datatr)
        tr.stats.sampling_rate = sampleratetr
        tr.stats.starttime = d 

        starttime = d.datetime
        endtime = (d+(tr.stats.npts/tr.stats.sampling_rate)).datetime

        # Work on a copy of the trace
        trace = tr.copy()

        # Detrend and apply cosine taper
        trace.detrend('linear')
        trace.detrend('demean')
        trace.taper(type='cosine', max_percentage=0.05)

        # Zero-phase Butterworth-lowpass filter at Nyquist frequency
        trace.filter("lowpass", freq=float(sampling_rate) / 2.000001, corners=2,zerophase=True)
        trace.decimate(factor=int(trace.stats.sampling_rate / sampling_rate), strict_length=False, no_filter=True)

        # Delete .wav file:
        os.remove('/home/diogoloc/dados_posdoc/Gliders_DATA/DATA_DOWNLOAD_GLIDER/tmp/'+filename)

        OUTPUT_TRACE = FOLDER_OUTPUT+'/MSEED/'+d.strftime("%Y")+'/'+d.strftime("%Y-%m-%d")+'/'
        os.makedirs(OUTPUT_TRACE,exist_ok=True)
        trace.write(OUTPUT_TRACE+filename.split('.')[0]+'.mseed', format='MSEED')

        return trace.stats.sampling_rate
       
    #except:
        #pass

#----------------------------

# =======
# Program
# =======

start_time = time.time()

wav_files = sorted(np.genfromtxt(wav_file_txt,delimiter=',',dtype='str'))

wav_files_filtered = []
for i in wav_files:
    if not any(ext in i for ext in ['dives_operacional','teste','defeito','BIN','AMAR','snipets']):
        wav_files_filtered.append(i)

wav_files_lst = sorted(wav_files_filtered)
files_datetime_input = sorted(wav_files_lst)

# =========================
# Converting wav into mseed
# =========================
sr_lst = []
# create and configure the process pool
with Pool(processes=6) as pool:
    # execute tasks
    for result in tqdm(pool.imap_unordered(download_downsampling_function, files_datetime_input),total=len(files_datetime_input), desc='Converting WAV --> MSEED'):
        sr_lst.append(result)
# process pool is closed automatically


print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print(str(len(files_datetime_input))+' waveforms processed!')
print('\n')