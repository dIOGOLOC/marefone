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
from multiprocessing import Pool
from tqdm.auto import tqdm
import time

# ======
# Config
# ======

FOLDER_INPUT = '/run/media/dIOGOLOC/FA03_ABR_2021/'

FOLDER_OUTPUT = '/run/media/dIOGOLOC/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/OUTPUT/'

NB_process = 4

# ========
# Function
# ========

def downsampling_function(file_wav, sampling_rate=100):

    try: 

        """
        Download wav file via SSH and decimate the trace to achieve the desired sampling rate, sr.

        NOTE:data will be detrended and a cosine taper applied before
            decimation, in order to avoid edge effects when applying the lowpass
            filter before decimating.

        SOURCE: https://quakemigrate.readthedocs.io/en/latest/_modules/quakemigrate/util.html#decimate

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
        #----------------------------
        #Collecting wav data (AMAR038.1.20191229T000001Z.wav)

        filename = file_wav.split('/')[-1].split('.wav')[0]

        # Retrieving header informations
        sensor_name = filename.split('.')[0]

        d = UTCDateTime(filename.split('.')[-1])

        #Check if the file exists:
        OUTPUT_TRACE = FOLDER_OUTPUT+'/MSEED_LFI/'+d.strftime("%Y")+'/'+d.strftime("%Y-%m-%d")+'/'

        if os.path.exists(OUTPUT_TRACE+filename.split('.')[0]+'.mseed') == False:
            
            os.makedirs(OUTPUT_TRACE,exist_ok=True)
    
            sampleratetr, datatr = wavfile.read(file_wav)

            tr = Trace(data=datatr)
            tr.stats.sampling_rate = sampleratetr
            tr.stats.starttime = d 

            starttime = d.datetime
            endtime = (d+(tr.stats.npts/tr.stats.sampling_rate)).datetime

            # Work on a copy of the trace
            trace = tr.copy()

            # Zero-phase Butterworth-lowpass filter at Nyquist frequency
            trace.filter("lowpass", freq=float(sampling_rate) / 2.000001, corners=2,zerophase=True)
            trace.decimate(factor=int(trace.stats.sampling_rate / sampling_rate), strict_length=False, no_filter=True)

            #Saving MSEED file
            trace.write(OUTPUT_TRACE+filename+'.mseed', format='MSEED')
    except:
        print('Problem: '+file_wav)

#----------------------------

# =======
# Program
# =======

start_time = time.time()

wav_files = []

for root, dirs, files in os.walk(FOLDER_INPUT):
    for file in files:
        if file.endswith('.wav'):
            wav_files.append(os.path.join(root, file))

wav_files_filtered = []
for i in wav_files:
    if not any(ext in i for ext in ['teste','defeito','BIN','snipets','sem_ruido',]):
        wav_files_filtered.append(i)

files_datetime_input = sorted(wav_files_filtered)

# =========================
# Converting wav into mseed
# =========================
sr_lst = []
# create and configure the process pool

with Pool(processes=NB_process) as p:
    max_ = len(files_datetime_input)
    with tqdm(total=max_, desc='Converting WAV --> MSEED') as pbar:
        for _ in p.imap_unordered(downsampling_function, files_datetime_input):
            pbar.update()
# process pool is closed automatically

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print(str(len(files_datetime_input))+' waveforms processed!')
print('\n')