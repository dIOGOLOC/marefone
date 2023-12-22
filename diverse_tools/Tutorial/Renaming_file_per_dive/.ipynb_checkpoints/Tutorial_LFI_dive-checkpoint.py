#!/usr/bin/env python
# coding: utf-8

# # Importando módulos 

import obspy

from multiprocessing import Pool
from obspy import read,UTCDateTime,Trace
from obspy.clients.fdsn import Client
import os
import glob
import numpy as np
from collections import defaultdict
import pandas as pd
from pandas import Timestamp, date_range

from tqdm import tqdm


import itertools

#para plotar as figuras
#import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition,inset_axes
import matplotlib
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from shapely.geometry.polygon import LinearRing

import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from IPython.display import HTML
from IPython import display

# ================
# Inputs e Outputs
# ================


FOLDER_OUTPUT = '/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/OUTPUT/'

MSEED_INPUT = "/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/OUTPUT/MSEED_LFI/"

EVENTs_DATA_FOLDER = '/home/sysop/Documents/codes_marefone/EVENT_MSEED/'

SELECT_EVENTs_DATA_FOLDER = '/home/sysop/Documents/codes_marefone/EVENT_MSEED_SEL/'

METADATA_OUTPUT = "/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/gliders_data/info_csv/metadados_glider_acustico_pmpas-bs.csv"

campanha_LFI_csv = '/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/Equipamentos/informacoes_sensores_lfi/Dados das aquisicoes_LFIs.xls'

campanha_LFI_posicoes_csv = '/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/Equipamentos/informacoes_sensores_lfi/LFIs_Posicoes_correto.xlsx'

# ----------------------------------------------------------------

# =======================================
# Lendo Dataframes (LFI e Lista ".mseed")
# =======================================

print('Lendo Dataframes (LFI e Lista ".mseed")')
df_LFI = pd.read_feather('/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/OUTPUT/NETTAB/df_LFI_nettab.feather')


dataframe_files = pd.read_feather('/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/OUTPUT/NETTAB/df_LFI_dataframe_files.feather')

# ================================
# Renomeando os arquivos ".mseed":
# ================================

dataframe_files_rows = [row for idx,row in dataframe_files.iterrows()]

# --------------------------------------------------------------------------------

def rename_files_LFI(row):
    try:
        df_LFI_info = df_LFI[(df_LFI['starttime'] < row['starttime']) & (df_LFI['endtime'] > row['starttime']) & (df_LFI['registrador'] == row['sensor'])]

        st = obspy.read(row['file'])
        st[0].stats.network = 'LF'
        st[0].stats.station = str(df_LFI_info['nettab_name'].values[0])
        st[0].stats.channel = 'HHH'
        
        #<SDSdir>/Year/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY
        #Saving in SDS structure        
        OUTPUT_STREAM = FOLDER_OUTPUT+'MSEED_LFI_SDS/'+st[0].stats.starttime.strftime("%Y")+'/'+st[0].stats.network+'/'+st[0].stats.station+'/HHH.D/'
        os.makedirs(OUTPUT_STREAM,exist_ok=True)
            
        station_name_file = st[0].stats.network+'.'+st[0].stats.station+'..'+st[0].stats.channel+'.D.'+st[0].stats.starttime.strftime("%Y.%j.%H.%M.%S")+'.mseed'
        if os.path.exists(OUTPUT_STREAM+station_name_file) == False:
            st.write(OUTPUT_STREAM+station_name_file, format='MSEED')    
    except:
        print(row)
# --------------------------------------------------------------------------------

with Pool(processes=4) as p:
    with tqdm(total=len(dataframe_files_rows)) as pbar:
        for result in p.imap_unordered(rename_files_LFI,dataframe_files_rows):
            pbar.update()

