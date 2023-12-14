#!/usr/bin/env python
# coding: utf-8

# # Importando módulos 

# In[2]:


import obspy
from obspy.taup import TauPyModel

from multiprocessing import Pool
from obspy import read,UTCDateTime,Trace
from obspy.clients.fdsn import Client
import os
import glob
import numpy as np
from collections import defaultdict
import pandas as pd

#para plotar as figuras
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition,inset_axes
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from datetime import datetime,timedelta,date
from tqdm import tqdm

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

MSEED_INPUT = "/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/OUTPUT/MSEED/"

METADATA_OUTPUT = "/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/OUTPUT/METADATA/"

filename_csv = '/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/gliders_data/info_csv/metadados_glider_acustico_pmpas-bs.csv'

campanha_csv = '/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/Equipamentos/informacoes_sensores_gliders/Dados_campanhas_Gliders.xls'

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# =====================
# Importando Dataframe
# =====================

dataframe_campanha_mseed_final = pd.read_feather('/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/OUTPUT/NETTAB/df_campanha_glider.feather')

dataframe_campanha_mseed_final_rows = [row for idx,row in dataframe_campanha_mseed_final.iterrows()]

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def rename_files(row):
    
    files_mseed_HD = row['filename_mseed_HD']

    
    files_path_mseed_HD = [glob.glob(MSEED_INPUT+'*/*/'+i+'.mseed') for i in files_mseed_HD]
    
    for y in files_path_mseed_HD:
        
        st = obspy.read(y[0])
        st[0].stats.network = 'GL'
        st[0].stats.station = row['station_nettab']
        st[0].stats.channel = 'HHH'

        #<SDSdir>/Year/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY
        #Saving in SDS structure        
        OUTPUT_STREAM = FOLDER_OUTPUT+'MSEED_GL_SDS/'+st[0].stats.starttime.strftime("%Y")+'/'+st[0].stats.network+'/'+st[0].stats.station+'/HHH.D/'
        os.makedirs(OUTPUT_STREAM,exist_ok=True)
        
        station_name_file = st[0].stats.network+'.'+st[0].stats.station+'..'+st[0].stats.channel+'.D.'+st[0].stats.starttime.strftime("%Y.%j.%H.%M.%S")+'.mseed'
        if os.path.exists(OUTPUT_STREAM+station_name_file) == False:
            st.write(OUTPUT_STREAM+station_name_file, format='MSEED')

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ======================
# Renomeando os arquivos
# ======================

with Pool(processes=4) as p:
    with tqdm(total=len(dataframe_campanha_mseed_final)) as pbar:
        for result in p.imap_unordered(rename_files,dataframe_campanha_mseed_final_rows):
            pbar.update()
