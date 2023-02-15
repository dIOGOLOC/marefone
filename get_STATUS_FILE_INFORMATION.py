#!/usr/bin/python -u

'''
--------------------------------------------------------------------------------
      Getting and Plotting information from status files
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 07/2022


Description:
This code will retrieve and plot information from status files recorded by a
SeaGlider.


More information in:
https://glidertools.readthedocs.io/en/latest/loading.html

'''
import os
import glob
import numpy as np
import time
from tqdm import tqdm
from multiprocessing import Pool
from obspy import read,read_inventory, UTCDateTime, Stream
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, DateFormatter
import matplotlib.dates as mdates
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,FixedLocator, FixedFormatter
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import glidertools as gt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter,LatitudeLocator,LongitudeLocator

import pyarrow.feather as feather
import xarray
import pandas as pd
from assessment_py.status_assessment import filelist,READ_PLOT_DATE_file

from parameters_py.config import (
					OUTPUT_FIGURE_DIR,DIR_STATUS_FILES,NUM_PROCESS,LABEL_LANG,OUTPUT_OUT_DIR
					)


# ============================
# Retrieving GCF status files
# ============================
if LABEL_LANG == 'br':
    print('Obtendo arquivos de funcionamento do Glider e salvando em OUTPUT/FEATHER_FILES:')
    print('\n')

else:
    print('\n')
    print('Getting Glider status files and saving in OUTPUT/FEATHER_FILES:')
    print('\n')

files_NC = filelist(DIR_STATUS_FILES)

def get_glider_data(f):
    # splitting subdir/basename
    subdir, filename = os.path.split(f)
    name_file = filename.split('.nc')[0]
    '''
    From the variable listing, one can choose multiple variables to load. 
    Note that one only needs the variable name to load the data. 
    Below, we've created a list of variables that we'll be using for this demo.

    The gt.load.seaglider_basestation_netCDFs function is used to load a list of variables. 
    It requires the filename string or list (as described above) and keys. 
    It may be that these variables are not sampled at the same frequency. 
    In this case, the loading function will load the sampling frequency dimensions separately. 
    The function will try to find a time variable for each sampling frequency/dimension.

    Coordinates and automatic time fetching
    All associated coordinate variables will also be loaded with the data if coordinates are documented.
    These may included latitude, longitude, depth and time (naming may vary). 
    If time cannot be found for a dimension, a time variable from a different dimension with the same 
    number of observations is used instead. 
    This insures that data can be merged based on the time of sampling.

    Merging data based on time
    If the return_merged is set to True, the function will merge the dimensions if the dimension has an 
    associated time variable.

    The function returns a dictionary of xarray.Datasets - a Python package that deals with coordinate 
    indexed multi-dimensional arrays. We recommend that you read the documentation 
    (http://xarray.pydata.org/en/stable/) as this package is used throughout GliderTools. 
    This allows the original metadata to be copied with the data. The dictionary keys are the names of 
    the dimensions. If return_merged is set to True an additional entry under the key merged will be 
    included.

    The structure of a dimension output is shown below. Note that the merged data will use the largest 
    dimension as the primary dataset and the other data will be merged onto that time index. 
    Data is linearly interpolated to the nearest time measurement of the primary index, 
    but only by one measurement to ensure transparancy.

    Metadata handling
    If the keyword arguement keep_global_attrs=True, the attributes from the original files 
    (for all that are the same) are passed on to the output Datasets from the original netCDF attributes. 
    The variable attributes (units, comments, axis...) are passed on by default, but can also be set to False if not wanted. 
    GliderTools functions will automatically pass on these attributes to function outputs if a xarray.
    DataArray with attributes is given. 
    All functions applied to data will also be recorded under the variable attribute processing.    
    '''
       
    names = [
            'ctd_depth',
            'ctd_time',
            'latitude',
            'longitude'
            'time',
            'detph',
            'audiolat',
            'audiolon'
            ]

    try: 

        #Check if file exists
        FEATHER_FOLDER = OUTPUT_OUT_DIR+'FEATHER_FILES/'
        file_feather_name = FEATHER_FOLDER+name_file+'_nc_data.feather'
        if os.path.isfile(file_feather_name):
            pass

        else:

            ds_dict = gt.load.seaglider_basestation_netCDFs(f, names,verbose=False,keep_global_attrs=True)
            # Creating a Pandas DataFrame:
            dat = ds_dict['sg_data_point'].to_dataframe()

            # ----------------------------------------------------------------------------------------------------
            # Convert from pandas to Arrow and saving in feather formart file
            os.makedirs(FEATHER_FOLDER,exist_ok=True)
            feather.write_feather(dat, file_feather_name)
            # ----------------------------------------------------------------------------------------------------

            return 0

    except:
        pass  

# ---------------------------------------------------------------------------

with Pool(NUM_PROCESS) as p:
    list(tqdm(p.imap(func=get_glider_data, iterable=files_NC), total=len(files_NC)))
print('\n')


# ---------------------------------------------------------------------------
if LABEL_LANG == 'br':
    print('Importantando  de arquivos de funcionamento do Glider(*.feather).')
    print('\n')

else:
    print('\n')
    print('Importing of Glider state of health files (*.feather).')
    print('\n')

nc_feather_files_lst = sorted(glob.glob(OUTPUT_OUT_DIR+'FEATHER_FILES/*'))

nc_feather_files = [pd.read_feather(i) for i in nc_feather_files_lst]

if LABEL_LANG == 'br':
    print('Número de arquivos de funcionamento do Glider(*.feather):',len(nc_feather_files))
    print('\n')

else:
    print('\n')
    print('Number of Glider state of health files (*.feather):',len(nc_feather_files))
    print('\n')

if LABEL_LANG == 'br':
    print('Concatenando os dataframes(*.feather).')
    print('\n')

else:
    print('\n')
    print('Dataframe concatenating (*.feather).')
    print('\n')

df = pd.concat(nc_feather_files, ignore_index=True)
df.sort_values('ctd_time_dt64')

print(df.info)
#-------------------------------------------

if LABEL_LANG == 'br':
    print('Plotando os dados *.feather:')
    print('\n')

else:
    print('\n')
    print('Plotting *.feather data:')
    print('\n')

fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(nrows=1, ncols=1)
#-------------------------------------------

crs = ccrs.NearsidePerspective(central_longitude=-40, central_latitude=-20)
map_loc = fig.add_subplot(gs[0],projection=crs)

LLCRNRLON_LARGE = -52
URCRNRLON_LARGE = -28
LLCRNRLAT_LARGE = -30
URCRNRLAT_LARGE = -12

#map_loc.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
map_loc.yaxis.set_ticks_position('both')
map_loc.xaxis.set_ticks_position('both')


map_loc.grid(True,which='major',color='gray',linewidth=1,linestyle='--')


# Create a Stamen Terrain instance.
stamen_terrain = cimgt.Stamen('terrain-background')

# Add the Stamen data at zoom level 8.
map_loc.add_image(stamen_terrain, 10)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cfeature.NaturalEarthFeature(
												category='cultural',
												name='admin_1_states_provinces_lines',
												scale='50m',
												facecolor='none'
												)

map_loc.add_feature(cfeature.LAND,)
map_loc.add_feature(cfeature.COASTLINE)
map_loc.add_feature(states_provinces, edgecolor='k',linewidth=0.5)

gl = map_loc.gridlines(color='gray',linewidth=0.5,linestyle='--',draw_labels=True)

gl.xlocator = LongitudeLocator(4)
gl.ylocator = LatitudeLocator(4)
gl.xformatter = LongitudeFormatter()
gl.yformatter = LatitudeFormatter()

gl.xlabel_style = {'size': 20, 'color': 'gray'}
gl.xlabel_style = {'color': 'black', 'weight': 'bold'}
gl.ylabel_style = {'size': 20, 'color': 'gray'}
gl.ylabel_style = {'color': 'black', 'weight': 'bold'}

map_loc.tick_params(labelbottom=True,labeltop=True,labelleft=True,labelright=True, labelsize=15)

smap = map_loc.scatter(df['longitude'], df['latitude'], marker='o',s=5,c=mdates.date2num(df['ctd_time_dt64']),cmap='cividis',vmin=mdates.date2num(df['ctd_time_dt64']).min(),vmax=mdates.date2num(df['ctd_time_dt64']).max(), transform=ccrs.PlateCarree())

cbar = plt.colorbar(smap,cmap='viridis', orientation='vertical',ticklocation='auto',format=DateFormatter('%b %Y'))

os.makedirs(OUTPUT_FIGURE_DIR,exist_ok=True)
fig.savefig(OUTPUT_FIGURE_DIR+'GLIDER_MAP_TRAJETORY.png',dpi=300)

#--------------------------------------------------------------------------------------------------------------------

a = READ_PLOT_DATE_file(FIG_FOLDER_OUTPUT=OUTPUT_FIGURE_DIR,TXT_FILE='/home/diogoloc/dados_posdoc/gliders_project/gliders_data/file_wav.txt')


result_m8 = []
pool = Pool(processes=NUM_PROCESS)
for result in tqdm(pool.imap_unordered(func=get_status_file_GURALP, iterable=daily_lst_m8), total=len(daily_lst_m8)):
    result_m8.append(result)
#--------------------------------------------------------------------------------------------------------------------

if LABEL_LANG == 'br':
    print('Canal: m9')

else:
    print('Channel: m9')
    print('\n')

result_m9 = []
pool = Pool(processes=NUM_PROCESS)
for result in tqdm(pool.imap_unordered(func=get_status_file_GURALP, iterable=daily_lst_m9), total=len(daily_lst_m9)):
    result_m9.append(result)
#--------------------------------------------------------------------------------------------------------------------

if LABEL_LANG == 'br':
    print('Canal: ma')

else:
    print('Channel: ma')
    print('\n')

result_ma = []
pool = Pool(processes=NUM_PROCESS)
for result in tqdm(pool.imap_unordered(func=get_status_file_GURALP, iterable=daily_lst_ma), total=len(daily_lst_ma)):
    result_ma.append(result)
#--------------------------------------------------------------------------------------------------------------------
if LABEL_LANG == 'br':
    print('Canal: me')

else:
    print('Channel: me')
    print('\n')

result_me = []
pool = Pool(processes=NUM_PROCESS)
for result in tqdm(pool.imap_unordered(func=get_status_file_GURALP, iterable=daily_lst_me), total=len(daily_lst_me)):
    result_me.append(result)
#--------------------------------------------------------------------------------------------------------------------

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')

# ====================
# Plotting status data
# ====================
print('Plotting status data')

#--------------------------------------------------------------------------------------------------------------------
days = DayLocator()   # every year
months = MonthLocator()  # every month
yearsFmt = DateFormatter('%Y-%m-%d')

fig, ax = plt.subplots(figsize=(15,7))
for i in result_m8:
	try:
		ax.plot(UTCDateTime(i[2]).matplotlib_date, i[3],'ok')
	except:
		pass
ax.set_ylabel('Z mass position (microvolts)')
y_lim = ax.get_ylim()
ax.set_ylim(y_lim[0],y_lim[1])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(days)

ax.grid(True)

# rotates and right aligns the x labels, and moves the bottom of the
# axes up to make room for them
fig.autofmt_xdate()
fig.savefig(OUTPUT_FIGURE_DIR+result_me[0][0]+'_Z_mass_position.png', dpi=300, facecolor='w')

#--------------------------------------------------------------------------------------------------------------------
days = DayLocator()   # every year
months = MonthLocator()  # every month
yearsFmt = DateFormatter('%Y-%m-%d')

fig, ax = plt.subplots(figsize=(15,7))
for i in result_m9:
	try:
		ax.plot(UTCDateTime(i[2]).matplotlib_date, i[3],'ok')
	except:
		pass
ax.set_ylabel('N mass position (microvolts)')
y_lim = ax.get_ylim()
ax.set_ylim(y_lim[0],y_lim[1])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(days)

ax.grid(True)

# rotates and right aligns the x labels, and moves the bottom of the
# axes up to make room for them
fig.autofmt_xdate()
fig.savefig(OUTPUT_FIGURE_DIR+result_me[0][0]+'_N_mass_position.png', dpi=300, facecolor='w')

#--------------------------------------------------------------------------------------------------------------------
days = DayLocator()   # every year
months = MonthLocator()  # every month
yearsFmt = DateFormatter('%Y-%m-%d')

fig, ax = plt.subplots(figsize=(15,7))
for i in result_ma:
	try:
		ax.plot(UTCDateTime(i[2]).matplotlib_date, i[3],'ok')
	except:
		pass
ax.set_ylabel('E mass position (microvolts)')
y_lim = ax.get_ylim()
ax.set_ylim(y_lim[0],y_lim[1])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(days)

ax.grid(True)

# rotates and right aligns the x labels, and moves the bottom of the
# axes up to make room for them
fig.autofmt_xdate()
fig.savefig(OUTPUT_FIGURE_DIR+result_me[0][0]+'_E_mass_position.png', dpi=300, facecolor='w')

#--------------------------------------------------------------------------------------------------------------------
days = DayLocator()   # every year
months = MonthLocator()  # every month
yearsFmt = DateFormatter('%Y-%m-%d')

fig, ax = plt.subplots(figsize=(15,7))
for i in result_me:
	try:
		ax.plot(UTCDateTime(i[2]).matplotlib_date, i[3],'ok')
	except:
		pass
ax.set_ylabel('Temperature (microvolts)')
y_lim = ax.get_ylim()
ax.set_ylim(y_lim[0],y_lim[1])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(days)

ax.grid(True)

# rotates and right aligns the x labels, and moves the bottom of the
# axes up to make room for them
fig.autofmt_xdate()
fig.savefig(OUTPUT_FIGURE_DIR+result_me[0][0]+'_temperature.png', dpi=300, facecolor='w')clear
