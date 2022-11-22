#!/usr/bin/env python
"""
This file makes transect map shown in SI

NOTE REQUIRES python 2.7 TO PLOT CORRECTLY. RESHAPE FEATURES ALTERED IN LATER VERSIONS

"""

import numpy as np
import matplotlib
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.lines as mlines
import cartopy.crs as ccrs
from math import *
import os
from netCDF4 import Dataset as nd
import iris
import iris.plot as iplt
import pylab
import pylab as pl
import cartopy.feature as cfeat
from collections import Counter
import numpy.ma as ma
from matplotlib.patches import Polygon
from matplotlib.backends.backend_pdf import PdfPages
import scipy.stats as ss
from scipy.stats import pearsonr

########################
## Read lat/lon from example file
########################

dir_base= '/gws/nopw/j04/acure/phase3/netcdf_output/'
file_name= 'Cloud_drop_number_concentration/monthly/jan/ACURE_P3_CDNC_jan_PD_220.nc'
dataframe= nd(dir_base+file_name,'r','NETCDF4')
lat_list=dataframe.variables['latitude'][:]
lon_list=dataframe.variables['longitude'][:]
nlat=lat_list.shape[0]
nlon=lon_list.shape[0]


######################################
## Read in MODIS f_c for map under-shading
## Jul, Nov
######################################

filename='MODIS_COSP_degraded_Cloud_Retrieval_Fraction_Liquid_jul2017_using_mask.nc'
dataframe= nd('/gws/nopw/j04/acure/lregayre/MODIS_data/'+filename,'r',format='NETCDF4')
CF_jul= dataframe.variables['Cloud_Retrieval_Fraction_Liquid'][:]

filename='MODIS_COSP_degraded_Cloud_Retrieval_Fraction_Liquid_nov2017_using_mask.nc'
dataframe= nd('/gws/nopw/j04/acure/lregayre/MODIS_data/'+filename,'r',format='NETCDF4')
CF_nov= dataframe.variables['Cloud_Retrieval_Fraction_Liquid'][:]

month_list=['dec2016','jan2017','feb2017','mar2017','apr2017','may2017','jun2017','jul2017','aug2017','sep2017','oct2017','nov2017']
dir_base= '/gws/nopw/j04/acure/lregayre/MODIS_data/'

def read_array(month):
    filename='MODIS_COSP_degraded_Cloud_Retrieval_Fraction_Liquid_'+month+'_using_mask.nc'
    dataframe= nd(dir_base+filename,'r',format='NETCDF4')
    var_array_local= dataframe.variables['Cloud_Retrieval_Fraction_Liquid'][:]
    return(var_array_local)

all_month_array= np.zeros((13,nlat,nlon))

for imonth in np.arange(12):
    month = month_list[imonth]
    print(month)
    all_month_array[imonth]=read_array(month)



for ilat in np.arange(nlat):
    for ilon in np.arange(nlon):
        month_vals=all_month_array[0:12,ilat,ilon]
        non_zero_months=month_vals[np.where(month_vals>=0.0)]
        all_month_array[12,ilat,ilon]= np.average(non_zero_months)


#####################
## Reshape for plotting
#####################

def reshape_array(local_array):
    array_return=np.zeros((local_array.shape))
    ilon_middle=np.int(local_array.shape[1]/2)
    array_return[:,0:ilon_middle]=local_array[:,ilon_middle:local_array.shape[1]]
    array_return[:,ilon_middle:local_array.shape[1]]=local_array[:,0:ilon_middle]
    return(array_return)


cmap_pos= plt.get_cmap('viridis') ## any -999 will get 0 color
levels_pos= np.linspace(0,1,11)
norm_pos=BoundaryNorm(levels_pos, ncolors=cmap_pos.N, clip=True)


#################################################################################################
## South Pacific 
#################################################################################################

lat_list[55]
lon_list[150]

lat_list[59]
lon_list[143]

SA_trajectory_lon_index=[150,149,148,147,146,145,144,143] 
SA_trajectory_lat_index=[55,56,57,57,58,58,59,59]



########################################################
## South Atlantic 
#########################################################

lat_list[62]
lon_list[190]

lat_list[62]
lon_list[184]

SAf_trajectory_lon_index=[190,189,188,187,186,185] 
SAf_trajectory_lat_index=[62,62,62,62,62,62]


#############################################
## North Pacific
#############################################

lat_list[96]
lon_list[123]

lat_list[87]
lon_list[121]

NA_trajectory_lon_index=[122,122,122,122,122,121,121,121,121,121]
NA_trajectory_lat_index=[96,95,94,93,92,91,90,89,88,87]


#################################################
## North Atlantic
#################################################


lat_list[116]
lon_list[179]

lat_list[108]
lon_list[176]

EU_trajectory_lon_index=[179,179,178,178,177,177,176,176]
EU_trajectory_lat_index=[115,114,113,112,111,110,109,108]



##########################################
## binary arry for transect locations
##########################################


binary_array=np.zeros((nlat,nlon))
for ilat in np.arange(len(EU_trajectory_lat_index)):
    binary_array[EU_trajectory_lat_index[ilat],EU_trajectory_lon_index[ilat]]=1.0

for ilat in np.arange(len(SA_trajectory_lat_index)):
    binary_array[SA_trajectory_lat_index[ilat],SA_trajectory_lon_index[ilat]]=1.0

for ilat in np.arange(len(SAf_trajectory_lat_index)):
    binary_array[SAf_trajectory_lat_index[ilat],SAf_trajectory_lon_index[ilat]]=1.0

for ilat in np.arange(len(NA_trajectory_lat_index)):
    binary_array[NA_trajectory_lat_index[ilat],NA_trajectory_lon_index[ilat]]=1.0




##############################
## Transform for plotting
## because default maps have longitude
## discrepancy with model output
##############################

binary_reframed= np.zeros((nlat,nlon))
binary_reframed[:,0:int(nlon/2)] = binary_array[:,int(nlon/2):nlon]
binary_reframed[:,int(nlon/2):nlon] = binary_array[:,0:int(nlon/2)]


###################################
## Combine jul and nov transects
###################################

binary_array_jul=np.zeros((nlat,nlon))
binary_array_nov=np.zeros((nlat,nlon))

for ilat in np.arange(len(EU_trajectory_lat_index)):
    binary_array_jul[EU_trajectory_lat_index[ilat],EU_trajectory_lon_index[ilat]]=1.0

for ilat in np.arange(len(SA_trajectory_lat_index)):
    binary_array_nov[SA_trajectory_lat_index[ilat],SA_trajectory_lon_index[ilat]]=1.0

for ilat in np.arange(len(SAf_trajectory_lat_index)):
    binary_array_nov[SAf_trajectory_lat_index[ilat],SAf_trajectory_lon_index[ilat]]=1.0

for ilat in np.arange(len(NA_trajectory_lat_index)):
    binary_array_jul[NA_trajectory_lat_index[ilat],NA_trajectory_lon_index[ilat]]=1.0


binary_reframed_jul= np.zeros((nlat,nlon))
binary_reframed_jul[:,0:int(nlon/2)] = binary_array_jul[:,int(nlon/2):nlon]
binary_reframed_jul[:,int(nlon/2):nlon] = binary_array_jul[:,0:int(nlon/2)]
binary_reframed_nov= np.zeros((nlat,nlon))
binary_reframed_nov[:,0:int(nlon/2)] = binary_array_nov[:,int(nlon/2):nlon]
binary_reframed_nov[:,int(nlon/2):nlon] = binary_array_nov[:,0:int(nlon/2)]



cmap_binary = mp.colors.ListedColormap(['w', 'r'])
levels_binary= ((0.0,0.5,1.0))
norm_binary= BoundaryNorm(levels_binary, ncolors=cmap_binary.N, clip=True)


### Plot
fig= plt.figure()
ax_a=plt.subplot(1,2,1,projection=ccrs.Mollweide())
CS_a= ax_a.imshow(reshape_array(CF_jul), cmap=cmap_pos, norm=norm_pos, transform=ccrs.PlateCarree())
CS_a1 = ax_a.imshow(binary_reframed_jul[:,:], cmap=cmap_binary, norm=norm_binary, transform=ccrs.PlateCarree(),alpha=0.6)
CB_a = plt.colorbar(CS_a, shrink= 0.8,orientation= 'horizontal', pad=0.05,use_gridspec=True)
ax_a.coastlines()
ax_a.annotate('a)', xy=(0.04,0.68), xycoords='figure fraction',horizontalalignment='center', verticalalignment='center',fontsize=12)
ax_b= plt.subplot(1,2,2,projection=ccrs.Mollweide())
CS_b= ax_b.imshow(reshape_array(CF_nov), cmap=cmap_pos, norm=norm_pos, transform=ccrs.PlateCarree())
CS_b1 = ax_b.imshow(binary_reframed_nov[:,:], cmap=cmap_binary, norm=norm_binary, transform=ccrs.PlateCarree(),alpha=0.6)
CB_b = plt.colorbar(CS_b, shrink= 0.8, orientation= 'horizontal', pad=0.05,use_gridspec=True)
ax_b.coastlines()
ax_b.annotate('b)', xy=(0.54,0.68), xycoords='figure fraction',horizontalalignment='center', verticalalignment='center',fontsize=12)
plt.tight_layout(pad=-0.3, w_pad=1.0, h_pad=0.1)
plt.savefig('/gws/nopw/j04/acure/lregayre/maps/transects/transects_w_Jul_Nov_MODIS_LCF_a_b.pdf',dpi=300,orientation='portrait',bbox_inches='tight',pad_inches=0)
plt.show()


