#!/usr/bin/env python
"""
Seasonal cycle of hemispheric difference in marine Nd

"""

import numpy as np
import matplotlib
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
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

iregion=0
region_list=['NH_SH_difference']
region_titles=['Hemispheric difference']
region_file_titles=['Hemispheric_difference']
month_list=[ 'dec','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','ann']
nmonths=len(month_list)

region= region_list[iregion]
region_title= region_titles[iregion]
region_file_title= region_file_titles[iregion]


var_list_obs= ['CDNC']
var_list_PPE=['CDNC']
var_list_def=['CDNC']
nvar= len(var_list_obs)
var_name_titles=['N$_{d}$']
units_list=['cm$^{-3}$']


dir_PPE='/gws/nopw/j04/acure/lregayre/data_for_emulation/'
dir_obs=dir_PPE+'observations/'
dir_def=dir_PPE+'default/'

PPE_array=np.zeros((nvar,nmonths,221))
obs_array=np.zeros((nvar,nmonths))
def_array=np.zeros((nvar,nmonths))

for ivar in np.arange(nvar):
    for imonth in np.arange(nmonths):
        print(imonth)
        PPE_array[ivar,imonth,:]=np.loadtxt(dir_PPE+month_list[imonth]+'/'+var_list_PPE[ivar]+'_PPE_'+region+'_'+month_list[imonth]+'.dat')



###########################################################
# Read in obs data, regridded and degraded to match model resolution
###########################################################

for ivar in (np.arange(nvar)):
    print(var_list_PPE[ivar])
    obs_array[ivar]= np.loadtxt(dir_obs+var_list_obs[ivar]+'_'+region+'_dec_to_ann_from_regridded_w_mask.dat')



#################################################
# Read in default data
## annual needs to be calculated here
#################################################

for ivar in np.arange(nvar):
    for imonth in np.arange(nmonths-1):
        month=month_list[imonth]
        def_array[ivar,imonth]= np.loadtxt(dir_def+month+'/'+var_list_def[ivar]+'_default_'+region+'_'+month+'.dat')
    def_array[ivar,12]= np.mean(def_array[ivar,0:12])



#####################################
## Calculate amplidues of seasonal cycles
#####################################

obs_amplitude=np.zeros((nvar))
def_amplitude=np.zeros((nvar))
PPE_amplitude=np.zeros((nvar,221))

for ivar in np.arange(nvar):
    obs_amplitude[ivar]=np.max(obs_array[ivar,0:12])-np.min(obs_array[ivar,0:12])
    def_amplitude[ivar]=np.max(def_array[ivar,0:12])-np.min(def_array[ivar,0:12])
    for ippe in np.arange(221):
        PPE_amplitude[ivar,ippe]=np.max(PPE_array[ivar,0:12,ippe])-np.min(PPE_array[ivar,0:12,ippe])

## write to file

for ivar in np.arange(nvar):
    outfile='/gws/nopw/j04/acure/lregayre/data_for_emulation/seasonal_amplitude/'+var_list_obs[ivar]+'_'+region+'_seasonal_amplitude_observed.dat'
    with open(outfile, 'w') as f:
        f.write(np.str(obs_amplitude[ivar]))
    outfile='/gws/nopw/j04/acure/lregayre/data_for_emulation/seasonal_amplitude/'+var_list_def[ivar]+'_'+region+'_seasonal_amplitude_default.dat'
    if (ivar!=4):
        with open(outfile, 'w') as f:
            f.write(np.str(def_amplitude[ivar]))
    outfile='/gws/nopw/j04/acure/lregayre/data_for_emulation/seasonal_amplitude/'+var_list_PPE[ivar]+'_'+region+'_seasonal_amplitude_PPE.dat'
    np.savetxt(outfile,PPE_amplitude[ivar,:])


############################
## Calculate relevant PPE stats
############################

PPE_means= np.zeros((nvar,nmonths))
PPE_min= np.zeros((nvar,nmonths))
PPE_5= np.zeros((nvar,nmonths))
PPE_33= np.zeros((nvar,nmonths))
PPE_67= np.zeros((nvar,nmonths))
PPE_95= np.zeros((nvar,nmonths))
PPE_max= np.zeros((nvar,nmonths))

for ivar in np.arange(nvar):
    for imonth in np.arange(nmonths):
        PPE_means[ivar,imonth]= np.mean(PPE_array[ivar,imonth,:])
        PPE_min[ivar,imonth]= np.min(PPE_array[ivar,imonth,:])
        PPE_33[ivar,imonth]= np.percentile(PPE_array[ivar,imonth,:],33)
        PPE_67[ivar,imonth]= np.percentile(PPE_array[ivar,imonth,:],67)
        PPE_95[ivar,imonth]= np.percentile(PPE_array[ivar,imonth,:],95)
        PPE_5[ivar,imonth]= np.percentile(PPE_array[ivar,imonth,:],5)
        PPE_max[ivar,imonth]= np.max(PPE_array[ivar,imonth,:])


#########
## Plot
#########

font_size=14
mp.rcParams['axes.labelsize'] = font_size*0.9
mp.rcParams['font.size'] = font_size*1.1
mp.rcParams['axes.linewidth'] = font_size*0.1
mp.rcParams['axes.titlesize'] = font_size*1.1
mp.rcParams['legend.fontsize'] = font_size*1.2
mp.rcParams['xtick.labelsize'] = font_size*0.8
mp.rcParams['ytick.labelsize'] = font_size*0.9
mp.rcParams['lines.markersize']=5.0
cmap_pos= plt.get_cmap('viridis')
cmap_pos_inv= plt.get_cmap('viridis_r')
cmap_pos_neg = plt.get_cmap('PiYG')
cmap_pos_neg_inv= plt.get_cmap('PiYG_r')


xvals=np.arange(nmonths+1)
xvals[-1]+=1
xvals[-2]+=1
xlabels=['D','J','F','M','A','M','J','J','A','S','O','N','Ann']
xticks=(0,1,2,3,4,5,6,7,8,9,10,11,13.5)


green_patch = mpatches.Patch(color='green', label='Observed')
red_patch = mpatches.Patch(color='red', label='Default')
black_patch = mpatches.Patch(color='k', label='PPE 66% CI')
Dgray_patch = mpatches.Patch(color='darkgray', label='PPE 90% CI')
Lgray_patch = mpatches.Patch(color='lightgray', label='PPE 100% CI')


def add_subplot(position_index,ivar):
    ax=fig.add_subplot(position_index)
    ax.fill_between(xvals[0:12],PPE_min[ivar,0:12],PPE_max[ivar,0:12],color='lightgray')
    ax.fill_between(xvals[0:12],PPE_5[ivar,0:12],PPE_95[ivar,0:12],color='darkgray')
    ax.fill_between(xvals[0:12],PPE_33[ivar,0:12],PPE_67[ivar,0:12],color='k')
    ax.fill_between(xvals[12:14],(PPE_min[ivar,12],PPE_min[ivar,12]),(PPE_max[ivar,12],PPE_max[ivar,12]),color='lightgray')
    ax.fill_between(xvals[12:14],(PPE_5[ivar,12],PPE_5[ivar,12]),(PPE_95[ivar,12],PPE_95[ivar,12]),color='darkgray')
    ax.fill_between(xvals[12:14],(PPE_33[ivar,12],PPE_33[ivar,12]),(PPE_67[ivar,12],PPE_67[ivar,12]),color='k')
    for imonth in np.arange(12):
        ax.plot(xvals[imonth],obs_array[ivar,imonth],'D',color='green')
        ax.plot(np.mean(xvals[12:14]),obs_array[ivar,-1],'D',color='green')
    for imonth in np.arange(12):
        ax.plot(xvals[imonth],def_array[ivar,imonth],'D',color='red')
        ax.plot(np.mean(xvals[12:14]),def_array[ivar,-1],'D',color='red')
    leg_centre=ax.legend(handles=[black_patch,Dgray_patch,Lgray_patch,green_patch,red_patch],ncol=2,loc='upper center',prop={'size': unique_legend_size})
    ax.set_xlabel('Month')
    ax.set_ylabel('H$_{d}$ ('+units_list[ivar]+')')
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_aspect(0.02)

dir_out='/gws/nopw/j04/acure/lregayre/seasonal_cycle/'
file_out='seasonal_cycle_'+region_file_title
unique_legend_size=10
fig=plt.figure()
fig.suptitle(region_title+' in '+var_name_titles[ivar])
add_subplot(111,0)


plt.savefig(dir_out+file_out+'.pdf',dpi=300,orientation='landscape')


