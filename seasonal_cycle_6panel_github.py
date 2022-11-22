#!/usr/bin/env python
"""
Plot seasonal cycle of any state variable, across regions

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

iregion=0  ## Ranges from 0-5, as includes global mean
region_list=['global','NE_pacific','SE_pacific','SE_atlantic','NN_atlantic','N_southern_ocean']
region_titles=['Global','North Pacific','South Pacific','South Atlantic','North Atlantic','Southern Ocean']
region_fname_titles=['global','NPacific','SPacific','SAtlantic','NAtlantic','SOcean']
month_list=[ 'dec','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','ann']
nmonths=len(month_list)

region= region_list[iregion]
region_title= region_titles[iregion]
region_fname_title= region_fname_titles[iregion]

greek_letters_lower=[chr(code) for code in range(945,970)]
TAU=greek_letters_lower[19]
greek_letters_upper=[chr(code) for code in range(913,938)]
DELTA=greek_letters_upper[3]

var_list_obs= ['CERES_SW_TOA','Cloud_Water_Path_Liquid','CDNC','Cloud_Retrieval_Fraction_Liquid','Cloud_Optical_Thickness_Liquid','Cloud_Particle_Size_Liquid']
var_list_PPE=['SW_TOA','LWP_MODIS','CDNC','CF_Liquid_MODIS','Tau_Liquid_MODIS','Reff_Liquid_MODIS']
var_list_PPE_out=['Fsw','LWP','Nd','Fc','Tau','Re']
var_list_def=['SW_TOA','LWP','CDNC','LCF','Tau_Liquid_MODIS','Reff_Liquid_MODIS']
nvar= len(var_list_obs)
var_name_titles=['$F_{SW}$','LWP','$N_{d}$','$f_{c}$',TAU+'$_{c}$','$r_{e}$']
units_list=['W $m^{-2}$','g m$^{-2}$','cm$^{-3}$','','','']

dir_PPE='/gws/nopw/j04/acure/lregayre/data_for_emulation/'
dir_obs=dir_PPE+'observations/'
dir_def=dir_PPE+'default/'

PPE_array=np.zeros((nvar,nmonths,221))
obs_array=np.zeros((nvar,nmonths))
def_array=np.zeros((nvar,nmonths))

for ivar in np.arange(nvar):
    for imonth in np.arange(nmonths):
        print(imonth)
        PPE_array[ivar,imonth,:]=np.loadtxt(dir_PPE+month_list[imonth]+'/'+var_list_PPE[ivar]+'_PPE_'+region+'_'+month_list[imonth]+'_revised_match_MODIS_CF.dat')

for ivar in np.arange(nvar):
    for ippe in np.arange(221):
        PPE_array[ivar,12,ippe]= np.mean(PPE_array[ivar,0:12,ippe])


############################################################
# Read in obs data, regridded and degraded to match model resolution
############################################################

obs_array[0]= np.loadtxt(dir_obs+var_list_obs[0]+'_'+region+'_dec_to_ann_revised_match_MODIS_CF.dat')

for ivar in (np.arange(nvar-1)+1):
    print(var_list_PPE[ivar])
    obs_array[ivar]= np.loadtxt(dir_obs+var_list_obs[ivar]+'_'+region+'_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')


#################################################
## Read in default data
## annual needs to be calculated here
#################################################


for ivar in np.arange(nvar):
    for imonth in np.arange(nmonths-1):
        month=month_list[imonth]
        def_array[ivar,imonth]= np.loadtxt(dir_def+month+'/'+var_list_def[ivar]+'_default_'+region+'_'+month+'.dat')
    def_array[ivar,12]= np.mean(def_array[ivar,0:12])


# LWP units conversion from kg/m2 to g/m2
PPE_array[1,:,:]= PPE_array[1,:,:]*1000

# Reff Liquid units conversion to rescale diagnostic by reverting the packing problem scaling
PPE_array[5,:,:] = PPE_array[5,:,:]*1000000


#####################################
## Calculate amplidues of seasonal cycles
## but for observations, rule out min 
## where observed value is  -999 !!
#####################################

obs_amplitude=np.zeros((nvar))
def_amplitude=np.zeros((nvar))
PPE_amplitude=np.zeros((nvar,221))

for ivar in np.arange(nvar):
    obs_subset=[]
    for imonth in np.arange(nmonths-1):
        if (obs_array[ivar,imonth]!=-999.0):
            obs_subset=np.append(obs_subset,obs_array[ivar,imonth])
    obs_amplitude[ivar]=np.max(obs_subset)-np.min(obs_subset)
    def_amplitude[ivar]=np.max(def_array[ivar,0:12])-np.min(def_array[ivar,0:12])
    for ippe in np.arange(221):
        PPE_amplitude[ivar,ippe]=np.max(PPE_array[ivar,0:12,ippe])-np.min(PPE_array[ivar,0:12,ippe])

## write to file

for ivar in np.arange(nvar):
    outfile='/gws/nopw/j04/acure/lregayre/data_for_emulation/seasonal_amplitude/'+var_list_obs[ivar]+'_'+region+'_seasonal_amplitude_observed_revised_match_MODIS_CF.dat'
    with open(outfile, 'w') as f:
        f.write(np.str(obs_amplitude[ivar]))
    outfile='/gws/nopw/j04/acure/lregayre/data_for_emulation/seasonal_amplitude/'+var_list_def[ivar]+'_'+region+'_seasonal_amplitude_default_revised_match_MODIS_CF.dat'
    with open(outfile, 'w') as f:
        f.write(np.str(def_amplitude[ivar]))
    outfile='/gws/nopw/j04/acure/lregayre/data_for_emulation/seasonal_amplitude/'+var_list_PPE_out[ivar]+'_'+region+'_seasonal_amplitude_PPE_revised_match_MODIS_CF.dat'
    np.savetxt(outfile,PPE_amplitude[ivar,:])


#####################
## Calculate PPE stats
#####################

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
        PPE_33[ivar,imonth]= np.percentile(PPE_array[ivar,imonth,:],(50-33))
        PPE_67[ivar,imonth]= np.percentile(PPE_array[ivar,imonth,:],(50+33))
        PPE_95[ivar,imonth]= np.percentile(PPE_array[ivar,imonth,:],95)
        PPE_5[ivar,imonth]= np.percentile(PPE_array[ivar,imonth,:],5)
        PPE_max[ivar,imonth]= np.max(PPE_array[ivar,imonth,:])


#########
## Plot 
#########

font_size=10
mp.rcParams['axes.labelsize'] = font_size*0.9
mp.rcParams['font.size'] = font_size*1
mp.rcParams['axes.linewidth'] = font_size*0.1
mp.rcParams['axes.titlesize'] = font_size*1.1
mp.rcParams['legend.fontsize'] = font_size*1.2
mp.rcParams['xtick.labelsize'] = font_size*0.7
mp.rcParams['ytick.labelsize'] = font_size*0.9
mp.rcParams['lines.markersize']=3.0
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
black_patch = mpatches.Patch(color='k', label='66% CI')
Dgray_patch = mpatches.Patch(color='darkgray', label='90% CI')
Lgray_patch = mpatches.Patch(color='lightgray', label='100% CI')


def add_subplot(position_index,ivar,legend_indicator,ivar_units):
    ax=fig.add_subplot(position_index)
    ax.fill_between(xvals[0:12],PPE_min[ivar,0:12],PPE_max[ivar,0:12],color='lightgray')
    ax.fill_between(xvals[0:12],PPE_5[ivar,0:12],PPE_95[ivar,0:12],color='darkgray')
    ax.fill_between(xvals[0:12],PPE_33[ivar,0:12],PPE_67[ivar,0:12],color='k')
    ax.fill_between(xvals[12:14],(PPE_min[ivar,12],PPE_min[ivar,12]),(PPE_max[ivar,12],PPE_max[ivar,12]),color='lightgray')
    ax.fill_between(xvals[12:14],(PPE_5[ivar,12],PPE_5[ivar,12]),(PPE_95[ivar,12],PPE_95[ivar,12]),color='darkgray')
    ax.fill_between(xvals[12:14],(PPE_33[ivar,12],PPE_33[ivar,12]),(PPE_67[ivar,12],PPE_67[ivar,12]),color='k')
    for imonth in np.arange(12):
        if (obs_array[ivar,imonth]>0.0):
            ax.plot(xvals[imonth],obs_array[ivar,imonth],'D',color='green')
    if (obs_array[ivar,12]>0.0):
            ax.plot(np.mean(xvals[12:14]),obs_array[ivar,-1],'D',color='green')
    for imonth in np.arange(12):
        if (def_array[ivar,imonth]>0.0):
            ax.plot(xvals[imonth],def_array[ivar,imonth],'D',color='red')
    if (def_array[ivar,12]>0.0):
            ax.plot(np.mean(xvals[12:14]),def_array[ivar,-1],'D',color='red')
    if (legend_indicator==1):
        leg=ax.legend(handles=[green_patch,red_patch],loc='upper right',prop={'size': unique_legend_size})
    if (legend_indicator==2):
        leg=ax.legend(handles=[black_patch,Dgray_patch,Lgray_patch],loc='upper right',prop={'size': unique_legend_size})
    ax.set_xlabel('Month')
    if ivar_units==0:
        ax.set_ylabel(var_name_titles[ivar])
    else:
        ax.set_ylabel(var_name_titles[ivar]+' / '+units_list[ivar])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)


dir_out='/gws/nopw/j04/acure/lregayre/seasonal_cycle/'
file_out='seasonal_cycle_'+region_fname_title+'_with_COD_CF_and_REVISED_default_revised_match_MODIS_CF'


unique_legend_size=6

fig=plt.figure()
add_subplot(231,0,1,1)
add_subplot(232,1,0,1)
add_subplot(233,2,2,1)
add_subplot(234,3,0,0)
add_subplot(235,4,0,0)
add_subplot(236,5,0,0)
fig.subplots_adjust(top=0.9)
fig.suptitle(region_title,y=0.95)
fig.subplots_adjust(wspace=0.5,hspace=0.25)
plt.savefig(dir_out+file_out+'.pdf',dpi=300,orientation='landscape')





