 ]] #!/usr/bin/env python
"""
Seasonal cycle of state variables in the North Atlantic, with additional sub panels to exemplify the effects of constraint and inconsistencies

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
from os.path import exists
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
region_list=['NN_atlantic']
nregion=len(region_list)
region_list_PPE=['NAtlantic']

icon_list=np.hstack((np.arange(105),106,np.arange(87)+108,np.arange(197)+197,np.arange(55)+395))

constraint_dir='/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/individual_NRMSE_indices/original_additional_and_Hd/'
integer_file_text='/all_1M_python_indices_retained_w_0pc_error_for_variable_'

icon=0 ## spans 0 to 215 (all regional, plus Nd hemispheric from preliminary constraint - additional constraint variables added at later stages)


sample_size=1000000

region_titles=['Northern_North_Atlantic']
month_list=[ 'dec','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','ann']
nmonths=len(month_list)

greek_letters_lower=[chr(code) for code in range(945,970)]
TAU=greek_letters_lower[19]
greek_letters_upper=[chr(code) for code in range(913,938)]
DELTA=greek_letters_upper[3]


var_list_obs= ['CERES_SW_TOA','CDNC','Cloud_Retrieval_Fraction_Liquid','Cloud_Water_Path_Liquid','Cloud_Optical_Thickness_Liquid','Cloud_Particle_Size_Liquid']
var_list_def=['SW_TOA','CDNC','LCF','LWP','Tau_Liquid_MODIS','Reff_Liquid_MODIS'] ## Note COT didn't exist in AR's run. Doesn't exist as a diagnostic, only for MODIS. DG provided.
var_list_PPE=['Fsw','Nd','Fc','LWP','Tau','Re']
var_list_PPE_221=['SW_TOA','CDNC','CF_Liquid_MODIS','LWP_MODIS','Tau_Liquid_MODIS','Reff_Liquid_MODIS']
nvar= len(var_list_obs)
var_name_titles=['$F_{SW}$','$N_{d}$','$f_{c}$','LWP',TAU+'$_{c}$','$r_{e}$']
units_list=['/ W $m^{-2}$','/ cm$^{-3}$','','/ g m$^{-2}$','','']
units_list_trim= ['W $m^{-2}$','cm$^{-3}$','','g m$^{-2}$','','']

dir_PPE='/gws/nopw/j04/acure/lregayre/data_post_emulation/'
dir_PPE_221='/gws/nopw/j04/acure/lregayre/data_for_emulation/'
dir_obs=dir_PPE_221+'observations/'
dir_def=dir_PPE_221+'default/'


##################
## Plot features
##################

font_size=10
mp.rcParams['axes.labelsize'] = font_size*1.2
mp.rcParams['font.size'] = font_size*1
mp.rcParams['axes.linewidth'] = font_size*0.1
mp.rcParams['axes.titlesize'] = font_size*1.0
mp.rcParams['legend.fontsize'] = font_size*1.2
mp.rcParams['xtick.labelsize'] = font_size*0.75
mp.rcParams['ytick.labelsize'] = font_size*1.0
mp.rcParams['lines.markersize']=3.0
inline_text_size=font_size*1.2
aspect_subplots=1.

unique_legend_size=8
ms_default=2

cmap_pos= plt.get_cmap('viridis') 
cmap_pos_inv= plt.get_cmap('viridis_r') 
cmap_pos_neg = plt.get_cmap('PiYG')
cmap_pos_neg_inv= plt.get_cmap('PiYG_r')

xvals=np.arange(nmonths+1)
xvals[-1]+=1
xvals[-2]+=1
xlabels=['D','J','F','M','A','M','J','J','A','S','O','N','Ann']
xticks=(0,1,2,3,4,5,6,7,8,9,10,11,13.5)

xvals2=xvals+16
xlabels1=['','J','','M','','M','','J','','S','','N','','D','','F','','A','','J','','A','','O','','Ann']
xlabels2=['D','','F','','A','','J','','A','','O','','Ann','','J','','M','','M','','J','','S','','N','']
xticks=(0,1,2,3,4,5,6,7,8,9,10,11,13.5,16,17,18,19,20,21,22,23,24,25,26,27,29.5)


green_patch = mpatches.Patch(color='green', label='Observed')
red_patch = mpatches.Patch(color='red', label='Default')
darkorange_patch = mpatches.Patch(color='darkorange', label='Constraint')
black_patch = mpatches.Patch(color='k', label='66% CI')
Dgray_patch = mpatches.Patch(color='darkgray', label='90% CI')
Lgray_patch = mpatches.Patch(color='lightgray', label='Range')

good_color='lightseagreen'
bad_color='orchid'
good_color_light='paleturquoise'
bad_color_light='plum'
good_color_dark='teal'
bad_color_dark='darkorchid'

good_patch= mpatches.Patch(color=good_color_light,label=(var_name_titles[1]+' constraint'))
bad_patch= mpatches.Patch(color=bad_color_light,label=(var_name_titles[3]+' constraint'))


region= region_list[iregion]
region_title= region_titles[iregion]
region_PPE=region_list_PPE[iregion]
print('iregion= '+str(iregion)+' '+region_PPE)
dir_out='/gws/nopw/j04/acure/lregayre/seasonal_cycle/sample_constrained/individual_regional_all_emulated/'
file_out_regional='/seasonal_cycle_'+region_PPE+'_constrained_all_individual_constraint_figures.pdf'
regional_pdf=matplotlib.backends.backend_pdf.PdfPages(dir_out+file_out_regional)
PPE_array_unconstrained=np.zeros((nvar,nmonths,sample_size))
obs_array=np.zeros((nvar,nmonths))
def_array=np.zeros((nvar,nmonths))
flag_221=np.zeros((nvar,nmonths))


for ivar in np.arange(nvar):
    for imonth in np.arange(nmonths):
        print(str(imonth)+' processing '+month_list[imonth]+' for '+var_list_PPE[ivar])
        file_in=dir_PPE+'emulated_mean_values_'+month_list[imonth]+'_'+var_list_PPE[ivar]+'_'+region_PPE+'_'+str(sample_size)+'_w_o_carb.dat'
        if exists(file_in):
            PPE_array_unconstrained[ivar,imonth,:]=np.loadtxt(file_in)
        else:
            file_in_221=dir_PPE_221+month_list[imonth]+'/'+var_list_PPE_221[ivar]+'_PPE_'+region+'_'+month_list[imonth]+'_revised_match_MODIS_CF.dat'
            PPE_array_unconstrained[ivar,imonth,0:221]=np.loadtxt(file_in_221)
            flag_221[ivar,imonth]=1.0
            print('Using 221 PPE data for '+month_list[imonth]+' '+var_list_PPE_221[ivar]+' for '+region)



# LWP units conversion from kg/m2 to g/m2
PPE_array_unconstrained[3,:,:]= PPE_array_unconstrained[3,:,:]*1000
# Reff Liquid units conversion to rescale diagnostic by reverting the packing problem scaling
PPE_array_unconstrained[5,:,:] = PPE_array_unconstrained[5,:,:]*1000000


# Read in obs data, regridded and degraded to match model resolution
obs_array[0]= np.loadtxt(dir_obs+var_list_obs[0]+'_'+region+'_dec_to_ann_revised_match_MODIS_CF.dat')
for ivar in (np.arange(nvar-1)+1):
    print(var_list_PPE_221[ivar])
    obs_array[ivar]= np.loadtxt(dir_obs+var_list_obs[ivar]+'_'+region+'_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')

for ivar in np.arange(nvar):
    for imonth in np.arange(nmonths-1):
        month=month_list[imonth]
        def_array[ivar,imonth]= np.loadtxt(dir_def+month+'/'+var_list_def[ivar]+'_default_'+region+'_'+month+'.dat')
    def_array[ivar,12]= np.mean(def_array[ivar,0:12])


##########################################
## Apply constraints to PPE_array_unconstrained
## for Nd Nov 
## and LWP Nov
## using indices from emulation process
##########################################

good_constraint_indices= np.loadtxt('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/individual_NRMSE_indices/original_additional_and_Hd/all_1M_python_indices_retained_w_0pc_error_for_variable_43.dat') ## Nd in North Atlantic Nov==43

PPE_array_good=np.zeros((PPE_array_unconstrained.shape[0],PPE_array_unconstrained.shape[1],good_constraint_indices.shape[0]))

for ivar in np.arange(PPE_array_unconstrained.shape[0]):
    print(ivar)
    for imonth in np.arange(PPE_array_unconstrained.shape[1]):
        for imember in np.arange(good_constraint_indices.shape[0]):
            PPE_array_good[ivar,imonth,imember]=PPE_array_unconstrained[ivar,imonth,int(good_constraint_indices[imember])]


bad_constraint_indices= np.loadtxt('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/individual_NRMSE_indices/original_additional_and_Hd/all_1M_python_indices_retained_w_0pc_error_for_variable_85.dat') ## LWP Nov==85

PPE_array_bad=np.zeros((PPE_array_unconstrained.shape[0],PPE_array_unconstrained.shape[1],bad_constraint_indices.shape[0]))

for ivar in np.arange(PPE_array_unconstrained.shape[0]):
    print(ivar)
    for imonth in np.arange(PPE_array_unconstrained.shape[1]):
        for imember in np.arange(bad_constraint_indices.shape[0]):
            PPE_array_bad[ivar,imonth,imember]=PPE_array_unconstrained[ivar,imonth,int(bad_constraint_indices[imember])]


PPE_means= np.zeros((nvar,nmonths))
PPE_min= np.zeros((nvar,nmonths))
PPE_5= np.zeros((nvar,nmonths))
PPE_17= np.zeros((nvar,nmonths))
PPE_83= np.zeros((nvar,nmonths))
PPE_95= np.zeros((nvar,nmonths))
PPE_max= np.zeros((nvar,nmonths))
PPE_good_5=np.zeros((nvar,nmonths))
PPE_good_95=np.zeros((nvar,nmonths))
PPE_bad_5=np.zeros((nvar,nmonths))
PPE_bad_95=np.zeros((nvar,nmonths))
PPE_good_min=np.zeros((nvar,nmonths))
PPE_good_max=np.zeros((nvar,nmonths))
PPE_bad_min=np.zeros((nvar,nmonths))
PPE_bad_max=np.zeros((nvar,nmonths))
PPE_good_17=np.zeros((nvar,nmonths))
PPE_good_83=np.zeros((nvar,nmonths))
PPE_bad_17=np.zeros((nvar,nmonths))
PPE_bad_83=np.zeros((nvar,nmonths))


for ivar in np.arange(nvar):
    for imonth in np.arange(nmonths):
        print(str(ivar)+' '+month_list[imonth])
        PPE_means[ivar,imonth]= np.mean(PPE_array_unconstrained[ivar,imonth,:])
        PPE_min[ivar,imonth]= np.min(PPE_array_unconstrained[ivar,imonth,:])
        PPE_17[ivar,imonth]= np.percentile(PPE_array_unconstrained[ivar,imonth,:],17)
        PPE_83[ivar,imonth]= np.percentile(PPE_array_unconstrained[ivar,imonth,:],83)
        PPE_95[ivar,imonth]= np.percentile(PPE_array_unconstrained[ivar,imonth,:],95)
        PPE_5[ivar,imonth]= np.percentile(PPE_array_unconstrained[ivar,imonth,:],5)
        PPE_max[ivar,imonth]= np.max(PPE_array_unconstrained[ivar,imonth,:])
        PPE_good_5[ivar,imonth]= np.percentile(PPE_array_good[ivar,imonth,:],5)
        PPE_good_95[ivar,imonth]= np.percentile(PPE_array_good[ivar,imonth,:],95)
        PPE_bad_5[ivar,imonth]= np.percentile(PPE_array_bad[ivar,imonth,:],5)
        PPE_bad_95[ivar,imonth]= np.percentile(PPE_array_bad[ivar,imonth,:],95)
        PPE_good_min[ivar,imonth]= np.min(PPE_array_good[ivar,imonth,:])
        PPE_good_max[ivar,imonth]= np.max(PPE_array_good[ivar,imonth,:])
        PPE_bad_min[ivar,imonth]= np.min(PPE_array_bad[ivar,imonth,:])
        PPE_bad_max[ivar,imonth]= np.max(PPE_array_bad[ivar,imonth,:])
        PPE_good_17[ivar,imonth]= np.percentile(PPE_array_good[ivar,imonth,:],17)
        PPE_good_83[ivar,imonth]= np.percentile(PPE_array_good[ivar,imonth,:],83)
        PPE_bad_17[ivar,imonth]= np.percentile(PPE_array_bad[ivar,imonth,:],17)
        PPE_bad_83[ivar,imonth]= np.percentile(PPE_array_bad[ivar,imonth,:],83)


def add_subplot(position_index1,position_index2,position_index3,ivar,legend_indicator,xaxis_label_indicator,inline_text,ypos_for_text,yaxis_label_right,arrow_index,iplot_constraint):
    ax=fig.add_subplot(position_index1,position_index2,position_index3)
    ax.fill_between(xvals[0:12],PPE_min[ivar,0:12],PPE_max[ivar,0:12],color='lightgray')
    ax.fill_between(xvals[0:12],PPE_5[ivar,0:12],PPE_95[ivar,0:12],color='darkgray')
    ax.fill_between(xvals[0:12],PPE_17[ivar,0:12],PPE_83[ivar,0:12],color='k')
    ax.fill_between(xvals[12:14],(PPE_min[ivar,12],PPE_min[ivar,12]),(PPE_max[ivar,12],PPE_max[ivar,12]),color='lightgray')
    ax.fill_between(xvals[12:14],(PPE_5[ivar,12],PPE_5[ivar,12]),(PPE_95[ivar,12],PPE_95[ivar,12]),color='darkgray')
    ax.fill_between(xvals[12:14],(PPE_17[ivar,12],PPE_17[ivar,12]),(PPE_83[ivar,12],PPE_83[ivar,12]),color='k')
    ax.fill_between(xvals2[0:12],PPE_good_17[ivar,0:12],PPE_good_83[ivar,0:12],color=good_color,alpha=con_alpha,edgecolor=good_color_dark)
    ax.fill_between(xvals2[0:12],PPE_bad_17[ivar,0:12],PPE_bad_83[ivar,0:12],color=bad_color,alpha=con_alpha,edgecolor=bad_color_dark)
    ax.fill_between(xvals2[12:14],(PPE_good_17[ivar,12],PPE_good_17[ivar,12]),(PPE_good_83[ivar,12],PPE_good_83[ivar,12]),color=good_color,alpha=con_alpha,edgecolor=good_color_dark)
    ax.fill_between(xvals2[12:14],(PPE_bad_17[ivar,12],PPE_bad_17[ivar,12]),(PPE_bad_83[ivar,12],PPE_bad_83[ivar,12]),color=bad_color,alpha=con_alpha,edgecolor=bad_color_dark)
    for imonth in np.arange(12):
        if (obs_array[ivar,imonth]>0.0):
            ax.plot(xvals[imonth],obs_array[ivar,imonth],'D',color='green',ms=ms_default)
            ax.plot(xvals2[imonth],obs_array[ivar,imonth],'D',color='green',ms=ms_default)
    if (obs_array[ivar,12]>0.0):
            ax.plot(np.mean(xvals[12:14]),obs_array[ivar,-1],'D',color='green',ms=ms_default)
            ax.plot(np.mean(xvals2[12:14]),obs_array[ivar,-1],'D',color='green',ms=ms_default)
    if (iplot_constraint==1):
            ax.plot(xvals2[11],obs_array[ivar,11],'D',color='darkorange',ms=ms_default)
    for imonth in np.arange(12):
        if (def_array[ivar,imonth]>0.0):
            ax.plot(xvals[imonth],def_array[ivar,imonth],'D',color='red',ms=ms_default)
    if (def_array[ivar,12]>0.0):
            ax.plot(np.mean(xvals[12:14]),def_array[ivar,-1],'D',color='red',ms=ms_default)
    if (legend_indicator==1):
        leg=ax.legend(handles=[green_patch,red_patch],loc='lower right',prop={'size': unique_legend_size},frameon=False)
    if (legend_indicator==2):
        leg=ax.legend(handles=[black_patch,Dgray_patch,Lgray_patch],loc='upper right',prop={'size': unique_legend_size},frameon=False)
    if (legend_indicator==3):
        leg=plt.legend(handles=[good_patch,bad_patch],loc='upper right',prop={'size': unique_legend_size},frameon=False)
        leg2=plt.legend(handles=[darkorange_patch],loc='lower right',prop={'size': unique_legend_size},frameon=False)
        ax.add_artist(leg)
        ax.add_artist(leg2)
    ax.set_xticks(xticks)
    ax.set_xticklabels([])
    if (xaxis_label_indicator==0):
        ax.set_xticklabels([])
        ax.set_xticks([])
    if (xaxis_label_indicator==2):
        ax.set_xlabel('Month')
        ax.set_xticklabels(xlabels1)
    if (xaxis_label_indicator==3):
        ax.set_xlabel('Month')
        ax.set_xticklabels(xlabels2)
    if (yaxis_label_right==1):
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        ax.set_ylabel(units_list_trim[ivar],rotation=text_rotation,labelpad=labelpad_rotation)
    else:
        ax.set_ylabel(units_list_trim[ivar])
    ax.set_ylim([ymin_list[ivar],ymax_list[ivar]])
    ax.annotate(inline_text, xy=((xstart,ypos_for_text)),  xycoords='data',horizontalalignment='left', verticalalignment='center',fontsize=inline_text_size)
    xvals_ax,yvals_ax = ax.axes.get_xlim(),ax.axes.get_ylim()
    xrange = xvals_ax[1]-xvals_ax[0]
    yrange = yvals_ax[1]-yvals_ax[0]


ypos_list=[180.,380.,0.64,245,42.,17.]
ymin_list=[40.,0.,0.1,40.,-5.,5.]
ymax_list=[200.,450.,0.72,280.,50.,19.0]

con_alpha=0.2
text_rotation=270
labelpad_rotation=15
labelpad_unrotated=2
xstart=-1.

fig=plt.figure()
add_subplot(3,2,1,0,0,0,'a) '+var_name_titles[0],ypos_list[0],0,0,0)
add_subplot(3,2,2,1,2,0,'b) '+var_name_titles[1],ypos_list[1],1,1,1)
add_subplot(3,2,3,2,0,0,'c) '+var_name_titles[2],ypos_list[2],0,0,0)
add_subplot(3,2,4,3,3,0,'d) '+var_name_titles[3],ypos_list[3],1,2,1)
add_subplot(3,2,5,4,0,2,'e) '+var_name_titles[4],ypos_list[4],0,0,0)
add_subplot(3,2,6,5,1,3,'f) '+var_name_titles[5],ypos_list[5],1,0,0)

fig.subplots_adjust(wspace=0.025,hspace=0.05)

plt.savefig('/gws/nopw/j04/acure/lregayre/seasonal_cycle/sample_constrained/seasonal_cycle_12panel_example_of_constraint_reordered_legends.pdf')

plt.show()

