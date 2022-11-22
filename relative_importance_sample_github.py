#!/usr/bin/env python
"""
caluclate and plot relative importance metrics from sample of 1 million model variants

REQUIRES Python 3.7 OR PINGOUIN DF.PCORR
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
import pandas as pd
import pingouin as pg ## NEEDS JASPY/3.7

nppe = 221
nsample=1000000
lower_limit=5. ## plot parameter contributions lower than this limit (%) are not plotted

region_list=['global']
region_titles=['Global']
iregion=0
region= region_list[iregion]
region_title= region_titles[iregion]

greek_letters_lower=[chr(code) for code in range(945,970)]
TAU=greek_letters_lower[19]
greek_letters_upper=[chr(code) for code in range(913,938)]
DELTA=greek_letters_upper[3]


var_list_PPE=['ERF','ACI','ARI']
nvar= len(var_list_PPE)
var_name_titles=[DELTA+' F$_{aer}$',DELTA+' F$_{aci}$',DELTA+' F$_{ari}$']
units_list=['W m$^{-2}$','W m$^{-2}$','W m$^{-2}$']
dir_sample='/gws/nopw/j04/acure/lregayre/data_post_emulation/'

PPE_array=np.zeros((nvar,nsample))

for ivar in np.arange(nvar):
    PPE_array[ivar,:]=np.loadtxt(dir_sample+'emulated_mean_values_'+var_list_PPE[ivar]+'_'+str(nsample)+'_w_o_carb.dat')



################################
################################
## Read in parameter design values
## UNIT SCALE
################################
################################

par_index_for_plotting_60=(0,1,2,23,24,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,47,48,49,50,51,52,53,54,55,56,57,58)
par_index_for_plotting_55=(0,1,2,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,43,44,45,46,47,48,49,50,51,52,53,54) ## 37 actually used
npar=55

sample_file_norm='/home/users/lregayre/R/ACURE/samples/constrained_multi_million_sample_first_million.dat'

with open(sample_file_norm,'r+') as f:
	LHC_stripped=f.read().splitlines()

LHC_all_norm=np.zeros((nsample,npar))
for isample in np.arange(nsample):
    LHC_all_norm[isample,:]=[float(x) for x in LHC_stripped[isample].split(' ')]


LHC_norm=LHC_all_norm

parameter_names_60= ['bl_nuc','ait_width','cloud_ph','carb_ff_ems','carb_ff_ems_eur','carb_ff_ems_nam','carb_ff_ems_chi','carb_ff_ems_asi','carb_ff_ems_mar','carb_ff_ems_r','carb_bb_ems','carb_bb_ems_sam','carb_bb_ems_naf','carb_bb_ems_saf','carb_bb_ems_bnh','carb_bb_ems_rnh','carb_bb_ems_rsh','carb_res_ems','carb_res_ems_chi','carb_res_ems_asi','carb_res_ems_afr','carb_res_ems_lat','carb_res_ems_r','carb_ff_diam','carb_bb_diam','carb_res_diam','prim_so4_diam','sea_spray','anth_so2','anth_so2_chi','anth_so2_asi','anth_so2_eur','anth_so2_nam','anth_so2_r','volc_so2','bvoc_soa','dms','prim_moc','dry_dep_ait','dry_dep_acc','dry_dep_so2','kappa_oc','sig_w','rain_frac','cloud_ice_thresh','conv_plume_scav','scav_diam','bc_ri','oxidants_oh','oxidants_o3','bparam','two_d_fsd_factor','c_r_correl','autoconv_exp_lwp','autoconv_exp_nd','dbsdtbs_turb_0','ai','m_ci','a_ent_1_rp','pcalc_index INDEX ONLY']

parameter_names_55= ['bl_nuc','ait_width','cloud_ph','carb_ff_ems_eur','carb_ff_ems_nam','carb_ff_ems_chi','carb_ff_ems_asi','carb_ff_ems_mar','carb_ff_ems_r','carb_bb_ems_sam','carb_bb_ems_naf','carb_bb_ems_saf','carb_bb_ems_bnh','carb_bb_ems_rnh','carb_bb_ems_rsh','carb_res_ems_chi','carb_res_ems_asi','carb_res_ems_afr','carb_res_ems_lat','carb_res_ems_r','carb_ff_diam','carb_bb_diam','carb_res_diam','prim_so4_diam','sea_spray','anth_so2_chi','anth_so2_asi','anth_so2_eur','anth_so2_nam','anth_so2_r','volc_so2','bvoc_soa','dms','prim_moc','dry_dep_ait','dry_dep_acc','dry_dep_so2','kappa_oc','sig_w','rain_frac','cloud_ice_thresh','conv_plume_scav','scav_diam','bc_ri','oxidants_oh','oxidants_o3','bparam','two_d_fsd_factor','c_r_correl','autoconv_exp_lwp','autoconv_exp_nd','dbsdtbs_turb_0','ai','m_ci','a_ent_1_rp']

color_list_55= ["black","red","goldenrod","darkorange","orange","moccasin","burlywood","tan","blanchedalmond","brown","firebrick","maroon","darkred","indianred","lightcoral","khaki","palegoldenrod","darkkhaki","olive","darkolivegreen","yellow","mediumaquamarine","hotpink","plum","thistle","blue","royalblue","cornflowerblue","mediumblue","darkblue","springgreen","palegreen","darkslategrey","forestgreen","rosybrown","mediumpurple","cyan","lightpink","lightseagreen","darkmagenta","slategrey","lightskyblue","cadetblue","deepskyblue","steelblue","magenta","salmon","sienna","deeppink","blueviolet","purple","palevioletred","silver","gray","orchid"]

parameter_names=parameter_names_55
color_list=color_list_55
par_index_for_plotting=par_index_for_plotting_55

###########################################
## Make patches for each paramter in case needed
## Put into a list for extraction along with other data
###########################################

patch_list=[]
for ipar in np.arange(len(parameter_names)):
    patch_list=np.append(patch_list,mpatches.Patch(color=color_list[ipar], label=parameter_names[ipar]))


##############################################################
## Make LHC_norm subset so that dataframes can be produced with ease
##############################################################

parameter_names_subset=[]
LHC_norm_subset=np.zeros((nsample,len(par_index_for_plotting)))
LHC_norm_subset_transpose=np.zeros((len(par_index_for_plotting),nsample))
for ipar in np.arange(len(par_index_for_plotting)):
    ipar_subset=par_index_for_plotting[ipar]
    LHC_norm_subset[:,ipar]=LHC_norm[:,ipar_subset]
    LHC_norm_subset_transpose[ipar,:]=LHC_norm[:,ipar_subset]
    parameter_names_subset.append(parameter_names[ipar_subset])



##########################################################################
## Calculate partial correlations
##########################################################################

par_cor_sign=np.zeros((nvar,npar))
par_cor_r2vals=np.zeros((nvar,npar))



for ivar in np.arange(nvar):
    print('ivar:'+str(ivar))
    forcing_for_hstack=np.zeros((nsample,1))
    forcing_for_hstack[:,0]=PPE_array[ivar,:]
    forcing_w_sample=np.hstack((forcing_for_hstack,LHC_norm_subset))
    df= pd.DataFrame(columns=np.hstack((['var'],parameter_names_subset)),data=forcing_w_sample)
    par_cor_array=df.pcorr()
    par_cor_w_var=np.array(par_cor_array)[0,1:len(parameter_names_subset)+1]
    for ipar in np.arange(len(par_index_for_plotting)): 
        print(str(ipar))
        ipar_subset=par_index_for_plotting[ipar]
        par_cor_r2vals[ivar,ipar_subset]=par_cor_w_var[ipar]**2
        if (par_cor_w_var[ipar]>0.0):
            par_cor_sign[ivar,ipar_subset]=1.0
        else:
            par_cor_sign[ivar,ipar_subset]=-1.0


################################
## Collate information to create figure
################################


parameter_name_var=[]
par_cor_sign_var=[]
par_cor_r2_vals_var=[]
var_list=[]
xvals_var=[]
bar_colors_var=[]
xlist=(0,1,2)
legend_names_var=[]
legend_colors_var=[]
legend_patch_var=[]

par_cor_pvals=np.zeros((nvar,npar))

for ivar in np.arange(nvar):
    var=var_list_PPE[ivar]
    xval=xlist[ivar]
    for ipar in par_index_for_plotting:
        if (par_cor_pvals[ivar,ipar]<1.0): ## redundant
            if parameter_names[ipar] not in parameter_name_var: ## want unique only
                legend_names_var=np.append(legend_names_var,parameter_names[ipar])
                legend_colors_var=np.append(legend_colors_var,color_list[ipar])
                legend_patch_var=np.append(legend_patch_var,patch_list[ipar])
            parameter_name_var=np.append(parameter_name_var,parameter_names[ipar])
            par_cor_r2_vals_var=np.append(par_cor_r2_vals_var,par_cor_r2vals[ivar,ipar])
            par_cor_sign_var=np.append(par_cor_sign_var,par_cor_sign[ivar,ipar])
            var_list= np.append(var_list,var)
            xvals_var=np.append(xvals_var,xval)
            bar_colors_var=np.append(bar_colors_var,color_list[ipar])


#########################################################
### Calculate percentages (Relative Importance metric) for each term
#########################################################

par_cor_r2_percent_var=np.zeros((len(var_list)))

for ivar in np.arange(nvar):
    var=var_list_PPE[ivar]
    par_cor_var_total=0.0
    for imember in np.arange(len(var_list)):
        if (var_list[imember]==var):
            par_cor_var_total+=par_cor_r2_vals_var[imember]
    for imember in np.arange(len(var_list)):
        if (var_list[imember]==var):
            par_cor_r2_percent_var[imember]=par_cor_r2_vals_var[imember]/par_cor_var_total*100


############################################################
## Restrict plotting to those with Relative Importance metric > lower_limit
############################################################


parameter_name_var_trim=[]
par_cor_sign_var_trim=[]
par_cor_r2_vals_var_trim=[]
var_list_trim=[]
xvals_var_trim=[]
bar_colors_var_trim=[]

par_cor_r2_percent_var_trim=[]

for imember in np.arange(len(var_list)):
    if (par_cor_r2_percent_var[imember]>lower_limit):
        parameter_name_var_trim=np.append(parameter_name_var_trim,parameter_name_var[imember])
        par_cor_sign_var_trim=np.append(par_cor_sign_var_trim,par_cor_sign_var[imember])
        par_cor_r2_vals_var_trim=np.append(par_cor_r2_vals_var_trim,par_cor_r2_vals_var[imember])
        var_list_trim=np.append(var_list_trim,var_list[imember])
        xvals_var_trim=np.append(xvals_var_trim,xvals_var[imember])
        bar_colors_var_trim=np.append(bar_colors_var_trim,bar_colors_var[imember])
        par_cor_r2_percent_var_trim=np.append(par_cor_r2_percent_var_trim,par_cor_r2_percent_var[imember])


legend_names_var_trim=[]
legend_colors_var_trim=[]
legend_patch_var_trim=[]


for imember in np.arange(len(var_list)):
    for ipar in par_index_for_plotting:
        if (np.logical_and((parameter_names[ipar] in parameter_name_var_trim),(parameter_names[ipar] not in legend_names_var_trim))):
            legend_names_var_trim=np.append(legend_names_var_trim,parameter_names[ipar])
            legend_colors_var_trim=np.append(legend_colors_var_trim,color_list[ipar])
            legend_patch_var_trim=np.append(legend_patch_var_trim,patch_list[ipar])



#######################################################################################################
## Reorganise order of legend components so that they will be identical across all figures, as per original order of parameters
#######################################################################################################


ordered_legend_names_var=[]
ordered_legend_colors_var=[]
ordered_legend_patches_var=[] ## redundant

for ipar in par_index_for_plotting:
    if (parameter_names[ipar] in legend_names_var_trim):
        ordered_legend_names_var=np.append(ordered_legend_names_var,parameter_names[ipar])
        ordered_legend_colors_var=np.append(ordered_legend_colors_var,color_list[ipar])
        ordered_legend_patches_var=np.append(ordered_legend_patches_var,patch_list[ipar])


nunique=len(ordered_legend_names_var)
print(ordered_legend_names_var)


#################################################
## Define patches
## maimum of 37 parameters to plot associated patches
#################################################

if nunique>36:
    ipar=36
    P36= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>35:
    ipar=35
    P35= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>34:
    ipar=34
    P34= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>33:
    ipar=33
    P33= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>32:
    ipar=32
    P32= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>31:
    ipar=31
    P31= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>30:
    ipar=30
    P30= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>29:
    ipar=29
    P29= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>28:
    ipar=28
    P28= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>27:
    ipar=27
    P27= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>26:
    ipar=26
    P26= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>25:
    ipar=25
    P25= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>24:
    ipar=24
    P24= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>23:
    ipar=23
    P23= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>22:
    ipar=22
    P22= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>21:
    ipar=21
    P21= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>20:
    ipar=20
    P20= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>19:
    ipar=19
    P19= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>18:
    ipar=18
    P18= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>17:
    ipar=17
    P17= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>16:
    ipar=16
    P16= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>15:
    ipar=15
    P15= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>14:
    ipar=14
    P14= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>13:
    ipar=13
    P13= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>12:
    ipar=12
    P12= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>11:
    ipar=11
    P11= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>10:
    ipar=10
    P10= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>9:
    ipar=9
    P9= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>8:
    ipar=8
    P8= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>7:
    ipar=7
    P7= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>6:
    ipar=6
    P6= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>5:
    ipar=5
    P5= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>4:
    ipar=4
    P4= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>3:
    ipar=3
    P3= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>2:
    ipar=2
    P2= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>1:
    ipar=1
    P1= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])

if nunique>0:
    ipar=0
    P0= mpatches.Patch(color=ordered_legend_colors_var[ipar],label=ordered_legend_names_var[ipar])



##################################################
##  make legend handles
##################################################

if nunique==0:
    print('Error - PROBLEM WITH NUMBER OF UNIQUE ELEMENTS TO BE PLOTTED. NO VIABLE DATA')
elif nunique==1:
    legend_handles_var=[P0]
elif nunique==2:
    legend_handles_var=[P0,P1]
elif nunique==3:
    legend_handles_var=[P0,P1,P2]
elif nunique==4:
    legend_handles_var=[P0,P1,P2,P3]
elif nunique==5:
    legend_handles_var=[P0,P1,P2,P3,P4]
elif nunique==6:
    legend_handles_var=[P0,P1,P2,P3,P4,P5]
elif nunique==7:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6]
elif nunique==8:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7]
elif nunique==9:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8]
elif nunique==10:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9]
elif nunique==11:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10]
elif nunique==12:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11]
elif nunique==13:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12]
elif nunique==14:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13]
elif nunique==15:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14]
elif nunique==16:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15]
elif nunique==17:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16]
elif nunique==18:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17]
elif nunique==19:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18]
elif nunique==20:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19]
elif nunique==21:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20]
elif nunique==22:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21]
elif nunique==23:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22]
elif nunique==24:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23]
elif nunique==25:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24]
elif nunique==26:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24,P25]
elif nunique==27:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24,P25,P26]
elif nunique==28:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24,P25,P26,P27]
elif nunique==29:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24,P25,P26,P27,P28]
elif nunique==30:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24,P25,P26,P27,P28,P29]
elif nunique==31:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24,P25,P26,P27,P28,P29,P30]
elif nunique==32:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24,P25,P26,P27,P28,P29,P30,P31]
elif nunique==33:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24,P25,P26,P27,P28,P29,P30,P31,P32]
elif nunique==34:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24,P25,P26,P27,P28,P29,P30,P31,P32,P33]
elif nunique==35:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24,P25,P26,P27,P28,P29,P30,P31,P32,P33,P34]
elif nunique==36:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24,P25,P26,P27,P28,P29,P30,P31,P32,P33,P34,P35]
elif nunique==37:
    legend_handles_var=[P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24,P25,P26,P27,P28,P29,P30,P31,P32,P33,P34,P35,P36]


########################
## Combine key elements of plot
## into pandas dataframe
########################

combined_sign_w_percentage_trim=par_cor_sign_var_trim*par_cor_r2_percent_var_trim
parameter_effects_var_trim = {'Parameter': parameter_name_var_trim, 'RI_metric': combined_sign_w_percentage_trim, 'Variable': var_list_trim, 'Position': xvals_var_trim,'Bar_colors':bar_colors_var_trim}

df= pd.DataFrame(parameter_effects_var_trim, columns=['Parameter','RI_metric','Variable','Position','Bar_colors'])


#############################
## Make bars for each
## month, including annual mean
## and seasonal amplitude
#############################


fig,ax=plt.subplots(1)
nmember=len(df.Position)
for ivar in np.arange(nvar):
    if ivar==0:
        first_bar=True
        previous_members=0 
    else:
        previous_members+=nmember_var ## carrying over position index from previous variable
    ipos=xlist[ivar] ## based on x position value, equivalent to df.Position[imember]
    nmember_var=0  ## number of members in this specific variable
    for imember in np.arange(nmember): ## all members of dataframe compared to extract correct number of indicies for variables (ordered)
        if df.Position[imember]==ipos:
            nmember_var+=1
    first_pos=True ## Starts true for each variable
    first_neg=True
    member_index_var=np.arange(nmember_var)+previous_members ## an index of indices
    for ipar in member_index_var: ## including 'previous_members' allow each variable to be accessed from single array
        if (np.logical_and((first_pos),(df.RI_metric[ipar]>0.0))): ## First bar of variable in positive direction
            if first_bar:
                ax=plt.bar(df.Position[ipar],df.RI_metric[ipar],color=df.Bar_colors[ipar])
                first_bar=False
                first_pos=False
            else:
                plt.bar(df.Position[ipar],df.RI_metric[ipar],color=df.Bar_colors[ipar])
                first_pos=False
        elif (np.logical_and((first_neg),(df.RI_metric[ipar]<0.0))): ## First bar of variable in negative direction
            if first_bar:
                ax=plt.bar(df.Position[ipar],df.RI_metric[ipar],color=df.Bar_colors[ipar])
                first_bar=False
                first_neg=False
            else:
                plt.bar(df.Position[ipar],df.RI_metric[ipar],color=df.Bar_colors[ipar])
                first_neg=False
        elif (df.RI_metric[ipar]>0): ## Another bar in positive direction, so start in correct place
            bottom_pos=0.0
            for ipar_sub in member_index_var:
                if (np.logical_and((df.RI_metric[ipar_sub]>0.0),(ipar_sub<ipar))):
                    bottom_pos+=df.RI_metric[ipar_sub]
            plt.bar(df.Position[ipar],df.RI_metric[ipar],bottom=bottom_pos,color=df.Bar_colors[ipar])
        else: ## Another bar in the negative direction, so start in correct place, as above
            bottom_neg=0.0
            for ipar_sub in member_index_var:
                if (np.logical_and((df.RI_metric[ipar_sub]<0.0),(ipar_sub<ipar))):
                    bottom_neg+=df.RI_metric[ipar_sub]
            plt.bar(df.Position[ipar],df.RI_metric[ipar],bottom=bottom_neg,color=df.Bar_colors[ipar])



#################################
## Adding a legend
## Requires patches for each parameter
#################################

ylim_default= plt.ylim()
xlim_default= plt.xlim()
plt.xticks(xlist,labels=var_name_titles)
plt.xlabel('')
plt.ylabel('Relative importance metric (%)')
plt.hlines(0,xmin=xlist[0]-1,xmax=xlist[-1]+1,linewidth=0.7)
plt.xlim(xlist[0]-1,xlist[-1]+1)

if nunique<10:
    leg_size=8
    leg_col=3
    plt.ylim(ylim_default[0],ylim_default[1]+20)
elif nunique<15:
    leg_size=8
    leg_col=3
    plt.ylim(ylim_default[0],ylim_default[1]+35)
elif nunique<20:
    leg_size=8
    leg_col=3
    plt.ylim(ylim_default[0],ylim_default[1]+70)
elif nunique<25:
    leg_size=7
    leg_col=3
    plt.ylim(ylim_default[0],ylim_default[1]+80)
elif nunique<30:
    leg_size=7
    leg_col=3
    plt.ylim(ylim_default[0],ylim_default[1]+120)
elif nunique<35:
    leg_size=7
    leg_col=3
    plt.ylim(ylim_default[0],ylim_default[1]+140)
elif nunique<40:
    leg_size=7
    leg_col=3
    plt.ylim(ylim_default[0],ylim_default[1]+200)
else:
    leg_size=7
    leg_col=3
    plt.ylim(ylim_default[0],ylim_default[1]+250)

leg_var=plt.legend(handles=legend_handles_var,loc='upper center',prop={'size': leg_size},ncol=leg_col)


outfile='/gws/nopw/j04/acure/lregayre/parameter_contributions/sample/ERF/parameter_relative_importance_metric_ERF_and_components_ann_no_scav_diam_partial_correlations_no_carbonaceous_'+str(np.int(lower_limit))+'pc_cutoff.pdf'
plt.savefig(outfile, dpi=300, bbox_inches='tight') 
