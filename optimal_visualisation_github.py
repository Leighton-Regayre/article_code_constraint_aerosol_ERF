##################################################################################################################
## Plot percentage reduction in 90% CI as variables are progressively added to constraint
## to find optimal constraint and evaluate effect of additional constraints
##################################################################################################################

from __future__ import division ## This allows all division to return floats
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from math import *
import os
from netCDF4 import Dataset as nd
import pylab
from matplotlib.backends.backend_pdf import PdfPages
import cf_units
import sys
import matplotlib.lines as mlines
import scipy
import pandas as pd
import time
import seaborn as sns
from pathlib import Path
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


#################################
# progressively add
#################################

## read file created during constraint process

constraint_summary_file_bottom_up='/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/finesse_reverse/all/range_reduction_number_of_constraints.txt'
constraint_array=np.loadtxt(constraint_summary_file_bottom_up,delimiter=',')
constraint_array.shape


## These indices of constraint variables are in order constraint variables were progressively added, but relate to order in which variables were emulated

names_file_bottom_up='/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/finesse_reverse/all/all_sample_indices_py_index_FINESSE_w_uncertainty_0pc_error_constrained_ERF_-1.4_to_0.04_with_51.5pc_reduction_from_225_constraints_INDICES_OF_RETAINED_w_o_Hd_regional.txt'
bottom_up_indices=np.loadtxt(names_file_bottom_up)
x_index=np.arange(constraint_array.shape[0])


greek_letters_lower=[chr(code) for code in range(945,970)]
TAU=greek_letters_lower[19]
greek_letters_upper=[chr(code) for code in range(913,938)]
DELTA=greek_letters_upper[3]

##############################################################################################
## Read in a list of indices and associated labels, for variables considered to be consistent w/  Nd
##############################################################################################

icon_consistent=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,102,103,104,106,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,192,197,198,199,201,202,203,204,205,206,207,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,229,230,231,232,233,234,235,237,238,281,282,283,284,285,286,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,395,396,397,398,399,400,401,402,403,404,405,406,407]

constraint_variable_text_list=['seasonal_amplitude_Nd_hemispheric_diff','jan_Hd','feb_Hd','mar_Hd','apr_Hd','may_Hd','jun_Hd','jul_Hd','aug_Hd','sep_Hd','oct_Hd','nov_Hd','dec_Hd','ann_Hd','jul_transect_Nd_NAtlantic','jul_transect_LWP_Nd_NAtlantic','jul_transect_Re_Nd_NAtlantic','jul_transect_Re_NAtlantic','seasonal_amplitude_Fsw_NAtlantic','jan_Fsw_NAtlantic','feb_Fsw_NAtlantic','mar_Fsw_NAtlantic','apr_Fsw_NAtlantic','may_Fsw_NAtlantic','jun_Fsw_NAtlantic','jul_Fsw_NAtlantic','aug_Fsw_NAtlantic','sep_Fsw_NAtlantic','oct_Fsw_NAtlantic','nov_Fsw_NAtlantic','dec_Fsw_NAtlantic','ann_Fsw_NAtlantic','seasonal_amplitude_Nd_NAtlantic','jan_Nd_NAtlantic','feb_Nd_NAtlantic','mar_Nd_NAtlantic','apr_Nd_NAtlantic','may_Nd_NAtlantic','jun_Nd_NAtlantic','jul_Nd_NAtlantic','aug_Nd_NAtlantic','sep_Nd_NAtlantic','oct_Nd_NAtlantic','nov_Nd_NAtlantic','dec_Nd_NAtlantic','ann_Nd_NAtlantic','seasonal_amplitude_Fc_NAtlantic','jan_Fc_NAtlantic','feb_Fc_NAtlantic','mar_Fc_NAtlantic','apr_Fc_NAtlantic','may_Fc_NAtlantic','jun_Fc_NAtlantic','jul_Fc_NAtlantic','aug_Fc_NAtlantic','sep_Fc_NAtlantic','oct_Fc_NAtlantic','nov_Fc_NAtlantic','dec_Fc_NAtlantic','ann_Fc_NAtlantic','seasonal_amplitude_Re_NAtlantic','jan_Re_NAtlantic','feb_Re_NAtlantic','mar_Re_NAtlantic','apr_Re_NAtlantic','may_Re_NAtlantic','jun_Re_NAtlantic','jul_Re_NAtlantic','aug_Re_NAtlantic','sep_Re_NAtlantic','oct_Re_NAtlantic','nov_Re_NAtlantic','dec_Re_NAtlantic','ann_Re_NAtlantic','seasonal_amplitude_LWP_NAtlantic',"jan_LWP_NAtlantic","feb_LWP_NAtlantic","mar_LWP_NAtlantic","apr_LWP_NAtlantic","may_LWP_NAtlantic","jun_LWP_NAtlantic","jul_LWP_NAtlantic","aug_LWP_NAtlantic","sep_LWP_NAtlantic","oct_LWP_NAtlantic","nov_LWP_NAtlantic","dec_LWP_NAtlantic","ann_LWP_NAtlantic",'seasonal_amplitude_Tau_NAtlantic',"jan_Tau_NAtlantic","feb_Tau_NAtlantic","mar_Tau_NAtlantic","apr_Tau_NAtlantic","may_Tau_NAtlantic","jun_Tau_NAtlantic","jul_Tau_NAtlantic","aug_Tau_NAtlantic","sep_Tau_NAtlantic","oct_Tau_NAtlantic","nov_Tau_NAtlantic","dec_Tau_NAtlantic","ann_Tau_NAtlantic",'jul_transect_Nd_NPacific','jul_transect_Cf_NPacific','jul_transect_LWP_NPacific','jul_transect_Cf_Nd_Pacific','jul_transect_LWP_Nd_NPacific','jul_transect_Re_Nd_NPacific','seasonal_amplitude_Fsw_NPacific','jan_Fsw_NPacific','feb_Fsw_NPacific','mar_Fsw_NPacific','apr_Fsw_NPacific','may_Fsw_NPacific','jun_Fsw_NPacific','jul_Fsw_NPacific','aug_Fsw_NPacific','sep_Fsw_NPacific','oct_Fsw_NPacific','nov_Fsw_NPacific','dec_Fsw_NPacific','ann_Fsw_NPacific',"seasonal_amplitude_Nd_NPacific","jan_Nd_NPacific","feb_Nd_NPacific","mar_Nd_NPacific","apr_Nd_NPacific","may_Nd_NPacific","jun_Nd_NPacific","jul_Nd_NPacific","aug_Nd_NPacific","sep_Nd_NPacific","oct_Nd_NPacific","nov_Nd_NPacific","dec_Nd_NPacific","ann_Nd_NPacific",'seasonal_amplitude_Fc_NPacific','jan_Fc_NPacific','feb_Fc_NPacific','mar_Fc_NPacific','apr_Fc_NPacific','may_Fc_NPacific','jun_Fc_NPacific','jul_Fc_NPacific','aug_Fc_NPacific','sep_Fc_NPacific','oct_Fc_NPacific','nov_Fc_NPacific','dec_Fc_NPacific','ann_Fc_NPacific',"seasonal_amplitude_Re_NPacific","jan_Re_NPacific","feb_Re_NPacific","mar_Re_NPacific","apr_Re_NPacific","may_Re_NPacific","jun_Re_NPacific","jul_Re_NPacific","aug_Re_NPacific","sep_Re_NPacific","oct_Re_NPacific","nov_Re_NPacific","dec_Re_NPacific","ann_Re_NPacific","seasonal_amplitude_LWP_NPacific","jan_LWP_NPacific","feb_LWP_NPacific","mar_LWP_NPacific","apr_LWP_NPacific","may_LWP_NPacific","jun_LWP_NPacific","jul_LWP_NPacific","aug_LWP_NPacific","sep_LWP_NPacific","oct_LWP_NPacific","nov_LWP_NPacific","dec_LWP_NPacific","ann_LWP_NPacific","seasonal_amplitude_Tau_NPacific","jan_Tau_NPacific","feb_Tau_NPacific","mar_Tau_NPacific","apr_Tau_NPacific","may_Tau_NPacific","jun_Tau_NPacific","jul_Tau_NPacific","aug_Tau_NPacific","sep_Tau_NPacific","oct_Tau_NPacific","nov_Tau_NPacific","dec_Tau_NPacific","ann_Tau_NPacific",'nov_transect_Nd_SAtlantic','nov_transect_LWP_Nd_SAtlantic','nov_transect_Fc_SAtlantic','nov_transect_AI_SAtlantic','nov_transect_Nd_AI_SAtlantic','seasonal_amplitude_Fsw_SAtlantic','jan_Fsw_SAtlantic','feb_Fsw_SAtlantic','mar_Fsw_SAtlantic','apr_Fsw_SAtlantic','may_Fsw_SAtlantic','jun_Fsw_SAtlantic','jul_Fsw_SAtlantic','aug_Fsw_SAtlantic','sep_Fsw_SAtlantic','oct_Fsw_SAtlantic','nov_Fsw_SAtlantic','dec_Fsw_SAtlantic','ann_Fsw_SAtlantic',"seasonal_amplitude_Nd_SAtlantic",'jan_Nd_SAtlantic','feb_Nd_SAtlantic','mar_Nd_SAtlantic','apr_Nd_SAtlantic','may_Nd_SAtlantic','jun_Nd_SAtlantic','jul_Nd_SAtlantic','aug_Nd_SAtlantic','sep_Nd_SAtlantic','oct_Nd_SAtlantic','nov_Nd_SAtlantic','dec_Nd_SAtlantic','ann_Nd_SAtlantic','seasonal_amplitude_Fc_SAtlantic','jan_Fc_SAtlantic','feb_Fc_SAtlantic','mar_Fc_SAtlantic','apr_Fc_SAtlantic','may_Fc_SAtlantic','jun_Fc_SAtlantic','jul_Fc_SAtlantic','aug_Fc_SAtlantic','sep_Fc_SAtlantic','oct_Fc_SAtlantic','nov_Fc_SAtlantic','dec_Fc_SAtlantic','ann_Fc_SAtlantic',"seasonal_amplitude_Re_SAtlantic","jan_Re_SAtlantic","feb_Re_SAtlantic","mar_Re_SAtlantic","apr_Re_SAtlantic","may_Re_SAtlantic","jun_Re_SAtlantic","jul_Re_SAtlantic","aug_Re_SAtlantic","sep_Re_SAtlantic","oct_Re_SAtlantic","nov_Re_SAtlantic","dec_Re_SAtlantic","ann_Re_SAtlantic","seasonal_amplitude_LWP_SAtlantic","jan_LWP_SAtlantic","feb_LWP_SAtlantic","mar_LWP_SAtlantic","apr_LWP_SAtlantic","may_LWP_SAtlantic","jun_LWP_SAtlantic","jul_LWP_SAtlantic","aug_LWP_SAtlantic","sep_LWP_SAtlantic","oct_LWP_SAtlantic","nov_LWP_SAtlantic","dec_LWP_SAtlantic","ann_LWP_SAtlantic","seasonal_amplitude_Tau_SAtlantic","jan_Tau_SAtlantic","feb_Tau_SAtlantic","mar_Tau_SAtlantic","apr_Tau_SAtlantic","may_Tau_SAtlantic","jun_Tau_SAtlantic","jul_Tau_SAtlantic","aug_Tau_SAtlantic","sep_Tau_SAtlantic","oct_Tau_SAtlantic","nov_Tau_SAtlantic","dec_Tau_SAtlantic","ann_Tau_SAtlantic",'nov_transect_Nd_SPacific','seasonal_amplitude_Fsw_SPacific','jan_Fsw_SPacific','feb_Fsw_SPacific','mar_Fsw_SPacific','apr_Fsw_SPacific','may_Fsw_SPacific','jun_Fsw_SPacific','jul_Fsw_SPacific','aug_Fsw_SPacific','sep_Fsw_SPacific','oct_Fsw_SPacific','nov_Fsw_SPacific','dec_Fsw_SPacific','ann_Fsw_SPacific','seasonal_amplitude_Nd_SPacific','jan_Nd_SPacific','feb_Nd_SPacific','mar_Nd_SPacific','apr_Nd_SPacific','may_Nd_SPacific','jun_Nd_SPacific','jul_Nd_SPacific','aug_Nd_SPacific','sep_Nd_SPacific','oct_Nd_SPacific','nov_Nd_SPacific','dec_Nd_SPacific','ann_Nd_SPacific','seasonal_amplitude_Fc_SPacific','jan_Fc_SPacific','feb_Fc_SPacific','mar_Fc_SPacific','apr_Fc_SPacific','may_Fc_SPacific','jun_Fc_SPacific','jul_Fc_SPacific','aug_Fc_SPacific','sep_Fc_SPacific','oct_Fc_SPacific','nov_Fc_SPacific','dec_Fc_SPacific','ann_Fc_SPacific',"seasonal_amplitude_Re_SPacific","jan_Re_SPacific","feb_Re_SPacific","mar_Re_SPacific","apr_Re_SPacific","may_Re_SPacific","jun_Re_SPacific","jul_Re_SPacific","aug_Re_SPacific","sep_Re_SPacific","oct_Re_SPacific","nov_Re_SPacific","dec_Re_SPacific","ann_Re_SPacific","seasonal_amplitude_LWP_SPacific","jan_LWP_SPacific","feb_LWP_SPacific","mar_LWP_SPacific","apr_LWP_SPacific","may_LWP_SPacific","jun_LWP_SPacific","jul_LWP_SPacific","aug_LWP_SPacific","sep_LWP_SPacific","oct_LWP_SPacific","nov_LWP_SPacific","dec_LWP_SPacific","ann_LWP_SPacific","seasonal_amplitude_Tau_SPacific","jan_Tau_SPacific","feb_Tau_SPacific","mar_Tau_SPacific","apr_Tau_SPacific","may_Tau_SPacific","jun_Tau_SPacific","jul_Tau_SPacific","aug_Tau_SPacific","sep_Tau_SPacific","oct_Tau_SPacific","nov_Tau_SPacific","dec_Tau_SPacific","ann_Tau_SPacific",'seasonal_amplitude_Fsw_SOcean','jan_Fsw_SOcean','feb_Fsw_SOcean','mar_Fsw_SOcean','apr_Fsw_SOcean','may_Fsw_SOcean','jun_Fsw_SOcean','jul_Fsw_SOcean','aug_Fsw_SOcean','sep_Fsw_SOcean','oct_Fsw_SOcean','nov_Fsw_SOcean','dec_Fsw_SOcean','ann_Fsw_SOcean','seasonal_amplitude_Nd_SOcean','jan_Nd_SOcean','feb_Nd_SOcean','mar_Nd_SOcean','apr_Nd_SOcean','may_Nd_SOcean','jun_Nd_SOcean','jul_Nd_SOcean','aug_Nd_SOcean','sep_Nd_SOcean','oct_Nd_SOcean','nov_Nd_SOcean','dec_Nd_SOcean','ann_Nd_SOcean',"seasonal_amplitude_Fc_SOcean",'jan_Fc_SOcean','feb_Fc_SOcean','mar_Fc_SOcean','apr_Fc_SOcean','may_Fc_SOcean','jun_Fc_SOcean','jul_Fc_SOcean','aug_Fc_SOcean','sep_Fc_SOcean','oct_Fc_SOcean','nov_Fc_SOcean','dec_Fc_SOcean','ann_Fc_SOcean',"seasonal_amplitude_Re_SOcean","jan_Re_SOcean","feb_Re_SOcean","mar_Re_SOcean","apr_Re_SOcean","may_Re_SOcean","jun_Re_SOcean","jul_Re_SOcean","aug_Re_SOcean","sep_Re_SOcean","oct_Re_SOcean","nov_Re_SOcean","dec_Re_SOcean","ann_Re_SOcean","seasonal_amplitude_LWP_SOcean","jan_LWP_SOcean","feb_LWP_SOcean","mar_LWP_SOcean","apr_LWP_SOcean","may_LWP_SOcean","jun_LWP_SOcean","jul_LWP_SOcean","aug_LWP_SOcean","sep_LWP_SOcean","oct_LWP_SOcean","nov_LWP_SOcean","dec_LWP_SOcean","ann_LWP_SOcean","seasonal_amplitude_Tau_SOcean","jan_Tau_SOcean","feb_Tau_SOcean","mar_Tau_SOcean","apr_Tau_SOcean","may_Tau_SOcean","jun_Tau_SOcean","jul_Tau_SOcean","aug_Tau_SOcean","sep_Tau_SOcean","oct_Tau_SOcean","nov_Tau_SOcean","dec_Tau_SOcean","ann_Tau_SOcean"]


################################
## Extract constraint variable names
## in order added
################################

x_labels_bottom_up=[]
for icon in np.arange(len(bottom_up_indices)):
    x_labels_bottom_up.append(constraint_variable_text_list[int(bottom_up_indices[icon])])


################
## Plot
################


bottom_up_color='goldenrod'
near_perfect_color='mediumpurple'
perfect_color='cornflowerblue'
inconsistent_color='k'

bottom_up_line= mpatches.Patch(color=bottom_up_color, label='Progressively added to constraint')
perfect_line= mpatches.Patch(color=perfect_color, label='No structural inadequacies')
near_perfect_line= mpatches.Patch(color=near_perfect_color, label='Few structural inadequacies')
inconsistent_line= mpatches.Patch(color=inconsistent_color,label='Added inconsistent constraint variables')


perfect_scale_factor=1.08
perfect_val=np.max(constraint_array[:,2])*perfect_scale_factor
ncon=len(constraint_array[:,2])

perfect_yvals=np.repeat(perfect_val,ncon)
perfect_yvals[0:12]=constraint_array[0:12,2]

split_pos=60
nsplit=split_pos-10
split_distance=perfect_val-perfect_yvals[9]

increment=split_distance/nsplit
for ival in np.arange(nsplit-2):
    perfect_yvals[ival+12]=perfect_yvals[ival-1+12]+increment

ymax=np.max(perfect_yvals)


perfect_smooth_scale_factor=1.08
perfect_smooth_val=perfect_val
perfect_smooth_yvals=np.repeat(perfect_smooth_val,ncon)
perfect_smooth_yvals[0:10]=constraint_array[0:10,2]

split_pos=60
nsplit=split_pos-10
split_distance=perfect_smooth_val-perfect_smooth_yvals[9]
increment_vals=split_distance - np.flip(np.geomspace(0.1,split_distance,nsplit))[1:nsplit+1]

for ival in np.arange(nsplit-1):
    perfect_smooth_yvals[ival-1+10]=perfect_smooth_yvals[9]+increment_vals[ival]


for ival in np.arange(ncon-split_pos+2)+split_pos-2:
    perfect_smooth_yvals[ival]=perfect_smooth_yvals[split_pos-3]


ymax=np.max(perfect_smooth_yvals)

##################################
## Plot 'near perfect', w/ smooth transition
##################################

near_perfect_smooth_val=ymax-3.
near_perfect_smooth_yvals=np.repeat(near_perfect_smooth_val,ncon)
near_perfect_smooth_yvals[0:10]=constraint_array[0:10,2]

split_pos=42
nsplit=split_pos-10
split_distance=near_perfect_smooth_val-near_perfect_smooth_yvals[9]
increment_vals=split_distance - np.flip(np.geomspace(0.1,split_distance,nsplit))[1:nsplit+1]

for ival in np.arange(nsplit-1):
    near_perfect_smooth_yvals[ival-1+10]=near_perfect_smooth_yvals[9]+increment_vals[ival]


near_perfect_smooth_val=near_perfect_smooth_yvals[split_pos-1]
min_val=65.5
nsplit=ncon-split_pos
split_distance=near_perfect_smooth_val-min_val
nsub=np.abs(np.random.normal((near_perfect_smooth_val-min_val)/nsplit,(near_perfect_smooth_val-min_val)/nsplit,nsplit))
np.sum(nsub)

for ival in np.arange(nsplit+2)+split_pos-2:
    near_perfect_smooth_yvals[ival]=near_perfect_smooth_yvals[ival-1]-nsub[ival-split_pos]



######################################
## Where perfect<near_perfect
## set near_perfect=perfect
######################################

for ival in np.arange(100):
    if (perfect_smooth_yvals[ival] < near_perfect_smooth_yvals[ival]):
        perfect_smooth_yvals[ival] = near_perfect_smooth_yvals[ival]



##############################################################
## Plot w/ indication of effect of inconsistent constraint variables
## reduced limit of variables shown as arrows
## only plot first XX to 'new_ylim'
##############################################################

new_lim=125
ymin=57
y_text_inc=1.2

font_size=12
mp.rcParams['axes.labelsize'] = 12
mp.rcParams['font.size'] = 12
mp.rcParams['axes.linewidth'] = 1
mp.rcParams['legend.fontsize'] = 12
mp.rcParams['xtick.labelsize'] = 12
mp.rcParams['ytick.labelsize'] = 12
mp.rcParams['lines.markersize']=5.0

fontsize_legend=11
fontsize_axis_label=13
fontsize_inline=12


fig,ax=plt.subplots(constrained_layout=True)
ax.scatter(x_index[0:new_lim],perfect_smooth_yvals[0:new_lim],color=perfect_color,marker='.')
ax.scatter(x_index[0:new_lim],near_perfect_smooth_yvals[0:new_lim],color=near_perfect_color,marker='.')
ax.scatter(x_index[0:new_lim],constraint_array[0:new_lim,2],color=bottom_up_color,marker='.')
plt.ylabel('Percentage reduction in '+DELTA+'F$_{aci}$ 90% CI',fontsize=fontsize_axis_label)
plt.xlabel('Number of constraint variables used',fontsize=fontsize_axis_label)
ax.set_ylim(ymin,int(ymax)+2)
ax.set_xlim(-3,new_lim+23)
legend1 = plt.legend(handles=[perfect_line,near_perfect_line,bottom_up_line,inconsistent_line], loc=8,prop={'size': fontsize_legend})
plt.gca().add_artist(legend1)
plt.arrow(1,ymin+2,0,-1.3,head_length=0.5, fc=bottom_up_color,head_width=3.5)
plt.text(2,ymin+y_text_inc,'44%',fontsize=fontsize_inline)
plt.arrow(new_lim+3,ymin+4,0,-1.5, head_length=0.5, fc=bottom_up_color,head_width=3.5)
plt.text(new_lim+4.5,ymin+2.1+y_text_inc,'52%',fontsize=fontsize_inline)
plt.arrow(new_lim+10,ymin+2,0,-1.5, head_length=0.5, fc=inconsistent_color,head_width=3.5)
plt.text(new_lim+11.5,ymin+y_text_inc,'37%',fontsize=fontsize_inline)
plt.savefig('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/finesse/all_points_progressively_add_w_hypothetical_best_case_flipped_no_shading_star_inconsistent_effect.pdf')
