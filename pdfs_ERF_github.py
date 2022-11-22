#!/usr/bin/env python
"""

Makes 1 million variant pdf from eulator output for aerosol ERF and components
Optimal constraint overlaid

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


greek_letters_lower=[chr(code) for code in range(945,970)]
TAU=greek_letters_lower[19]
greek_letters_upper=[chr(code) for code in range(913,938)]
DELTA=greek_letters_upper[3]
import seaborn as sns

## Requires python 3.7 for sns

iregion=0
region_list=['global']
region_titles=['Global mean']
region= region_list[iregion]
region_title= region_titles[iregion]

nppe=221
nsample=1000000
var_list_PPE = ['ERF','ACI','ARI']
var_list_titles=[DELTA+' F$_{aer}$',DELTA+' F$_{aci}$',DELTA+' $F_{ari}$']
units_list=['W m$^{-2}$','W m$^{-2}$','W m$^{-2}$']

nregion=len(region_list)
nvar=len(var_list_PPE)
dir_in= '/gws/nopw/j04/acure/lregayre/data_post_emulation/'


############################
## Set up read of constrained data
############################

constraint_dir='/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/finesse_reverse/all/'
constraint_file='all_sample_indices_py_index_FINESSE_w_uncertainty_0pc_error_constrained_ERF_-1.26_to_-0.13_with_69.44pc_reduction_from_13_constraints_w_o_Hd_regional.dat'
outdir='/gws/nopw/j04/acure/lregayre/pdfs/ERF/FINESSE/'
file_out_pdf='/global_mean_aerosol_ERF_ACI_ARI_constraint_from_All_Regions_FINESSE_BOTTOM_UP_peak_w_B2020.pdf'
constraint_text=' All regions finesse bottom-up peak'


########################################
## Plot features
########################################

font_size=15
mp.rcParams['axes.labelsize'] = font_size*0.8
mp.rcParams['font.size'] = font_size*1.2
mp.rcParams['axes.linewidth'] = font_size*0.05
mp.rcParams['axes.titlesize'] = font_size*1.4
mp.rcParams['legend.fontsize'] = font_size*1.4
mp.rcParams['xtick.labelsize'] = font_size*1
mp.rcParams['ytick.labelsize'] = font_size*1.4



#####################################
## Read in PPE data for each variable
## all 1 million
#####################################

PPE_array=np.zeros((nvar,nsample))

for ivar in np.arange(nvar):
    file_in=dir_in+'/emulated_mean_values_'+var_list_PPE[ivar]+'_1000000_w_o_carb.dat'
    PPE_array[ivar,:]=np.loadtxt(file_in)



##########################
## Calculate percentiles
##########################

unconstrained_box_values=np.zeros((nvar,5))

for ivar in np.arange(nvar):
    unconstrained_box_values[ivar,0]=np.percentile(PPE_array[ivar,:],5)
    unconstrained_box_values[ivar,1]=np.percentile(PPE_array[ivar,:],25)
    unconstrained_box_values[ivar,2]=np.percentile(PPE_array[ivar,:],50)
    unconstrained_box_values[ivar,3]=np.percentile(PPE_array[ivar,:],75)
    unconstrained_box_values[ivar,4]=np.percentile(PPE_array[ivar,:],95)


## Bellouin2020
B2020_box=np.zeros((nvar,2))
B2020_box[0,0]=-3.15
B2020_box[0,1]=-0.35
B2020_box[1,0]=-2.65
B2020_box[1,1]=-0.07
B2020_box[2,0]=-0.71
B2020_box[2,1]=-0.14


###############################
## Constrain PPE using indices
## of 1M from optimal constraint
###############################

integers_raw=np.loadtxt(constraint_dir+constraint_file)
nsample=len(integers_raw)
PPE_ERF_constrained=np.zeros((nsample))
PPE_ACI_constrained=np.zeros((nsample))
PPE_ARI_constrained=np.zeros((nsample))
for isample in np.arange(nsample):
    PPE_ERF_constrained[isample]=PPE_array[0,int(integers_raw[isample])]
    PPE_ACI_constrained[isample]=PPE_array[1,int(integers_raw[isample])]
    PPE_ARI_constrained[isample]=PPE_array[2,int(integers_raw[isample])]


box_values=np.zeros((nvar,5))
box_values[0,0]=np.percentile(PPE_ERF_constrained,5)
box_values[0,1]=np.percentile(PPE_ERF_constrained,25)
box_values[0,2]=np.percentile(PPE_ERF_constrained,50)
box_values[0,3]=np.percentile(PPE_ERF_constrained,75)
box_values[0,4]=np.percentile(PPE_ERF_constrained,95)
box_values[1,0]=np.percentile(PPE_ACI_constrained,5)
box_values[1,1]=np.percentile(PPE_ACI_constrained,25)
box_values[1,2]=np.percentile(PPE_ACI_constrained,50)
box_values[1,3]=np.percentile(PPE_ACI_constrained,75)
box_values[1,4]=np.percentile(PPE_ACI_constrained,95)
box_values[2,0]=np.percentile(PPE_ARI_constrained,5)
box_values[2,1]=np.percentile(PPE_ARI_constrained,25)
box_values[2,2]=np.percentile(PPE_ARI_constrained,50)
box_values[2,3]=np.percentile(PPE_ARI_constrained,75)
box_values[2,4]=np.percentile(PPE_ARI_constrained,95)


print('box_values for constraint:'+constraint_text)
print(box_values)


outdir='/gws/nopw/j04/acure/lregayre/pdfs/ERF/'


##########
## Plot
##########


fig= plt.figure()

fig.set_figheight(15)
fig.set_figwidth(7)
fontsize_labels=12
fontsize_sublabels=9
fontsize_axes= 18
fontsize_legend= 18
linewidth_pdf=3
linewidth_box=1.0 #1.2

ivar=0
ax_a = plt.subplot(1,3,1, adjustable='box')
ax_a.set_xlabel(var_list_titles[ivar]+' / W m$^{-2}$',fontsize=fontsize_axes)
t1=sns.distplot(PPE_array[ivar,:],color='r',hist_kws={"alpha":0.2})
t2=sns.distplot(PPE_ERF_constrained,color='Orange',hist_kws={"alpha":0.2})

ylim_starter= ax_a.get_ylim()[1]
ylim_upper=ylim_starter*1.6
ylim_top_space=ylim_upper - ylim_starter
ax_a.set_ylim((0,ylim_upper))

yval_box_middle_B2020=ylim_starter+ylim_top_space/4.*2.5
yval_box_upper_B2020= yval_box_middle_B2020+(ylim_top_space/3.)*0.16
yval_box_lower_B2020= yval_box_middle_B2020-(ylim_top_space/3.)*0.13
yvals_line_B2020=np.repeat(yval_box_middle_B2020, 30)
xvals_line_B2020=np.linspace(B2020_box[ivar,0],B2020_box[ivar,1],30)
yvals_vertical_B2020=np.linspace(yval_box_lower_B2020,yval_box_upper_B2020,30)
xvals_left_line_B2020=np.repeat(B2020_box[ivar,0],30)
xvals_right_line_B2020=np.repeat(B2020_box[ivar,1],30)

plt.plot(xvals_line_B2020,yvals_line_B2020,'--',color='indigo',linewidth=linewidth_box)
plt.plot(xvals_left_line_B2020,yvals_vertical_B2020,'-',color='indigo',linewidth=linewidth_box)
plt.plot(xvals_right_line_B2020,yvals_vertical_B2020,'-',color='indigo',linewidth=linewidth_box)


## Constrained
yval_box_middle=ylim_starter+ylim_top_space/4.*0.37
yval_box_upper= yval_box_middle+(ylim_top_space/3.*2.)*0.147
yval_box_lower= yval_box_middle-(ylim_top_space/3.*2.)*0.14
yvals_line=np.repeat(yval_box_middle, 30)
xvals_line=np.linspace(box_values[ivar,0],box_values[ivar,4],30)
box_line_left_x=np.repeat(box_values[ivar,1], 30)
box_line_right_x=np.repeat(box_values[ivar,3], 30)
box_xvals= np.linspace(box_values[ivar,1],box_values[ivar,3],30)
box_line_top_y=np.repeat(yval_box_upper,30)
box_line_bottom_y=np.repeat(yval_box_lower,30)
yvals_vertical=np.linspace(yval_box_lower,yval_box_upper,30)
xvals_left_line=np.repeat(box_values[ivar,0],30)
xvals_right_line=np.repeat(box_values[ivar,4],30)
xvals_middle_line=np.repeat(box_values[ivar,2],30)
plt.plot(xvals_line,yvals_line,'--',color='Orange',linewidth=linewidth_box)
plt.plot(xvals_left_line,yvals_vertical,'-',color='Orange',linewidth=linewidth_box)
plt.plot(xvals_right_line,yvals_vertical,'-',color='Orange',linewidth=linewidth_box)
box_patch=plt.Rectangle((box_values[ivar,1],yval_box_lower),box_values[ivar,3]-box_values[ivar,1],yval_box_upper-yval_box_lower,color='Navajowhite')
ax_a.add_patch(box_patch)
plt.plot(xvals_middle_line,yvals_vertical,'-',color='Orange',linewidth=linewidth_box+0.2)
plt.plot(box_line_left_x,yvals_vertical,'-',color='Orange',linewidth=linewidth_box)
plt.plot(box_line_right_x,yvals_vertical,'-',color='Orange',linewidth=linewidth_box)
plt.plot(box_xvals,box_line_top_y,'-',color='Orange',linewidth=linewidth_box)
plt.plot(box_xvals,box_line_bottom_y,'-',color='Orange',linewidth=linewidth_box)


## Unconstrained
unconstrained_yval_box_middle=ylim_starter+ylim_top_space/4.*1.45
unconstrained_yval_box_upper= unconstrained_yval_box_middle+(ylim_top_space/3.*2.)*0.151
unconstrained_yval_box_lower= unconstrained_yval_box_middle-(ylim_top_space/3.*2.)*0.149
unconstrained_yvals_line=np.repeat(unconstrained_yval_box_middle, 30)
unconstrained_xvals_line=np.linspace(unconstrained_box_values[ivar,0],unconstrained_box_values[ivar,4],30)
unconstrained_box_line_left_x=np.repeat(unconstrained_box_values[ivar,1], 30)
unconstrained_box_line_right_x=np.repeat(unconstrained_box_values[ivar,3], 30)
unconstrained_box_xvals= np.linspace(unconstrained_box_values[ivar,1],unconstrained_box_values[ivar,3],30)
unconstrained_box_line_top_y=np.repeat(unconstrained_yval_box_upper,30)
unconstrained_box_line_bottom_y=np.repeat(unconstrained_yval_box_lower,30)
unconstrained_yvals_vertical=np.linspace(unconstrained_yval_box_lower,unconstrained_yval_box_upper,30)
unconstrained_xvals_left_line=np.repeat(unconstrained_box_values[ivar,0],30)
unconstrained_xvals_right_line=np.repeat(unconstrained_box_values[ivar,4],30)
unconstrained_xvals_middle_line=np.repeat(unconstrained_box_values[ivar,2],30)
plt.plot(unconstrained_xvals_line,unconstrained_yvals_line,'--',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_xvals_left_line,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_xvals_right_line,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box)
unconstrained_box_patch=plt.Rectangle((unconstrained_box_values[ivar,1],unconstrained_yval_box_lower),unconstrained_box_values[ivar,3]-unconstrained_box_values[ivar,1],unconstrained_yval_box_upper-unconstrained_yval_box_lower,color='lightcoral')
ax_a.add_patch(unconstrained_box_patch)
plt.plot(unconstrained_xvals_middle_line,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box+0.2)
plt.plot(unconstrained_box_line_left_x,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_box_line_right_x,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_box_xvals,unconstrained_box_line_top_y,'-',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_box_xvals,unconstrained_box_line_bottom_y,'-',color='r',linewidth=linewidth_box)

ax_a.set_xticks((-3,-2,-1,0.0,1.0,2.0))
ax_a.set_xticklabels(('-3','-2','-1','0','1','2'),fontsize=fontsize_axes)
ax_a.set_xlim(-3.5,3.2)
ax_a.set_xlabel(var_list_titles[ivar]+' / '+units_list[ivar],size=14)
ax_a.set_ylabel('Density', size= 14)
ax_a.annotate('a) ', xy=(-2.8,yval_box_upper_B2020*1.06),  xycoords='data',horizontalalignment='center', verticalalignment='center',fontsize=fontsize_labels)
ax_a.annotate('Bellouin et al.', xy=(B2020_box[ivar,1]+0.2,yval_box_upper_B2020*1.05),  xycoords='data',horizontalalignment='left', verticalalignment='center',fontsize=fontsize_sublabels)
ax_a.annotate('(2020)', xy=(B2020_box[ivar,1]+0.2,yval_box_lower_B2020),  xycoords='data',horizontalalignment='left', verticalalignment='center',fontsize=fontsize_sublabels)

ax_a.annotate('Constrained', xy=(box_values[ivar,4]+0.2,(yval_box_lower+yval_box_middle)/2.),  xycoords='data',horizontalalignment='left', verticalalignment='center',fontsize=fontsize_sublabels)
ax_a.annotate('Original', xy=(unconstrained_box_values[ivar,4]+0.2,(unconstrained_yval_box_lower+unconstrained_yval_box_middle)/2.),  xycoords='data',horizontalalignment='left', verticalalignment='center',fontsize=fontsize_sublabels)

ax_a.tick_params(which='both',left='off',right='off')
ax_a.set_yticks([])
for tick in ax_a.xaxis.get_major_ticks():
    tick.label.set_fontsize(10)


## ACI

ivar=1
ax_b = plt.subplot(1,3,2, adjustable='box')
ax_b.set_xlabel(var_list_titles[ivar]+' / W m$^{-2}$',fontsize=fontsize_axes)
t1=sns.distplot(PPE_array[ivar,:],color='r',hist_kws={"alpha":0.2})
t2=sns.distplot(PPE_ACI_constrained,color='Orange',hist_kws={"alpha":0.2})

ylim_starter= ax_b.get_ylim()[1]
ylim_upper=ylim_starter*1.6
ylim_top_space=ylim_upper - ylim_starter
ax_b.set_ylim((0,ylim_upper))

yval_box_middle_B2020=ylim_starter+ylim_top_space/4.*2.5
yval_box_upper_B2020= yval_box_middle_B2020+(ylim_top_space/3.)*0.16
yval_box_lower_B2020= yval_box_middle_B2020-(ylim_top_space/3.)*0.13
yvals_line_B2020=np.repeat(yval_box_middle_B2020, 30)
xvals_line_B2020=np.linspace(B2020_box[ivar,0],B2020_box[ivar,1],30)
yvals_vertical_B2020=np.linspace(yval_box_lower_B2020,yval_box_upper_B2020,30)
xvals_left_line_B2020=np.repeat(B2020_box[ivar,0],30)
xvals_right_line_B2020=np.repeat(B2020_box[ivar,1],30)

plt.plot(xvals_line_B2020,yvals_line_B2020,'--',color='indigo',linewidth=linewidth_box)
plt.plot(xvals_left_line_B2020,yvals_vertical_B2020,'-',color='indigo',linewidth=linewidth_box)
plt.plot(xvals_right_line_B2020,yvals_vertical_B2020,'-',color='indigo',linewidth=linewidth_box)


## Constrained
yval_box_middle=ylim_starter+ylim_top_space/4.*0.37
yval_box_upper= yval_box_middle+(ylim_top_space/3.*2.)*0.147
yval_box_lower= yval_box_middle-(ylim_top_space/3.*2.)*0.14
yvals_line=np.repeat(yval_box_middle, 30)
xvals_line=np.linspace(box_values[ivar,0],box_values[ivar,4],30)
box_line_left_x=np.repeat(box_values[ivar,1], 30)
box_line_right_x=np.repeat(box_values[ivar,3], 30)
box_xvals= np.linspace(box_values[ivar,1],box_values[ivar,3],30)
box_line_top_y=np.repeat(yval_box_upper,30)
box_line_bottom_y=np.repeat(yval_box_lower,30)
yvals_vertical=np.linspace(yval_box_lower,yval_box_upper,30)
xvals_left_line=np.repeat(box_values[ivar,0],30)
xvals_right_line=np.repeat(box_values[ivar,4],30)
xvals_middle_line=np.repeat(box_values[ivar,2],30)
plt.plot(xvals_line,yvals_line,'--',color='Orange',linewidth=linewidth_box)
plt.plot(xvals_left_line,yvals_vertical,'-',color='Orange',linewidth=linewidth_box)
plt.plot(xvals_right_line,yvals_vertical,'-',color='Orange',linewidth=linewidth_box)
box_patch=plt.Rectangle((box_values[ivar,1],yval_box_lower),box_values[ivar,3]-box_values[ivar,1],yval_box_upper-yval_box_lower,color='Navajowhite')
ax_b.add_patch(box_patch)
plt.plot(xvals_middle_line,yvals_vertical,'-',color='Orange',linewidth=linewidth_box+0.2)
plt.plot(box_line_left_x,yvals_vertical,'-',color='Orange',linewidth=linewidth_box)
plt.plot(box_line_right_x,yvals_vertical,'-',color='Orange',linewidth=linewidth_box)
plt.plot(box_xvals,box_line_top_y,'-',color='Orange',linewidth=linewidth_box)
plt.plot(box_xvals,box_line_bottom_y,'-',color='Orange',linewidth=linewidth_box)


## Unconstrained
unconstrained_yval_box_middle=ylim_starter+ylim_top_space/4.*1.45
unconstrained_yval_box_upper= unconstrained_yval_box_middle+(ylim_top_space/3.*2.)*0.151
unconstrained_yval_box_lower= unconstrained_yval_box_middle-(ylim_top_space/3.*2.)*0.149
unconstrained_yvals_line=np.repeat(unconstrained_yval_box_middle, 30)
unconstrained_xvals_line=np.linspace(unconstrained_box_values[ivar,0],unconstrained_box_values[ivar,4],30)
unconstrained_box_line_left_x=np.repeat(unconstrained_box_values[ivar,1], 30)
unconstrained_box_line_right_x=np.repeat(unconstrained_box_values[ivar,3], 30)
unconstrained_box_xvals= np.linspace(unconstrained_box_values[ivar,1],unconstrained_box_values[ivar,3],30)
unconstrained_box_line_top_y=np.repeat(unconstrained_yval_box_upper,30)
unconstrained_box_line_bottom_y=np.repeat(unconstrained_yval_box_lower,30)
unconstrained_yvals_vertical=np.linspace(unconstrained_yval_box_lower,unconstrained_yval_box_upper,30)
unconstrained_xvals_left_line=np.repeat(unconstrained_box_values[ivar,0],30)
unconstrained_xvals_right_line=np.repeat(unconstrained_box_values[ivar,4],30)
unconstrained_xvals_middle_line=np.repeat(unconstrained_box_values[ivar,2],30)
plt.plot(unconstrained_xvals_line,unconstrained_yvals_line,'--',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_xvals_left_line,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_xvals_right_line,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box)
unconstrained_box_patch=plt.Rectangle((unconstrained_box_values[ivar,1],unconstrained_yval_box_lower),unconstrained_box_values[ivar,3]-unconstrained_box_values[ivar,1],unconstrained_yval_box_upper-unconstrained_yval_box_lower,color='lightcoral')
ax_b.add_patch(unconstrained_box_patch)
plt.plot(unconstrained_xvals_middle_line,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box+0.2)
plt.plot(unconstrained_box_line_left_x,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_box_line_right_x,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_box_xvals,unconstrained_box_line_top_y,'-',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_box_xvals,unconstrained_box_line_bottom_y,'-',color='r',linewidth=linewidth_box)

ax_b.set_xticks((-3,-2,-1,0.0,1.0,2.0))
ax_b.set_xticklabels(('-3','-2','-1','0','1','2'),fontsize=fontsize_axes)
ax_b.set_xlim(-3.5,3.2)
ax_b.set_xlabel(var_list_titles[ivar]+' / '+units_list[ivar],size=14)
ax_b.set_ylabel('', size= 14)
ax_b.annotate('b) ', xy=(-2.8,yval_box_upper_B2020*1.06),  xycoords='data',horizontalalignment='center', verticalalignment='center',fontsize=fontsize_labels)

ax_b.tick_params(which='both',left='off',right='off')
ax_b.set_yticks([])
for tick in ax_b.xaxis.get_major_ticks():
    tick.label.set_fontsize(10)







## ARI

ivar=2
ax_c = plt.subplot(1,3,3, adjustable='box')
ax_c.set_xlabel(var_list_titles[ivar]+' / W m$^{-2}$',fontsize=fontsize_axes)
t1=sns.distplot(PPE_array[ivar,:],color='r',hist_kws={"alpha":0.2})
t2=sns.distplot(PPE_ARI_constrained,color='Orange',hist_kws={"alpha":0.2})

ylim_starter= ax_c.get_ylim()[1]
ylim_upper=ylim_starter*1.6
ylim_top_space=ylim_upper - ylim_starter
ax_c.set_ylim((0,ylim_upper))

yval_box_middle_B2020=ylim_starter+ylim_top_space/4.*2.5
yval_box_upper_B2020= yval_box_middle_B2020+(ylim_top_space/3.)*0.16
yval_box_lower_B2020= yval_box_middle_B2020-(ylim_top_space/3.)*0.13
yvals_line_B2020=np.repeat(yval_box_middle_B2020, 30)
xvals_line_B2020=np.linspace(B2020_box[ivar,0],B2020_box[ivar,1],30)
yvals_vertical_B2020=np.linspace(yval_box_lower_B2020,yval_box_upper_B2020,30)
xvals_left_line_B2020=np.repeat(B2020_box[ivar,0],30)
xvals_right_line_B2020=np.repeat(B2020_box[ivar,1],30)

plt.plot(xvals_line_B2020,yvals_line_B2020,'--',color='indigo',linewidth=linewidth_box)
plt.plot(xvals_left_line_B2020,yvals_vertical_B2020,'-',color='indigo',linewidth=linewidth_box)
plt.plot(xvals_right_line_B2020,yvals_vertical_B2020,'-',color='indigo',linewidth=linewidth_box)


## Constrained
yval_box_middle=ylim_starter+ylim_top_space/4.*0.37
yval_box_upper= yval_box_middle+(ylim_top_space/3.*2.)*0.147
yval_box_lower= yval_box_middle-(ylim_top_space/3.*2.)*0.14
yvals_line=np.repeat(yval_box_middle, 30)
xvals_line=np.linspace(box_values[ivar,0],box_values[ivar,4],30)
box_line_left_x=np.repeat(box_values[ivar,1], 30)
box_line_right_x=np.repeat(box_values[ivar,3], 30)
box_xvals= np.linspace(box_values[ivar,1],box_values[ivar,3],30)
box_line_top_y=np.repeat(yval_box_upper,30)
box_line_bottom_y=np.repeat(yval_box_lower,30)
yvals_vertical=np.linspace(yval_box_lower,yval_box_upper,30)
xvals_left_line=np.repeat(box_values[ivar,0],30)
xvals_right_line=np.repeat(box_values[ivar,4],30)
xvals_middle_line=np.repeat(box_values[ivar,2],30)
plt.plot(xvals_line,yvals_line,'--',color='Orange',linewidth=linewidth_box)
plt.plot(xvals_left_line,yvals_vertical,'-',color='Orange',linewidth=linewidth_box)
plt.plot(xvals_right_line,yvals_vertical,'-',color='Orange',linewidth=linewidth_box)
box_patch=plt.Rectangle((box_values[ivar,1],yval_box_lower),box_values[ivar,3]-box_values[ivar,1],yval_box_upper-yval_box_lower,color='Navajowhite')
ax_c.add_patch(box_patch)
plt.plot(xvals_middle_line,yvals_vertical,'-',color='Orange',linewidth=linewidth_box+0.2)
plt.plot(box_line_left_x,yvals_vertical,'-',color='Orange',linewidth=linewidth_box)
plt.plot(box_line_right_x,yvals_vertical,'-',color='Orange',linewidth=linewidth_box)
plt.plot(box_xvals,box_line_top_y,'-',color='Orange',linewidth=linewidth_box)
plt.plot(box_xvals,box_line_bottom_y,'-',color='Orange',linewidth=linewidth_box)


## Unconstrained
unconstrained_yval_box_middle=ylim_starter+ylim_top_space/4.*1.45
unconstrained_yval_box_upper= unconstrained_yval_box_middle+(ylim_top_space/3.*2.)*0.151
unconstrained_yval_box_lower= unconstrained_yval_box_middle-(ylim_top_space/3.*2.)*0.149
unconstrained_yvals_line=np.repeat(unconstrained_yval_box_middle, 30)
unconstrained_xvals_line=np.linspace(unconstrained_box_values[ivar,0],unconstrained_box_values[ivar,4],30)
unconstrained_box_line_left_x=np.repeat(unconstrained_box_values[ivar,1], 30)
unconstrained_box_line_right_x=np.repeat(unconstrained_box_values[ivar,3], 30)
unconstrained_box_xvals= np.linspace(unconstrained_box_values[ivar,1],unconstrained_box_values[ivar,3],30)
unconstrained_box_line_top_y=np.repeat(unconstrained_yval_box_upper,30)
unconstrained_box_line_bottom_y=np.repeat(unconstrained_yval_box_lower,30)
unconstrained_yvals_vertical=np.linspace(unconstrained_yval_box_lower,unconstrained_yval_box_upper,30)
unconstrained_xvals_left_line=np.repeat(unconstrained_box_values[ivar,0],30)
unconstrained_xvals_right_line=np.repeat(unconstrained_box_values[ivar,4],30)
unconstrained_xvals_middle_line=np.repeat(unconstrained_box_values[ivar,2],30)
plt.plot(unconstrained_xvals_line,unconstrained_yvals_line,'--',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_xvals_left_line,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_xvals_right_line,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box)
unconstrained_box_patch=plt.Rectangle((unconstrained_box_values[ivar,1],unconstrained_yval_box_lower),unconstrained_box_values[ivar,3]-unconstrained_box_values[ivar,1],unconstrained_yval_box_upper-unconstrained_yval_box_lower,color='lightcoral')
ax_c.add_patch(unconstrained_box_patch)
plt.plot(unconstrained_xvals_middle_line,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box+0.2)
plt.plot(unconstrained_box_line_left_x,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_box_line_right_x,unconstrained_yvals_vertical,'-',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_box_xvals,unconstrained_box_line_top_y,'-',color='r',linewidth=linewidth_box)
plt.plot(unconstrained_box_xvals,unconstrained_box_line_bottom_y,'-',color='r',linewidth=linewidth_box)

ax_c.set_xticks((-3,-2,-1,0.0,1.0,2.0))
ax_c.set_xticklabels(('-3','-2','-1','0','1','2'),fontsize=fontsize_axes)
ax_c.set_xlim(-3.5,3.2)
ax_c.set_xlabel(var_list_titles[ivar]+' / '+units_list[ivar],size=14)
ax_c.set_ylabel('', size= 14)
ax_c.annotate('c) ', xy=(-2.8,yval_box_upper_B2020*1.06),  xycoords='data',horizontalalignment='center', verticalalignment='center',fontsize=fontsize_labels)

ax_c.tick_params(which='both',left='off',right='off')
ax_c.set_yticks([])
for tick in ax_c.xaxis.get_major_ticks():
    tick.label.set_fontsize(10)


ax_a.set_box_aspect(1.0)
ax_b.set_box_aspect(1.0)
ax_c.set_box_aspect(1.0)
plt.tight_layout()
plt.subplots_adjust(wspace=0.1, hspace=0.1)
fname_out= 'pdfs_ERF_ACI_ARI_global_ann_seaborne_w_B2020_from_sample_and_constraint'
plt.savefig(outdir+fname_out+'.pdf', format='pdf', dpi=300,pad_inches = 0.2,bbox_inches='tight')



