#################################################
## Pairwise plot of average NRMSE effects 
## with subpanel to exemplify effect of individual constraint
#################################################


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

nvar=450
sample_size=1000000

region_list_unsorted=np.hstack((np.repeat('NAtlantic',54),np.repeat('NPacific',34),np.repeat('SPacific',43),np.repeat('SAtlantic',43),np.repeat('SOcean',41),'hemispheric_diff',np.repeat('SOcean',43),np.repeat('SAtlantic',46),np.repeat('SPacific',42),np.repeat('NPacific',56),np.repeat('NAtlantic',34)))

greek_letters_lower=[chr(code) for code in range(945,970)]
TAU=greek_letters_lower[19]
greek_letters_upper=[chr(code) for code in range(913,938)]
DELTA=greek_letters_upper[3]


#############################################
## Read in array of cross-variable constraint effects
#############################################

NRMSE_change_array=np.loadtxt('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/NRMSE_values_entire_sample_w_Hd/cross_constraint_effects_resorted.csv', delimiter=',')


#####################
## Plotting features
#####################

pc_change_labels_values=[-90,-70,-50,-30,-10,0,10,30,50,70,90]
pc_change_labels=['<-80','-80 to -60','-60 to -40','-40 to -20','-20 to -5','-5 to 5','5 to 20','20 to 40','40 to 60','60 to 80','80<']
pc_change_labels_spaced=['<-80','-80  to  -60','-60  to  -40','-40  to  -20','-20  to  -5','-5  to  5','5  to  20','20  to  40','40  to  60','60  to  80','80<']
from matplotlib import patches as mpatches
from matplotlib import colors

cmap_divergent = colors.ListedColormap(['teal','lightseagreen','turquoise','paleturquoise','lightcyan','white','thistle','plum','violet','orchid','darkorchid'])##11
bounds_divergent=[-100,-80,-60,-40,-20,-5,5,20,40,60,80,250]
norm_divergent = colors.BoundaryNorm(bounds_divergent, cmap_divergent.N)


#####################################
## Plot NAtlantic effect on NAtlantic
## w/ hemispheric diff seasonal amplitude
## as test of plotting features
#####################################

NRMSE_change_array_subset= NRMSE_change_array[0:102,0:102]

fig,ax=plt.subplots(1,1)
im=ax.imshow(NRMSE_change_array_subset,cmap=cmap_divergent,norm=norm_divergent)
plt.ylabel('Variable used for constraint')
plt.xlabel('Constrained variable')
colors = [ im.cmap(im.norm(value)) for value in pc_change_labels_values]
patches = [ mpatches.Patch(color=colors[i], label="{l}".format(l=pc_change_labels[i]) ) for i in range(len(pc_change_labels)) ]
plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
x_label_list = ['$H_{d}$', 'T', '$F_{SW}$', '$N_{d}$','$f_{c}$','$r_{e}$','LWP',TAU+'$_{c}$']
ax.set_xticks([0.0,14.0,18.0,32.0,46.0,60.0,74.0,88.0])
ax.set_xticklabels(x_label_list,rotation='90')
y_label_list = ['$H_{d}$', 'T', '$F_{SW}$', '$N_{d}$','$f_{c}$','$r_{e}$','LWP',TAU+'$_{c}$']
ax.set_yticks([0.0,14.0,18.0,32.0,46.0,60.0,74.0,88.0])
ax.set_yticklabels(x_label_list,rotation='0')
fig.tight_layout()
plt.savefig('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/NRMSE_values_entire_sample_w_Hd/cross_variable_NRMSE_effects/cross_var_NRMSE_effects_NAtlantic_ON_NAtlantic.png')
plt.savefig('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/NRMSE_values_entire_sample_w_Hd/cross_variable_NRMSE_effects/cross_var_NRMSE_effects_NAtlantic_ON_NAtlantic.pdf')

NRMSE_change_array_subset= NRMSE_change_array[0:102,0:102]


##########################################################################################
## JULY Fsw constrained by other individual constraint variables
## also in jul (b/c cross-variable figure indicates clear effects)
## Inidices here are according to order of statistical emulation
##########################################################################################

Fsw_obs=np.loadtxt('/gws/nopw/j04/acure/lregayre/data_for_emulation/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIC_CF.dat')[7]
Fsw_unconstrained=np.loadtxt('/gws/nopw/j04/acure/lregayre/data_post_emulation/emulated_mean_values_jul_Fsw_NAtlantic_1000000_w_o_carb.dat')
Fsw_constrained_Fsw_indices=np.loadtxt('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/individual_NRMSE_indices/original_additional_and_Hd/all_1M_python_indices_retained_w_0pc_error_for_variable_14.dat')
Fsw_constrained_LWP_indices=np.loadtxt('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/individual_NRMSE_indices/original_additional_and_Hd/all_1M_python_indices_retained_w_0pc_error_for_variable_417.dat')
Fsw_constrained_Nd_indices=np.loadtxt('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/individual_NRMSE_indices/original_additional_and_Hd/all_1M_python_indices_retained_w_0pc_error_for_variable_26.dat')


Fsw_constrained_Fsw=np.zeros((len(Fsw_constrained_Fsw_indices)))
Fsw_constrained_LWP=np.zeros((len(Fsw_constrained_LWP_indices)))
Fsw_constrained_Nd=np.zeros((len(Fsw_constrained_Nd_indices)))


for imember in np.arange(len(Fsw_constrained_Fsw_indices)):
    Fsw_constrained_Fsw[imember]=Fsw_unconstrained[int(Fsw_constrained_Fsw_indices[imember])]

for imember in np.arange(len(Fsw_constrained_LWP_indices)):
    Fsw_constrained_LWP[imember]=Fsw_unconstrained[int(Fsw_constrained_LWP_indices[imember])]

for imember in np.arange(len(Fsw_constrained_Nd_indices)):
    Fsw_constrained_Nd[imember]=Fsw_unconstrained[int(Fsw_constrained_Nd_indices[imember])]


line_obs= mpatches.Patch(color='darkgray', label='Observed value')
line_unconstrained= mpatches.Patch(color='k', label='Unconstrained (*)')
line_con_by_Nd= mpatches.Patch(color='lightseagreen', label='Constrained by $N_{d}$')
line_con_by_LWP= mpatches.Patch(color='orchid', label='Constrained by LWP')

plt.rcParams['font.size'] = '8' ## for subplot title, which is otherwise inaccessible


fig=plt.figure()
gs=fig.add_gridspec(20,24)
gs.update(wspace=0.025, hspace=0.025)
ax2=fig.add_subplot(gs[1:8,18:20])
ax2.legend(handles=patches,loc=2,fontsize=8)
ax2.text(0.1,1.05,'NRMSE change (%)', fontsize=8.5)
ax2.axis('off')
ax1=fig.add_subplot(gs[0:19,0:18])
im=ax1.imshow(NRMSE_change_array_subset,cmap=cmap_divergent,norm=norm_divergent)
plt.ylabel('Variable used for constraint')
plt.xlabel('Constrained variable')
colors = [ im.cmap(im.norm(value)) for value in pc_change_labels_values]
patches = [ mpatches.Patch(color=colors[i], label="{l}".format(l=pc_change_labels[i]) ) for i in range(len(pc_change_labels)) ]
x_label_list = ['$H_{d}$', 'T', '$F_{SW}$', '*', '$N_{d}$','$f_{c}$','$r_{e}$','LWP',TAU+'$_{c}$']
ax1.set_xticks([0.0,14.0,18.0,25.0,32.0,46.0,60.0,74.0,88.0])
ax1.set_xticklabels(x_label_list,rotation='90')
y_label_list = ['$H_{d}$', 'T', '$F_{SW}$', '$N_{d}$', '$f_{c}$','$r_{e}$','LWP', TAU+'$_{c}$']
ax1.set_yticks([0.0,14.0,32.0,40.0,60.0,74.0,82.0,88.0])
ax1.set_yticklabels(y_label_list,rotation='0')

ax3=fig.add_subplot(gs[12:19,18:])
ax3.set_yticklabels([])
ax3.set_yticks([])

ax3.set_ylim((0,0.12))
t1=sns.distplot(Fsw_unconstrained,color='k',hist_kws={"alpha":0.2},kde_kws={'linewidth':0.8})
t2=sns.distplot(Fsw_constrained_Nd,color='lightseagreen',hist_kws={"alpha":0.2},kde_kws={'linewidth':0.8})
t3=sns.distplot(Fsw_constrained_LWP,color='orchid',hist_kws={"alpha":0.2},kde_kws={'linewidth':0.8})
plt.axvline(x=Fsw_obs,color='darkgray',linestyle='dashed',linewidth=1.)
ax3.set_xlim((75,175))
leg = plt.legend(handles=[line_obs,line_unconstrained,line_con_by_Nd,line_con_by_LWP], loc=1, prop={'size':6})
for line in leg.get_lines():
    line.set_linewidth(0.4)



plt.xlabel('F$_{SW}$ / W m$^{-2}$')
fig.tight_layout()
plt.savefig('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/NRMSE_values_entire_sample_w_Hd/cross_variable_NRMSE_effects/cross_var_NRMSE_effects_NAtlantic_ON_NAtlantic_w_subpanel_Fsw_version.png')
plt.savefig('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/NRMSE_values_entire_sample_w_Hd/cross_variable_NRMSE_effects/cross_var_NRMSE_effects_NAtlantic_ON_NAtlantic_w_subpanel_Fsw_version.pdf')


