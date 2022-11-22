#################################################
## Use NRMSE values from 1 million member sample
## to plot the pairwise-consistency
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

greek_letters_lower=[chr(code) for code in range(945,970)]
TAU=greek_letters_lower[19]
greek_letters_upper=[chr(code) for code in range(913,938)]
DELTA=greek_letters_upper[3]


############################################
## Read in array of cross-variable constraint effects
############################################


NRMSE_change_array=np.loadtxt('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/NRMSE_values_entire_sample_w_Hd/cross_constraint_effects_resorted.csv', delimiter=',')


#########################
## Plotting features
## consistent across figures
##########################

pc_change_labels_values=[-90,-70,-50,-30,-10,0,10,30,50,70,90]
pc_change_labels=['<-80','-80 to -60','-60 to -40','-40 to -20','-20 to -5','-5 to 5','5 to 20','20 to 40','40 to 60','60 to 80','80<']
from matplotlib import patches as mpatches
from matplotlib import colors

cmap_divergent = colors.ListedColormap(['teal','lightseagreen','turquoise','paleturquoise','lightcyan','white','thistle','plum','violet','orchid','darkorchid'])
bounds_divergent=[-100,-80,-60,-40,-20,-5,5,20,40,60,80,250]
norm_divergent = colors.BoundaryNorm(bounds_divergent, cmap_divergent.N)


#####################################
## Plot NAtlantic effect on NAtlantic
## w/ hemispheric diff seasonal amplitude
## as a test of ploitting routine
#####################################

NRMSE_change_array_subset= NRMSE_change_array[0:102,0:102]

fig,ax=plt.subplots(1,1)
im=ax.imshow(NRMSE_change_array_subset,cmap=cmap_divergent,norm=norm_divergent)
plt.ylabel('Variable used for constraint',fontsize=12)
plt.xlabel('Constrained variable',fontsize=12)
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
plt.show()




x_label_list_SPacific = ['$H_{d}$', 'T', '\n $F_{SW}$', '$N_{d}$','$f_{c}$','$r_{e}$','LWP',TAU+'$_{c}$']
y_label_list_SPacific = ['$H_{d}$', 'T', '\n $F_{SW}$', '$N_{d}$','$f_{c}$','$r_{e}$','LWP',TAU+'$_{c}$']
y_label_list_SOcean= ['$H_{d}$', '$F_{SW}$', '$N_{d}$','$f_{c}$','$r_{e}$','LWP',TAU+'$_{c}$']
x_label_list_SOcean= ['$H_{d}$', '$F_{SW}$', '$N_{d}$','$f_{c}$','$r_{e}$','LWP',TAU+'$_{c}$']

def make_NRMSE_plot(x_indices,y_indices,xtick_pos,ytick_pos,x_label,y_label,fname):
    array_len_x=x_indices[1]-x_indices[0]+14
    array_len_y=y_indices[1]-y_indices[0]+14
    NRMSE_subset=np.zeros((array_len_y,array_len_x))
    NRMSE_subset[0:14,0:14]=NRMSE_change_array[0:14,0:14]
    NRMSE_subset[0:14,14:array_len_x]=NRMSE_change_array[0:14,x_indices[0]:x_indices[1]]
    NRMSE_subset[14:array_len_y,0:14]=NRMSE_change_array[y_indices[0]:y_indices[1],0:14]
    NRMSE_subset[14:array_len_y,14:array_len_x]=NRMSE_change_array[y_indices[0]:y_indices[1],x_indices[0]:x_indices[1]]
    fig,ax=plt.subplots(1,1)
    im=ax.imshow(NRMSE_subset,cmap=cmap_divergent,norm=norm_divergent)
    plt.ylabel('Variable used for constraint',fontsize=12)
    plt.xlabel('Constrained variable',fontsize=12)
    colors = [ im.cmap(im.norm(value)) for value in pc_change_labels_values]
    patches = [ mpatches.Patch(color=colors[i], label="{l}".format(l=pc_change_labels[i]) ) for i in range(len(pc_change_labels)) ]
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
    ax.set_xticks(xtick_pos)
    if (x_label=='Southern Ocean'):
        ax.set_xticklabels(x_label_list_SOcean,rotation='90')
    elif (x_label=='South Pacific'):
        ax.set_xticklabels(x_label_list_SPacific,rotation='90')
    else:
        ax.set_xticklabels(x_label_list,rotation='90')
    ax.set_yticks(ytick_pos)
    if (y_label=='Southern Ocean'):
        ax.set_yticklabels(y_label_list_SOcean,rotation='0')
    elif (y_label=='South Pacific'):
        ax.set_yticklabels(y_label_list_SPacific,rotation='0')
    else:
        ax.set_yticklabels(y_label_list,rotation='0')
    fig.tight_layout()
    plt.savefig('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/NRMSE_values_entire_sample_w_Hd/cross_variable_NRMSE_effects/cross_var_NRMSE_effects_'+fname+'.png')
    plt.savefig('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/NRMSE_values_entire_sample_w_Hd/cross_variable_NRMSE_effects/cross_var_NRMSE_effects_'+fname+'.pdf')




NAtlantic_ticks=[-0.25,13.75,17.75,31.75,45.75,59.75,74.25,88.25]
NPacific_ticks=[0.0,14.0,20.0,34.0,48.0,62.0,76.0,90.0]
SAtlantic_ticks=[0.0,14.0,19.0,33.0,47.0,61.0,75.0,89.0]
SPacific_ticks=[0.0,13.75,14.75,28.75,42.75,57.0,71.0,85.0]
SOcean_ticks=[0.0,14.0,28.0,41.75,56.0,70.0,84.0]


NAtlantic_indices=[14,102]
NPacific_indices=[102,192]
SAtlantic_indices=[192,281]
SPacific_indices=[281,366]
SOcean_indices=[366,450]

y_indices=NAtlantic_indices ## constraint
x_indices=NAtlantic_indices ## constrained

make_NRMSE_plot(NAtlantic_indices,NAtlantic_indices,NAtlantic_ticks,NAtlantic_ticks,'North Atlantic','North Atlantic','NAtlantic_ON_NAtlantic') ## y-axis ON x-axis

#########################################################################################################
## Combine constraint regions and constrained regions
## y-axis ON x-axis
#########################################################################################################

# Natlantic on others

make_NRMSE_plot(NPacific_indices,NAtlantic_indices,NPacific_ticks,NAtlantic_ticks,'North Pacific','North Atlantic','NAtlantic_ON_NPacific') ## y-axis ON x-axis
make_NRMSE_plot(SPacific_indices,NAtlantic_indices,SPacific_ticks,NAtlantic_ticks,'South Pacific','North Atlantic','NAtlantic_ON_SPacific') ## y-axis ON x-axis
make_NRMSE_plot(SAtlantic_indices,NAtlantic_indices,SAtlantic_ticks,NAtlantic_ticks,'South Atlantic','North Atlantic','NAtlantic_ON_SAtlantic') ## y-axis ON x-axis
make_NRMSE_plot(SOcean_indices,NAtlantic_indices,SOcean_ticks,NAtlantic_ticks,'Southern Ocean','North Atlantic','NAtlantic_ON_SOcean') ## y-axis ON x-axis

# NPacific on self and others

make_NRMSE_plot(NAtlantic_indices,NPacific_indices,NAtlantic_ticks,NPacific_ticks,'North Atlantic','North Pacific','NPacific_ON_NAtlantic') ## y-axis ON x-axis
make_NRMSE_plot(NPacific_indices,NPacific_indices,NPacific_ticks,NPacific_ticks,'North Pacific','North Pacific','NPacific_ON_NPacific') ## y-axis ON x-axis
make_NRMSE_plot(SPacific_indices,NPacific_indices,SPacific_ticks,NPacific_ticks,'South Pacific','North Pacific','NPacific_ON_SPacific') ## y-axis ON x-axis
make_NRMSE_plot(SAtlantic_indices,NPacific_indices,SAtlantic_ticks,NPacific_ticks,'South Atlantic','North Pacific','NPacific_ON_SAtlantic') ## y-axis ON x-axis
make_NRMSE_plot(SOcean_indices,NPacific_indices,SOcean_ticks,NPacific_ticks,'Southern Ocean','North Pacific','NPacific_ON_SOcean') ## y-axis ON x-axis

# SPacific on self and others

make_NRMSE_plot(NAtlantic_indices,SPacific_indices,NAtlantic_ticks,SPacific_ticks,'North Atlantic','South Pacific','SPacific_ON_NAtlantic') ## y-axis ON x-axis
make_NRMSE_plot(NPacific_indices,SPacific_indices,NPacific_ticks,SPacific_ticks,'North Pacific','South Pacific','SPacific_ON_NPacific') ## y-axis ON x-axis
make_NRMSE_plot(SPacific_indices,SPacific_indices,SPacific_ticks,SPacific_ticks,'South Pacific','South Pacific','SPacific_ON_SPacific') ## y-axis ON x-axis
make_NRMSE_plot(SAtlantic_indices,SPacific_indices,SAtlantic_ticks,SPacific_ticks,'South Atlantic','South Pacific','SPacific_ON_SAtlantic') ## y-axis ON x-axis
make_NRMSE_plot(SOcean_indices,SPacific_indices,SOcean_ticks,SPacific_ticks,'Southern Ocean','South Pacific','SPacific_ON_SOcean') ## y-axis ON x-axis

# SAtlantic on self and others

make_NRMSE_plot(NAtlantic_indices,SAtlantic_indices,NAtlantic_ticks,SAtlantic_ticks,'North Atlantic','South Atlantic','SAtlantic_ON_NAtlantic') ## y-axis ON x-axis
make_NRMSE_plot(NPacific_indices,SAtlantic_indices,NPacific_ticks,SAtlantic_ticks,'North Pacific','South Atlantic','SAtlantic_ON_NPacific') ## y-axis ON x-axis
make_NRMSE_plot(SPacific_indices,SAtlantic_indices,SPacific_ticks,SAtlantic_ticks,'South Pacific','South Atlantic','SAtlantic_ON_SPacific') ## y-axis ON x-axis
make_NRMSE_plot(SAtlantic_indices,SAtlantic_indices,SAtlantic_ticks,SAtlantic_ticks,'South Atlantic','South Atlantic','SAtlantic_ON_SAtlantic') ## y-axis ON x-axis
make_NRMSE_plot(SOcean_indices,SAtlantic_indices,SOcean_ticks,SAtlantic_ticks,'Southern Ocean','South Atlantic','SAtlantic_ON_SOcean') ## y-axis ON x-axis

# SOcean on self and others

make_NRMSE_plot(NAtlantic_indices,SOcean_indices,NAtlantic_ticks,SOcean_ticks,'North Atlantic','Southern Ocean','SOcean_ON_NAtlantic') ## y-axis ON x-axis
make_NRMSE_plot(NPacific_indices,SOcean_indices,NPacific_ticks,SOcean_ticks,'North Pacific','Southern Ocean','SOcean_ON_NPacific') ## y-axis ON x-axis
make_NRMSE_plot(SPacific_indices,SOcean_indices,SPacific_ticks,SOcean_ticks,'South Pacific','Southern Ocean','SOcean_ON_SPacific') ## y-axis ON x-axis
make_NRMSE_plot(SAtlantic_indices,SOcean_indices,SAtlantic_ticks,SOcean_ticks,'South Atlantic','Southern Ocean','SOcean_ON_SAtlantic') ## y-axis ON x-axis
make_NRMSE_plot(SOcean_indices,SOcean_indices,SOcean_ticks,SOcean_ticks,'Southern Ocean','Southern Ocean','SOcean_ON_SOcean') ## y-axis ON x-axis


