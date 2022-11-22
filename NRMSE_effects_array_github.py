#################################################
## Evaluate NRMSE values from 1 million member sample
## and indices of individual constraint 
## to quantify pairwise constraint effects
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

NRMSE_path_all='/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/NRMSE_values_entire_sample/'
NRMSE_path_w_Hd='/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/NRMSE_values_entire_sample_w_Hd/'

indices_path='/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/individual_NRMSE_indices/additional/'
indices_path_w_Hd='/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/individual_NRMSE_indices/only_Hd/'

region_list_unsorted=np.hstack((np.repeat('NAtlantic',54),np.repeat('NPacific',34),np.repeat('SPacific',43),np.repeat('SAtlantic',43),np.repeat('SOcean',41),'hemispheric_diff',np.repeat('SOcean',43),np.repeat('SAtlantic',46),np.repeat('SPacific',42),np.repeat('NPacific',56),np.repeat('NAtlantic',34)))

constraint_variable_text_list_unsorted=['jul_transect_Nd_NAtlantic','jul_transect_LWP_Nd_NAtlantic','jul_transect_Re_Nd_NAtlantic','jul_transect_Re_NAtlantic','seasonal_amplitude_Nd_NAtlantic','seasonal_amplitude_Tau_NAtlantic','seasonal_amplitude_LWP_NAtlantic','seasonal_amplitude_Fsw_NAtlantic','jan_Fsw_NAtlantic','feb_Fsw_NAtlantic','mar_Fsw_NAtlantic','apr_Fsw_NAtlantic','may_Fsw_NAtlantic','jun_Fsw_NAtlantic','jul_Fsw_NAtlantic','aug_Fsw_NAtlantic','sep_Fsw_NAtlantic','oct_Fsw_NAtlantic','nov_Fsw_NAtlantic','dec_Fsw_NAtlantic','ann_Fsw_NAtlantic','feb_Nd_NAtlantic','mar_Nd_NAtlantic','apr_Nd_NAtlantic','may_Nd_NAtlantic','jun_Nd_NAtlantic','jul_Nd_NAtlantic','aug_Nd_NAtlantic','sep_Nd_NAtlantic','oct_Nd_NAtlantic','ann_Nd_NAtlantic','jan_Fc_NAtlantic','feb_Fc_NAtlantic','mar_Fc_NAtlantic','apr_Fc_NAtlantic','may_Fc_NAtlantic','jun_Fc_NAtlantic','jul_Fc_NAtlantic','aug_Fc_NAtlantic','sep_Fc_NAtlantic','oct_Fc_NAtlantic','nov_Fc_NAtlantic','dec_Fc_NAtlantic','ann_Fc_NAtlantic','jan_Re_NAtlantic','apr_Re_NAtlantic','may_Re_NAtlantic','jun_Re_NAtlantic','jul_Re_NAtlantic','aug_Re_NAtlantic','sep_Re_NAtlantic','oct_Re_NAtlantic','dec_Re_NAtlantic','ann_Re_NAtlantic','jul_transect_Nd_NPacific','jul_transect_Cf_NPacific','jul_transect_LWP_NPacific','jul_transect_Cf_Nd_NPacific','jul_transect_LWP_Nd_NPacific','jul_transect_Re_Nd_NPacific','seasonal_amplitude_Fc_NPacific','seasonal_amplitude_Fsw_NPacific','jan_Fsw_NPacific','feb_Fsw_NPacific','mar_Fsw_NPacific','apr_Fsw_NPacific','may_Fsw_NPacific','jun_Fsw_NPacific','jul_Fsw_NPacific','aug_Fsw_NPacific','sep_Fsw_NPacific','oct_Fsw_NPacific','nov_Fsw_NPacific','dec_Fsw_NPacific','ann_Fsw_NPacific','jan_Fc_NPacific','feb_Fc_NPacific','mar_Fc_NPacific','apr_Fc_NPacific','may_Fc_NPacific','jun_Fc_NPacific','jul_Fc_NPacific','aug_Fc_NPacific','sep_Fc_NPacific','oct_Fc_NPacific','nov_Fc_NPacific','dec_Fc_NPacific','ann_Fc_NPacific','nov_transect_Nd_SPacific','seasonal_amplitude_Nd_SPacific','seasonal_amplitude_Fc_SPacific','seasonal_amplitude_Fsw_SPacific','jan_Fsw_SPacific','feb_Fsw_SPacific','mar_Fsw_SPacific','apr_Fsw_SPacific','may_Fsw_SPacific','jun_Fsw_SPacific','jul_Fsw_SPacific','aug_Fsw_SPacific','sep_Fsw_SPacific','oct_Fsw_SPacific','nov_Fsw_SPacific','dec_Fsw_SPacific','ann_Fsw_SPacific','jan_Nd_SPacific','feb_Nd_SPacific','mar_Nd_SPacific','apr_Nd_SPacific','may_Nd_SPacific','jun_Nd_SPacific','jul_Nd_SPacific','aug_Nd_SPacific','sep_Nd_SPacific','oct_Nd_SPacific','nov_Nd_SPacific','dec_Nd_SPacific','ann_Nd_SPacific','jan_Fc_SPacific','feb_Fc_SPacific','mar_Fc_SPacific','apr_Fc_SPacific','may_Fc_SPacific','jun_Fc_SPacific','jul_Fc_SPacific','aug_Fc_SPacific','sep_Fc_SPacific','oct_Fc_SPacific','nov_Fc_SPacific','dec_Fc_SPacific','ann_Fc_SPacific','nov_transect_Nd_SAtlantic','nov_transect_LWP_Nd_SAtlantic','nov_transect_Fc_SAtlantic','nov_transect_AI_SAtlantic','nov_transect_Nd_AI_SAtlantic','seasonal_amplitude_Fc_SAtlantic','seasonal_amplitude_Fsw_SAtlantic','jan_Fsw_SAtlantic','feb_Fsw_SAtlantic','mar_Fsw_SAtlantic','apr_Fsw_SAtlantic','may_Fsw_SAtlantic','jun_Fsw_SAtlantic','jul_Fsw_SAtlantic','aug_Fsw_SAtlantic','sep_Fsw_SAtlantic','oct_Fsw_SAtlantic','nov_Fsw_SAtlantic','dec_Fsw_SAtlantic','ann_Fsw_SAtlantic','jan_Nd_SAtlantic','feb_Nd_SAtlantic','mar_Nd_SAtlantic','apr_Nd_SAtlantic','may_Nd_SAtlantic','jun_Nd_SAtlantic','jul_Nd_SAtlantic','aug_Nd_SAtlantic','nov_Nd_SAtlantic','dec_Nd_SAtlantic','ann_Nd_SAtlantic','jan_Fc_SAtlantic','feb_Fc_SAtlantic','mar_Fc_SAtlantic','apr_Fc_SAtlantic','may_Fc_SAtlantic','jun_Fc_SAtlantic','jul_Fc_SAtlantic','aug_Fc_SAtlantic','sep_Fc_SAtlantic','oct_Fc_SAtlantic','nov_Fc_SAtlantic','ann_Fc_SAtlantic','seasonal_amplitude_Fsw_SOcean','seasonal_amplitude_Nd_SOcean','jan_Fsw_SOcean','feb_Fsw_SOcean','mar_Fsw_SOcean','apr_Fsw_SOcean','may_Fsw_SOcean','jun_Fsw_SOcean','jul_Fsw_SOcean','aug_Fsw_SOcean','sep_Fsw_SOcean','oct_Fsw_SOcean','nov_Fsw_SOcean','dec_Fsw_SOcean','ann_Fsw_SOcean','jan_Nd_SOcean','feb_Nd_SOcean','mar_Nd_SOcean','apr_Nd_SOcean','may_Nd_SOcean','jun_Nd_SOcean','jul_Nd_SOcean','aug_Nd_SOcean','sep_Nd_SOcean','oct_Nd_SOcean','nov_Nd_SOcean','dec_Nd_SOcean','ann_Nd_SOcean','jan_Fc_SOcean','feb_Fc_SOcean','mar_Fc_SOcean','apr_Fc_SOcean','may_Fc_SOcean','jun_Fc_SOcean','jul_Fc_SOcean','aug_Fc_SOcean','sep_Fc_SOcean','oct_Fc_SOcean','nov_Fc_SOcean','dec_Fc_SOcean','ann_Fc_SOcean','seasonal_amplitude_Nd_hemispheric_diff',"seasonal_amplitude_LWP_SOcean","seasonal_amplitude_Fc_SOcean","seasonal_amplitude_Re_SOcean","seasonal_amplitude_Tau_SOcean","jan_LWP_SOcean","feb_LWP_SOcean","mar_LWP_SOcean","apr_LWP_SOcean","may_LWP_SOcean","jun_LWP_SOcean","jul_LWP_SOcean","aug_LWP_SOcean","sep_LWP_SOcean","oct_LWP_SOcean","nov_LWP_SOcean","dec_LWP_SOcean","ann_LWP_SOcean","jan_Re_SOcean","feb_Re_SOcean","mar_Re_SOcean","apr_Re_SOcean","may_Re_SOcean","jun_Re_SOcean","jul_Re_SOcean","aug_Re_SOcean","sep_Re_SOcean","oct_Re_SOcean","nov_Re_SOcean","dec_Re_SOcean","ann_Re_SOcean","jan_Tau_SOcean","feb_Tau_SOcean","mar_Tau_SOcean","apr_Tau_SOcean","may_Tau_SOcean","jun_Tau_SOcean","jul_Tau_SOcean","aug_Tau_SOcean","sep_Tau_SOcean","oct_Tau_SOcean","nov_Tau_SOcean","dec_Tau_SOcean","ann_Tau_SOcean","seasonal_amplitude_Nd_SAtlantic","seasonal_amplitude_Re_SAtlantic","seasonal_amplitude_LWP_SAtlantic","seasonal_amplitude_Tau_SAtlantic","sep_Nd_SAtlantic","oct_Nd_SAtlantic","dec_Nd_SAtlantic","jan_LWP_SAtlantic","feb_LWP_SAtlantic","mar_LWP_SAtlantic","apr_LWP_SAtlantic","may_LWP_SAtlantic","jun_LWP_SAtlantic","jul_LWP_SAtlantic","aug_LWP_SAtlantic","sep_LWP_SAtlantic","oct_LWP_SAtlantic","nov_LWP_SAtlantic","dec_LWP_SAtlantic","ann_LWP_SAtlantic","jan_Re_SAtlantic","feb_Re_SAtlantic","mar_Re_SAtlantic","apr_Re_SAtlantic","may_Re_SAtlantic","jun_Re_SAtlantic","jul_Re_SAtlantic","aug_Re_SAtlantic","sep_Re_SAtlantic","oct_Re_SAtlantic","nov_Re_SAtlantic","dec_Re_SAtlantic","ann_Re_SAtlantic","jan_Tau_SAtlantic","feb_Tau_SAtlantic","mar_Tau_SAtlantic","apr_Tau_SAtlantic","may_Tau_SAtlantic","jun_Tau_SAtlantic","jul_Tau_SAtlantic","aug_Tau_SAtlantic","sep_Tau_SAtlantic","oct_Tau_SAtlantic","nov_Tau_SAtlantic","dec_Tau_SAtlantic","ann_Tau_SAtlantic","seasonal_amplitude_LWP_SPacific","seasonal_amplitude_Re_SPacific","seasonal_amplitude_Tau_SPacific","jan_LWP_SPacific","feb_LWP_SPacific","mar_LWP_SPacific","apr_LWP_SPacific","may_LWP_SPacific","jun_LWP_SPacific","jul_LWP_SPacific","aug_LWP_SPacific","sep_LWP_SPacific","oct_LWP_SPacific","nov_LWP_SPacific","dec_LWP_SPacific","ann_LWP_SPacific","jan_Re_SPacific","feb_Re_SPacific","mar_Re_SPacific","apr_Re_SPacific","may_Re_SPacific","jun_Re_SPacific","jul_Re_SPacific","aug_Re_SPacific","sep_Re_SPacific","oct_Re_SPacific","nov_Re_SPacific","dec_Re_SPacific","ann_Re_SPacific","jan_Tau_SPacific","feb_Tau_SPacific","mar_Tau_SPacific","apr_Tau_SPacific","may_Tau_SPacific","jun_Tau_SPacific","jul_Tau_SPacific","aug_Tau_SPacific","sep_Tau_SPacific","oct_Tau_SPacific","nov_Tau_SPacific","dec_Tau_SPacific","ann_Tau_SPacific","seasonal_amplitude_LWP_NPacific","seasonal_amplitude_Nd_NPacific","seasonal_amplitude_Re_NPacific","seasonal_amplitude_Tau_NPacific","jan_LWP_NPacific","feb_LWP_NPacific","mar_LWP_NPacific","apr_LWP_NPacific","may_LWP_NPacific","jun_LWP_NPacific","jul_LWP_NPacific","aug_LWP_NPacific","sep_LWP_NPacific","oct_LWP_NPacific","nov_LWP_NPacific","dec_LWP_NPacific","ann_LWP_NPacific","jan_Re_NPacific","feb_Re_NPacific","mar_Re_NPacific","apr_Re_NPacific","may_Re_NPacific","jun_Re_NPacific","jul_Re_NPacific","aug_Re_NPacific","sep_Re_NPacific","oct_Re_NPacific","nov_Re_NPacific","dec_Re_NPacific","ann_Re_NPacific","jan_Tau_NPacific","feb_Tau_NPacific","mar_Tau_NPacific","apr_Tau_NPacific","may_Tau_NPacific","jun_Tau_NPacific","jul_Tau_NPacific","aug_Tau_NPacific","sep_Tau_NPacific","oct_Tau_NPacific","nov_Tau_NPacific","dec_Tau_NPacific","ann_Tau_NPacific","jan_Nd_NPacific","feb_Nd_NPacific","mar_Nd_NPacific","apr_Nd_NPacific","may_Nd_NPacific","jun_Nd_NPacific","jul_Nd_NPacific","aug_Nd_NPacific","sep_Nd_NPacific","oct_Nd_NPacific","nov_Nd_NPacific","dec_Nd_NPacific","ann_Nd_NPacific","seasonal_amplitude_Fc_NAtlantic","seasonal_amplitude_Re_NAtlantic","jan_Nd_NAtlantic","nov_Nd_NAtlantic","dec_Nd_NAtlantic","feb_Re_NAtlantic","mar_Re_NAtlantic","nov_Re_NAtlantic","jan_LWP_NAtlantic","feb_LWP_NAtlantic","mar_LWP_NAtlantic","apr_LWP_NAtlantic","may_LWP_NAtlantic","jun_LWP_NAtlantic","jul_LWP_NAtlantic","aug_LWP_NAtlantic","sep_LWP_NAtlantic","oct_LWP_NAtlantic","nov_LWP_NAtlantic","dec_LWP_NAtlantic","ann_LWP_NAtlantic","jan_Tau_NAtlantic","feb_Tau_NAtlantic","mar_Tau_NAtlantic","apr_Tau_NAtlantic","may_Tau_NAtlantic","jun_Tau_NAtlantic","jul_Tau_NAtlantic","aug_Tau_NAtlantic","sep_Tau_NAtlantic","oct_Tau_NAtlantic","nov_Tau_NAtlantic","dec_Tau_NAtlantic","ann_Tau_NAtlantic"]

constraint_variable_list_Hd=['jan_Hd','feb_Hd','mar_Hd','apr_Hd','may_Hd','jun_Hd','jul_Hd','aug_Hd','sep_Hd','oct_Hd','nov_Hd','dec_Hd','ann_Hd']

######################################################
## Revise array order so that each region can be treated distinctly
## Put hemispheric difference first
######################################################

region_based_order=np.hstack((215,0,1,2,3,7,8,9,10,11,12,13,14,15,16,17,18,19,20,4,405,21,22,23,24,25,26,27,28,29,406,407,30,403,31,32,33,34,35,36,37,38,39,40,41,42,43,404,44,408,409,45,46,47,48,49,50,51,410,52,53,6,411,412,413,414,415,416,417,418,419,420,421,422,423,5,424,425,426,427,428,429,430,431,432,433,434,435,436,54,55,56,57,58,59,61,62,63,64,65,66,67,68,69,70,71,72,73,74,348,390,391,392,393,394,395,396,397,398,399,400,401,402,60,75,76,77,78,79,80,81,82,83,84,85,86,87,349,364,365,366,367,368,369,370,371,372,373,374,375,376,347,351,352,353,354,355,356,357,358,359,360,361,362,363,350,377,378,379,380,381,382,383,384,385,386,387,388,389,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,259,151,152,153,154,155,156,157,158,263,264,159,160,161,136,162,163,164,165,166,167,168,169,170,171,172,265,173,260,279,280,281,282,283,284,285,286,287,288,289,290,291,261,266,267,268,269,270,271,272,273,274,275,276,277,278,262,292,293,294,295,296,297,298,299,300,301,302,303,304,88,91,92,93,94,95,96,97,98,99,100,101,102,103,104,89,105,106,107,108,109,110,111,112,113,114,115,116,117,90,118,119,120,121,122,123,124,125,126,127,128,129,130,306,321,322,323,324,325,326,327,328,329,330,331,332,333,305,308,309,310,311,312,313,314,315,316,317,318,319,320,307,334,335,336,337,338,339,340,341,342,343,344,345,346,174,176,177,178,179,180,181,182,183,184,185,186,187,188,175,189,190,191,192,193,194,195,196,197,198,199,200,201,217,202,203,204,205,206,207,208,209,210,211,212,213,214,218,233,234,235,236,237,238,239,240,241,242,243,244,245,216,220,221,222,223,224,225,226,227,228,229,230,231,232,219,246,247,248,249,250,251,252,253,254,255,256,257,258))



region_list_sorted=[]
constraint_variable_text_list_sorted=[]

for ivar in np.arange(nvar-13):## accounting for 13 Hd months
    region_list_sorted.append(region_list_unsorted[region_based_order[ivar]])
    constraint_variable_text_list_sorted.append(constraint_variable_text_list_unsorted[region_based_order[ivar]])


unconstrained_NRMSE_array=np.zeros((nvar,sample_size))

for ivar in np.arange(nvar):
    if ivar==0:
        ivar_resorted=region_based_order[ivar]
        print(np.str(ivar)+' '+np.str(ivar_resorted)+': NRMSE_value_full_sample_no_zero_'+np.str(ivar_resorted)+'_additional.dat '+region_list_sorted[ivar]+' '+constraint_variable_text_list_sorted[ivar])
        unconstrained_NRMSE_array[ivar,:]=np.loadtxt(NRMSE_path_all+'NRMSE_value_full_sample_no_zero_'+np.str(ivar_resorted)+'_additional.dat')
    elif ivar<=13:
        print(np.str(ivar)+': NRMSE_value_full_sample_no_zero_'+np.str(ivar)+'_additional.dat for '+constraint_variable_list_Hd[ivar-1])
        unconstrained_NRMSE_array[ivar,:]=np.loadtxt(NRMSE_path_w_Hd+'NRMSE_value_full_sample_no_zero_'+np.str(ivar)+'_additional.dat')
    else:
        ivar_resorted=region_based_order[ivar-13]## accounting for inserted Hd months
        print(np.str(ivar)+' '+np.str(ivar_resorted)+': NRMSE_value_full_sample_no_zero_'+np.str(ivar_resorted)+'_additional.dat '+region_list_sorted[ivar-13]+' '+constraint_variable_text_list_sorted[ivar-13])
        unconstrained_NRMSE_array[ivar,:]=np.loadtxt(NRMSE_path_all+'NRMSE_value_full_sample_no_zero_'+np.str(ivar_resorted)+'_additional.dat')

unconstrained_NRMSE_means=unconstrained_NRMSE_array.mean(axis=1)

NRMSE_change_array=np.zeros((nvar,nvar))



constraint_integers_as_real_array=np.empty(nvar,dtype=object) ## completely empty array of objects

for ivar in np.arange(nvar):
    if ivar==0:
        ivar_resorted=region_based_order[ivar]
        print(np.str(ivar)+' '+np.str(ivar_resorted)+': '+region_list_sorted[ivar]+'_1M_python_indices_retained_w_0pc_error_for_variable_'+np.str(ivar_resorted)+'_additional.dat '+constraint_variable_text_list_sorted[ivar])
        constraint_integers_as_real_array[ivar]=np.loadtxt(indices_path+region_list_sorted[ivar]+'_1M_python_indices_retained_w_0pc_error_for_variable_'+np.str(ivar_resorted)+'_additional.dat')
    elif ivar<=13:
        print(np.str(ivar)+': '+'Hd_1M_python_indices_retained_w_0pc_error_for_variable_'+np.str(ivar)+'_additional.dat '+constraint_variable_list_Hd[ivar-1])
        constraint_integers_as_real_array[ivar]=np.loadtxt(indices_path_w_Hd+'Hd_1M_python_indices_retained_w_0pc_error_for_variable_'+np.str(ivar)+'_additional.dat')
    else:
        ivar_resorted=region_based_order[ivar-13]
        print(np.str(ivar)+' '+np.str(ivar_resorted)+': '+region_list_sorted[ivar-13]+'_1M_python_indices_retained_w_0pc_error_for_variable_'+np.str(ivar_resorted)+'_additional.dat '+constraint_variable_text_list_sorted[ivar-13])
        constraint_integers_as_real_array[ivar]=np.loadtxt(indices_path+region_list_sorted[ivar-13]+'_1M_python_indices_retained_w_0pc_error_for_variable_'+np.str(ivar_resorted)+'_additional.dat')


pc_change_array=np.zeros((nvar,nvar))
for ivar_dim1 in np.arange(nvar): ## The constrainer
    for ivar_dim2 in np.arange(nvar): ## the constrained
        print(np.str(ivar_dim1)+', '+np.str(ivar_dim2))
        integers_for_constraint=constraint_integers_as_real_array[ivar_dim1].astype(int)
        constrained_NRMSE_local= [unconstrained_NRMSE_array[ivar_dim2,i] for i in integers_for_constraint[:]]
        constrained_mean_local=np.mean(constrained_NRMSE_local)
        pc_change=(constrained_mean_local/unconstrained_NRMSE_means[ivar]-1)*100
        NRMSE_change_array[ivar_dim1,ivar_dim2]=pc_change



### Write to file for later use

np.savetxt('/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/NRMSE_values_entire_sample_w_Hd/cross_constraint_effects_resorted.csv', NRMSE_change_array, delimiter=',')
