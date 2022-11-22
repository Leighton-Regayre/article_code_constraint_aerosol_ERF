#########################################################################################
## Progressively adding constraint variables to achieve optimal constraint
###########################################################################################

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

restart_partial=False
region='all' ## only option here

if (restart_partial):
    region='all'

## restart folder an option for restarting the process part-way through if dropped from HPC. Replace 'XXX' with appropriate values

restart_folder='/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/finesse_reverse/all/'

#only_one=True ## case where only 1 index rejected so far
#pc_reduction_90_updated=XXX
#restart_indices_rejected_file=restart_folder+'all_sample_indices_py_index_FINESSE_w_uncertainty_0pc_error_constrained_ERF_XXX_to_XXX_with_XXXpc_reduction_from_XXX_constraints_INDICES_OF_REJECTED_w_o_Hd_regional.txt'
#restart_indices_retained_file=restart_folder+'all_sample_indices_py_index_FINESSE_w_uncertainty_0pc_error_constrained_ERF_XXX_to_XXX_with_XXXpc_reduction_from_XXX_constraints_INDICES_OF_RETAINED_w_o_Hd_regional.txt'

#only_one=False ## case where multiple indices already rejected
#pc_reduction_90_updated=XXX
#restart_indices_rejected_file=restart_folder+'all_sample_indices_py_index_FINESSE_w_uncertainty_0pc_error_constrained_ERF_XXX_to_XXX_with_XXXpc_reduction_from_XXX_constraints_INDICES_OF_REJECTED_w_o_Hd_regional.txt'
#restart_indices_retained_file=restart_folder+'all_sample_indices_py_index_FINESSE_w_uncertainty_0pc_error_constrained_ERF_XXX_to_XXX_with_XXXpc_reduction_from_XXX_constraints_INDICES_OF_RETAINED_w_o_Hd_regional.txt'


if restart_partial:
    if only_one:
        icon_reject_updated_decimal=np.zeros((1))
        icon_reject_updated_decimal[0]=np.loadtxt(restart_indices_rejected_file)
    else:
        icon_reject_updated_decimal=np.loadtxt(restart_indices_rejected_file)
    icon_keep_updated_decimal=np.loadtxt(restart_indices_retained_file)
    constraint_TF=True
    icon_reject_updated=[]
    icon_keep_updated=[]
    for icon in np.arange(len(icon_reject_updated_decimal)):
        icon_reject_updated.append(np.int(icon_reject_updated_decimal[icon]))
    for icon in np.arange(len(icon_keep_updated_decimal)):
        icon_keep_updated.append(np.int(icon_keep_updated_decimal[icon]))


representation_error=0 ## spatial and temporal could be included if well-defined
obs_error=0 ## low/unquantified b/c using MODIS-COSP diagnostics
total_error_percent=representation_error+obs_error
sample_size=1000000

sample_type_text='1000000_w_o_carb.dat'
constraint_percent=0.5 ## retained percentage of original 1 million. Defined to match Johnson et al., 2020 and Regayre et al. 2020 values from constraints w/ well-defined errors
optimal_path_out='/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/finesse_reverse/'+region+'/'
obs_dir_base='/gws/nopw/j04/acure/lregayre/data_for_emulation/'


## constraints in revised order for this process

constraint_variable_text_list=['seasonal_amplitude_Nd_hemispheric_diff','jan_Hd','feb_Hd','mar_Hd','apr_Hd','may_Hd','jun_Hd','jul_Hd','aug_Hd','sep_Hd','oct_Hd','nov_Hd','dec_Hd','ann_Hd','jul_transect_Nd_NAtlantic','jul_transect_LWP_Nd_NAtlantic','jul_transect_Re_Nd_NAtlantic','jul_transect_Re_NAtlantic','seasonal_amplitude_Fsw_NAtlantic','jan_Fsw_NAtlantic','feb_Fsw_NAtlantic','mar_Fsw_NAtlantic','apr_Fsw_NAtlantic','may_Fsw_NAtlantic','jun_Fsw_NAtlantic','jul_Fsw_NAtlantic','aug_Fsw_NAtlantic','sep_Fsw_NAtlantic','oct_Fsw_NAtlantic','nov_Fsw_NAtlantic','dec_Fsw_NAtlantic','ann_Fsw_NAtlantic','seasonal_amplitude_Nd_NAtlantic','jan_Nd_NAtlantic','feb_Nd_NAtlantic','mar_Nd_NAtlantic','apr_Nd_NAtlantic','may_Nd_NAtlantic','jun_Nd_NAtlantic','jul_Nd_NAtlantic','aug_Nd_NAtlantic','sep_Nd_NAtlantic','oct_Nd_NAtlantic','nov_Nd_NAtlantic','dec_Nd_NAtlantic','ann_Nd_NAtlantic','seasonal_amplitude_Fc_NAtlantic','jan_Fc_NAtlantic','feb_Fc_NAtlantic','mar_Fc_NAtlantic','apr_Fc_NAtlantic','may_Fc_NAtlantic','jun_Fc_NAtlantic','jul_Fc_NAtlantic','aug_Fc_NAtlantic','sep_Fc_NAtlantic','oct_Fc_NAtlantic','nov_Fc_NAtlantic','dec_Fc_NAtlantic','ann_Fc_NAtlantic','seasonal_amplitude_Re_NAtlantic','jan_Re_NAtlantic','feb_Re_NAtlantic','mar_Re_NAtlantic','apr_Re_NAtlantic','may_Re_NAtlantic','jun_Re_NAtlantic','jul_Re_NAtlantic','aug_Re_NAtlantic','sep_Re_NAtlantic','oct_Re_NAtlantic','nov_Re_NAtlantic','dec_Re_NAtlantic','ann_Re_NAtlantic','seasonal_amplitude_LWP_NAtlantic',"jan_LWP_NAtlantic","feb_LWP_NAtlantic","mar_LWP_NAtlantic","apr_LWP_NAtlantic","may_LWP_NAtlantic","jun_LWP_NAtlantic","jul_LWP_NAtlantic","aug_LWP_NAtlantic","sep_LWP_NAtlantic","oct_LWP_NAtlantic","nov_LWP_NAtlantic","dec_LWP_NAtlantic","ann_LWP_NAtlantic",'seasonal_amplitude_Tau_NAtlantic',"jan_Tau_NAtlantic","feb_Tau_NAtlantic","mar_Tau_NAtlantic","apr_Tau_NAtlantic","may_Tau_NAtlantic","jun_Tau_NAtlantic","jul_Tau_NAtlantic","aug_Tau_NAtlantic","sep_Tau_NAtlantic","oct_Tau_NAtlantic","nov_Tau_NAtlantic","dec_Tau_NAtlantic","ann_Tau_NAtlantic",'jul_transect_Nd_NPacific','jul_transect_Cf_NPacific','jul_transect_LWP_NPacific','jul_transect_Cf_Nd_NPacific','jul_transect_LWP_Nd_NPacific','jul_transect_Re_Nd_NPacific','seasonal_amplitude_Fsw_NPacific','jan_Fsw_NPacific','feb_Fsw_NPacific','mar_Fsw_NPacific','apr_Fsw_NPacific','may_Fsw_NPacific','jun_Fsw_NPacific','jul_Fsw_NPacific','aug_Fsw_NPacific','sep_Fsw_NPacific','oct_Fsw_NPacific','nov_Fsw_NPacific','dec_Fsw_NPacific','ann_Fsw_NPacific',"seasonal_amplitude_Nd_NPacific","jan_Nd_NPacific","feb_Nd_NPacific","mar_Nd_NPacific","apr_Nd_NPacific","may_Nd_NPacific","jun_Nd_NPacific","jul_Nd_NPacific","aug_Nd_NPacific","sep_Nd_NPacific","oct_Nd_NPacific","nov_Nd_NPacific","dec_Nd_NPacific","ann_Nd_NPacific",'seasonal_amplitude_Fc_NPacific','jan_Fc_NPacific','feb_Fc_NPacific','mar_Fc_NPacific','apr_Fc_NPacific','may_Fc_NPacific','jun_Fc_NPacific','jul_Fc_NPacific','aug_Fc_NPacific','sep_Fc_NPacific','oct_Fc_NPacific','nov_Fc_NPacific','dec_Fc_NPacific','ann_Fc_NPacific',"seasonal_amplitude_Re_NPacific","jan_Re_NPacific","feb_Re_NPacific","mar_Re_NPacific","apr_Re_NPacific","may_Re_NPacific","jun_Re_NPacific","jul_Re_NPacific","aug_Re_NPacific","sep_Re_NPacific","oct_Re_NPacific","nov_Re_NPacific","dec_Re_NPacific","ann_Re_NPacific","seasonal_amplitude_LWP_NPacific","jan_LWP_NPacific","feb_LWP_NPacific","mar_LWP_NPacific","apr_LWP_NPacific","may_LWP_NPacific","jun_LWP_NPacific","jul_LWP_NPacific","aug_LWP_NPacific","sep_LWP_NPacific","oct_LWP_NPacific","nov_LWP_NPacific","dec_LWP_NPacific","ann_LWP_NPacific","seasonal_amplitude_Tau_NPacific","jan_Tau_NPacific","feb_Tau_NPacific","mar_Tau_NPacific","apr_Tau_NPacific","may_Tau_NPacific","jun_Tau_NPacific","jul_Tau_NPacific","aug_Tau_NPacific","sep_Tau_NPacific","oct_Tau_NPacific","nov_Tau_NPacific","dec_Tau_NPacific","ann_Tau_NPacific",'nov_transect_Nd_SAtlantic','nov_transect_LWP_Nd_SAtlantic','nov_transect_Fc_SAtlantic','nov_transect_AI_SAtlantic','nov_transect_Nd_AI_SAtlantic','seasonal_amplitude_Fsw_SAtlantic','jan_Fsw_SAtlantic','feb_Fsw_SAtlantic','mar_Fsw_SAtlantic','apr_Fsw_SAtlantic','may_Fsw_SAtlantic','jun_Fsw_SAtlantic','jul_Fsw_SAtlantic','aug_Fsw_SAtlantic','sep_Fsw_SAtlantic','oct_Fsw_SAtlantic','nov_Fsw_SAtlantic','dec_Fsw_SAtlantic','ann_Fsw_SAtlantic',"seasonal_amplitude_Nd_SAtlantic",'jan_Nd_SAtlantic','feb_Nd_SAtlantic','mar_Nd_SAtlantic','apr_Nd_SAtlantic','may_Nd_SAtlantic','jun_Nd_SAtlantic','jul_Nd_SAtlantic','aug_Nd_SAtlantic','sep_Nd_SAtlantic','oct_Nd_SAtlantic','nov_Nd_SAtlantic','dec_Nd_SAtlantic','ann_Nd_SAtlantic','seasonal_amplitude_Fc_SAtlantic','jan_Fc_SAtlantic','feb_Fc_SAtlantic','mar_Fc_SAtlantic','apr_Fc_SAtlantic','may_Fc_SAtlantic','jun_Fc_SAtlantic','jul_Fc_SAtlantic','aug_Fc_SAtlantic','sep_Fc_SAtlantic','oct_Fc_SAtlantic','nov_Fc_SAtlantic','dec_Fc_SAtlantic','ann_Fc_SAtlantic',"seasonal_amplitude_Re_SAtlantic","jan_Re_SAtlantic","feb_Re_SAtlantic","mar_Re_SAtlantic","apr_Re_SAtlantic","may_Re_SAtlantic","jun_Re_SAtlantic","jul_Re_SAtlantic","aug_Re_SAtlantic","sep_Re_SAtlantic","oct_Re_SAtlantic","nov_Re_SAtlantic","dec_Re_SAtlantic","ann_Re_SAtlantic","seasonal_amplitude_LWP_SAtlantic","jan_LWP_SAtlantic","feb_LWP_SAtlantic","mar_LWP_SAtlantic","apr_LWP_SAtlantic","may_LWP_SAtlantic","jun_LWP_SAtlantic","jul_LWP_SAtlantic","aug_LWP_SAtlantic","sep_LWP_SAtlantic","oct_LWP_SAtlantic","nov_LWP_SAtlantic","dec_LWP_SAtlantic","ann_LWP_SAtlantic","seasonal_amplitude_Tau_SAtlantic","jan_Tau_SAtlantic","feb_Tau_SAtlantic","mar_Tau_SAtlantic","apr_Tau_SAtlantic","may_Tau_SAtlantic","jun_Tau_SAtlantic","jul_Tau_SAtlantic","aug_Tau_SAtlantic","sep_Tau_SAtlantic","oct_Tau_SAtlantic","nov_Tau_SAtlantic","dec_Tau_SAtlantic","ann_Tau_SAtlantic",'nov_transect_Nd_SPacific','seasonal_amplitude_Fsw_SPacific','jan_Fsw_SPacific','feb_Fsw_SPacific','mar_Fsw_SPacific','apr_Fsw_SPacific','may_Fsw_SPacific','jun_Fsw_SPacific','jul_Fsw_SPacific','aug_Fsw_SPacific','sep_Fsw_SPacific','oct_Fsw_SPacific','nov_Fsw_SPacific','dec_Fsw_SPacific','ann_Fsw_SPacific','seasonal_amplitude_Nd_SPacific','jan_Nd_SPacific','feb_Nd_SPacific','mar_Nd_SPacific','apr_Nd_SPacific','may_Nd_SPacific','jun_Nd_SPacific','jul_Nd_SPacific','aug_Nd_SPacific','sep_Nd_SPacific','oct_Nd_SPacific','nov_Nd_SPacific','dec_Nd_SPacific','ann_Nd_SPacific','seasonal_amplitude_Fc_SPacific','jan_Fc_SPacific','feb_Fc_SPacific','mar_Fc_SPacific','apr_Fc_SPacific','may_Fc_SPacific','jun_Fc_SPacific','jul_Fc_SPacific','aug_Fc_SPacific','sep_Fc_SPacific','oct_Fc_SPacific','nov_Fc_SPacific','dec_Fc_SPacific','ann_Fc_SPacific',"seasonal_amplitude_Re_SPacific","jan_Re_SPacific","feb_Re_SPacific","mar_Re_SPacific","apr_Re_SPacific","may_Re_SPacific","jun_Re_SPacific","jul_Re_SPacific","aug_Re_SPacific","sep_Re_SPacific","oct_Re_SPacific","nov_Re_SPacific","dec_Re_SPacific","ann_Re_SPacific","seasonal_amplitude_LWP_SPacific","jan_LWP_SPacific","feb_LWP_SPacific","mar_LWP_SPacific","apr_LWP_SPacific","may_LWP_SPacific","jun_LWP_SPacific","jul_LWP_SPacific","aug_LWP_SPacific","sep_LWP_SPacific","oct_LWP_SPacific","nov_LWP_SPacific","dec_LWP_SPacific","ann_LWP_SPacific","seasonal_amplitude_Tau_SPacific","jan_Tau_SPacific","feb_Tau_SPacific","mar_Tau_SPacific","apr_Tau_SPacific","may_Tau_SPacific","jun_Tau_SPacific","jul_Tau_SPacific","aug_Tau_SPacific","sep_Tau_SPacific","oct_Tau_SPacific","nov_Tau_SPacific","dec_Tau_SPacific","ann_Tau_SPacific",'seasonal_amplitude_Fsw_SOcean','jan_Fsw_SOcean','feb_Fsw_SOcean','mar_Fsw_SOcean','apr_Fsw_SOcean','may_Fsw_SOcean','jun_Fsw_SOcean','jul_Fsw_SOcean','aug_Fsw_SOcean','sep_Fsw_SOcean','oct_Fsw_SOcean','nov_Fsw_SOcean','dec_Fsw_SOcean','ann_Fsw_SOcean','seasonal_amplitude_Nd_SOcean','jan_Nd_SOcean','feb_Nd_SOcean','mar_Nd_SOcean','apr_Nd_SOcean','may_Nd_SOcean','jun_Nd_SOcean','jul_Nd_SOcean','aug_Nd_SOcean','sep_Nd_SOcean','oct_Nd_SOcean','nov_Nd_SOcean','dec_Nd_SOcean','ann_Nd_SOcean',"seasonal_amplitude_Fc_SOcean",'jan_Fc_SOcean','feb_Fc_SOcean','mar_Fc_SOcean','apr_Fc_SOcean','may_Fc_SOcean','jun_Fc_SOcean','jul_Fc_SOcean','aug_Fc_SOcean','sep_Fc_SOcean','oct_Fc_SOcean','nov_Fc_SOcean','dec_Fc_SOcean','ann_Fc_SOcean',"seasonal_amplitude_Re_SOcean","jan_Re_SOcean","feb_Re_SOcean","mar_Re_SOcean","apr_Re_SOcean","may_Re_SOcean","jun_Re_SOcean","jul_Re_SOcean","aug_Re_SOcean","sep_Re_SOcean","oct_Re_SOcean","nov_Re_SOcean","dec_Re_SOcean","ann_Re_SOcean","seasonal_amplitude_LWP_SOcean","jan_LWP_SOcean","feb_LWP_SOcean","mar_LWP_SOcean","apr_LWP_SOcean","may_LWP_SOcean","jun_LWP_SOcean","jul_LWP_SOcean","aug_LWP_SOcean","sep_LWP_SOcean","oct_LWP_SOcean","nov_LWP_SOcean","dec_LWP_SOcean","ann_LWP_SOcean","seasonal_amplitude_Tau_SOcean","jan_Tau_SOcean","feb_Tau_SOcean","mar_Tau_SOcean","apr_Tau_SOcean","may_Tau_SOcean","jun_Tau_SOcean","jul_Tau_SOcean","aug_Tau_SOcean","sep_Tau_SOcean","oct_Tau_SOcean","nov_Tau_SOcean","dec_Tau_SOcean","ann_Tau_SOcean"]


ncon_full=len(constraint_variable_text_list)

greek_letters_lower=[chr(code) for code in range(945,970)]
TAU=greek_letters_lower[19]
greek_letters_upper=[chr(code) for code in range(913,938)]
DELTA=greek_letters_upper[3]


###################################################
## observation data read in same order as constraint names
###################################################


obs_set=np.zeros((ncon_full))

## Hd
obs_set[0]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/CDNC_NH_SH_difference_seasonal_amplitude_observed.dat')
obs_set[1]=np.loadtxt(obs_dir_base+'/observations/CDNC_NH_SH_difference_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1] ## jan
obs_set[2]=np.loadtxt(obs_dir_base+'/observations/CDNC_NH_SH_difference_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2] ## feb
obs_set[3]=np.loadtxt(obs_dir_base+'/observations/CDNC_NH_SH_difference_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3] ## mar
obs_set[4]=np.loadtxt(obs_dir_base+'/observations/CDNC_NH_SH_difference_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4] ## apr
obs_set[5]=np.loadtxt(obs_dir_base+'/observations/CDNC_NH_SH_difference_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5] ## may
obs_set[6]=np.loadtxt(obs_dir_base+'/observations/CDNC_NH_SH_difference_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6] ## jun
obs_set[7]=np.loadtxt(obs_dir_base+'/observations/CDNC_NH_SH_difference_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7] ## jul
obs_set[8]=np.loadtxt(obs_dir_base+'/observations/CDNC_NH_SH_difference_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8] ## aug
obs_set[9]=np.loadtxt(obs_dir_base+'/observations/CDNC_NH_SH_difference_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9] ## sep
obs_set[10]=np.loadtxt(obs_dir_base+'/observations/CDNC_NH_SH_difference_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10] ## oct
obs_set[11]=np.loadtxt(obs_dir_base+'/observations/CDNC_NH_SH_difference_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11] ## nov
obs_set[12]=np.loadtxt(obs_dir_base+'/observations/CDNC_NH_SH_difference_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0] ## dec
obs_set[13]=np.loadtxt(obs_dir_base+'/observations/CDNC_NH_SH_difference_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12] ## ann

###############
## NAtlantic
###############

## NAtlantic transects
obs_set[14]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/jul2017/CDNC_European_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[15]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/jul2017/LWP_CDNC_European_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[16]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/jul2017/Reff_CDNC_European_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[17]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/jul2017/Reff_European_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')

# NAtlantic Fsw
obs_set[18]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/CERES_SW_TOA_NN_atlantic_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[19]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[1] ## jan w/ 0=dec
obs_set[20]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[2] ## feb w/ 0=dec
obs_set[21]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[3] ## mar w/ 0=dec
obs_set[22]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[4] ## apr w/ 0=dec
obs_set[23]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[5] ## may w/ 0=dec
obs_set[24]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[6] ## jun w/ 0=dec
obs_set[25]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[7] ## jul w/ 0=dec
obs_set[26]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[8] ## aug w/ 0=dec
obs_set[27]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[9] ## sep w/ 0=dec
obs_set[28]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[10] ## oct w/ 0=dec
obs_set[29]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[11] ## nov w/ 0=dec
obs_set[30]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[0] ## dec w/ 0=dec
obs_set[31]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NN_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[12] ## ann w/ 0=dec

# NAtlantic Nd
obs_set[32]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/CDNC_NN_atlantic_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[33]=np.loadtxt(obs_dir_base+'/observations/CDNC_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1] # Jan
obs_set[34]=np.loadtxt(obs_dir_base+'/observations/CDNC_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2] ## feb w/ 0=dec
obs_set[35]=np.loadtxt(obs_dir_base+'/observations/CDNC_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3] ##  w/ 0=dec
obs_set[36]=np.loadtxt(obs_dir_base+'/observations/CDNC_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4] ##  w/ 0=dec
obs_set[37]=np.loadtxt(obs_dir_base+'/observations/CDNC_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5] ##  w/ 0=dec
obs_set[38]=np.loadtxt(obs_dir_base+'/observations/CDNC_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6] ##  w/ 0=dec
obs_set[39]=np.loadtxt(obs_dir_base+'/observations/CDNC_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7] ##  w/ 0=dec
obs_set[40]=np.loadtxt(obs_dir_base+'/observations/CDNC_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8] ##  w/ 0=dec
obs_set[41]=np.loadtxt(obs_dir_base+'/observations/CDNC_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9] ##  w/ 0=dec
obs_set[42]=np.loadtxt(obs_dir_base+'/observations/CDNC_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10] ##  w/ 0=dec
obs_set[43]=np.loadtxt(obs_dir_base+'/observations/CDNC_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11] #nov
obs_set[44]=np.loadtxt(obs_dir_base+'/observations/CDNC_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0] #dec
obs_set[45]=np.loadtxt(obs_dir_base+'/observations/CDNC_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12] ##  w/ 0=dec

# NAtlantic Fc
obs_set[46]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[47]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]
obs_set[48]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]
obs_set[49]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]
obs_set[50]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]
obs_set[51]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]
obs_set[52]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]
obs_set[53]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]
obs_set[54]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]
obs_set[55]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]
obs_set[56]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]
obs_set[57]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]
obs_set[58]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0] # dec
obs_set[59]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12] #ann

# NAtlantic Re
obs_set[60]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Particle_Size_Liquid_NN_atlantic_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[61]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]/1000000 #jan
obs_set[62]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]/1000000 
obs_set[63]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]/1000000 
obs_set[64]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]/1000000 #
obs_set[65]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]/1000000 #
obs_set[66]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]/1000000 #
obs_set[67]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]/1000000 #
obs_set[68]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]/1000000 #
obs_set[69]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]/1000000 
obs_set[70]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]/1000000 
obs_set[71]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]/1000000 
obs_set[72]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]/1000000 #dec
obs_set[73]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]/1000000 #ann


# NAtlantic LWP
obs_set[74]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Water_Path_Liquid_NN_atlantic_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[75]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]/1000
obs_set[76]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]/1000
obs_set[77]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]/1000
obs_set[78]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]/1000
obs_set[79]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]/1000
obs_set[80]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]/1000
obs_set[81]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]/1000
obs_set[82]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]/1000
obs_set[83]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]/1000
obs_set[84]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]/1000
obs_set[85]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]/1000
obs_set[86]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]/1000
obs_set[87]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]/1000


# NAtlantic Tc
obs_set[88]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Optical_Thickness_Liquid_NN_atlantic_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[89]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1] 
obs_set[90]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2] 
obs_set[91]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3] 
obs_set[92]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4] 
obs_set[93]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5] 
obs_set[94]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6] 
obs_set[95]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7] 
obs_set[96]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8] 
obs_set[97]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9] 
obs_set[98]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10] 
obs_set[99]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11] 
obs_set[100]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0] 
obs_set[101]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NN_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12] 


#################
## NPacific
#################

## NPacific transects
obs_set[102]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/jul2017/CDNC_NorthAmerican_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[103]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/jul2017/LCF_NorthAmerican_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[104]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/jul2017/LWP_NorthAmerican_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[105]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/jul2017/LCF_CDNC_NorthAmerican_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[106]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/jul2017/LWP_CDNC_NorthAmerican_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[107]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/jul2017/Reff_CDNC_NorthAmerican_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[108]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/CERES_SW_TOA_NE_pacific_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[109]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[1] ## jan w/ 0=dec
obs_set[110]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[2]
obs_set[111]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[3]
obs_set[112]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[4]
obs_set[113]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[5]
obs_set[114]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[6]
obs_set[115]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[7]
obs_set[116]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[8]
obs_set[117]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[9]
obs_set[118]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[10]
obs_set[119]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[11]
obs_set[120]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[0] #dec
obs_set[121]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_NE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[12] # ann


## NPacific Nd
obs_set[122]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/CDNC_NE_pacific_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[123]=np.loadtxt(obs_dir_base+'/observations/CDNC_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]
obs_set[124]=np.loadtxt(obs_dir_base+'/observations/CDNC_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]
obs_set[125]=np.loadtxt(obs_dir_base+'/observations/CDNC_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]
obs_set[126]=np.loadtxt(obs_dir_base+'/observations/CDNC_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]
obs_set[127]=np.loadtxt(obs_dir_base+'/observations/CDNC_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]
obs_set[128]=np.loadtxt(obs_dir_base+'/observations/CDNC_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]
obs_set[129]=np.loadtxt(obs_dir_base+'/observations/CDNC_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]
obs_set[130]=np.loadtxt(obs_dir_base+'/observations/CDNC_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]
obs_set[131]=np.loadtxt(obs_dir_base+'/observations/CDNC_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]
obs_set[132]=np.loadtxt(obs_dir_base+'/observations/CDNC_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]
obs_set[133]=np.loadtxt(obs_dir_base+'/observations/CDNC_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]
obs_set[134]=np.loadtxt(obs_dir_base+'/observations/CDNC_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]
obs_set[135]=np.loadtxt(obs_dir_base+'/observations/CDNC_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]


## NPacific Fc
obs_set[136]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Retrieval_Fraction_Liquid_NE_pacific_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[137]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1] # jan
obs_set[138]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]
obs_set[139]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]
obs_set[140]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]
obs_set[141]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]
obs_set[142]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]
obs_set[143]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]
obs_set[144]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]
obs_set[145]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]
obs_set[146]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]
obs_set[147]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]
obs_set[148]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]
obs_set[149]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12] # ann


## NPacific Re
obs_set[150]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Particle_Size_Liquid_NE_pacific_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[151]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]/1000000  # jan
obs_set[152]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]/1000000 
obs_set[153]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]/1000000 
obs_set[154]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]/1000000 
obs_set[155]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]/1000000 
obs_set[156]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]/1000000 
obs_set[157]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]/1000000 
obs_set[158]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]/1000000 
obs_set[159]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]/1000000 
obs_set[160]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]/1000000 
obs_set[161]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]/1000000 
obs_set[162]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]/1000000 
obs_set[163]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]/1000000 


## NPacific LWP
obs_set[164]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Water_Path_Liquid_NE_pacific_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[165]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]/1000 # jan
obs_set[166]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]/1000
obs_set[167]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]/1000
obs_set[168]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]/1000
obs_set[169]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]/1000
obs_set[170]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]/1000
obs_set[171]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]/1000
obs_set[172]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]/1000
obs_set[173]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]/1000
obs_set[174]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]/1000
obs_set[175]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]/1000
obs_set[176]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]/1000
obs_set[177]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]/1000



## NPacific Tc
obs_set[178]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Optical_Thickness_Liquid_NE_pacific_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[179]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]
obs_set[180]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]
obs_set[181]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]
obs_set[182]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]
obs_set[183]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]
obs_set[184]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]
obs_set[185]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]
obs_set[186]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]
obs_set[187]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]
obs_set[188]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]
obs_set[189]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]
obs_set[190]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]
obs_set[191]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_NE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]


##################
## SAtlantic
##################

## SAtlantic transects
obs_set[192]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/nov2017/CDNC_Namibian_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[193]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/nov2017/LWP_CDNC_Namibian_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[194]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/nov2017/LCF_Namibian_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[195]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/nov2017/AI_Namibian_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[196]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/nov2017/CDNC_AI_Namibian_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')


## SAtlantic Fsw
obs_set[197]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/CERES_SW_TOA_SE_atlantic_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[198]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[1] ## jan
obs_set[199]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[2] 
obs_set[200]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[3] 
obs_set[201]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[4] 
obs_set[202]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[5] 
obs_set[203]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[6] 
obs_set[204]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[7] 
obs_set[205]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[8] 
obs_set[206]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[9] 
obs_set[207]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[10] 
obs_set[208]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[11] 
obs_set[209]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[0] 
obs_set[210]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_atlantic_dec_to_ann_revised_match_MODIS_CF.dat')[12] ## ann

## SAtlantic Nd
obs_set[211]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/CDNC_SE_atlantic_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[212]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1] ## jan
obs_set[213]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]
obs_set[214]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]
obs_set[215]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]
obs_set[216]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]
obs_set[217]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]
obs_set[218]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]
obs_set[219]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8] ## aug
obs_set[220]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9] # sep
obs_set[221]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10] 
obs_set[222]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11] #nov
obs_set[223]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]
obs_set[224]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12] # ann

## SAtlantic Fc
obs_set[225]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[226]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1] ## Jan
obs_set[227]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2] #feb
obs_set[228]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3] #mar
obs_set[229]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4] #apr
obs_set[230]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]
obs_set[231]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]
obs_set[232]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]
obs_set[233]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]
obs_set[234]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]
obs_set[235]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]
obs_set[236]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11] #nov
obs_set[237]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0] # dec
obs_set[238]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12] # ann

## SAtlantic Re
obs_set[239]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Particle_Size_Liquid_SE_atlantic_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[240]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]/1000000 # jan
obs_set[241]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]/1000000
obs_set[242]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]/1000000
obs_set[243]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]/1000000
obs_set[244]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]/1000000
obs_set[245]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]/1000000
obs_set[246]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]/1000000
obs_set[247]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]/1000000
obs_set[248]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]/1000000
obs_set[249]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]/1000000
obs_set[250]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]/1000000
obs_set[251]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]/1000000 # dec
obs_set[252]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]/1000000 # ann

## SAtlantic LWP
obs_set[253]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Water_Path_Liquid_SE_atlantic_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[254]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]/1000 # jan
obs_set[255]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]/1000
obs_set[256]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]/1000 
obs_set[257]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]/1000
obs_set[258]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]/1000
obs_set[259]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]/1000
obs_set[260]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]/1000
obs_set[261]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]/1000
obs_set[262]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]/1000
obs_set[263]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]/1000
obs_set[264]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]/1000
obs_set[265]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]/1000 #dec
obs_set[266]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]/1000 # ann

## SAtlantic Tc
obs_set[267]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Optical_Thickness_Liquid_SE_atlantic_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[268]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1] # jan
obs_set[269]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]
obs_set[270]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]
obs_set[271]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]
obs_set[272]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]
obs_set[273]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]
obs_set[274]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]
obs_set[275]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]
obs_set[276]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]
obs_set[277]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]
obs_set[278]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]
obs_set[279]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]
obs_set[280]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_atlantic_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]


#############
## SPacific
#############

## SPacific transect and Fsw
obs_set[281]=np.loadtxt(obs_dir_base+'/observations/transect_gradients/nov2017/CDNC_SouthAmerican_trajectory_gradients_observed_meters_revised_MODIS_CF.dat')
obs_set[282]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/CERES_SW_TOA_SE_pacific_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[283]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[1] ## jan
obs_set[284]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[2] 
obs_set[285]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[3] 
obs_set[286]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[4] 
obs_set[287]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[5] 
obs_set[288]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[6] 
obs_set[289]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[7] 
obs_set[290]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[8] 
obs_set[291]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[9] 
obs_set[292]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[10] 
obs_set[293]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[11] 
obs_set[294]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[0] 
obs_set[295]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_SE_pacific_dec_to_ann_revised_match_MODIS_CF.dat')[12] 


## SPacific Nd
obs_set[296]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/CDNC_SE_pacific_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[297]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1] # jan
obs_set[298]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2] 
obs_set[299]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3] 
obs_set[300]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4] 
obs_set[301]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5] 
obs_set[302]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6] 
obs_set[303]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7] 
obs_set[304]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8] ## aug, w/0=dec
obs_set[305]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9] 
obs_set[306]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10] 
obs_set[307]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11] 
obs_set[308]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0] 
obs_set[309]=np.loadtxt(obs_dir_base+'/observations/CDNC_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12] 

## SPacific Fc
obs_set[310]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Retrieval_Fraction_Liquid_SE_pacific_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[311]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1] #J
obs_set[312]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]
obs_set[313]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]
obs_set[314]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]
obs_set[315]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]
obs_set[316]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]
obs_set[317]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]
obs_set[318]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]
obs_set[319]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]
obs_set[320]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]
obs_set[321]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]
obs_set[322]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]
obs_set[323]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]

## SPacific Re
obs_set[324]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Particle_Size_Liquid_SE_pacific_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[325]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]/1000000   ## jan
obs_set[326]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]/1000000 
obs_set[327]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]/1000000 
obs_set[328]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]/1000000 
obs_set[329]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]/1000000 
obs_set[330]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]/1000000 
obs_set[331]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]/1000000 
obs_set[332]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]/1000000 
obs_set[333]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]/1000000 
obs_set[334]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]/1000000 
obs_set[335]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]/1000000 
obs_set[336]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]/1000000 
obs_set[337]=np.loadtxt(obs_dir_base+'/observations/Cloud_Particle_Size_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]/1000000 

## SPacific LWP
obs_set[338]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Water_Path_Liquid_SE_pacific_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[339]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]/1000 # jan
obs_set[340]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]/1000
obs_set[341]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]/1000
obs_set[342]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]/1000
obs_set[343]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]/1000
obs_set[344]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]/1000
obs_set[345]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]/1000
obs_set[346]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]/1000
obs_set[347]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]/1000
obs_set[348]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]/1000
obs_set[349]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]/1000
obs_set[350]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]/1000
obs_set[351]=np.loadtxt(obs_dir_base+'/observations/Cloud_Water_Path_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]/1000

## SPacific Tc
obs_set[352]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Optical_Thickness_Liquid_SE_pacific_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[353]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]  ## jan
obs_set[354]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2] 
obs_set[355]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3] 
obs_set[356]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4] 
obs_set[357]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5] 
obs_set[358]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6] 
obs_set[359]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7] 
obs_set[360]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8] 
obs_set[361]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9] 
obs_set[362]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10] 
obs_set[363]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11] 
obs_set[364]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0] 
obs_set[365]=np.loadtxt(obs_dir_base+'/observations/Cloud_Optical_Thickness_Liquid_SE_pacific_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]



##########
## SOcean
##########

## SOcean Fsw
obs_set[366]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/CERES_SW_TOA_N_southern_ocean_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[367]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_N_southern_ocean_dec_to_ann_revised_match_MODIS_CF.dat')[1] ## jan
obs_set[368]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_N_southern_ocean_dec_to_ann_revised_match_MODIS_CF.dat')[2]
obs_set[369]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_N_southern_ocean_dec_to_ann_revised_match_MODIS_CF.dat')[3]
obs_set[370]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_N_southern_ocean_dec_to_ann_revised_match_MODIS_CF.dat')[4]
obs_set[371]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_N_southern_ocean_dec_to_ann_revised_match_MODIS_CF.dat')[5]
obs_set[372]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_N_southern_ocean_dec_to_ann_revised_match_MODIS_CF.dat')[6]
obs_set[373]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_N_southern_ocean_dec_to_ann_revised_match_MODIS_CF.dat')[7]
obs_set[374]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_N_southern_ocean_dec_to_ann_revised_match_MODIS_CF.dat')[8]
obs_set[375]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_N_southern_ocean_dec_to_ann_revised_match_MODIS_CF.dat')[9]
obs_set[376]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_N_southern_ocean_dec_to_ann_revised_match_MODIS_CF.dat')[10]
obs_set[377]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_N_southern_ocean_dec_to_ann_revised_match_MODIS_CF.dat')[11]
obs_set[378]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_N_southern_ocean_dec_to_ann_revised_match_MODIS_CF.dat')[0]
obs_set[379]=np.loadtxt(obs_dir_base+'/observations/CERES_SW_TOA_N_southern_ocean_dec_to_ann_revised_match_MODIS_CF.dat')[12] ## ann

## SOcean Nd
obs_set[380]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/CDNC_N_southern_ocean_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[381]=np.loadtxt(obs_dir_base+'/observations/CDNC_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1] ## jan
obs_set[382]=np.loadtxt(obs_dir_base+'/observations/CDNC_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]
obs_set[383]=np.loadtxt(obs_dir_base+'/observations/CDNC_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]
obs_set[384]=np.loadtxt(obs_dir_base+'/observations/CDNC_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]
obs_set[385]=np.loadtxt(obs_dir_base+'/observations/CDNC_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]
obs_set[386]=np.loadtxt(obs_dir_base+'/observations/CDNC_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]
obs_set[387]=np.loadtxt(obs_dir_base+'/observations/CDNC_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]
obs_set[388]=np.loadtxt(obs_dir_base+'/observations/CDNC_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]
obs_set[389]=np.loadtxt(obs_dir_base+'/observations/CDNC_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]
obs_set[390]=np.loadtxt(obs_dir_base+'/observations/CDNC_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]
obs_set[391]=np.loadtxt(obs_dir_base+'/observations/CDNC_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]
obs_set[392]=np.loadtxt(obs_dir_base+'/observations/CDNC_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0] #dec
obs_set[393]=np.loadtxt(obs_dir_base+'/observations/CDNC_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12] # ann

## SOcean Fc
obs_set[394]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[395]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]
obs_set[396]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]
obs_set[397]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]
obs_set[398]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]
obs_set[399]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]
obs_set[400]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]
obs_set[401]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]
obs_set[402]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]
obs_set[403]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]
obs_set[404]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]
obs_set[405]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]
obs_set[406]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]
obs_set[407]=np.loadtxt(obs_dir_base+'/observations/Cloud_Retrieval_Fraction_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]

## SOcean Re
obs_set[408]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Particle_Size_Liquid_N_southern_ocean_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[409]=np.loadtxt(obs_dir_base+'observations/Cloud_Particle_Size_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]/1000000 
obs_set[410]=np.loadtxt(obs_dir_base+'observations/Cloud_Particle_Size_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]/1000000 
obs_set[411]=np.loadtxt(obs_dir_base+'observations/Cloud_Particle_Size_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]/1000000 
obs_set[412]=np.loadtxt(obs_dir_base+'observations/Cloud_Particle_Size_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]/1000000 
obs_set[413]=np.loadtxt(obs_dir_base+'observations/Cloud_Particle_Size_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]/1000000 
obs_set[414]=np.loadtxt(obs_dir_base+'observations/Cloud_Particle_Size_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]/1000000 
obs_set[415]=np.loadtxt(obs_dir_base+'observations/Cloud_Particle_Size_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]/1000000 
obs_set[416]=np.loadtxt(obs_dir_base+'observations/Cloud_Particle_Size_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]/1000000 
obs_set[417]=np.loadtxt(obs_dir_base+'observations/Cloud_Particle_Size_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]/1000000 
obs_set[418]=np.loadtxt(obs_dir_base+'observations/Cloud_Particle_Size_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]/1000000 
obs_set[419]=np.loadtxt(obs_dir_base+'observations/Cloud_Particle_Size_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]/1000000 
obs_set[420]=np.loadtxt(obs_dir_base+'observations/Cloud_Particle_Size_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]/1000000 
obs_set[421]=np.loadtxt(obs_dir_base+'observations/Cloud_Particle_Size_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]/1000000 

## SOcean LWP
obs_set[422]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Water_Path_Liquid_N_southern_ocean_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[423]=np.loadtxt(obs_dir_base+'observations/Cloud_Water_Path_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask.dat')[1]/1000
obs_set[424]=np.loadtxt(obs_dir_base+'observations/Cloud_Water_Path_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask.dat')[2]/1000
obs_set[425]=np.loadtxt(obs_dir_base+'observations/Cloud_Water_Path_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask.dat')[3]/1000
obs_set[426]=np.loadtxt(obs_dir_base+'observations/Cloud_Water_Path_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask.dat')[4]/1000
obs_set[427]=np.loadtxt(obs_dir_base+'observations/Cloud_Water_Path_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask.dat')[5]/1000
obs_set[428]=np.loadtxt(obs_dir_base+'observations/Cloud_Water_Path_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask.dat')[6]/1000
obs_set[429]=np.loadtxt(obs_dir_base+'observations/Cloud_Water_Path_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask.dat')[7]/1000
obs_set[430]=np.loadtxt(obs_dir_base+'observations/Cloud_Water_Path_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask.dat')[8]/1000
obs_set[431]=np.loadtxt(obs_dir_base+'observations/Cloud_Water_Path_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask.dat')[9]/1000
obs_set[432]=np.loadtxt(obs_dir_base+'observations/Cloud_Water_Path_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask.dat')[10]/1000
obs_set[433]=np.loadtxt(obs_dir_base+'observations/Cloud_Water_Path_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask.dat')[11]/1000
obs_set[434]=np.loadtxt(obs_dir_base+'observations/Cloud_Water_Path_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask.dat')[0]/1000
obs_set[435]=np.loadtxt(obs_dir_base+'observations/Cloud_Water_Path_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask.dat')[12]/1000

## SOcean Tc
obs_set[436]=np.loadtxt(obs_dir_base+'/seasonal_amplitude/Cloud_Optical_Thickness_Liquid_N_southern_ocean_seasonal_amplitude_observed_revised_match_MODIS_CF.dat')
obs_set[437]=np.loadtxt(obs_dir_base+'observations/Cloud_Optical_Thickness_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[1]
obs_set[438]=np.loadtxt(obs_dir_base+'observations/Cloud_Optical_Thickness_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[2]
obs_set[439]=np.loadtxt(obs_dir_base+'observations/Cloud_Optical_Thickness_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[3]
obs_set[440]=np.loadtxt(obs_dir_base+'observations/Cloud_Optical_Thickness_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[4]
obs_set[441]=np.loadtxt(obs_dir_base+'observations/Cloud_Optical_Thickness_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[5]
obs_set[442]=np.loadtxt(obs_dir_base+'observations/Cloud_Optical_Thickness_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[6]
obs_set[443]=np.loadtxt(obs_dir_base+'observations/Cloud_Optical_Thickness_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[7]
obs_set[444]=np.loadtxt(obs_dir_base+'observations/Cloud_Optical_Thickness_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[8]
obs_set[445]=np.loadtxt(obs_dir_base+'observations/Cloud_Optical_Thickness_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[9]
obs_set[446]=np.loadtxt(obs_dir_base+'observations/Cloud_Optical_Thickness_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[10]
obs_set[447]=np.loadtxt(obs_dir_base+'observations/Cloud_Optical_Thickness_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[11]
obs_set[448]=np.loadtxt(obs_dir_base+'observations/Cloud_Optical_Thickness_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[0]
obs_set[449]=np.loadtxt(obs_dir_base+'observations/Cloud_Optical_Thickness_Liquid_N_southern_ocean_dec_to_ann_from_regridded_w_mask_revised_match_MODIS_CF.dat')[12]


################
## Sample data
## from emulation
################

icon=0

dir_sample= '/gws/nopw/j04/acure/lregayre/data_post_emulation/'
sample_mean=np.loadtxt(dir_sample+'emulated_mean_values_'+constraint_variable_text_list[icon]+'_'+sample_type_text)
sample_sd=np.loadtxt(dir_sample+'emulated_sd_values_'+constraint_variable_text_list[icon]+'_'+sample_type_text)


if region=='NAtlantic':
    icon_out=np.arange(102-14)+14 ## -14 to remove  hemispheric difference constraints regionally
elif region=='NPacific':
    icon_out=np.arange(90)+102 ## add on NAtlantic
elif region=='SPacific':
    icon_out=np.arange(85)+281
elif region=='SAtlantic':
    icon_out=np.arange(89)+192
elif region=='SOcean':
    icon_out=np.arange(84)+366
else:
    icon_out=np.arange(ncon_full)

ncon=len(icon_out)


#####################################
## Read in PPE data for each variable
## global and regional
#####################################

var_list_PPE = ['ERF','ACI','ARI']
var_list_titles=[DELTA+' F',DELTA+' F$_{aci}$',DELTA+' F$_{ari}$']
units_list=['W m$^{-2}$','W m$^{-2}$','W m$^{-2}$']

nvar=len(var_list_PPE)
sample_array=np.zeros((nvar,sample_size))

for ivar in np.arange(nvar):
    file_in=dir_sample+'emulated_mean_values_'+var_list_PPE[ivar]+'_'+sample_type_text
    sample_array[ivar,:]=np.loadtxt(file_in)


##############################
## Calculate NRMSE and order
##############################

## example first, which is partly used for length of sample later

icon=0
obs_data=obs_set[icon]
diff_list=obs_data-sample_mean
abs_diff_list=np.abs(diff_list)
adjusted_abs_diff_list=np.abs(abs_diff_list-sample_sd)
adjusted_abs_diff_list[np.where(abs_diff_list<sample_sd)]=0.0 ## reset to zero if actual difference is less than emulator error

max_val=np.max(adjusted_abs_diff_list)
min_val=np.min(adjusted_abs_diff_list)

norm_abs_diff_list=(1- (max_val-adjusted_abs_diff_list)/(max_val-min_val))

sample_index=np.arange(sample_mean.shape[0])

sorted_norm_abs_diff_list,sorted_sample_index= zip(*sorted(zip(norm_abs_diff_list,sample_index)))

sample_size=len(sorted_norm_abs_diff_list)


def calc_NRMSE_w_all_uncertainty_sources(icon):
    sample_mean=np.loadtxt(dir_sample+'emulated_mean_values_'+constraint_variable_text_list[icon]+'_'+sample_type_text)
    sample_sd=np.loadtxt(dir_sample+'emulated_sd_values_'+constraint_variable_text_list[icon]+'_'+sample_type_text)
    diff_list=obs_set[icon]-sample_mean
    print('     ')
    print(constraint_variable_text_list[icon])
    print('Observed value: '+np.str(obs_set[icon]))
    print('Sample mean: '+np.str(np.mean(sample_mean)))
    print('    ')
    abs_diff_list=np.abs(diff_list)
    abs_obs_uncertainty=np.abs(obs_set[icon]*total_error_percent/100)
    uncertainty_adjustment=sample_sd+abs_obs_uncertainty
    adjusted_abs_diff_list=np.abs(abs_diff_list-uncertainty_adjustment)
    max_val=np.max(adjusted_abs_diff_list)
    min_val=np.min(adjusted_abs_diff_list)
    max_val
    min_val
    norm_abs_diff_list=(1- (max_val-adjusted_abs_diff_list)/(max_val-min_val))
    norm_abs_diff_list[np.where(abs_diff_list<uncertainty_adjustment)]=0.0
    return(norm_abs_diff_list,sample_mean,sample_sd,abs_diff_list)




NRMSE_path_all='/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/individual_NRMSE_indices/original_additional_and_Hd/'


#######################################
## Write out indices for individual constraints
## using constraint_variable_text_list 
#######################################

file_out_indices=NRMSE_path_all+'/all/Indices_only_for_variables_w_decent_emulators_and_constraint_potential_all.txt'
indices_as_decimals=np.loadtxt(file_out_indices)
icon_decent_emulator_some_constraint_all_regions=[]
for icon in np.arange(len(indices_as_decimals)):
    icon_decent_emulator_some_constraint_all_regions.append(np.int(indices_as_decimals[icon]))


icon_decent_emulator_some_constraint=[]
for icon in icon_out:
    if icon in (icon_decent_emulator_some_constraint_all_regions):
        icon_decent_emulator_some_constraint.append(icon)



###############################################
## Remake NRMSE_decent_emulator_some_constraint
## - a duplicate for case where first_pass==FALSE
###############################################

file_out='/gws/nopw/j04/acure/lregayre/constraint/combined_constraint_sample/individual_NRMSE_indices/original_additional_and_Hd/all/NRMSE_including_discarded.dat'

## Read in NRMSE_including_discarded
print('About to read large file. Hold tight')
NRMSE_including_discarded=np.loadtxt(file_out)
for icon in icon_decent_emulator_some_constraint:
    print(icon)
    try: NRMSE_decent_emulator_some_constraint
    except NameError: NRMSE_decent_emulator_some_constraint=NRMSE_including_discarded[icon,:]
    else: NRMSE_decent_emulator_some_constraint=np.vstack((NRMSE_decent_emulator_some_constraint,NRMSE_including_discarded[icon,:]))




#############################################
## Quantify constraint using ALL POSSIBLE viable variables
## (NRMSE_decent_emulator_some_constraint)
#############################################


average_error_w_index=np.zeros((2,sample_size))
average_error_w_index[0,:]=np.copy(sample_index) ## UNSORTED
average_error_w_index[1,:]=NRMSE_decent_emulator_some_constraint.mean(axis=0)
(unique, counts) = np.unique(average_error_w_index[1,:], return_counts=True)
frequencies = np.asarray((unique, counts)).T
percentage_zero=frequencies[0,1]/sample_size*100
if np.logical_or((percentage_zero<=constraint_percent),(unique[0]>0.0)):
    retain_number=np.int(constraint_percent*sample_size/100)
else:
    retain_number=np.int(percentage_zero*sample_size/100)

sorted_average_error=np.zeros((2,sample_size))
sorted_average_error[1,:],sorted_average_error[0,:]= zip(*sorted(zip(average_error_w_index[1,:],average_error_w_index[0,:])))
constraint_subset_vals=sorted_average_error[0,0:np.int(retain_number)]
constraint_integers=[]
for ival in np.arange(constraint_subset_vals.shape[0]):
    constraint_integers.append(np.int(constraint_subset_vals[ival]))

sample_array_constrained=np.zeros((nvar,np.int(retain_number))) ## in sample number order
for isubset in np.arange(retain_number):
    sample_array_constrained[:,np.int(isubset)]= sample_array[:,constraint_integers[np.int(isubset)]]


constrained_box_values=np.zeros((nvar,5))
for ivar in np.arange(nvar):
    constrained_box_values[ivar,0]=np.percentile(sample_array_constrained[ivar,:],5)
    constrained_box_values[ivar,1]=np.percentile(sample_array_constrained[ivar,:],25)
    constrained_box_values[ivar,2]=np.percentile(sample_array_constrained[ivar,:],50)
    constrained_box_values[ivar,3]=np.percentile(sample_array_constrained[ivar,:],75)
    constrained_box_values[ivar,4]=np.percentile(sample_array_constrained[ivar,:],95)

constrained_90pc=constrained_box_values[1,4]-constrained_box_values[1,0]
constrained_66pc=np.percentile(sample_array[1,:],66)-np.percentile(sample_array[1,:],33)



#####################################################
## Now apply only constraints consistent w/ Nd in target region
#####################################################


icon_consistent=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,102,103,104,106,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,192,197,198,199,201,202,203,204,205,206,207,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,229,230,231,232,233,234,235,237,238,281,282,283,284,285,286,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,395,396,397,398,399,400,401,402,403,404,405,406,407] ## includes hemispheric differences


#######################################
## Keep icon_out if
## emulates well and has constraint potential
## ONLY USE icon_out_subset 
#######################################

icon_out_subset=[]
for icon in icon_consistent:
    if icon in icon_decent_emulator_some_constraint:
        icon_out_subset.append(icon)


## remake NRMSE array for Nd pairwise consistent variables

NRMSE_consistent=NRMSE_including_discarded[icon_out_subset[0],:]
for icon in icon_out_subset[1:]:
    print(icon)
    print(np.mean(NRMSE_including_discarded[icon,:]))
    NRMSE_consistent=np.vstack((NRMSE_consistent,NRMSE_including_discarded[icon,:]))


average_error_w_index=np.zeros((2,sample_size))
average_error_w_index[0,:]=np.copy(sample_index) 
average_error_w_index[1,:]=NRMSE_consistent.mean(axis=0)
(unique, counts) = np.unique(average_error_w_index[1,:], return_counts=True)
frequencies = np.asarray((unique, counts)).T
percentage_zero=frequencies[0,1]/sample_size*100
if np.logical_or((percentage_zero<=constraint_percent),(unique[0]>0.0)):
    retain_number=np.int(constraint_percent*sample_size/100)
else:
    retain_number=np.int(percentage_zero*sample_size/100)


sorted_average_error=np.zeros((2,sample_size))
sorted_average_error[1,:],sorted_average_error[0,:]= zip(*sorted(zip(average_error_w_index[1,:],average_error_w_index[0,:])))
constraint_subset_vals=sorted_average_error[0,0:np.int(retain_number)]
constraint_integers=[]
for ival in np.arange(constraint_subset_vals.shape[0]):
    constraint_integers.append(np.int(constraint_subset_vals[ival]))

sample_array_constrained=np.zeros((nvar,np.int(retain_number))) ## in sample number order
for isubset in np.arange(retain_number):
    sample_array_constrained[:,np.int(isubset)]= sample_array[:,constraint_integers[np.int(isubset)]]


constrained_box_values=np.zeros((nvar,5))
for ivar in np.arange(nvar):
    constrained_box_values[ivar,0]=np.percentile(sample_array_constrained[ivar,:],5)
    constrained_box_values[ivar,1]=np.percentile(sample_array_constrained[ivar,:],25)
    constrained_box_values[ivar,2]=np.percentile(sample_array_constrained[ivar,:],50)
    constrained_box_values[ivar,3]=np.percentile(sample_array_constrained[ivar,:],75)
    constrained_box_values[ivar,4]=np.percentile(sample_array_constrained[ivar,:],95)

constrained_90pc=constrained_box_values[1,4]-constrained_box_values[1,0]
constrained_66pc=np.percentile(sample_array[1,:],66)-np.percentile(sample_array[1,:],33)


####################################
## 90% credible interval box
## to match CMIP, Bellouin et al., 2020, etc.
####################################

unconstrained_box_values=np.zeros((nvar,5))
for ivar in np.arange(nvar):
    unconstrained_box_values[ivar,0]=np.percentile(sample_array[ivar,:],5)
    unconstrained_box_values[ivar,1]=np.percentile(sample_array[ivar,:],25)
    unconstrained_box_values[ivar,2]=np.percentile(sample_array[ivar,:],50)
    unconstrained_box_values[ivar,3]=np.percentile(sample_array[ivar,:],75)
    unconstrained_box_values[ivar,4]=np.percentile(sample_array[ivar,:],95)

unconstrained_90pc=unconstrained_box_values[1,4]-unconstrained_box_values[1,0]
unconstrained_66pc=np.percentile(sample_array[1,:],66)-np.percentile(sample_array[1,:],33)


##########################################
## Method searches ALL possible constraints to 
## remove effect of order constraints applied
## 1/ search through ALL possible constraints
## to find strongest 
## 2/ Repeat until no viable constraints
## 3/ Test removing and replacing
##########################################


def add_strongest_constraint(index_list_used,index_list_additional,NRMSE_array,constrained_target,pc_diff):
    icon_logic=False
    icon_keep=[]
    icon_reject=[]
    pc_reduction_tracker_90=[]
    pc_reduction_tracker_66=[]
    constrained_so_far_90pc_reduction=constrained_target
    strongest_so_far=-100.0 ## percentage constraint increase over previous constraints - set to negative value here b/c always want to add at least 1
    print('Initial constraint percentage: '+np.str(constrained_so_far_90pc_reduction))
    ncon_so_far=len(index_list_used)
    print('ncon_so_far='+np.str(ncon_so_far))
    NRMSE_keep_so_far=np.array([])
    if ncon_so_far>0:
        for icon in index_list_used:
            NRMSE_keep_so_far= np.vstack((NRMSE_keep_so_far, NRMSE_array[icon,:])) if NRMSE_keep_so_far.size else NRMSE_array[icon,:]
            icon_keep.append(icon)
            print(np.str(icon_keep))
            print(np.str(NRMSE_keep_so_far.shape))
    for icon in index_list_additional:
        average_error_w_index=np.zeros((2,sample_size)) ## keeps index imember from 1st icon, 3rd column
        average_error_w_index[0,:]=np.copy(sample_index) ## UNSORTED
        average_error_w_index[1,:]=np.vstack((NRMSE_keep_so_far,NRMSE_array[icon,:])).mean(axis=0) if NRMSE_keep_so_far.size else NRMSE_array[icon,:]
        (unique, counts) = np.unique(average_error_w_index[1,:], return_counts=True)
        frequencies = np.asarray((unique, counts)).T
        percentage_zero=frequencies[0,1]/sample_size*100
        if np.logical_or((percentage_zero<=constraint_percent),(unique[0]>0.0)):
            retain_number=np.int(constraint_percent*sample_size/100)
        else:
            retain_number=np.int(percentage_zero*sample_size/100)
        sorted_average_error=np.zeros((2,sample_size))
        sorted_average_error[1,:],sorted_average_error[0,:]= zip(*sorted(zip(average_error_w_index[1,:],average_error_w_index[0,:])))
        constraint_subset_vals=sorted_average_error[0,0:np.int(retain_number)]
        constraint_integers=[]
        for ival in np.arange(constraint_subset_vals.shape[0]):
            constraint_integers.append(np.int(constraint_subset_vals[ival]))
        sample_array_constrained=np.zeros((nvar,np.int(retain_number))) ## in sample number order
        for isubset in np.arange(retain_number):
            sample_array_constrained[:,np.int(isubset)]= sample_array[:,constraint_integers[np.int(isubset)]]
        constrained_box_values=np.zeros((nvar,5))
        for ivar in np.arange(nvar):
            constrained_box_values[ivar,0]=np.percentile(sample_array_constrained[ivar,:],5)
            constrained_box_values[ivar,1]=np.percentile(sample_array_constrained[ivar,:],25)
            constrained_box_values[ivar,2]=np.percentile(sample_array_constrained[ivar,:],50)
            constrained_box_values[ivar,3]=np.percentile(sample_array_constrained[ivar,:],75)
            constrained_box_values[ivar,4]=np.percentile(sample_array_constrained[ivar,:],95)
            constrained_1sigma_neg=np.percentile(sample_array_constrained[1,:],33)
            constrained_1sigma_plus=np.percentile(sample_array_constrained[1,:],66)
        print(unconstrained_box_values)
        print(constrained_box_values)
        unconstrained_90pc=unconstrained_box_values[1,4]-unconstrained_box_values[1,0]
        constrained_90pc=constrained_box_values[1,4]-constrained_box_values[1,0]
        constrained_90pc_reduction=(unconstrained_90pc-constrained_90pc)/unconstrained_90pc*100
        constrained_66pc=constrained_1sigma_plus-constrained_1sigma_neg
        constrained_66pc_reduction=(unconstrained_66pc-constrained_66pc)/unconstrained_66pc*100
        print('90% CI reduced by: ')
        print(np.str(constrained_90pc_reduction))
        print('1sigma CI reduced by: ')
        print(np.str(constrained_66pc_reduction))
        if ((constrained_90pc_reduction - constrained_so_far_90pc_reduction)>strongest_so_far):
            icon_logic=True
            print('STRONGEST SO FAR CHANGED TO icon='+np.str(icon)+' '+constraint_variable_text_list[icon]+' WITH DIFFERENCE IN CONSTRAINT OF :'+np.str(constrained_90pc_reduction - constrained_so_far_90pc_reduction))
            strongest_so_far=(constrained_90pc_reduction - constrained_so_far_90pc_reduction)
            if 'icon_strongest' in locals(): ## exists after 1st rejected var
                icon_reject.append(icon_strongest) ## previous strongest constraint
            icon_strongest=icon
            pc_reduction_tracker_90=constrained_90pc_reduction
            pc_reduction_tracker_66=constrained_66pc_reduction
        else:
            icon_reject.append(icon)
            print('REJECTING icon='+np.str(icon)+' '+constraint_variable_text_list[icon]+' FOR NOW AS NOT STRONGEST CONSTRAINT')
            print('Current constrained 90% CI difference='+np.str(strongest_so_far)+'%')
        if 'icon_strongest' in locals():
            print('STRONGEST SO FAR IS icon='+np.str(icon_strongest)+' '+constraint_variable_text_list[icon_strongest]+' WITH CONSTRAINT OF :'+np.str(pc_reduction_tracker_90))
        else:
            print('No additional constraints detected yet')
        print('                         ')
    ## identify strongest
    if 'icon_strongest' in locals():
        icon_keep.append(icon_strongest)
        NRMSE_keep_so_far=np.vstack((NRMSE_keep_so_far,NRMSE_array[icon_strongest,:])) if NRMSE_keep_so_far.size else NRMSE_array[icon_strongest,:]
    pc_reduction_tracker_90
    pc_reduction_tracker_66
    return(icon_keep,icon_reject,pc_reduction_tracker_90,NRMSE_keep_so_far,icon_logic)



#################################
## First constraint step, or restart skip
#################################


if restart_partial:
    print('skipping top down initialisation step of: remove weakest')
else:
    icon_keep_updated,icon_reject_updated,pc_reduction_90_updated,NRMSE_keep_updated,constraint_TF=add_strongest_constraint([],icon_consistent,NRMSE_including_discarded,0.0,90.0)

    
##############################
## First write to file, or restart skip
##############################

#################################
## Calculate constraint so far
## Write indices to file
## Write constraint_variable_text to file
## Put 90pc reduction in title
## Also write out icon_reject_updated 
## so process can be restarted
#################################


def write_indices_to_file(NRMSE_array,icon_list_reject,icon_list_retain,fname_start):
    average_error_w_index=np.zeros((2,sample_size)) ## keeps index imember from 1st icon, or 3rd column
    average_error_w_index[0,:]=np.copy(sample_index) ## UNSORTED
    if (NRMSE_array.shape[0]==sample_size): ## only one constraint so far
        average_error_w_index[1,:]=NRMSE_array
    else:
        average_error_w_index[1,:]=NRMSE_array.mean(axis=0)
    (unique, counts) = np.unique(average_error_w_index[1,:], return_counts=True)
    frequencies = np.asarray((unique, counts)).T
    percentage_zero=frequencies[0,1]/sample_size*100
    if np.logical_or((percentage_zero<=constraint_percent),(unique[0]>0.0)):
        retain_number=np.int(constraint_percent*sample_size/100)
    else:
        retain_number=np.int(percentage_zero*sample_size/100)
    sorted_average_error=np.zeros((2,sample_size))
    sorted_average_error[1,:],sorted_average_error[0,:]= zip(*sorted(zip(average_error_w_index[1,:],average_error_w_index[0,:])))
    constraint_subset_vals=sorted_average_error[0,0:np.int(retain_number)]
    constraint_integers=[]
    for ival in np.arange(constraint_subset_vals.shape[0]):
        constraint_integers.append(np.int(constraint_subset_vals[ival]))
    sample_array_constrained=np.zeros((nvar,np.int(retain_number))) ## in sample number order
    for isubset in np.arange(retain_number):
        sample_array_constrained[:,np.int(isubset)]= sample_array[:,constraint_integers[np.int(isubset)]]
    for ivar in np.arange(nvar):
        constrained_box_values[ivar,0]=np.percentile(sample_array_constrained[ivar,:],5)
        constrained_box_values[ivar,1]=np.percentile(sample_array_constrained[ivar,:],25)
        constrained_box_values[ivar,2]=np.percentile(sample_array_constrained[ivar,:],50)
        constrained_box_values[ivar,3]=np.percentile(sample_array_constrained[ivar,:],75)
        constrained_box_values[ivar,4]=np.percentile(sample_array_constrained[ivar,:],95)
        constrained_1sigma_neg=np.percentile(sample_array_constrained[1,:],33)
        constrained_1sigma_plus=np.percentile(sample_array_constrained[1,:],66)
    print(unconstrained_box_values)
    print(constrained_box_values)
    unconstrained_90pc=unconstrained_box_values[1,4]-unconstrained_box_values[1,0]
    constrained_90pc=constrained_box_values[1,4]-constrained_box_values[1,0]
    constrained_90pc_reduction=(unconstrained_90pc-constrained_90pc)/unconstrained_90pc*100
    if NRMSE_array.shape[0]==1000000: ## first constraint
        ncon=1
    else:
        ncon=NRMSE_array.shape[0]
    print('Constrained percentage: '+np.str(constrained_90pc_reduction))
    file_out_labels=optimal_path_out+'/'+fname_start+'_'+np.str(total_error_percent)+'pc_error_constrained_ERF_'+np.str(np.round(constrained_box_values[0,0],2))+'_to_'+np.str(np.round(constrained_box_values[0,4],2))+'_with_'+np.str(np.round(constrained_90pc_reduction,2))+'pc_reduction_from_'+np.str(ncon)+'_constraints_CONSTRAINT_NAMES_OF_REMOVED_w_o_Hd_regional.txt'
    file_out_reject_indices=optimal_path_out+'/'+fname_start+'_'+np.str(total_error_percent)+'pc_error_constrained_ERF_'+np.str(np.round(constrained_box_values[0,0],2))+'_to_'+np.str(np.round(constrained_box_values[0,4],2))+'_with_'+np.str(np.round(constrained_90pc_reduction,2))+'pc_reduction_from_'+np.str(ncon)+'_constraints_INDICES_OF_REJECTED_w_o_Hd_regional.txt'
    file_out_retain_indices=optimal_path_out+'/'+fname_start+'_'+np.str(total_error_percent)+'pc_error_constrained_ERF_'+np.str(np.round(constrained_box_values[0,0],2))+'_to_'+np.str(np.round(constrained_box_values[0,4],2))+'_with_'+np.str(np.round(constrained_90pc_reduction,2))+'pc_reduction_from_'+np.str(ncon)+'_constraints_INDICES_OF_RETAINED_w_o_Hd_regional.txt'
    if os.path.exists(file_out_labels):
        os.remove(file_out_labels)
    for icon in icon_list_reject:
        with open(file_out_labels, 'a', encoding='utf-8') as file:
            file.write(constraint_variable_text_list[icon]+'\n')
    for icon in np.arange(len(icon_list_reject)):
         with open(file_out_reject_indices, 'a', encoding='utf-8') as file:
            file.write(np.str(icon_list_reject[icon])+'\n')
    for icon in np.arange(len(icon_list_retain)):
         with open(file_out_retain_indices, 'a', encoding='utf-8') as file:
            file.write(np.str(icon_list_retain[icon])+'\n')


if restart_partial:
    print('Skipping initialisation write to file')
else:
    write_indices_to_file(NRMSE_keep_updated,icon_reject_updated,icon_keep_updated,region+'_sample_indices_py_index_FINESSE_w_uncertainty')


## By now, each of icon_keep_updated,icon_reject_updatedmNRMSE_including_discarded,pc_reduction_90_updated
## have been calculated, or read in for a restart


###################################
## Loop to progressively add constraints
###################################


while (len(icon_reject_updated)>1):
    icon_keep_updated,icon_reject_updated,pc_reduction_90_updated,NRMSE_keep_updated,constraint_TF=add_strongest_constraint(icon_keep_updated,icon_reject_updated,NRMSE_including_discarded,pc_reduction_90_updated,90.0)
    write_indices_to_file(NRMSE_keep_updated,icon_reject_updated,icon_keep_updated,region+'_sample_indices_py_index_FINESSE_w_uncertainty')


