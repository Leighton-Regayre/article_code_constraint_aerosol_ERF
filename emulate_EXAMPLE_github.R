## Leighton Regayre. Cope adapted from original shared by Lindsay Lee
## Create emulators of monthly mean and/or transect PPE output
## 220 or 221 members depending on month (as stored in filename)
##
## Predict from emulator, using large sample for later constraint

# Load the packages for the emulation and sensitivity analysis. ## may need to 'install.packages("DiceKriging")' first - note pop-up window for CRAN mirror
library(DiceKriging)
library(sensitivity)
library(truncnorm)
library(lhs)
library(trapezoid)
library(readr) ## for writing large data files

iemulate=0 ## this index ranges from 0-453 spanning all constraint variables, including those that do not emulate well and were not included in the 450 described in the main article
reduce_1M=1 ## binary operator to temporarily reduce the retained set of variants to the first million members
reduce_size=1000000
sample_size=2 ## 1 or 2 ## index for file type

## Set working directory

setwd("/gws/nopw/j04/acure/lregayre/emulation/")


## Constraint features

constraint_itype_list=c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1)+1 ## +1 for R index starting at 1, not 0 as python, but retained structure to match python code

## 0+1=monthly mean; 2+1= transect; 1+1=hemispheric difference

ivar_list=c(2,3,6,1,2,3,4,5,6,7,8,0,2,4,5,0,6,7,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,5,5,5,5,5,5,5,5,5,5,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,5,5,5,5,5,5,5,5,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,5,5,5,5,5,5,5,5,5,5,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,2,4,2,3,4,2,4,1,0,1,4,5,1,2,2,2,2,2,2,2,2,0,0,0,4,4,4,4,4,4,4,4,4,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)+1
##var_list=c('CDNC','CF_Liquid_MODIS','LWP_MODIS','Reff_Liquid_MODIS','Tau_Liquid_MODIS','SW_TOA')


imonth_list=c(0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,3,4,5,6,7,8,9,10,11,13,0,1,2,3,4,6,7,8,9,10,11,12,13,0,1,2,7,8,0,1,2,5,6,7,8,9,10,11,13,0,1,2,7,8,0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,1,2,3,4,5,6,7,8,9,12,13,0,2,3,5,6,7,8,9,10,11,12,13,0,1,2,4,12,0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,1,2,3,4,0,1,2,3,4,5,6,7,8,9,10,11,12,13,1,2,3,4,5,6,7,8,9,10,11,12,0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,5,6,7,8,9,12,0,2,3,0,6,7,8,0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,0,0,0,0,0,0,5,7,1,1,1,4,4,5,6,9,10,11,12,13,1,2,12,3,4,5,6,9,10,11,12,13,3,4,12,1,2,3,4,5,6,7,8,9,10,11,12,13,1,2,3,4,5,6,7,8,9,10,11,12,13,1,2,3,4,5,6,7,8,9,10,11,12,13,1,2,3,4,10,11,13,1,2,3,4,5,9,10,11,12,13,1,4,5,6,7,8,9,10,11,12,13,1,2,3,4,5,6,7,8,9,10,11,12,13,10,11,1,2,3,4,5,6,7,8,9,10,11,12,13,3,4,5,6,7,8,9,10,13,1,2,3,4,5,6,7,8,9,10,11,12,13,1,2,3,4,5,6,7,8,9,10,11,12,13,5,6,7,8,9,10,11,12,13,1,0,1,2,3,4,5,6,7,8,9,10,11,12,13)+1

## 99 for annual mean, ignored; 0=amp; 1=dec; 2=jan; 3=feb; 4=mar; 5=apr; 6=may; 7=jun; 8=jul; 9=aug; 10=sep; 11=oct; 12=nov; 13=ann


iregion_list=c(0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,1,1,2,2,2,3,3,0,3,3,3,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,1,5,5,5,5,5,5,5,5,5,5,5,5,5,5)+1

## region_list=c('NN_atlantic','SE_atlantic','NE_pacific','SE_pacific','N_southern_ocean') added hemispheric difference

## order of constraints was determined organically over multiple research stages

constraint_type_text_1=c('jul_transect_LWP_NAtlantic','jul_transect_Re_NAtlantic','jul_transect_Cf_Nd_NAtlantic','jul_transect_Cf_NPacific','jul_transect_LWP_NPacific','jul_transect_Re_NPacific','jul_transect_AI_NPacific','jul_transect_Nd_AI_NPacific','jul_transect_Cf_Nd_NPacific','jul_transect_LWP_Nd_NPacific','jul_transect_Re_Nd_NPacific','nov_transect_Nd_SPacific','nov_transect_LWP_SPacific','nov_transect_AI_SPacific','nov_transect_Nd_AI_SPacific','nov_transect_Nd_SAtlantic','nov_transect_Cf_Nd_SAtlantic','nov_transect_LWP_Nd_SAtlantic',"seasonal_amplitude_Nd_NAtlantic","feb_Nd_NAtlantic","mar_Nd_NAtlantic","apr_Nd_NAtlantic","may_Nd_NAtlantic","jun_Nd_NAtlantic","jul_Nd_NAtlantic","aug_Nd_NAtlantic","sep_Nd_NAtlantic","oct_Nd_NAtlantic","ann_Nd_NAtlantic","seasonal_amplitude_Fc_NAtlantic","dec_Fc_NAtlantic","jan_Fc_NAtlantic","feb_Fc_NAtlantic","mar_Fc_NAtlantic","may_Fc_NAtlantic","jun_Fc_NAtlantic","jul_Fc_NAtlantic","aug_Fc_NAtlantic","sep_Fc_NAtlantic","oct_Fc_NAtlantic","nov_Fc_NAtlantic","ann_Fc_NAtlantic","seasonal_amplitude_LWP_NAtlantic","dec_LWP_NAtlantic","jan_LWP_NAtlantic","jun_LWP_NAtlantic","jul_LWP_NAtlantic","seasonal_amplitude_Re_NAtlantic","dec_Re_NAtlantic","jan_Re_NAtlantic","apr_Re_NAtlantic","may_Re_NAtlantic","jun_Re_NAtlantic","jul_Re_NAtlantic","aug_Re_NAtlantic","sep_Re_NAtlantic","oct_Re_NAtlantic","ann_Re_NAtlantic","seasonal_amplitude_Tau_NAtlantic","dec_Tau_NAtlantic","jan_Tau_NAtlantic","jun_Tau_NAtlantic","jul_Tau_NAtlantic","seasonal_amplitude_Fsw_NAtlantic","dec_Fsw_NAtlantic","jan_Fsw_NAtlantic","feb_Fsw_NAtlantic","mar_Fsw_NAtlantic","apr_Fsw_NAtlantic","may_Fsw_NAtlantic","jun_Fsw_NAtlantic","jul_Fsw_NAtlantic","aug_Fsw_NAtlantic","sep_Fsw_NAtlantic","oct_Fsw_NAtlantic","nov_Fsw_NAtlantic","ann_Fsw_NAtlantic","seasonal_amplitude_Nd_SAtlantic","dec_Nd_SAtlantic","jan_Nd_SAtlantic","feb_Nd_SAtlantic","mar_Nd_SAtlantic","apr_Nd_SAtlantic","may_Nd_SAtlantic","jun_Nd_SAtlantic","jul_Nd_SAtlantic","aug_Nd_SAtlantic","nov_Nd_SAtlantic","ann_Nd_SAtlantic","seasonal_amplitude_Fc_SAtlantic","jan_Fc_SAtlantic","feb_Fc_SAtlantic","apr_Fc_SAtlantic","may_Fc_SAtlantic","jun_Fc_SAtlantic","jul_Fc_SAtlantic","aug_Fc_SAtlantic","sep_Fc_SAtlantic","oct_Fc_SAtlantic","nov_Fc_SAtlantic","ann_Fc_SAtlantic","seasonal_amplitude_Re_SAtlantic","dec_Re_SAtlantic","jan_Re_SAtlantic","mar_Re_SAtlantic","nov_Re_SAtlantic","seasonal_amplitude_Fsw_SAtlantic","dec_Fsw_SAtlantic","jan_Fsw_SAtlantic","feb_Fsw_SAtlantic","mar_Fsw_SAtlantic","apr_Fsw_SAtlantic","may_Fsw_SAtlantic","jun_Fsw_SAtlantic","jul_Fsw_SAtlantic","aug_Fsw_SAtlantic","sep_Fsw_SAtlantic","oct_Fsw_SAtlantic","nov_Fsw_SAtlantic","ann_Fsw_SAtlantic")

constraint_type_text_2=c("seasonal_amplitude_Nd_NPacific","dec_Nd_NPacific","jan_Nd_NPacific","feb_Nd_NPacific","mar_Nd_NPacific","apr_Nd_NPacific","may_Nd_NPacific","jun_Nd_NPacific","jul_Nd_NPacific","aug_Nd_NPacific","sep_Nd_NPacific","oct_Nd_NPacific","nov_Nd_NPacific","ann_Nd_NPacific","seasonal_amplitude_Fc_NPacific","dec_Fc_NPacific","jan_Fc_NPacific","feb_Fc_NPacific","mar_Fc_NPacific","apr_Fc_NPacific","may_Fc_NPacific","jun_Fc_NPacific","jul_Fc_NPacific","aug_Fc_NPacific","sep_Fc_NPacific","oct_Fc_NPacific","nov_Fc_NPacific","ann_Fc_NPacific","seasonal_amplitude_Fsw_NPacific","dec_Fsw_NPacific","jan_Fsw_NPacific","feb_Fsw_NPacific","mar_Fsw_NPacific","apr_Fsw_NPacific","may_Fsw_NPacific","jun_Fsw_NPacific","jul_Fsw_NPacific","aug_Fsw_NPacific","sep_Fsw_NPacific","oct_Fsw_NPacific","nov_Fsw_NPacific","ann_Fsw_NPacific","seasonal_amplitude_Nd_SPacific","dec_Nd_SPacific","jan_Nd_SPacific","feb_Nd_SPacific","mar_Nd_SPacific","apr_Nd_SPacific","may_Nd_SPacific","jun_Nd_SPacific","jul_Nd_SPacific","aug_Nd_SPacific","sep_Nd_SPacific","oct_Nd_SPacific","nov_Nd_SPacific","ann_Nd_SPacific","seasonal_amplitude_Fc_SPacific","dec_Fc_SPacific","jan_Fc_SPacific","feb_Fc_SPacific","mar_Fc_SPacific","apr_Fc_SPacific","may_Fc_SPacific","jun_Fc_SPacific","jul_Fc_SPacific","aug_Fc_SPacific","sep_Fc_SPacific","oct_Fc_SPacific","nov_Fc_SPacific","ann_Fc_SPacific","seasonal_amplitude_Re_SPacific","dec_Re_SPacific","jan_Re_SPacific","feb_Re_SPacific","mar_Re_SPacific","seasonal_amplitude_Fsw_SPacific","dec_Fsw_SPacific","jan_Fsw_SPacific","feb_Fsw_SPacific","mar_Fsw_SPacific","apr_Fsw_SPacific","may_Fsw_SPacific","jun_Fsw_SPacific","jul_Fsw_SPacific","aug_Fsw_SPacific","sep_Fsw_SPacific","oct_Fsw_SPacific","nov_Fsw_SPacific","ann_Fsw_SPacific","dec_Nd_SOcean","jan_Nd_SOcean","feb_Nd_SOcean","mar_Nd_SOcean","apr_Nd_SOcean","may_Nd_SOcean","jun_Nd_SOcean","jul_Nd_SOcean","aug_Nd_SOcean","sep_Nd_SOcean","oct_Nd_SOcean","nov_Nd_SOcean","seasonal_amplitude_Fc_SOcean","dec_Fc_SOcean","jan_Fc_SOcean","feb_Fc_SOcean","mar_Fc_SOcean","apr_Fc_SOcean","may_Fc_SOcean","jun_Fc_SOcean","jul_Fc_SOcean","aug_Fc_SOcean","sep_Fc_SOcean","oct_Fc_SOcean","nov_Fc_SOcean","ann_Fc_SOcean","seasonal_amplitude_LWP_SOcean","apr_LWP_SOcean","may_LWP_SOcean","jun_LWP_SOcean","jul_LWP_SOcean","aug_LWP_SOcean","nov_LWP_SOcean","seasonal_amplitude_Re_SOcean","jan_Re_SOcean","feb_Re_SOcean","seasonal_amplitude_Tau_SOcean","may_Tau_SOcean","jun_Tau_SOcean","jul_Tau_SOcean","seasonal_amplitude_Fsw_SOcean","dec_Fsw_SOcean","jan_Fsw_SOcean","feb_Fsw_SOcean","mar_Fsw_SOcean","apr_Fsw_SOcean","may_Fsw_SOcean","jun_Fsw_SOcean","jul_Fsw_SOcean","aug_Fsw_SOcean","sep_Fsw_SOcean","oct_Fsw_SOcean","nov_Fsw_SOcean","ann_Fsw_SOcean","seasonal_amplitude_LWP_SAtlantic","seasonal_amplitude_Tau_SAtlantic","seasonal_amplitude_LWP_NPacific","seasonal_amplitude_Re_NPacific","seasonal_amplitude_Tau_NPacific","seasonal_amplitude_LWP_SPacific","seasonal_amplitude_Tau_SPacific")

constraint_type_text_3=c("apr_Fc_NAtlantic","jun_Nd_SPacific","nov_transect_Fc_SAtlantic","nov_transect_AI_SAtlantic","nov_transect_Nd_AI_SAtlantic","mar_Fc_SAtlantic")

constraint_type_text_4=c("mar_LWP_NAtlantic","apr_LWP_NAtlantic","may_LWP_NAtlantic","aug_LWP_NAtlantic","sep_LWP_NAtlantic","oct_LWP_NAtlantic","nov_LWP_NAtlantic","ann_LWP_NAtlantic","dec_Nd_NAtlantic","jan_Nd_NAtlantic","nov_Nd_NAtlantic","feb_Tau_NAtlantic","mar_Tau_NAtlantic","apr_Tau_NAtlantic","may_Tau_NAtlantic","aug_Tau_NAtlantic","sep_Tau_NAtlantic","oct_Tau_NAtlantic","nov_Tau_NAtlantic","ann_Tau_NAtlantic","feb_Re_NAtlantic","mar_Re_NAtlantic","nov_Re_NAtlantic","dec_LWP_NPacific","jan_LWP_NPacific","feb_LWP_NPacific","mar_LWP_NPacific","apr_LWP_NPacific","may_LWP_NPacific","jun_LWP_NPacific","jul_LWP_NPacific","aug_LWP_NPacific","sep_LWP_NPacific","oct_LWP_NPacific","nov_LWP_NPacific","ann_LWP_NPacific","dec_Tau_NPacific","jan_Tau_NPacific","feb_Tau_NPacific","mar_Tau_NPacific","apr_Tau_NPacific","may_Tau_NPacific","jun_Tau_NPacific","jul_Tau_NPacific","aug_Tau_NPacific","sep_Tau_NPacific","oct_Tau_NPacific","nov_Tau_NPacific","ann_Tau_NPacific","dec_Re_NPacific","jan_Re_NPacific","feb_Re_NPacific","mar_Re_NPacific","apr_Re_NPacific","may_Re_NPacific","jun_Re_NPacific","jul_Re_NPacific","aug_Re_NPacific","sep_Re_NPacific","oct_Re_NPacific","nov_Re_NPacific","ann_Re_NPacific","dec_LWP_SOcean","jan_LWP_SOcean","feb_LWP_SOcean","mar_LWP_SOcean","sep_LWP_SOcean","oct_LWP_SOcean","ann_LWP_SOcean","dec_Tau_SOcean","jan_Tau_SOcean","feb_Tau_SOcean","mar_Tau_SOcean","apr_Tau_SOcean","aug_Tau_SOcean","sep_Tau_SOcean","oct_Tau_SOcean","nov_Tau_SOcean","ann_Tau_SOcean","dec_Re_SOcean","mar_Re_SOcean","apr_Re_SOcean","may_Re_SOcean","jun_Re_SOcean","jul_Re_SOcean","aug_Re_SOcean","sep_Re_SOcean","oct_Re_SOcean","nov_Re_SOcean","ann_Re_SOcean","dec_LWP_SAtlantic","jan_LWP_SAtlantic","feb_LWP_SAtlantic","mar_LWP_SAtlantic","apr_LWP_SAtlantic","may_LWP_SAtlantic","jun_LWP_SAtlantic","jul_LWP_SAtlantic","aug_LWP_SAtlantic","sep_LWP_SAtlantic","oct_LWP_SAtlantic","nov_LWP_SAtlantic","ann_LWP_SAtlantic","sep_Nd_SAtlantic","oct_Nd_SAtlantic","dec_Tau_SAtlantic","jan_Tau_SAtlantic","feb_Tau_SAtlantic","mar_Tau_SAtlantic","apr_Tau_SAtlantic","may_Tau_SAtlantic","jun_Tau_SAtlantic","jul_Tau_SAtlantic","aug_Tau_SAtlantic","sep_Tau_SAtlantic","oct_Tau_SAtlantic","nov_Tau_SAtlantic","ann_Tau_SAtlantic","feb_Re_SAtlantic","apr_Re_SAtlantic","may_Re_SAtlantic","jun_Re_SAtlantic","jul_Re_SAtlantic","aug_Re_SAtlantic","sep_Re_SAtlantic","oct_Re_SAtlantic","ann_Re_SAtlantic")

constraint_type_text_5=c("dec_LWP_SPacific","jan_LWP_SPacific","feb_LWP_SPacific","mar_LWP_SPacific","apr_LWP_SPacific","may_LWP_SPacific","jun_LWP_SPacific","jul_LWP_SPacific","aug_LWP_SPacific","sep_LWP_SPacific","oct_LWP_SPacific","nov_LWP_SPacific","ann_LWP_SPacific","dec_Tau_SPacific","jan_Tau_SPacific","feb_Tau_SPacific","mar_Tau_SPacific","apr_Tau_SPacific","may_Tau_SPacific","jun_Tau_SPacific","jul_Tau_SPacific","aug_Tau_SPacific","sep_Tau_SPacific","oct_Tau_SPacific","nov_Tau_SPacific","ann_Tau_SPacific","apr_Re_SPacific","may_Re_SPacific","jun_Re_SPacific","jul_Re_SPacific","aug_Re_SPacific","sep_Re_SPacific","oct_Re_SPacific","nov_Re_SPacific","ann_Re_SPacific","dec_Fc_SAtlantic","seasonal_amplitude_Hd","dec_Hd","jan_Hd","feb_Hd","mar_Hd","apr_Hd","may_Hd","jun_Hd","jul_Hd","aug_Hd","sep_Hd","oct_Hd","nov_Hd","ann_Hd")

constraint_type_text=c(constraint_type_text_1,constraint_type_text_2,constraint_type_text_3,constraint_type_text_4,constraint_type_text_5)


constraint_type_text[iemulate]

################################
## Automated process below here
################################

path_base='/gws/nopw/j04/acure/lregayre/data_for_emulation/'
month_list=c('seasonal_amplitude','dec','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','ann')
region_list=c('NN_atlantic','SE_atlantic','NE_pacific','SE_pacific','N_southern_ocean','hemispheric_difference')
var_list=c('CDNC','CF_Liquid_MODIS','LWP_MODIS','Reff_Liquid_MODIS','Tau_Liquid_MODIS','SW_TOA')
transect_var_list=c('CDNC','LCF','LWP','Reff','AI','CDNC_AI','LCF_CDNC','LWP_CDNC','Reff_CDNC')
transect_month_list=c('jul','nov')
transect_region_list=c('European','NorthAmerican','SouthAmerican','Namibian')
global_mean_var_list=c('ERF','ACI','ARI')

##########################################
## Read in training data needed to create emulator
##########################################


if (constraint_itype_list[iemulate]==4) { ## ERF and components
    datafile= paste(path_base,month_list[imonth_list[iemulate]],"/",global_mean_var_list[ivar_list[iemulate]],"_PPE_global_ann.dat",sep='')
} else if (constraint_itype_list[iemulate]==1) { ## monthly mean
    if (month_list[imonth_list[iemulate]]=='seasonal_amplitude') { ## seasonal amplitude
        datafile= paste(path_base,month_list[imonth_list[iemulate]],"/",var_list[ivar_list[iemulate]],"_",region_list[iregion_list[iemulate]],"_",month_list[imonth_list[iemulate]],"_PPE_revised_match_MODIS_CF.dat",sep='')
    } else {
        datafile= paste(path_base,month_list[imonth_list[iemulate]],"/",var_list[ivar_list[iemulate]],"_PPE_",region_list[iregion_list[iemulate]],"_",month_list[imonth_list[iemulate]],"_revised_match_MODIS_CF.dat",sep='')
    }
} else if (constraint_itype_list[iemulate]==2) { ## Hemispheric difference
    if (month_list[imonth_list[iemulate]]=='seasonal_amplitude') { ## seasonal amplitude
        datafile= paste(path_base,month_list[imonth_list[iemulate]],"/",var_list[ivar_list[iemulate]],"_NH_SH_difference_",month_list[imonth_list[iemulate]],"_PPE.dat",sep='')
    } else {
        datafile= paste(path_base,month_list[imonth_list[iemulate]],"/",var_list[ivar_list[iemulate]],"_PPE_NH_SH_difference_",month_list[imonth_list[iemulate]],"_revised_match_MODIS_CF.dat",sep='')
   }
} else if (constraint_itype_list[iemulate]==3) { ## transect
    datafile= paste(path_base,transect_month_list[imonth_list[iemulate]],"/transect_gradients/",transect_var_list[ivar_list[iemulate]],"_",transect_region_list[iregion_list[iemulate]],"_trajectory_gradients_meters_revised_MODIS_CF.dat",sep='')
}


## Scan goes through the rows first so when you put them back into the matrix goes across the rows too.

data_to_emulate= scan(datafile)

###########################################
## Read in uniformly distributed parameter values
## subset of parameter values to match 
## final set perturbed, not those originally considered
###########################################

par_names_original= c('bl_nuc','ait_width','cloud_ph','carb_ff_ems_eur','carb_ff_ems_nam','carb_ff_ems_chi','carb_ff_ems_asi','carb_ff_ems_mar','carb_ff_ems_r','carb_bb_ems_sam','carb_bb_ems_naf','carb_bb_ems_saf','carb_bb_ems_bnh','carb_bb_ems_rnh','carb_bb_ems_rsh','carb_res_ems_chi','carb_res_ems_asi','carb_res_ems_afr','carb_res_ems_lat','carb_res_ems_r','carb_ff_diam','carb_bb_diam','carb_res_diam','prim_so4_diam','sea_spray','anth_so2_chi','anth_so2_asi','anth_so2_eur','anth_so2_nam','anth_so2_r','volc_so2','bvoc_soa','dms','prim_moc','dry_dep_ait','dry_dep_acc','dry_dep_so2','kappa_oc','sig_w','rain_frac','cloud_ice_thresh','convective_plume_scavenging','scav_diam','bc_ri','oxidants_oh','oxidants_o3','bparam','two_d_fsd_factor','c_r_correl','autoconv_exp_lwp','autoconv_exp_nd','dbsdtbs_turb_0','ai','m_ci','a_ent_1_rp')

par_names= c('bl_nuc','ait_width','cloud_ph','carb_ff_ems_eur','carb_ff_ems_nam','carb_ff_ems_chi','carb_ff_ems_asi','carb_ff_ems_mar','carb_ff_ems_r','carb_bb_ems_sam','carb_bb_ems_naf','carb_bb_ems_saf','carb_bb_ems_bnh','carb_bb_ems_rnh','carb_bb_ems_rsh','carb_res_ems_chi','carb_res_ems_asi','carb_res_ems_afr','carb_res_ems_lat','carb_res_ems_r','carb_ff_diam','carb_bb_diam','carb_res_diam','prim_so4_diam','sea_spray','anth_so2_chi','anth_so2_asi','anth_so2_eur','anth_so2_nam','anth_so2_r','volc_so2','bvoc_soa','dms','prim_moc','dry_dep_ait','dry_dep_acc','dry_dep_so2','kappa_oc','sig_w','rain_frac','cloud_ice_thresh','convective_plume_scavenging','bc_ri','oxidants_oh','oxidants_o3','bparam','two_d_fsd_factor','c_r_correl','autoconv_exp_lwp','autoconv_exp_nd','dbsdtbs_turb_0','ai','m_ci','a_ent_1_rp')

par_names_w_o_carb=c('bl_nuc','ait_width','cloud_ph','carb_ff_diam','carb_bb_diam','carb_res_diam','prim_so4_diam','sea_spray','anth_so2_chi','anth_so2_asi','anth_so2_eur','anth_so2_nam','anth_so2_r','volc_so2','bvoc_soa','dms','prim_moc','dry_dep_ait','dry_dep_acc','dry_dep_so2','kappa_oc','sig_w','rain_frac','cloud_ice_thresh','convective_plume_scavenging','bc_ri','oxidants_oh','oxidants_o3','bparam','two_d_fsd_factor','c_r_correl','autoconv_exp_lwp','autoconv_exp_nd','dbsdtbs_turb_0','ai','m_ci','a_ent_1_rp')

desfile= paste("/gws/nopw/j04/acure/lregayre/emulation/UKESM_PPE_Unit_emulation_ready.dat",sep='')
design_values_original=as.matrix(read.table(desfile))

nppe_design=dim(design_values_original)[1]
npar_design=dim(design_values_original)[2]

#############################
## remove header and scav_diam
#############################

design_values_no_header = matrix(design_values_original, ncol = npar_design, dimnames = NULL)

design_values_strings= cbind(design_values_no_header[2:222,1:42],design_values_no_header[2:222,44:55])
design_values=matrix(0,dim(design_values_strings)[1],dim(design_values_strings)[2])

design_values_strings_w_o_carb= cbind(design_values_no_header[2:222,1:3],design_values_no_header[2:222,21:42],design_values_no_header[2:222,44:55])
design_values_w_o_carb=matrix(0,dim(design_values_strings_w_o_carb)[1],dim(design_values_strings_w_o_carb)[2])

for (i in 1:dim(design_values)[1]) {
  for (j in 1:dim(design_values)[2]) {
    design_values[i,j]=as.numeric(design_values_strings[i,j])
  }
}


for (i in 1:dim(design_values_w_o_carb)[1]) {
  for (j in 1:dim(design_values_w_o_carb)[2]) {
    design_values_w_o_carb[i,j]=as.numeric(design_values_strings_w_o_carb[i,j])
  }
}



###################################
## Read in large sample parameter values
###################################


sample_path= '/gws/nopw/j04/acure/lregayre/emulation/samples/'
if (sample_size==1) {
  sample_file='constrained_near_million_sample.RData'
} else if (sample_size==2) {
  sample_file='constrained_multi_million_sample.RData'
} else {
  sample_file='constrained_sample.RData'
}

load(paste(sample_path,'/',sample_file,sep=''))

if (sample_size==1) {
  sample_values_original=near_million_sample_constrained
} else if (sample_size==2) {
    sample_values_original=multi_million_sample_constrained
} else {
   sample_values_original=sample_constrained
}

nsample=dim(sample_values_original)[1]

if (reduce_1M==0) {
  sample_values= cbind(sample_values_original[1:nsample,1:42],sample_values_original[1:nsample,44:55])
  sample_values_w_o_carb=cbind(sample_values_original[1:nsample,1:3],sample_values_original[1:nsample,21:42],sample_values_original[1:nsample,44:55])
} else {
  sample_values= cbind(sample_values_original[1:reduce_size,1:42],sample_values_original[1:reduce_size,44:55])
  sample_values_w_o_carb=cbind(sample_values_original[1:reduce_size,1:3],sample_values_original[1:reduce_size,21:42],sample_values_original[1:reduce_size,44:55])
}


############################
## Emulate and test emulator
## & save to file
############################

R_file=paste('/gws/nopw/j04/acure/lregayre/emulation/emulators_RData/emulator_for_satellite_constraint_number_',as.character(iemulate),'_w_and_w_o_carb.RData',sep='')

if (file.exists(R_file)) {
  load(R_file)
} else {
  m_w_o_carb=km(~.,design=data.frame(design_values_w_o_carb),response=data_to_emulate,covtype='matern5_2',control=list(maxit=500))
  pdf(paste('/gws/nopw/j04/acure/lregayre/emulation/LOOut/emulator_LOOut_',as.character(iemulate),'_w_o_carb.pdf',sep=''))
  plot(m_w_o_carb)
  dev.off()
  save(m_w_o_carb, file=paste('/gws/nopw/j04/acure/lregayre/emulation/emulators_RData/emulator_for_satellite_constraint_number_',as.character(iemulate),'_w_and_w_o_carb.RData',sep=''))
}



############################
## Predict from emulators using sample
############################


if (sample_size==1) {
  p_w_o_carb=predict.km(m_w_o_carb,newdata=data.frame(sample_values_w_o_carb),'UK',checknames=F)
} else if (sample_size==0) {
  p_w_o_carb=predict.km(m_w_o_carb,newdata=data.frame(sample_values_w_o_carb),'UK',checknames=F)
} else { ## multi-million, so read in large chunks of 100k
  source('/gws/nopw/j04/acure/lregayre/emulation/JJ_code/JJCode_EmPred_LargeSample.r') ## This functionality of this code can be provided on request. Alternatively, predict at sample points using HPC
  p_w_o_carb=JJCode_PredictFromEm_UsingLargeSample(EmModIn=m_w_o_carb,LargeSampInputCombs=sample_values_w_o_carb,nPredBreakVal=10000,PredMean=TRUE,PredSD=TRUE,Pred95CIs=FALSE)
}

mean_file_w_o_carb=paste('/gws/nopw/j04/acure/lregayre/data_post_emulation/emulated_mean_values_',constraint_type_text[iemulate],'_',as.character(length(p_w_o_carb$mean)),'_w_o_carb.dat',sep='')
sd_file_w_o_carb=paste('/gws/nopw/j04/acure/lregayre/data_post_emulation/emulated_sd_values_',constraint_type_text[iemulate],'_',as.character(length(p_w_o_carb$mean)),'_w_o_carb.dat',sep='')
write_lines(p_w_o_carb$mean,mean_file_w_o_carb)
write_lines(p_w_o_carb$sd,sd_file_w_o_carb)


