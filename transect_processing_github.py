#!/usr/bin/env python
"""
Extract data on specified transects for analysis
for ONE of multiple transect variables for a SINGLE MONTH
from monthly mean and annual mean .nc files.

"""

import numpy as np
import matplotlib
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.lines as mlines
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

ivar=0 # choice of 0-4 transect variables
imonth=7  ## 11=Nov, 7=Jul

var_folder_list=['Cloud_drop_number_concentration','Reff','AOD','LWP','Cloud_Fraction'] ## AOD extracts 440 and 550 below to calculate AI
var_source_list=['CDNC','Reff','AOD','LWP','CF']
var_name_list=['liquid_cloud_droplet_number_concentration_at_cloud_top_weighted_by_1_299','cosp_modis_weighted_reff_liquid','AOD','cosp_modis_weighted_liquid_water_path','cosp_modis_weighted_liquid_cloud_fraction']

var_name_out_list=['CDNC','Reff_Liquid_MODIS','AI','LWP','LCF']
divisor_flag_list=[1,1,0,1,3]
var_name_divisor_list=['weight_for_cdnc_at_cloud_top','cosp_modis_weighted_liquid_cloud_fraction','','cosp_modis_weighted_liquid_cloud_fraction','cosp_cloud_weights']

var_name_numerator_list_3=['','','','','cosp:_modis_liquid_cloud_fraction___']
var_name_divisor_folder_list_3=['','','','','m01s02i330_COSP:_ISCCPdivMISRdivMODIS_CLOUD_WEIGHTS']
var_name_divisor_list_3=['','','','','cosp:_isccpdivmisrdivmodis_cloud_weights']
var_mask_name_list=['','Cloud_Particle_Size_Liquid','','Cloud_Water_Path_Liquid','Cloud_Retrieval_Fraction_Liquid']
mask_flag_list=[0,1,0,1,1]


var_folder= var_folder_list[ivar]
var_source= var_source_list[ivar]
var_name= var_name_list[ivar]
var_name_out= var_name_out_list[ivar]
divisor_flag= divisor_flag_list[ivar]
var_name_divisor=var_name_divisor_list[ivar]

var_name_divisor_3= var_name_divisor_list_3[ivar]
var_name_numerator_3= var_name_numerator_list_3[ivar]
var_folder_divisor_3=var_name_divisor_folder_list_3[ivar]
var_mask=var_mask_name_list[ivar]
mask_flag=mask_flag_list[ivar]


month_list=['dec','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','ann']
nmonths=len(month_list)
month= month_list[imonth]
dir_base= '/gws/nopw/j04/acure/phase3/netcdf_output/'
dir_base_cloud_weights='/gws/nopw/j04/acure/lregayre/DG_job_output/'

nlat=144
nlon=192


month_list_default=['2016dec','2017jan','2017feb','2017mar','2017apr','2017may','2017jun','2017jul','2017aug','2017sep','2017oct','2017nov']##,'ann']
nmonths_default= len(month_list_default)
month_default=month_list_default[imonth]
default_extension_list=['m01s01i298','m01s01i245_product_of_effective_radius_of_cloud_liquid_water_particle_and_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask','m01s02i295_atmosphere_optical_thickness_due_to_fossil_fuel_organic_carbon_ambient_aerosol','m01s30i405_TOTAL_COLUMN_QCL__RHO_GRID__________','m01s02i452_COSP:_MODIS_LIQUID_CLOUD_FRACTION___']
default_divisor_list=['m01s01i299','m01s02i452_COSP:_MODIS_LIQUID_CLOUD_FRACTION___','','m01s02i452_COSP:_MODIS_LIQUID_CLOUD_FRACTION___','m01s02i330_COSP:_ISCCPdivMISRdivMODIS_CLOUD_WEIGHTS']


#######################
## Default simulation values
#######################

default_file_LR='/gws/nopw/j04/acure/lregayre/default_job_output/2017'+month+'/2017'+month+'_'+default_extension_list[ivar]+'.nc'
default_file_DG='/gws/nopw/j04/acure/lregayre/DG_job_output/'+month_default+'/'+month_default+'_'+default_extension_list[ivar]+'.nc'
default_file_index=[1,1,2,1,1]

if default_file_index[ivar]==1:
    default_divisor_file='/gws/nopw/j04/acure/lregayre/DG_job_output/'+month_default+'/'+month_default+'_'+default_divisor_list[ivar]+'.nc'


var_default=np.zeros((nlat,nlon))
if (default_file_index[ivar]==2):
    default_dataframe=nd(default_file_LR,'r','NETCDF4')
else:
    default_dataframe=nd(default_file_DG,'r','NETCDF4')
    default_divisor_dataframe=nd(default_divisor_file,'r','NETCDF4')


## Default 0 (Nd) 
if (ivar==0):
    var_default_numerator=default_dataframe.variables['unknown'][:]
    var_default_denominator=default_divisor_dataframe.variables['unknown'][:]
    var_default_numerator[np.where(var_default_denominator!=0.0)] = var_default_numerator[np.where(var_default_denominator!=0.0)] / var_default_denominator[np.where(var_default_denominator!=0.0)]
    var_default=var_default_numerator


## default 1 (Re) 
if (ivar==1):
    var_default_numerator=default_dataframe.variables['product_of_effective_radius_of_cloud_liquid_water_particle_and_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask'][:]
    var_default_denominator=default_divisor_dataframe.variables['cosp:_modis_liquid_cloud_fraction___'][:]
    var_default_numerator[np.where(var_default_denominator!=0.0)] = var_default_numerator[np.where(var_default_denominator!=0.0)] / var_default_denominator[np.where(var_default_denominator!=0.0)]
    var_default=var_default_numerator

## AOD -> AI
if (ivar==2):
    var_default_1_440=default_dataframe.variables['atmosphere_optical_thickness_due_to_fossil_fuel_organic_carbon_ambient_aerosol'][0,:,:]
    var_default_1_550=default_dataframe.variables['atmosphere_optical_thickness_due_to_fossil_fuel_organic_carbon_ambient_aerosol'][1,:,:]
    default_file_LR_2='/gws/nopw/j04/acure/lregayre/default_job_output/2017'+month+'/2017'+month+'_m01s02i285_atmosphere_optical_thickness_due_to_mineral_dust_aerosol_in_radiation.nc'
    default_dataframe_2=nd(default_file_LR_2,'r','NETCDF4')
    var_default_2_440=default_dataframe_2.variables['oddust_diag'][0,:,:]
    var_default_2_550=default_dataframe_2.variables['oddust_diag'][1,:,:]
    default_file_LR_3='/gws/nopw/j04/acure/lregayre/default_job_output/2017'+month+'/2017'+month+'_m01s02i303_atmosphere_optical_thickness_due_to_insoluble_aitken_mode_sulphate_aerosol.nc'
    default_dataframe_3=nd(default_file_LR_3,'r','NETCDF4')
    var_default_3_440=default_dataframe_3.variables['odaitins'][0,:,:]
    var_default_3_550=default_dataframe_3.variables['odaitins'][1,:,:]
    default_file_LR_4='/gws/nopw/j04/acure/lregayre/default_job_output/2017'+month+'/2017'+month+'_m01s02i302_atmosphere_optical_thickness_due_to_soluble_coarse_mode_sulphate_aerosol.nc'
    default_dataframe_4=nd(default_file_LR_4,'r','NETCDF4')
    var_default_4_440=default_dataframe_4.variables['odcorsol'][0,:,:]
    var_default_4_550=default_dataframe_4.variables['odcorsol'][1,:,:]
    default_file_LR_5='/gws/nopw/j04/acure/lregayre/default_job_output/2017'+month+'/2017'+month+'_m01s02i301_atmosphere_optical_thickness_due_to_soluble_accumulation_mode_sulphate_aerosol.nc'
    default_dataframe_5=nd(default_file_LR_5,'r','NETCDF4')
    var_default_5_440=default_dataframe_5.variables['odaccsol'][0,:,:]
    var_default_5_550=default_dataframe_5.variables['odaccsol'][1,:,:]
    default_file_LR_6='/gws/nopw/j04/acure/lregayre/default_job_output/2017'+month+'/2017'+month+'_m01s02i300_atmosphere_optical_thickness_due_to_soluble_aitken_mode_sulphate_aerosol.nc'
    default_dataframe_6=nd(default_file_LR_6,'r','NETCDF4')
    var_default_6_440=default_dataframe_6.variables['odaitsol'][0,:,:]
    var_default_6_550=default_dataframe_6.variables['odaitsol'][1,:,:]


### Combine AOD elements and calculate AI
if (ivar==2):
    var_default_550=var_default_1_550+var_default_2_550+var_default_3_550+var_default_4_550+var_default_5_550+var_default_6_550
    var_default_440=var_default_1_440+var_default_2_440+var_default_3_440+var_default_4_440+var_default_5_440+var_default_6_440
    var_default_angstrom=np.zeros((nlat,nlon))
    var_default_AI=np.zeros((nlat,nlon))
    for ilat in np.arange(nlat):
        for ilon in np.arange(nlon):
            if (var_default_550[ilat,ilon]==0.0):
                var_default_angstrom[ilat,ilon]=0.0
            elif (var_default_440[ilat,ilon]==0.0):
                var_default_angstrom[ilat,ilon]=0.0
            else:
                var_default_angstrom[ilat,ilon]=-np.log(var_default_550[ilat,ilon]/var_default_440[ilat,ilon])/np.log(550/440) ## Angstrom exponent
            var_default_AI[ilat,ilon]= var_default_550[ilat,ilon]*var_default_angstrom[ilat,ilon]
    var_default=var_default_AI


## default 3 (LWP)
if (ivar==3):
    var_default_numerator=default_dataframe.variables['atmosphere_cloud_liquid_water_content'][:]*1000 ##to match LWP measurements
    var_default_denominator=default_divisor_dataframe.variables['cosp:_modis_liquid_cloud_fraction___'][:]
    var_default_numerator[np.where(var_default_denominator!=0.0)] = var_default_numerator[np.where(var_default_denominator!=0.0)] / var_default_denominator[np.where(var_default_denominator!=0.0)]
    var_default=var_default_numerator


## Default 4 (fc)
if (ivar==4):
    var_default_numerator=default_dataframe.variables['cosp:_modis_liquid_cloud_fraction___'][:]
    var_default_denominator=default_divisor_dataframe.variables['cosp:_isccpdivmisrdivmodis_cloud_weights'][:]
    var_default_numerator[np.where(var_default_denominator!=0.0)] = var_default_numerator[np.where(var_default_denominator!=0.0)] / var_default_denominator[np.where(var_default_denominator!=0.0)]
    var_default=var_default_numerator


###############
## PPE members
###############


print(month)
if np.logical_or(np.logical_or((month=='dec'),month=='jan'),month=='feb'):
    nPPE='220'
    nppe=220
else:
    nPPE='221'
    nppe=221


file_name= var_folder+'/monthly/'+month+'/ACURE_P3_'+var_source+'_'+month+'_PD_'+nPPE+'.nc'
dataframe= nd(dir_base+file_name,'r','NETCDF4')
if (ivar!=2):
    var_array=dataframe.variables[var_name][:]
    if (ivar==3):
        file_name_divisor=var_folder_list[ivar+1]+'/monthly/'+month+'/ACURE_P3_'+var_source_list[ivar+1]+'_'+month+'_PD_'+nPPE+'.nc'
        dataframe_divisor=nd(dir_base+file_name_divisor,'r','NETCDF4')
        var_weights_array=dataframe_divisor[var_name_divisor][:]
        var_array[np.where(var_weights_array!=0.0)] = var_array[np.where(var_weights_array!=0.0)] / var_weights_array[np.where(var_weights_array!=0.0)]
    elif (ivar==4):
        dataframe_divisor= nd(dir_base_cloud_weights+'/'+month_default+'/'+month_default+'_m01s02i330_COSP:_ISCCPdivMISRdivMODIS_CLOUD_WEIGHTS.nc','r','NETCDF4')
        var_weights_array=dataframe_divisor['cosp:_isccpdivmisrdivmodis_cloud_weights'][:]
        for ippe in np.arange(nppe):
            var_array_local=var_array[ippe,:]
            print(ippe)
            for ilat in np.arange(nlat):
                for ilon in np.arange(nlon):
                    if (var_weights_array[ilat,ilon]>0.0):
                        var_array[ippe,ilat,ilon]= var_array[ippe,ilat,ilon] / var_weights_array[ilat,ilon]
                    else:
                        var_array[ippe,ilat,ilon]=0.0
    else:
        print('divisor flag stage 1')
        var_weights_array= dataframe.variables[var_name_divisor][:]
        var_array[np.where(var_weights_array!=0.0)] = var_array[np.where(var_weights_array!=0.0)] / var_weights_array[np.where(var_weights_array!=0.0)]
else: ##AOD
    var_array_1_440=dataframe.variables['atmosphere_optical_thickness_due_to_soluble_accumulation_mode_ambient_aerosol'][:,1,:,:]
    var_array_2_440=dataframe.variables['atmosphere_optical_thickness_due_to_soluble_coarse_mode_ambient_aerosol'][:,1,:,:]
    var_array_3_440=dataframe.variables['atmosphere_optical_thickness_due_to_insoluble_aitken_mode_ambient_aerosol'][:,1,:,:]
    var_array_4_440=dataframe.variables['atmosphere_optical_thickness_due_to_insoluble_accumulation_mode_ambient_aerosol'][:,1,:,:]
    var_array_5_440=dataframe.variables['atmosphere_optical_thickness_due_to_insoluble_coarse_mode_ambient_aerosol'][:,1,:,:]
    var_array_6_440=dataframe.variables['atmosphere_optical_thickness_due_to_dust_ambient_aerosol'][:,1,:,:]
    var_array_1_550=dataframe.variables['atmosphere_optical_thickness_due_to_soluble_accumulation_mode_ambient_aerosol'][:,2,:,:]
    var_array_2_550=dataframe.variables['atmosphere_optical_thickness_due_to_soluble_coarse_mode_ambient_aerosol'][:,2,:,:]
    var_array_3_550=dataframe.variables['atmosphere_optical_thickness_due_to_insoluble_aitken_mode_ambient_aerosol'][:,2,:,:]
    var_array_4_550=dataframe.variables['atmosphere_optical_thickness_due_to_insoluble_accumulation_mode_ambient_aerosol'][:,2,:,:]
    var_array_5_550=dataframe.variables['atmosphere_optical_thickness_due_to_insoluble_coarse_mode_ambient_aerosol'][:,2,:,:]
    var_array_6_550=dataframe.variables['atmosphere_optical_thickness_due_to_dust_ambient_aerosol'][:,2,:,:]


lat_list=dataframe.variables['latitude'][:]
nlat=lat_list.shape[0]
lon_list=dataframe.variables['longitude'][:]
nlon=lon_list.shape[0]



if (ivar==2):
    var_array_440=np.zeros((nppe,nlat,nlon))
    var_array_550=np.zeros((nppe,nlat,nlon))
    var_array_angstrom=np.zeros((nppe,nlat,nlon))
    var_array=np.zeros((nppe,nlat,nlon))
    for ippe in np.arange(nppe):
        print(ippe)
        for ilat in np.arange(nlat):
            for ilon in np.arange(nlon):
                var_array_440[ippe,ilat,ilon]=var_array_1_440[ippe,ilat,ilon]+var_array_2_440[ippe,ilat,ilon]+var_array_3_440[ippe,ilat,ilon]+var_array_4_440[ippe,ilat,ilon]+var_array_5_440[ippe,ilat,ilon]+var_array_6_440[ippe,ilat,ilon]
                var_array_550[ippe,ilat,ilon]=var_array_1_550[ippe,ilat,ilon]+var_array_2_550[ippe,ilat,ilon]+var_array_3_550[ippe,ilat,ilon]+var_array_4_550[ippe,ilat,ilon]+var_array_5_550[ippe,ilat,ilon]+var_array_6_550[ippe,ilat,ilon]
                if (var_array_550[ippe,ilat,ilon]==0.0):
                    var_array_angstrom[ippe,ilat,ilon]=0.0
                elif (var_array_440[ippe,ilat,ilon]==0.0):
                    var_array_angstrom[ippe,ilat,ilon]=0.0
                else:
                    var_array_angstrom[ippe,ilat,ilon]=-np.log(var_array_550[ippe,ilat,ilon]/var_array_440[ippe,ilat,ilon])/np.log(550/440)
                var_array[ippe,ilat,ilon]= var_array_550[ippe,ilat,ilon]*var_array_angstrom[ippe,ilat,ilon]


if (nPPE=='220'):
    var_array[215,:,:]= np.mean(var_array,axis=0)
    if (ivar==2):
        var_array_550[215,:,:]=np.mean(var_array_550,axis=0)


if (var_name=='cosp_modis_weighted_reff_liquid'):
    var_array=var_array*1000000 ## Rescaling by 10^6, since model scaled diagnostic to avoid packing problems. Now comparable to measurements
elif (var_name=='cosp_modis_weighted_liquid_water_path'):
    var_array=var_array*1000 ## converting from kg/m2 to g/m2 to match measurements



#############################################################
## masks created in regional_mean_MODIS_observation_data_github.py
## read in for consistent comparison
#############################################################


month_list_in=['dec2016', 'jan2017', 'feb2017', 'mar2017', 'apr2017', 'may2017', 'jun2017', 'jul2017', 'aug2017', 'sep2017', 'oct2017', 'nov2017']
nmonths_in= len(month_list_in)
if mask_flag==1:
    print(imonth)
    month_default=month_list_in[imonth]
    file_in= '/gws/nopw/j04/acure/lregayre/masks/'+var_mask+'/'+var_mask+'_monthly_mask_N96_same_lat_lon_'+month_default+'.dat'
    var_mask_array=np.loadtxt(file_in)



######################################
## South Pacific transect (South American)
######################################

SA_trajectory_lon_index=[150,149,148,147,146,145,144,143] ## specific gridbox indices from lat_list and lon_list
SA_trajectory_lat_index=[55,56,57,57,58,58,59,59]

if mask_flag==1:
    SA_mask= var_mask_array[SA_trajectory_lat_index,SA_trajectory_lon_index]
    print(SA_mask)

SA_trajectory_var= var_array[:,SA_trajectory_lat_index,SA_trajectory_lon_index]
SA_trajectory_default= var_default[SA_trajectory_lat_index,SA_trajectory_lon_index]
if (ivar==2):
    SA_trajectory_var550= var_array_550[:,SA_trajectory_lat_index,SA_trajectory_lon_index]
    SA_trajectory_default550= var_default_550[SA_trajectory_lat_index,SA_trajectory_lon_index]


standard_cf = 111000
SA_transect_lat_differences=[]
SA_transect_lon_differences=[]
SA_transect_hypotenuse=[0] ## in meters
for ilat in np.arange(len(SA_trajectory_lat_index)-1): ## ilat/ilon are paired here, so only incrementing one of them
    local_lat_diff=np.abs(lat_list[SA_trajectory_lat_index[ilat]]-lat_list[SA_trajectory_lat_index[ilat+1]])*standard_cf
    local_lon_diff=np.abs(lon_list[SA_trajectory_lon_index[ilat]]-lon_list[SA_trajectory_lon_index[ilat+1]])*standard_cf
    SA_transect_hypotenuse= np.append(SA_transect_hypotenuse,np.sqrt(local_lat_diff**2+local_lon_diff**2)+SA_transect_hypotenuse[-1])


#############
## Write to file
#############

outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/'+var_name_out+'_SouthAmerican_trajectory_'+month+'_revised_MODIS_CF.dat'
np.savetxt(outfile,SA_trajectory_var)
outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/'+var_name_out+'_SouthAmerican_trajectory_'+month+'_default_revised_MODIS_CF.dat'
np.savetxt(outfile,SA_trajectory_default)
if (ivar==2):
    outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/AOD_550_SouthAmerican_trajectory_'+month+'_revised_MODIS_CF.dat'
    np.savetxt(outfile,SA_trajectory_var550)
    outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/AOD_550_SouthAmerican_trajectory_'+month+'_default_revised_MODIS_CF.dat'
    np.savetxt(outfile,SA_trajectory_default550)

if (ivar==0):
    outfile = '/gws/nopw/j04/acure/lregayre/transects/distance_in_meters_along_SouthAmerican_transect_revised_MODIS_CF.dat'
    np.savetxt(outfile,SA_transect_hypotenuse)


#################################
## South Atlantic (Namibian) transect
#################################


SAf_trajectory_lon_index=[190,189,188,187,186,185,184] 
SAf_trajectory_lat_index=[62,62,62,62,62,62,62]


if mask_flag==1:
    SAf_mask= var_mask_array[SAf_trajectory_lat_index,SAf_trajectory_lon_index]
    print(SAf_mask)

SAf_trajectory_var= var_array[:,SAf_trajectory_lat_index,SAf_trajectory_lon_index]
SAf_trajectory_default= var_default[SAf_trajectory_lat_index,SAf_trajectory_lon_index]
if (ivar==2):
    SAf_trajectory_var550= var_array_550[:,SAf_trajectory_lat_index,SAf_trajectory_lon_index]
    SAf_trajectory_default550= var_default_550[SAf_trajectory_lat_index,SAf_trajectory_lon_index]


lat_list[SAf_trajectory_lat_index]

SAf_transect_lat_differences=[]
SAf_transect_lon_differences=[]
SAf_transect_hypotenuse=[0] ## in meters
for ilat in np.arange(len(SAf_trajectory_lat_index)-1): ## ilat/ilon are paired here, so only incrementing one of them
    local_lat_diff=np.abs(lat_list[SAf_trajectory_lat_index[ilat]]-lat_list[SAf_trajectory_lat_index[ilat+1]])*standard_cf
    local_lon_diff=np.abs(lon_list[SAf_trajectory_lon_index[ilat]]-lon_list[SAf_trajectory_lon_index[ilat+1]])*standard_cf
    print(np.str(local_lon_diff),np.str(local_lat_diff))
    SAf_transect_hypotenuse= np.append(SAf_transect_hypotenuse,np.sqrt(local_lat_diff**2+local_lon_diff**2)+SAf_transect_hypotenuse[-1])


############
## Write to file
############

outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/'+var_name_out+'_Namibian_trajectory_'+month+'_revised_MODIS_CF.dat'
np.savetxt(outfile,SAf_trajectory_var)
outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/'+var_name_out+'_Namibian_trajectory_'+month+'_default_revised_MODIS_CF.dat'
np.savetxt(outfile,SAf_trajectory_default)
if (ivar==2):
    outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/AOD_550_Namibian_trajectory_'+month+'_revised_MODIS_CF.dat'
    np.savetxt(outfile,SAf_trajectory_var550)
    outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/AOD_550_Namibian_trajectory_'+month+'_default_revised_MODIS_CF.dat'
    np.savetxt(outfile,SAf_trajectory_default550)



if (ivar==0):
    outfile = '/gws/nopw/j04/acure/lregayre/transects/distance_in_meters_along_Namibian_transect_revised_MODIS_CF.dat'
    np.savetxt(outfile,SAf_transect_hypotenuse)



#############################################
## North Pacific (North American) transect
#############################################

NA_trajectory_lon_index=[122,122,122,122,122,121,121,121,121,121]
NA_trajectory_lat_index=[96,95,94,93,92,91,90,89,88,87]

if mask_flag==1:
    NA_mask= var_mask_array[NA_trajectory_lat_index,NA_trajectory_lon_index]
    print(NA_mask)


NA_trajectory_var= var_array[:,NA_trajectory_lat_index,NA_trajectory_lon_index]
NA_trajectory_default= var_default[NA_trajectory_lat_index,NA_trajectory_lon_index]
if (ivar==2):
    NA_trajectory_var550= var_array_550[:,NA_trajectory_lat_index,NA_trajectory_lon_index]
    NA_trajectory_default550= var_default_550[NA_trajectory_lat_index,NA_trajectory_lon_index]



lat_list[NA_trajectory_lat_index]

adjusted_cf= 110.773 ## based on 25 degrees N (http://www.csgnetwork.com/degreelenllavcalc.html)
NA_transect_lat_differences=[]
NA_transect_lon_differences=[]
NA_transect_hypotenuse=[0] ## in meters
for ilat in np.arange(len(NA_trajectory_lat_index)-1): ## ilat/ilon are paired here, so only incrementing one of them
    local_lat_diff=np.abs(lat_list[NA_trajectory_lat_index[ilat]]-lat_list[NA_trajectory_lat_index[ilat+1]])*standard_cf
    local_lon_diff=np.abs(lon_list[NA_trajectory_lon_index[ilat]]-lon_list[NA_trajectory_lon_index[ilat+1]])*adjusted_cf
    print(np.str(local_lon_diff),np.str(local_lat_diff))
    NA_transect_hypotenuse= np.append(NA_transect_hypotenuse,np.sqrt(local_lat_diff**2+local_lon_diff**2)+NA_transect_hypotenuse[-1])


############
## Write to file
############

outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/'+var_name_out+'_NorthAmerican_trajectory_'+month+'_revised_MODIS_CF.dat'
np.savetxt(outfile,NA_trajectory_var)
outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/'+var_name_out+'_NorthAmerican_trajectory_'+month+'_default_revised_MODIS_CF.dat'
np.savetxt(outfile,NA_trajectory_default)
if (ivar==2):
    outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/AOD_550_NorthAmerican_trajectory_'+month+'_revised_MODIS_CF.dat'
    np.savetxt(outfile,NA_trajectory_var550)
    outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/AOD_550_NorthAmerican_trajectory_'+month+'_default_revised_MODIS_CF.dat'
    np.savetxt(outfile,NA_trajectory_default550)



if (ivar==0):
    outfile = '/gws/nopw/j04/acure/lregayre/transects/distance_in_meters_along_NorthAmerican_transect_revised_MODIS_CF.dat'
    np.savetxt(outfile,NA_transect_hypotenuse)



#################################################
## North Atlantic (European/EU) transect
#################################################

EU_trajectory_lon_index=[179,179,178,178,177,177,176,176]
EU_trajectory_lat_index=[115,114,113,112,111,110,109,108]

if mask_flag==1:
    EU_mask= var_mask_array[EU_trajectory_lat_index,EU_trajectory_lon_index] 
    print(EU_mask)

EU_trajectory_var= var_array[:,EU_trajectory_lat_index,EU_trajectory_lon_index]
EU_trajectory_default= var_default[EU_trajectory_lat_index,EU_trajectory_lon_index]
if (ivar==2):
    EU_trajectory_var550= var_array_550[:,EU_trajectory_lat_index,EU_trajectory_lon_index]
    EU_trajectory_default550= var_default_550[EU_trajectory_lat_index,EU_trajectory_lon_index]

EU_trajectory_default


lat_list[EU_trajectory_lat_index]

adjusted_cf= 110.773 ## based on 25 degrees N (http://www.csgnetwork.com/degreelenllavcalc.html)
EU_transect_lat_differences=[]
EU_transect_lon_differences=[]
EU_transect_hypotenuse=[0] ## in meters;
for ilat in np.arange(len(EU_trajectory_lat_index)-1): ## ilat/ilon are paired here, so only incrementing one of them
    local_lat_diff=np.abs(lat_list[EU_trajectory_lat_index[ilat]]-lat_list[EU_trajectory_lat_index[ilat+1]])*standard_cf
    local_lon_diff=np.abs(lon_list[EU_trajectory_lon_index[ilat]]-lon_list[EU_trajectory_lon_index[ilat+1]])*adjusted_cf
    print(np.str(local_lon_diff),np.str(local_lat_diff))
    EU_transect_hypotenuse= np.append(EU_transect_hypotenuse,np.sqrt(local_lat_diff**2+local_lon_diff**2)+EU_transect_hypotenuse[-1])


############
## Write to file
############

outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/'+var_name_out+'_European_trajectory_'+month+'_revised_MODIS_CF.dat'
np.savetxt(outfile,EU_trajectory_var)
outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/'+var_name_out+'_European_trajectory_'+month+'_default_revised_MODIS_CF.dat'
np.savetxt(outfile,EU_trajectory_default)
if (ivar==2):
    outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/AOD_550_European_trajectory_'+month+'_revised_MODIS_CF.dat'
    np.savetxt(outfile,EU_trajectory_var550)
    outfile = '/gws/nopw/j04/acure/lregayre/data_for_emulation/'+month+'/trajectories/AOD_550_European_trajectory_'+month+'_default_revised_MODIS_CF.dat'
    np.savetxt(outfile,EU_trajectory_default550)


if (ivar==0):
    outfile = '/gws/nopw/j04/acure/lregayre/transects/distance_in_meters_along_European_transect_revised_MODIS_CF.dat'
    np.savetxt(outfile,EU_transect_hypotenuse)


