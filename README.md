# article_code_constraint_aerosol_ERF
Data and code for the analysis and figure creation of the article "Identifying climate model structural inconsistencies allows for tight constraint of aerosol radiative forcing" by Regayre et al.

Monthly mean data from the perturbed parameter ensemble (PPE) is available at https://catalogue.ceda.ac.uk/uuid/b735718d66c1403fbf6b93ba3bd3b1a9 and observational data sources are described in the article.

All code here is written in python, with the exception of statistical emulation code, which is in R.

The pertrubed parameter design table is labelled:
*ppe_dataframe.csv*
The values in this file are used by UKESM1 suite u-bs714 and others to produce the PPE members.


Model and observation MODIS data were processed to create regional mean values for analysis using:
*regional_mean_MODIS_observation_data_github.py*
This file exemplifies our method for calculating regional means of observation data. Other (non-MODIS) state variables were calculated in an identical manner, with minor code modifications to match the observational data, as outlined in the article.
Regional means of state variables, for each PPE member, were calculated using:
*regional_mean_model_github.py*
Note that cloud weights for MODIS output were derived from a distinct simulation due to a diagnotstic processing error on ARCHER HPC. All cloud weight files are stored here in the folder:
*cloud_weights*


A map showing transects was created using:
*transect_map_github.py*
and transect variables were calculated using:
*transect_processing_github.py*

Statistical emulation of PPE output was performed using:
*emulate_EXAMPLE_github.R* - note this is R code, not Python code
Indices in this code can be changed to emulate any of the constraint variables described in the article

The emulation code makes use of a large sample of parameter combinations to extend the 221 PPE members to 1 million model variants. This sample is several Tb so can be shared directly on request. The emulation code also makes use of functions that increase the efficiency of creation of covariance matrcies that are the intellectual property of colleagues and can be provided upon request. This code can be run without relying on the underlying routines on a suitable HPC machine.


Probability distributions of aerosol radiative forcings were created using:
*pdfs_ERF_github.py*


Parameter relative importance values were plotted from either the full 1 million member sample (for aerosol ERF), or the PPE members, using:
*relative_importance_sample_github.py*
and either
*relative_importance_state_variables_github.py* or *relative_importance_transect_variables_github.py*
respectively

We plotted the seasonal cycles of state variables using:
*seasonal_12panel_github.py* for the figure in the main article, which has subpanels that exemplify the effects of 2 distinct constraints, and
*seasonal_cycle_6panel_github.py* and *seasonal_cycle_hemispheric_github.py*
for other figures.


Our method for progressively adding constraint variables is exemplified by:
*progressively_add_optimal_github.py*
The percentage retained can be adjusted within this file to evaluate sensitivity to this subjective choice.


Pairwise-consistency of multiple constraint variables in each region were created using:
*pairwise_general_github.py*
and an adaptation to include subpanels that exemplify the effects of constraint to specific constraint variables is:
*pairwise_w_subpanel_github.py*
Both files rely on an NRSME effects array created using:
*NRMSE_effects_array_github.py*


The effect of progressivley adding additional constraint variables, to identify an optimal constraint with this model and observational data set are visualised using:
*optimal_visualisation_github.py*





