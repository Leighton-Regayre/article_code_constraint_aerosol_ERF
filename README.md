# article_code_constraint_aerosol_ERF
Data and code for the analysis and figure creation of the article "Identifying climate model structural inconsistencies allows for tight constraint of aerosol radiative forcing" by Regayre et al.

Monthly mean data from the perturbed parameter ensemble (PPE) is available at https://catalogue.ceda.ac.uk/uuid/b735718d66c1403fbf6b93ba3bd3b1a9 and observational data sources are described in the article.

All code here is written in python, with the exception of statistical emulation code, which is in R.

The pertrubed parameter design table is labelled:

*ppe_dataframe.csv*

The values in this file are used by UKESM1 suite u-bs714 and others to produce the PPE members.


Model and observation data were processed to create regional mean values for analysis using:

XXX

A map showing transects was created using:

XXX

and transect variables were calculated using:

XXX


Statistical emulation of PPE output was performed using:

XXX

Indices in this code can be changed to emulate distinct constraint variables

The emulation code makes use of a large sample of parameter combinations to extend the 221 PPE members to 1 million model variants. This sample is labelled:

XXX

The emulation code also makes use of functions that increase the efficiency of creation of covariance matrcies that are the intellectual property of colleagues and can be provided upon request. This code can be run without relying on the underlying routines on a suitable HPC machine.

Probability distributions of aerosol radiative forcings were created using:

XXX


Parameter relative importance values were plotted using:

XXX

XXX

XXX

These 3 files were used to create all relative importance figures in the article


We plotted the seasonal cycles of state variables using:

XXX

XXX

This enhanced version of this code exemplifies the effect of constraint to distinct North Atlantic state variables:

XXX


Our method for progressively adding constraint variables is here:

XXX

The percentage retained can be adjusted to evaluate sensitivity to this subjective choice.

Pairwise-consistency of multiple constraint variables in each region were created using:

XXX

Regions can be specified using indices in the file
This enhanced version for North Atlantic variables includes a subpanel with distinct constraints of F_SW as shown in the main article:

XXX


The effect of progressivley adding additional constraint variables, to identify an optimal constraint with this model and observational data set are visualised using:

XXX






