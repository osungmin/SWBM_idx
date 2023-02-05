# SWBM_idx

This is a collection of example python scripts to run the rainfall-runoff model SWBM at individual pixels.


See more details about the SWBM.

More studies where SWBM is used. 
 

## preparation

 1. prep_runoff.py: to prepare the list of lat/lon at a desire resolution (e.g. 0.1 or 0.25 deg) and runoff observation at each grid pixel
 2. extract_vars.py: to extract tempertaure and precipitation data from E-OBS (or other netcdf files) at the pixels defined at 1.
 
## model run

 1. run_swbm.py: a python script to create random parameters (for model calibration), prepare run scripts, and excute the run scripts.
 2. swbm.r: the SWBM model itself, written in R.
