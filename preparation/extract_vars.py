#!/usr/bin/env python
"""
extract forcing data at a lat/lon grid point
"""
import os
import pandas as pd
import numpy as np
import datetime
from netCDF4 import Dataset

def load_idxs(): #you create this list from prep_list.py
    path = 'your_path_to_runoff_lat_lon_list/obs_d25.lis'
    load = pd.read_csv(path, header=0, index_col=0, sep=',')
    return (load.index)

def create_daterange(date0, date1, freq):
    start = datetime.datetime.strptime(date0, '%Y-%m-%d')
    end = datetime.datetime.strptime(date1, '%Y-%m-%d')
    daterange366 = pd.date_range(start, end, freq=freq) # D: daily 
    daterange = daterange366[~((daterange366.month == 2) & (daterange366.day == 29))]
    return(daterange, daterange366)


####################
idxs = load_idxs()
lats = [float(x.split("_")[0]) for x in idxs]
lons = [float(x.split("_")[1]) for x in idxs]
####################
daterange = create_daterange("1980-01-01", "2014-12-31", "D")[1] #adjust time range as you want
vars = ["rr", "tg", "snr"] #mm/day, kelven, net radiation, MJ m2 #three inputs for the model
###################
in_path  ="your_path_to_input_netcdfs"
out_path ="your_path_to_save_output"
###################
print()
print("***** extracting forcing data")
print(" > input:", in_path)
print(" > output:",out_path)
print()
print(len(idxs))
wait=input()


for v in vars:
    print()
    print ("  working on variable >>>", v)

    #file name for each varaible
    if v=="rr":  var_path ="rr/rr_ens_mean_0.25deg_reg_v20.0e."
    if v=="tg": var_path ="tg/tg_ens_mean_0.25deg_reg_v20.0e."
    if v=="snr": var_path ="snr/snr.daily.calc.era5.1440.720."

    nc_files=[]

    for yr in np.unique(daterange.year):

        #open netcdf file for each year
        f1 = Dataset(in_path+var_path+str(yr)+".nc")
        nc = (f1.variables[v][:,:,:]).filled(np.nan)
        print(yr, nc.shape)

        if yr==daterange.year[-1]:
           datalons = f1.variables['longitude'][:]
           datalats = f1.variables['latitude'][:]

        nc_files.append(nc)
        f1.close()

    print ("done.:", len(nc_files))
    print (". . . appending")
    append  = np.concatenate(nc_files, axis=0)

    print (". . . netcdf2array ready to extract idxs:", append.shape)
    print ()

    #extracing values
    for idx in range(len(idxs)):

       target={"lat":lats[idx], "lon":lons[idx]}
       print (" idx #", idx, idxs[idx], target, " out of", len(idxs))

       if os.path.exists(out_path+v+"_"+idxs[idx]+".dat"):
          print (" alreay prepared!!!")
          continue

       #output dataframe; i will save extracted variables in this dataframe
       df=pd.DataFrame(index=daterange, columns=[v]) 

       lat_idx = (np.abs(datalats - target["lat"])).argmin()
       lon_idx = (np.abs(datalons - target["lon"])).argmin()
       print ("** check lon/lat in nc:", lat_idx, lon_idx, datalats[lat_idx],datalons[lon_idx])

       df[v]  = append[:,lat_idx,lon_idx]

       if v=="tg": df[v]=df[v]-273.15 #deg C
        
       df = df[~((df.index.month == 2) & (df.index.day == 29))] #to remove feb 29

       #precipitation data in a netcdf file could have negative values due to format issue. 
       if v in ["tp","rr"]:
          df.loc[df[v] < 0.0001, v] =0.
            
       #saving extracted values
       df.to_csv(out_path+v+"/"+v+"_"+idxs[idx]+".dat", 
                 float_format='%.5f', 
                 index=True, header=True, 
                 sep=",", na_rep=-9999) 

       #just to to check
       if (df.count()[0] != 365*len(np.unique(daterange.year)):
           print ("missing data included !!??", idxs[idx], v)
           print (df.describe())
           wait=input()
       
    print ("done.", v)

print("end", out_path)

