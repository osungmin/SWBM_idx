#!/usr/bin/env python
"""
 extract idx list for locations where we have obs
 ewa and camel runoff data is considered
"""
import os
import pandas as pd
import numpy as np
from netCDF4 import Dataset

print ("imported modules")



class Cali_mit_obs(object):

    def __init__(self):
   
       self.camel_lis = self._load_camel()
       self.ewa_lis  = self._load_ewa()


    def _load_camel(self):
        area_path='/Net/Groups/BGI/scratch/suo/dataset/camel/basin_area.lis'
        attri_path='/Net/Groups/BGI/scratch/suo/dataset/camel/camels_attributes_v2.0/camels_topo.txt'

        #with open(basin_path) as fp:
        #    basins = fp.readlines()
        #basins = [basin.strip() for basin in basins]

        area_load = pd.read_csv(area_path, header=0, dtype={'basin':str}, index_col=None, sep=',', na_values=-999.)
        area_load.set_index('basin', inplace=True)
        area_load.area = area_load.area/1000000  #are in km2

        load = pd.read_csv(attri_path, header=0, dtype={'gauge_id':str}, index_col=None, sep=';')
        load.set_index('gauge_id', inplace=True)
        load = load[load.index.isin(area_load.index)].copy()

        load = pd.concat([load, area_load], axis=1)
        load = load.dropna() # bad data (too many zeros) has nan for area

        print("___ camel:", load.shape)
        
        return(load)


    def _load_ewa(self):

        path = '/Net/Groups/BGI/people/sungmino/runoff_data/EWA_metadata_v1.csv'
        load = pd.read_csv(path, header=0, index_col=0, sep=',', encoding='cp1252')
        load = load[load["Quality"]==1].copy()  #only best quality catchments
        return(load)


    def get_list(self):
        return self.camel_lis, self.ewa_lis


    def gridded(self, idx_lats, idx_lons):
        #to find idxs nearest to lon/lat points
        dlats=np.arange(89.875, -90, -0.25)
        dlons=np.arange(-179.875, 180, 0.25)        

        lats, lons = [], []
        for i in range(len(idx_lats)):
           ilat = np.abs(dlats - idx_lats[i]).argmin()
           ilon = np.abs(dlons - idx_lons[i]).argmin()
           lats.append(dlats[ilat])
           lons.append(dlons[ilon])
        return(lats, lons)


def _load_runoff(source, idx):
    if source=="camel": path='/Net/Groups/BGI/scratch/suo/dataset/camel/runoff/'+idx+".dat"
    if source=="ewa":  path='/Net/Groups/BGI/scratch/suo/models/inputs/runoff/runoff_'+idx+".dat"

    load = pd.read_csv(path, header=0, index_col=0, sep=',', na_values=-999)

    return load


def _load_clim():
   dirpath='/Net/Groups/BGI/scratch/suo/dataset/era5_0d25/static/'
   arid_path = dirpath+'aridity.era5.1440.720.1980-2014.nc'
   temp_path = dirpath+'t2m.average.era5.1440.720.1980-2014.nc'

   print("")
   print("reading aridity")
   f0  = Dataset(arid_path, 'r')#
   arid = f0.variables['snr'][:,:,:].filled(np.nan)
   print(f0.variables['latitude'][:3])
   print(arid.shape)

   print("")
   print("reading t2m_ave")
   f1  = Dataset(temp_path, 'r')#
   temp = f1.variables['t2m'][:,:,:].filled(np.nan)
   print(f1.variables['latitude'][:3])
   print(temp.shape)

   datalons = f1.variables['longitude'][:]
   datalats = f1.variables['latitude'][:]

   return(arid, temp, datalons, datalats)

"""
if __name__ == '__main__':
 
     cali = Cali_mit_obs()
     
     #station info
     camel_lis, ewa_lis = cali.get_list()
     print("*** loaded idx for camel and ewa", len(camel_lis), len(ewa_lis))
     print("click")
     wait=input()

     #output dataframe
     df=pd.DataFrame(columns=["lat_d25","lon_d25","st_id","arid","temp"])
     outpath='/Net/Groups/BGI/scratch/suo/ml_paras/obs_runoff/'

     #list of grid lat/lon 
     lats0,lons0 = cali.gridded(camel_lis.gauge_lat.values, camel_lis.gauge_lon.values)
     lats1,lons1 = cali.gridded(ewa_lis.Latitude.values, ewa_lis.Longitude.values) 

     #print(camel_lis[::50])

     #arid and temp ncdf
     arid_nc, temp_nc, datalons, datalats = _load_clim()

     #to check unique grid lat/lon: combine list + list
     df["lat_d25"]= lats0+lats1
     df["lon_d25"]= lons0+lons1
     df.index     = df["lat_d25"].astype(str)+"_"+df["lon_d25"].astype(str)

     df["st_id"]  = list(camel_lis.index)+list(ewa_lis.index)
     df["source"] = ["camel"]*len(camel_lis)+["ewa"]*len(ewa_lis)
     df["area"]   = list(camel_lis.area_geospa_fabric)+list(ewa_lis.Basin_area)
     
     #one of multiple stats in a grid pixel: flag will be assigned 
     df["cal_id"] = [np.nan]*len(df)
    
     print()
     print("*** output dataframe: ewa camel stations")
     print(df.shape)
     print("*** unique index of grid pixel (=where multiple stats)")
     print(len(np.unique(df.index)))
     print("")
     print(df[:10])
     print("")
    
     for i in np.unique(df.index): #by grid idx

       print("")
       idx=df[df.index==i]
       print(">>>", i)


       if len(idx)>1: 
          print("...averaged")
 
          check=0
          total_area= np.sum(idx["area"])

          if len((np.unique(idx["source"]))) != 1:
             print(" data source is mixed up?")
             wait=input()

          for j in range(len(idx)): #by stat id
              st_idx = str(idx.iloc[j]["st_id"])
              source = idx.iloc[j]["source"]
              area   = idx.iloc[j]["area"]

              #print ("**")
              #print (j, st_idx, source, area, round(area/total_area,3))

              # runoff x (area/total_area)
              runoff= _load_runoff(source, st_idx) * (area/total_area)

              if j==0: concat= runoff.copy()
              else:    concat= pd.concat([concat.copy(), runoff.copy()], axis=1)
         
              ###print(_load_runoff(source, st_idx).describe())

              df.loc[df.st_id==st_idx, "cal_id"]= i+"_avg"

              #extract clim
              lat_idx = (np.abs(datalats - df[df.st_id==st_idx]['lat_d25'][0])).argmin()
              lon_idx = (np.abs(datalons - df[df.st_id==st_idx]['lon_d25'][0])).argmin()

              df.loc[df.st_id==st_idx, "arid"]= np.round(arid_nc[:,lat_idx,lon_idx],3)
              df.loc[df.st_id==st_idx, "temp"]= np.round(temp_nc[:,lat_idx,lon_idx]-273.15,3)

              check+=area/total_area

          #sum of weight should be 1
          if round(check,3) != 1:
             print(" error in weight", check)
             wait=input()

          #i will take when all stats have values (to avoid single vs multiple runoff)
          concat=concat.dropna()
          concat['mean'] = concat.sum(axis=1) # use 'sum' here as runoff x area weight is used. 

          concat[["mean"]].to_csv(outpath+i+".dat" , index=True, header=['QObs'])
          concat=np.nan

       else: 

          st_idx = str(idx.iloc[0]["st_id"])
          source = idx.iloc[0]["source"]

          concat= _load_runoff(source, st_idx)

          df.loc[df.index==i, "cal_id"]=i  

          #extract clim
          lat_idx = (np.abs(datalats - df[df.index==i]['lat_d25'][0])).argmin()
          lon_idx = (np.abs(datalons - df[df.index==i]['lon_d25'][0])).argmin()

          df.loc[df.index==i, "arid"]= np.round(arid_nc[:,lat_idx,lon_idx],3)
          df.loc[df.index==i, "temp"]= np.round(temp_nc[:,lat_idx,lon_idx]-273.15,3)
          
          concat.to_csv(outpath+i+".dat" , index=True, header=['QObs'])
          concat=np.nan

 
     print("done")
     print(df.describe())
   
     df.to_csv('/Net/Groups/BGI/scratch/suo/ml_paras/meta/obs_d25.lis')
     print("")
     print(">>> unique grid number:", len(np.unique(df['cal_id'])))
     print("end.")

"""


path='/Net/Groups/BGI/scratch/suo/ml_paras/meta/runoffdata_summary.df'
load = pd.read_csv(path, header=0, index_col=0, sep=',', na_values='na')
load = load.dropna(subset=['grid'])
#load = load[(load.ndays>=365*10)&(load.carea>10)&(load.carea<1000)&(load.forest_gain<0.2)&(load.forest_loss<0.2)&(load.irrig<2)&(load.urban<2)&(load.nonreference!=1)].copy()
#print(load)

obs_idxs=pd.read_csv('/Net/Groups/BGI/scratch/suo/ml_paras/meta/obs_d25.lis', header=0, index_col=0)


print(len(np.unique(load.grid)))
print(len(np.unique(obs_idxs.index)))

idxs=list(load.grid)+list(obs_idxs.index)
print(len(np.unique(idxs)))

print( [x for x in obs_idxs.index if x not in load.grid])


