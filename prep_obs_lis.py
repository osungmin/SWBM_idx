#!/usr/bin/env python
"""
 extract a list of lat/lon where we have camel runoff observation
"""
import os
import pandas as pd
import numpy as np
print ("imported modules")


class Cali_mit_obs(object):

    def __init__(self):
   
       self.camel_lis = self._load_camel()
  
  
    def _load_camel(self):
        area_path='your_path_to_camel/camel/basin_area.lis'
        attri_path='your_path_to_camel/camels_attributes_v2.0/camels_topo.txt'

        area_load = pd.read_csv(area_path, header=0, dtype={'basin':str}, index_col=None, sep=',', na_values=-999.)
        area_load.set_index('basin', inplace=True)
        area_load.area = area_load.area/1000000  #area in km2

        load = pd.read_csv(attri_path, header=0, dtype={'gauge_id':str}, index_col=None, sep=';')
        load.set_index('gauge_id', inplace=True)
        load = load[load.index.isin(area_load.index)].copy()

        load = pd.concat([load, area_load], axis=1)
        load = load.dropna() # bad data (too many zeros) has nan for area

        print("___ camel:", load.shape)
        
        return(load)

       
    def get_list(self):
        return self.camel_lis


    def gridded(self, idx_lats, idx_lons):
        #to find idxs nearest to lon/lat points (example of 0.25 deg)
        dlats=np.arange(89.875, -90, -0.25) 
        dlons=np.arange(-179.875, 180, 0.25)        

        lats, lons = [], []
        for i in range(len(idx_lats)):
           ilat = np.abs(dlats - idx_lats[i]).argmin()
           ilon = np.abs(dlons - idx_lons[i]).argmin()
           lats.append(dlats[ilat])
           lons.append(dlons[ilon])
        return(lats, lons)


def _load_runoff(idx):
 
    path='your_path_to_camel/camel/runoff/'+idx+".dat"
    load = pd.read_csv(path, header=0, index_col=0, sep=',', na_values=-999)

    return load


if __name__ == '__main__':
 
     cali = Cali_mit_obs()
     
     #station info
     camel_lis = cali.get_list()
     print("*** loaded idx for camel", len(camel_lis))
     print("click")
     wait=input()

     #output dataframe
     df=pd.DataFrame(columns=["lat_d25","lon_d25","st_id","area","cal_id"])
     outpath='your_path_to_save_output_list'

     #list of grid lat/lon 
     lats0, lons0 = cali.gridded(camel_lis.gauge_lat.values, camel_lis.gauge_lon.values)
     #print(camel_lis[::50])

     #to save idx info to the dataframe 
     df["lat_d25"]= lats0
     df["lon_d25"]= lons0
     df.index     = df["lat_d25"].astype(str)+"_"+df["lon_d25"].astype(str)

     df["st_id"]  = list(camel_lis.index)
     df["area"]   = list(camel_lis.area_geospa_fabric)
     
     #to mark if there are multiple stations in a grid pixel: flag will be assigned later
     df["cal_id"] = [np.nan]*len(df)
    
     print()
     print("*** output dataframe: camel stations")
     print(df.shape)
     print("*** unique index of grid pixel (=where multiple stats)")
     print(len(np.unique(df.index)))
     print()
     print(df[:10])
     print()
    
     for i in np.unique(df.index): #by grid idx

       print()
       idx=df[df.index==i]
       print(">>>", i)

       if len(idx)>1: 
          #there are more than one stations at a grid pixel 
          #I will take averaege of runoff observation
          
          print("... averaging")
 
          check=0
          total_area= np.sum(idx["area"])

          for j in range(len(idx)): #by stat id
              st_idx = str(idx.iloc[j]["st_id"])
              area   = idx.iloc[j]["area"]

              #print ("**")
              #print (j, st_idx, area, round(area/total_area,3))

              # runoff x (area/total_area); compute average runoff with wighting based on basin area.
              runoff= _load_runoff(st_idx) * (area/total_area)

              if j==0: concat= runoff.copy()
              else:    concat= pd.concat([concat.copy(), runoff.copy()], axis=1)
         
              ###print(_load_runoff(st_idx).describe())

              df.loc[df.st_id==st_idx, "cal_id"]= i+"_avg"

              #just to check if sum of weight is 1
              check+=area/total_area

          #sum of weight should be 1
          if round(check,3) != 1:
             print(" ?? error in weight", check)
             wait=input()

          #i will take time series when all stats have values (to avoid single vs multiple runoff)
          concat=concat.dropna()
          concat['mean'] = concat.sum(axis=1) # use 'sum' here as runoff x area weight is computed. 

          concat[["mean"]].to_csv(outpath+i+".dat" , index=True, header=['QObs'])
          concat=np.nan

       else: 

          st_idx = str(idx.iloc[0]["st_id"])
         
          concat= _load_runoff(st_idx)

          df.loc[df.index==i, "cal_id"]=i  

          concat.to_csv(outpath+i+".dat" , index=True, header=['QObs'])
          concat=np.nan

 
     print("done")
     print(df.describe())
   
     df.to_csv('your_path_to_save_output_list/obs_d25.lis')
     print("")
     print(">>> unique grid number:", len(np.unique(df['cal_id'])))
     print("end.")

