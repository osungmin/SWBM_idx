#!/usr/bin/env python
import pandas as pd
import numpy as np
import datetime as dt
print("modules imported")


def _load_attributes(opt=None):
	f0='your_path_to_caravan/caravan/attributes/'+opt+'/attributes_caravan_'+opt+'.csv'
	load0=pd.read_csv(f0, header=0, index_col=0)
	f1='your_path_to_caravan/caravan/attributes/'+opt+'/attributes_hydroatlas_'+opt+'.csv'
	load1=pd.read_csv(f1, header=0, index_col=0)
	load=pd.concat([load0, load1], axis=1).dropna()
	print(" to check", load0.shape, load1.shape, load.shape)
	print(" __ loading attributes", opt)
	return(load)


def _load_runoff(opt=None, gauge_id=None):
	f='your_path_to_caravan/caravan/timeseries/csv/'+opt+'/'+gauge_id+'.csv'
	load=pd.read_csv(f, header=0, index_col=0, parse_dates=True, sep=',')
	load=load['streamflow'].copy()
		
	print(" __ loading time series", opt, gauge_id)
	#time period can be adjusted
	load=load[(load.index>=dt.datetime(1990,1,1))&(load.index<=dt.datetime(2019,12,31))].copy()
	return(load)


def compute_mean_runoff(df):
	#i will take mean of runoff (area-weighted) if there are multiple runoff observations at a grid pixel
	check= 0
	len_check= []

	total_area= np.sum(df['area'])

	for i in range(len(df)):

		# runoff x (area/total_area)
		area= df['area'][i]
		
		runoff=_load_runoff(df['source'][i], df.index[i]) * (area/total_area)
		if i==0: concat= runoff.copy()
		else:    concat= pd.concat([concat.copy(), runoff.copy()], axis=1)

		check+=area/total_area
		len_check.append(len(runoff.dropna()))

	#sum of weight should be 1
	if round(check,3) != 1:
		print(" error in weight", check)
		wait=input()
		
	#i will take when all idxs have values (to avoid single vs multiple runoff)
	concat=concat.dropna()
	concat['mean'] = concat.sum(axis=1) # use 'sum' here as runoff x area weight is used. 

	return (concat, len_check)


def gridded(idx_lats, idx_lons):
	#to find idxs nearest to lon/lat points
	#resolution can be adjusted
	dlats=np.arange(89.875, -90, -0.25)
	dlons=np.arange(-179.875, 180, 0.25) 
	
	grids = []
	for i in range(len(idx_lats)):
		ilat = np.abs(dlats - idx_lats[i]).argmin()
		ilon = np.abs(dlons - idx_lons[i]).argmin()
		
		if idx_lats[i] != idx_lats[i] : grids.append(np.nan)
		else:  grids.append(str(dlats[ilat])+"_"+str(dlons[ilon]))
		
	return(grids)


def _prep_meta(opts):	
	#  
	keys=['idx','gauge_lat','gauge_lon','area','source']
	outs={key: list() for key in keys}
	
	#
	for opt in opts:
		# variable names for each data
		xlat='Latitude' if opt=='ewa' else 'gauge_lat'
		xlon='Longitude' if opt=='ewa' else 'gauge_lon'
		xarea='Basin_area' if opt=='ewa' else 'basin_area'

        #
		print(" >>", opt)
		meta= _load_attributes(opt=opt)
	 	
		outs['idx'].extend([x for x in meta.index.values])
		outs['gauge_lat'].extend([x for x in meta[xlat].values])
		outs['gauge_lon'].extend([x for x in meta[xlon].values])
		outs['area'].extend([x for x in meta[xarea].values])
		outs['source'].extend([opt]*len(meta))   

	#
	print("done")
	meta = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in outs.items() ]))
	meta = meta[keys]
	
	#
	meta['grid']=gridded(meta['gauge_lat'], meta['gauge_lon'])
	meta=meta[(meta['area']>=100)&(meta['area']<=3000)].copy()     

	print(" > saving")
	meta.to_csv("your_path_to_save_list/meta_runoff.df", index=False, na_rep=-9999.)
	return(meta)


def _save_outs(outs, info):
	outs['idx'].append(info[0])
	outs['source'].append(info[1])
	outs['nday'].append(info[2])
	return(outs)


def collect_runoff():
	#
	meta=pd.read_csv("your_path_to_save_list/meta_runoff.df", header=0, index_col=0)
	outpath="your_path_to_save_runoff/caravan_runoff/"

	##to create meta data for collected runoff
	keys=['idx','source','avg','nday']
	outs={key: list() for key in keys}

	ii=0
	for idx in np.unique(meta.grid):
		#
		_grids = meta[meta.grid==idx].copy()
		#
		print(" ***", ii)

		#
		_runoff_selected=False

		if len(_grids)==1: # only one observatino per a grid pixel => i will take it. 

			runoff=_load_runoff(_grids['source'][0], _grids.index[0])
			if len(runoff.dropna())>=10*365:  #i will take data with >=10 yr record
				if np.nanmin(runoff)<0:
					runoff[runoff<0]=np.nan
					print(" error in runoff 1 ",idx)
					
				runoff.to_csv(outpath+"grid_"+idx+".dat", index=True, header=['QObs'], na_rep='-9999')
				_runoff_selected=True	 
				
		else: # there are multiple observations per a grid pixel  
			
			# i will take an average of the multiple observations 
			runoff, len_check = compute_mean_runoff(_grids)
			
			if len(runoff.dropna())>=10*365:
				print("mean", idx)
				if np.nanmin(runoff["mean"])<0:
					runoff.loc[runoff['mean'] <0, 'mean'] = np.nan
					print(" error in runoff 2 ",idx)

				runoff[["mean"]].to_csv(outpath+"grid_"+idx+".dat", index=True, header=['QObs'], na_rep='-9999')
				_runoff_selected=True

			else: #if averaged runoff is short, i will take a single longest data
				longest_idx= np.argmax(len_check)
				runoff=_load_runoff(_grids['source'][longest_idx], _grids.index[longest_idx])
				if len(runoff.dropna())>=10*365: #2. longest global data selected
					if np.nanmin(runoff)<0:
						runoff[runoff<0]=np.nan
						print(" error in runoff 3 ",idx)

					runoff.to_csv(outpath+"grid_"+idx+".dat", index=True, header=['QObs'], na_rep='-9999')
					_runoff_selected=True
				else: print("___ too short even the longest one", _grids['source'], len(runoff.dropna()))

		#
		if _runoff_selected: _save_outs(outs, [idx, _grids['source'][0], len(runoff.dropna())])
		ii+=1
	
	#
	print("done")
	meta=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in outs.items()]))
	meta.to_csv("your_path_to_save_list/meta_runoff_selected.df", index=False, na_rep=-9999.)


#####
locs=["camels","camelsaus","camelsbr","camelscl","camelsgb","hysets","lamah"]

#####
print("***** prepare meta *****")
meta = _prep_meta(locs)
print("meta.describe")
print(meta.describe())
print("len of unique grids")
print(len(np.unique(meta['grid'])))
print()
#####


######
print("***** collect runoff *****")
collect_runoff()
#####




