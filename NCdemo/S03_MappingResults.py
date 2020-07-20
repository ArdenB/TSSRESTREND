"""
Script goal, 

Open the results from TSS-RESTREND attribution analysis and convering that 
back to a spatial data format.  


"""

# ==============================================================================

__title__   = "csv to netcdf and maps"
__author__  = "Arden Burrell"
__version__ = "v1.0(23.06.2020)"
__email__   = "aburrell@whrc.org"

# ==============================================================================

# ========== Import packages ==========
import os
import sys
import numpy as np
import pandas as pd
# import datetime as dt
import warnings as warn
import xarray as xr
import dask
import bottleneck as bn
from collections import OrderedDict, defaultdict

import matplotlib as mpl
import matplotlib.pyplot as plt
# import seaborn as sns

# ==============================================================================
def main():
	# ========== Open the csv results file ==========
	fn = "./results/AttributionResults.csv"
	df = pd.read_csv(fn, index_col=0)
	
	# ========== Fix the indexing and the error messaging ==========
	lonlat = np.array([_unpackLatLon(index) for index in df.index])
	df["longitude"] = lonlat[:, 0]
	df["latitude" ] = lonlat[:, 1]
	df = df.reset_index(drop=True).set_index(["longitude","latitude" ])
	df["errors"] = df["errors"].apply(_df_error)

	# ========== Open a reference dataset ==========
	# This is a hack so i don't have to harrd code all the infomation 
	fnNDVI  = "./data/AUSdemo_GIMMS_ndvi.nc"
	ds_ref  =  xr.open_dataset(fnNDVI)
	# ========== Convert back to xr dataset ==========
	ds = xr.Dataset.from_dataframe(df)	
	
	# ++++++++++ Fix the time ++++++++++
	# add a time dim, this case the end of the original timeseris. 
	# Netcdf files need a time dimenstion
	ds = ds.assign_coords(time=pd.Timestamp("2015-12-31")).expand_dims("time")
	ds = ds.transpose('time', 'latitude', 'longitude')
	
	# ++++++++++ Fix the Global Attributes ++++++++++
	ds = GlobalAttributes(ds, ds_ref=ds_ref)

	# ++++++++++ Setup the netcdf file encoding and add  ++++++++++
	encoding = OrderedDict()
	for va in df.columns:
		# ========== Check and see if all th e values are nan, if yes drop them ==========
		if bn.allnan(ds[va]):
			ds = ds.drop_vars(va)
		else:
			encoding[va] = ({'shuffle':True, 
				'zlib':True,
				'complevel':5.
				# 'chunksizes':[1, ensinfo.lats.shape[0], 100],
				})
	warn.warn("TO BE IMPLEMENTED: The Attributes of the data variables need to be added here")
	# To ADD, Fix variable attributes

	# ========== Write the dataset to a netcdf file ==========
	print("Starting write of data at:", pd.Timestamp.now())
	ds.to_netcdf('./results/TSSRattribution_Results.nc', 
		format         = 'NETCDF4', 
		encoding       = encoding,
		unlimited_dims = ["time"])




	breakpoint()


# ==============================================================================
# FIX THE METADATA OF THE XR DATASET
# ==============================================================================

def GlobalAttributes(ds, ds_ref=None):
	"""
	Creates the global attributes for the netcdf file that is being written
	these attributes come from :
	https://www.unidata.ucar.edu/software/thredds/current/netcdf-java/metadata/DataDiscoveryAttConvention.html
	args
		ds: xarray ds
			Dataset containing the infomation im intepereting
		ds_ref: xarray ds
			Dataset that contains the original data to copy metadata from. Defualt is None so it can be skipped
	returns:
		attributes 	Ordered Dictionary cantaining the attribute infomation
	"""
	# ========== Create the ordered dictionary ==========
	attr = ds.attrs

	# fetch the references for my publications
	# pubs = puplications()
	
	# ========== Fill the Dictionary ==========

	# ++++++++++ Highly recomended ++++++++++ 
	attr["title"]               = "TSS-RESTREND Attribution Results" 
	attr["summary"]             = "The restacked results of a TSSRattribution function applied to GIMMS and TERRCLIMATE DATA"
	attr["Conventions"]         = "CF-1.7"
	
	# ++++++++++ Data Provinance ++++++++++ 
	attr["history"]             = "%s: Netcdf file created using %s (%s):%s by %s" % (
		str(pd.Timestamp.now()), __title__, __file__, __version__, __author__)
	# attr["history"]            += ds.history

	attr["creator_name"]        = __author__
	attr["creator_url"]         = "ardenburrell.com"
	attr["creator_email"]       = __email__
	attr["institution"]         = "WHRC"
	attr["date_created"]        = str(pd.Timestamp.now())
	
	# ++++++++++ Netcdf Summary infomation ++++++++++ 
	attr["time_coverage_start"] = str(pd.Timestamp("2015-12-31"))
	attr["time_coverage_end"]   = str(pd.Timestamp("1982-01-01"))
	if not ds_ref is None:
		ds.longitude.attrs = ds_ref.longitude.attrs
		ds.latitude.attrs  = ds_ref.latitude.attrs
		ds.time.attrs      = ds_ref.time.attrs
		# breakpoint()
		# ds.time.attrs["calendar"] = 'standard'
		# ds.time.attrs["units"]    = 'days since 1900-01-01 00:00'
	return ds

# ==============================================================================
# FUNCTIONS FOR FIXING THE INPUT DATAFRAME
# ==============================================================================

def _unpackLatLon(index):
	"""
	Fixes the string formatting of the indexs 
	args:
		index:	str
			string of format "(longitude, latitude)""
	return:
		longitude, latitude: Floats
			"""
	# fix the latitude and longitude 
	lonST, latST = index.split(", ")
	return float(lonST[1:]), float(latST[:-1])

def _df_error(error):
	""" 
	Takes in the error codes which are catogorical and encodes them into integers
	"""
	# ========== Keys for error messages ==========
	keys = ({
		"NANinCTSR.VI":1,
		"NANinCTSR.RF":2,
		"NANinCTSR.TM":3,
		"AttributionFailed":4
		})

	# breakpoint()
	if error in keys.keys():
		return keys[error]
	elif np.isnan(error):
		return 0 
	else:
		# unknow errror
		return 5





# ==============================================================================
if __name__ == '__main__':
	main()