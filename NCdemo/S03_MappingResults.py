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
import sys
import os
# Aappend the Current path and load custom plotting modules 

# ========== Import packages ==========
import numpy as np
import pandas as pd
import argparse
# import datetime as dt
import warnings as warn
import xarray as xr
import dask
import bottleneck as bn
from collections import OrderedDict, defaultdict
import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.stats.multitest as smsM

# ==============================================================================
def main(args):
		# =========== Read the metadata file in ==========
	infofile = './data/infomation.json'
	with open(infofile, "r+") as f:
		info = json.load(f)
		# ========== fix the timedelta ==========
		if type(info["ComputeTime"]) == float:
			info["ComputeTime"] = pd.Timedelta(info["ComputeTime"], unit="sec")


	# ========== Process system arguments ==========
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
	ds = GlobalAttributes(ds, info, ds_ref=ds_ref)

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

	# ========== Make a series of maps ==========
	if args.plots:
		for var in ds.variables:
			# check its a variable to be skipped or mapped
			if (var in ['longitude', 'latitude', "OtherFactorsValid", "errors", "time"] or var.endswith(".Pvalue")):
				continue

			# ========== Perform significance correction ==========
			if args.sig:
				signif, s_method = FDRSignificanceCorrection(ds, var, args.method)
			else:
				s_method = ""
				signif   = None
			
			print(var)
			breakpoint()





	breakpoint()
# ==============================================================================
# BUILD THE MAPS
# ==============================================================================
def MapMaker(ensinfo, results, va, va_col, maininfo, gitinfo, pshow, mask, photo, fdr_col=None, mask_col=None):
	"""
	Builds the maps of the change attribution
	args:
		ensinfo
		results 	npA : with the atttribution results
		va  		str : the name of the change variable
		va_col 		int : the column number
		maininfo	str : plot string
		gitinfo 	str : summary of the git header
		pshow		bool: show the plots or just save
		fdr_col		int : The colum of the zeros and ones that are used for fdr significance checking
		mask_col	int : The colum of the Nan's and ones that are the dryland mask
	""" 
	# ========== Create the mapdet object ==========
	mapdet = pf.mapclass(ensinfo, pshow)
	mapdet.dpi = 130
	
	# ========== Remove the non dryland areas ==========
	if not (mask_col is None):
		# ========== account for dryland mask ========== 
		results[1:-1, va_col] *= results[1:-1, mask_col]
	
	# ========== Set the colormap ==========
	tks = np.array([-0.30, -0.24, -0.12, -0.06, -0.02, 0, 0.02, 0.06, 0.12, 0.24, 0.30])
	
	if va == "Climate":
		cmapHex     = palettable.colorbrewer.diverging.BrBG_10.hex_colors
		mapdet.cmin = -0.3#-0.1 	# the min of the colormap
		mapdet.cmax =  0.3#0.1	# the max of the colormap
		# colour      = "#00bfff"
		colour      = "#005f7f"
	elif "CO2_" in va:
				# cmapHex = palettable.colorbrewer.diverging.PuOr_10.hex_colors
		cmapHex = palettable.colorbrewer.sequential.Purples_7.hex_colors
		mapdet.cmin =  0.0#-0.1 	# the min of the colormap
		mapdet.cmax =  0.3#0.1	# the max of the colormap
		tks         = np.array([0, 0.001, 0.02, 0.03, 0.04, 0.06, 0.24, 0.30])
		colour      = "#006400"
	elif va == "LandUse":
		cmapHex     = palettable.colorbrewer.diverging.PiYG_10.hex_colors
		mapdet.cmin = -0.3#-0.1 	# the min of the colormap
		mapdet.cmax =  0.3#0.1	# the max of the colormap
		# colour      = "#ffc125"
		colour       = "#eaa700"
	elif va == "OtherFactors":
		cmapHex = palettable.colorbrewer.diverging.RdBu_10.hex_colors
		mapdet.cmin = -0.3#-0.1 	# the min of the colormap
		mapdet.cmax =  0.3#0.1	# the max of the colormap
		colour = cmapHex[-2]
		
	
	# ========== Set the variables ==========
	mapdet.var   =  va + ensinfo.desc 
	mapdet.column=  va_col	# the column to be mapped
	mapdet.cZero =  0.5 	# the zero point of the colormap
	if not(mask is None):
		mapdet.mask	 = np.load(mask)
	
	# ========== Add the Font info ========== 
	mapdet.gridalp  = 0.05


	if not (fdr_col is None):
		# ========== account for significance in the results ========== 
		if isinstance(fdr_col, (list,)):
			for fdr_c in fdr_col:
				results[1:-1, va_col] *= results[1:-1, fdr_c]	
		else:
			results[1:-1, va_col] *= results[1:-1, fdr_col]
	
	# ========== Setup the cmap ==========
	if not ("CO2_Change" in va):
		cmap, norm, ticks, cbounds, spacing = pf.ReplaceHexColor(
			cmapHex, mapdet.cmin, mapdet.cmax, ticks=tks, zeroR = 0.001)
	else:
		spacing = 'uniform'
		cbounds = tks
		ticks   = tks
		# ========== create the colormap ==========
		cmap   = mpl.colors.ListedColormap(cmapHex)

		cmap.set_over(cmapHex[-1])
		
		norm   = mpl.colors.BoundaryNorm(cbounds, cmap.N)
		mapdet.extend  = "max"

	# cmap, norm, ticks, cbounds, spacing = pf.ReplaceHexColor(
	# 	cmapHex, mapdet.cmin, mapdet.cmax, ticks=tks, zeroR = 0.001)

	cmap.set_bad('dimgrey',1.)
	
	mapdet.cmap    = cmap 
	mapdet.norm    = norm
	mapdet.ticks   = ticks	
	mapdet.cbounds = cbounds
	mapdet.spacing = spacing

	# ========== Make a map ========== 
	plotinfo, fname = pf.mapmaker(results, ensinfo, mapdet) 

	# ========== Save the metadata ========== 
	if fname:
		infomation = [maininfo, histinfo, hist_fname, plotinfo, fname, gitinfo]
		cf.writemetadata(fname, infomation)
	# ipdb.set_trace()
	return colour

def FDRSignificanceCorrection(ds, var, FDRmethod, alpha = 0.10, ):
	"""
	Takes the results of an existing trend detection aproach and modifies them to
	account for multiple comparisons.  
	args
		ds: xr dataset
			list of numpy arrays containing results of trend analysis
		kys: list 
			list of what is in results
		years:
			years of accumulation 
	
	"""
	if FDRmethod == "fdr_by":
		t_method = "Adjusting for multiple comparisons using Benjamini/Yekutieli"
	elif FDRmethod == "fdr_bh":
		t_method = "Adjusting for multiple comparisons using Benjamini/Hochberg"
	else:
		raise ValueError("unknown MultipleComparisons method: %s, must be either fdr_by or fdr_bh " % FDRmethod)

	# ========== Work out the pvalue name =========== 
	if var == "ObservedChange":
		pnm = "Obs.Pvalue"
	elif var == "OtherFactors":
		t_method = "OtherFactors is Zero masked in areas where 1 or more components could not be calculated"
		print(t_method)
		return ds["OtherFactorsValid"].astype(float), t_method
	else:
		pnm = var + ".Pvalues"


	# ========== Locate the p values and reshape them into a 1d array ==========
	# ++++++++++ Find the pvalues ++++++++++
	ODimDa     = ds[pnm].stack(loc = ('time', 'latitude', 'longitude')).copy()
	isnan      = np.isnan(ODimDa)
	
	# ++++++++++ pull out the non nan pvalus ++++++++++
	pvalue1d   = ODimDa[~isnan]
	
	# =========== Perform the MC correction ===========
	pvalue_adj = smsM.multipletests(pvalue1d, method=FDRmethod, alpha=alpha)

	# =========== Get the significance in bool ===========
	sigfloat   = pvalue_adj[0].astype(float)

	# make an empty dataarray
	re         = xr.zeros_like(ODimDa)
	re[~isnan] = sigfloat

	print(t_method)
	return re.unstack(), t_method
# ==============================================================================
# FIX THE METADATA OF THE XR DATASET
# ==============================================================================

def GlobalAttributes(ds, info, ds_ref=None):
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
	info["history"]             = "%s: Netcdf file created using %s (%s):%s. Script originally developed by %s" % (
		str(pd.Timestamp.now()), __title__, __file__, __version__, __author__) + info["history"]
	# attr["history"]            += ds.history
	attr["history"]             = info["history"]

	attr["creator_name"]        = __author__
	attr["creator_url"]         = "ardenburrell.com"
	attr["creator_email"]       = __email__
	attr["institution"]         = "WHRC"
	attr["date_created"]        = str(pd.Timestamp.now())
	
	# ++++++++++ Netcdf Summary infomation ++++++++++ 
	attr["time_coverage_start"] = str(pd.Timestamp("2015-12-31"))
	attr["time_coverage_end"]   = str(pd.Timestamp("1982-01-01"))
	# ++++++++++ TSS.RESTREND INFOMATION ++++++++++ 
	attr["package_version"]     = "TSSRESTREND version: " + info["TSSRESTREND.version"]
	attr["package_url"]         = "https://cran.r-project.org/web/packages/TSS.RESTREND/index.html"
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
	description='Arguments that can be passed to the saving and plotting script'
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument(
		"-p", "--plots", action="store_false", 
		help="Quiet plots: if selected plots will not be displayed ")
	parser.add_argument(
		"-s", "--sig", action="store_false", 
		help="Significance: Apply a zero mask using FDR adjustment and the Benjamini/Hochberg method")
	parser.add_argument(
		"--method", type=str, default="fdr_bh", help="The method used to adjust for False Discovery Rate. must be fdr_bh or fdr_by")
	args = parser.parse_args() 
	main(args)