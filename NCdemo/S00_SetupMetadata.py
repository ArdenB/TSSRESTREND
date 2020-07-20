"""
Script goal, 

Script contains the preprocessing steps and run creattion


"""

# ==============================================================================

__title__   = "Dataset preprocessing and run creattion"
__author__  = "Arden Burrell"
__version__ = "v1.0(23.06.2020)"
__email__   = "aburrell@whrc.org"

# ==============================================================================
import xarray as xr 
import dask 
import bottleneck as bn
import pandas as pd
import os
import argparse
from collections import OrderedDict
import json

def main(args):
	fnV    = "./data/vod_AUS_GIMMS_monthly_1982_2015.nc"
	fnoutV = "./data/AUSdemo_GIMMS_ndvi.nc"
	fnP    = "./data/TerraClimate_ppt_AUS_remapped.nc"
	fnoutP = "./data/AUSdemo_TERRACLIMATE_ppt.nc"
	fnT    = "./data/TerraClimate_tmean_AUS_remapped.nc"
	fnoutT = "./data/AUSdemo_TERRACLIMATE_tmean.nc"

	# ========== This will determine by hoe much to carsen the data ======= 
	# pass 0 to get original resolution, higher numbers make coarser data which if faster to process
	coarsen = args.coarsen

	# ========== Add a loop to coarsen the data ==========
	if coarsen > 0:
		for va, fn in zip(["ndvi", "ppt", "tmean"],[fnoutV, fnoutP, fnoutT]):
			xr.set_options(keep_attrs=True)
			dsin = xr.open_dataset(fn) 
			# ========== Coarsen the dataset ==========
			dsout = dsin.coarsen(dim={"longitude":coarsen, "latitude":coarsen}, boundary="pad").mean()
			hist2 = "%s: Coarsend using Xarray coarsen with a window of size %d" % (pd.Timestamp.now(), coarsen)
			dsout.attrs = dsin.attrs
			dsout.attrs["history"]  = hist2 + dsout.attrs["history"]
			fnout = fn[:-3] +"_xrcoarsen_%dwin.nc" % coarsen
				# write out
			dsout.to_netcdf(fnout, 
				format         = 'NETCDF4', 
				encoding       = {va:encoding},
				unlimited_dims = ["time"])
			warn.warn("Need to implement an approximate timer here")


	# ========== This is where i plan to implement some tracking and auto formatiing ==========
	Metadata = OrderedDict()

	Metadata["history"] = "%s: Datasets built and run setup with GIMMS NDVI and TERRACLIMATE data. Downscaled to approximatly %dkm." % (str(pd.Timestamp.now()), (25*(coarsen-1)+25) )
	Metadata["coarsen"] = coarsen
	
	#  TO DO HERE
	# 	1. Add the dataset metadata
	# 	2. Add a funtion to set change units, this can be passed sys argas
	# 	3. Write to a Json file so it can be passed from file to file
	
	# Metadata["NDVI"]  = {GIMMS, }
	# Metadata["ppt"]   = {TERRACLIAMTE}
	# Metadata["tmean"] = {}

	warn.warn("Implement the creation of a dataset and run metadata file")
	breakpoint()

def _internalsaves(fnV, fnoutV, fnP, fnoutP, fnT, fnoutT):
	"""This function is not available to the user as it depends on data too large to put in a 
	git repo, It lef here so i know how to redo the data if needs be."""

	encoding = ({'shuffle':True, 
		'zlib':True,
		'complevel':5})
	hist = "%s: Regridded and modified version of GIMMS3gv1 NDVI data. It was produced by AB to demonstrate the TSS-RESTREND method" % (pd.Timestamp.now())
	if not os.path.isfile(fnoutV):
		# ========== NDVI ==========
		# drop percentile and rename lon and lat
		dsV = xr.open_dataset(fnV).rename({"lon":"longitude", "lat":"latitude"}).drop("percentile")
		# resscale and remove nan
		dsV["ndvi"] /= 10000
		dsV = dsV.where(dsV["ndvi"]>= 0)

		# Add to the history
		dsV.attrs["history"]  = hist + dsV.attrs["history"]
		dsV.attrs["FileName"] = fnoutV

		# write out
		dsV.to_netcdf(fnoutV, 
			format         = 'NETCDF4', 
			encoding       = {'ndvi':encoding},
			unlimited_dims = ["time"])

	if not os.path.isfile(fnoutP):
		# ========== ppt ==========
		dsP = xr.open_dataset(fnP).rename({"lon":"longitude", "lat":"latitude"})
		# Add to the history
		hist = "%s: Regridded and modified version of TERRACLIMATE ppt data. It was produced by AB to demonstrate the TSS-RESTREND method" % (pd.Timestamp.now())
		dsP.attrs["history"]  = hist + dsP.attrs["history"]
		dsP.attrs["FileName"] = fnoutP
		dsP["ppt"] = dsP["ppt"].astype("float32")

		# write out
		dsP.to_netcdf(fnoutP, 
			format         = 'NETCDF4', 
			encoding       = {'ppt':encoding},
			unlimited_dims = ["time"])

	if not os.path.isfile(fnoutT):
		# ========== tmean ==========
		dsT = xr.open_dataset(fnT).rename({"lon":"longitude", "lat":"latitude"})
		# Add to the history
		hist = "%s: Regridded and modified version of TERRACLIMATE temperature data. It was produced by AB to demonstrate the TSS-RESTREND method" % (pd.Timestamp.now())
		dsT.attrs["history"]  = hist + dsT.attrs["history"]
		dsT.attrs["FileName"] = fnoutT
		dsT["tmean"] = dsT["tmean"].astype("float32")

		# write out
		dsT.to_netcdf(fnoutT, 
			format         = 'NETCDF4', 
			encoding       = {'tmean':encoding},
			unlimited_dims = ["time"])
#==============================================================================
if __name__ == '__main__':
	# ========== Set the args Description ==========
	description='Build the standard plots for each multi-run ensemble'
	parser = argparse.ArgumentParser(description=description)
	
	# ========== Add additional arguments ==========
	parser.add_argument(
		"--coarsen", type=int, default=0, help="The size of the box used to downscale data, Defualt is zeros")
	args = parser.parse_args() 
	
	# ========== Call the main function ==========
	main(args)
