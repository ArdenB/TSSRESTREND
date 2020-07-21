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
import argparse
from collections import OrderedDict
import warnings as warn
import json
# ========== Load my custom functions ==========
import os
import sys
sys.path.append(os.getcwd())
import CustomFunctions as cf 


def main(args):
	# ========== Make Folders ========== 
	folders = ([
		"./data/",
		"./results/",
		"./results/plots/",
		])
	for fol in folders:
		cf.pymkdir(fol)

	# ========== Setup the data paths ========== 
	fnV    = "./data/vod_AUS_GIMMS_monthly_1982_2015.nc"
	fnoutV = "./data/AUSdemo_GIMMS_ndvi.nc"
	fnP    = "./data/TerraClimate_ppt_AUS_remapped.nc"
	fnoutP = "./data/AUSdemo_TERRACLIMATE_ppt.nc"
	fnT    = "./data/TerraClimate_tmean_AUS_remapped.nc"
	fnoutT = "./data/AUSdemo_TERRACLIMATE_tmean.nc"

	encoding = ({
		'shuffle':True, 
		'zlib':True,
		'complevel':5})

	# ========== Create the  ========== 
	dssize = _internalsaves(fnV, fnoutV, fnP, fnoutP, fnT, fnoutT, encoding)

	# ========== This will determine by hoe much to carsen the data ======= 
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
			dssize = dsout.latitude.size * dsout.longitude.size
	
	print("\n A run at this resolution will take aproximatly:", (dssize * pd.Timedelta(5, unit="sec")), 
		"to complete on a single core. \n For shorter run times use a larger coarsen value to shring the resoluntion or use parellel processing. \n")

	# ========== This is where i plan to implement some tracking and auto formatiing ==========
	Metadata = OrderedDict()
	Metadata["Title"]   = "TSS.RESTREND change estimate"
	Metadata["history"] = "%s: Datasets built and run setup with GIMMS NDVI and TERRACLIMATE data. Downscaled to approximatly %dkm." % (str(pd.Timestamp.now()), (25*(coarsen-1)+25) )
	Metadata["coarsen"] = coarsen
	Metadata["NDVI"]   = ({
		"name"     :"GIMMS", 
		"var"      :"ndvi", 
		"reference":"Pinzon, Jorge E., and Compton J. Tucker. A Non-Stationary 1981-2012 AVHRR NDVI3g Time Series. Remote Sensing 6, no. 8 (25 July 2014): 6929--60. https://doi.org/10.3390/rs6086929."
		})
	Metadata["precipitation"] = ({
		"name"     :"TERRACLIMATE", 
		"var"      :"ppt", 
		"reference":"Abatzoglou, John T., Solomon Z. Dobrowski, Sean A. Parks, and Katherine C. Hegewisch. TerraClimate, a High-Resolution Global Dataset of Monthly Climate and Climatic Water Balance from 1958-2015. Scientific Data 5 (9 January 2018): 170191. https://doi.org/10.1038/sdata.2017.191."})
	Metadata["temperature"]   = ({
		"name"     :"TERRACLIMATE", 
		"var"      :"tmean", 
		"reference":"Abatzoglou, John T., Solomon Z. Dobrowski, Sean A. Parks, and Katherine C. Hegewisch. TerraClimate, a High-Resolution Global Dataset of Monthly Climate and Climatic Water Balance from 1958-2015. Scientific Data 5 (9 January 2018): 170191. https://doi.org/10.1038/sdata.2017.191."})
	# ========== Set the max ops and acp ==========
	Metadata["maxacp"] = args.maxacp
	Metadata["maxosp"] = args.maxosp
	Metadata["annual"] = args.annual
	Metadata["pixelN"] = dssize
	if args.annual:
		Metadata["units" ] =  r"$NDVI_{max}$ per year"
	else:
		Metadata["units" ] =   r"$\Delta NDVI_{max}$"

	Metadata["DataProcessed"]       = ""
	Metadata["TSSRESTREND.version"] = ""
	Metadata["ComputeTime"]         = ""
	Metadata["ResultsProcessed"]    = ""


	# ========== writing JSON object ==========
	with open('./data/infomation.json', 'w') as f:
		json.dump(Metadata, f, indent=4)
	print("Run Setup Complete")




def _internalsaves(fnV, fnoutV, fnP, fnoutP, fnT, fnoutT, encoding):
	"""This function is not available to the user as it depends on data too large to put in a 
	git repo, It lef here so i know how to redo the data if needs be."""


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
	else:
		# ========== Read in the file and get the size ==========
		dsV = xr.open_dataset(fnoutV)

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
	return dsV.latitude.size * dsV.longitude.size
#==============================================================================
if __name__ == '__main__':
	# ========== Set the args Description ==========
	description='Build the standard plots for each multi-run ensemble'
	parser = argparse.ArgumentParser(description=description)
	
	# ========== Add additional arguments ==========
	parser.add_argument(
		"--coarsen", type=int, default=0, help="The size of the box used to downscale data, Defualt is zeros")
	parser.add_argument(
		"-a","--annual", action="store_true", help="The size of the box used to downscale data, Defualt is zeros")
	parser.add_argument(
		"--maxacp", type=int, default=12, help="The maximum accumulation period")
	parser.add_argument(
		"--maxosp", type=int, default=4, help="The maximim ofset period")
	args = parser.parse_args() 
	
	# ========== Call the main function ==========
	main(args)
