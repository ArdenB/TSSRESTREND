"""
Script goal, 

Open the ndvi, precipitation and temperature netcdf then convert them to a csv 
that R can read in.  


"""

# ==============================================================================

__title__ = "Process Netcdf into csv"
__author__ = "Arden Burrell"
__version__ = "v1.0(23.06.2020)"
__email__ = "aburrell@whrc.org"

# ==============================================================================

# ========== Import packages ==========
import os
import sys
import numpy as np
import pandas as pd
# import datetime as dt
# import warnings as warn
import xarray as xr
# import dask
# import bottleneck as bn
from collections import OrderedDict, defaultdict
import json

# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import seaborn as sns

# ==============================================================================
def main():
	# =========== Read the metadata file in ==========
	infofile = './data/infomation.json'
	with open(infofile, "r+") as f:
		info = json.load(f)
	coarsen = info["coarsen"]
	# ========== Set the filenames ==========
	if coarsen == 0:
		fnNDVI  = "./data/AUSdemo_GIMMS_ndvi.nc"
		fnPPT   = "./data/AUSdemo_TERRACLIMATE_ppt.nc"
		fnTMEAN = "./data/AUSdemo_TERRACLIMATE_tmean.nc"
	else:
		fnNDVI  = "./data/AUSdemo_GIMMS_ndvi_xrcoarsen_%dwin.nc" % coarsen
		fnPPT   = "./data/AUSdemo_TERRACLIMATE_ppt_xrcoarsen_%dwin.nc" % coarsen
		fnTMEAN = "./data/AUSdemo_TERRACLIMATE_tmean_xrcoarsen_%dwin.nc" % coarsen

	# ========== Loop over the three datasets ==========
	for dsname, dsdesc in zip([fnNDVI, fnPPT, fnTMEAN], ["ndvi", "ppt", "tmean"]):

		# ========== Read the dataset in ==========
		dsin = xr.open_dataset(dsname)

		# ========== Stack the dataset ==========
		ds_stack = dsin.stack(cord=('longitude', 'latitude'))

		# ========== Convert to a pandas datadrame ==========
		df_out =  pd.DataFrame(ds_stack[dsdesc].values.T, index=ds_stack.cord.values, columns=ds_stack.time.values)

		# ========== Create a file name out ==========
		fnout = "./data/demo_dataframe_%s.csv" % (dsdesc)

		# ========== Save the file out ==========
		df_out.to_csv(fnout)

		print(dsdesc, " CSV exported")
	
	# ========== writing JSON object ==========
	info["history"] = "%s: Datasets converted to CSV (%s). " % (str(pd.Timestamp.now()), __file__) + info["history"]
	info["DataProcessed"] = str(pd.Timestamp.now())
	os.remove(infofile)
	with open(infofile, 'w') as f:
		json.dump(info, f, indent=4)
	print("Data converted into csv and info file updated at:", pd.Timestamp.now())

# ==============================================================================
if __name__ == '__main__':
	main()