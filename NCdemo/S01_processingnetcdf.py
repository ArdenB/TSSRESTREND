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
import argparse
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
def main(args):
	# =========== Read the metadata file in ==========
	if args.use_archived is None:
		infofile = './results/infomation.json'
	else:
		'./results/archive/infomation_%02d.json' % args.use_archived

	with open(infofile, "r+") as f:
		info = json.load(f)
	coarsen = info["coarsen"]

	# ========== Set the filenames ==========
	fnNDVI   = info["NDVI"]["fname"]
	fnPPT    = info["precipitation"]["fname"]
	fnTMEAN  = info["temperature"]["fname"]
	fnC4frac = info["C4fraction"]["fname"]

	# ========== Loop over the three datasets ==========
	for dsname, dsdesc in zip([fnNDVI, fnPPT, fnTMEAN, fnC4frac], ["ndvi", "ppt", "tmean", "C4frac"]):

		# ========== Read the dataset in ==========
		dsin = xr.open_dataset(dsname)

		# ========== deal with C3 and C4 ==========
		if dsdesc == "C4frac":
			# ========== Change the c4 fraction ==========
			if info["photo"] == "C4":
				dsin = xr.ones_like(dsin)
			elif info["photo"] == "C3":
				dsin = xr.zeros_like(dsin)
				
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
	# ========== Set the args Description ==========
	description='Convert the dataset to CSV'
	parser = argparse.ArgumentParser(description=description)
	
	# ========== Add additional arguments ==========
	parser.add_argument(
		"--use_archived", type=int, default=None, help="Use this argument to redo archived infomation.json files")
	args = parser.parse_args() 
	
	# ========== Call the main function ==========
	main(args)
