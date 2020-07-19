"""
Script goal, 

Open the results from TSS-RESTREND attribution analysis and convering that 
back to a spatial data format.  


"""

# ==============================================================================

__title__ = "csv to netcdf"
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

# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import seaborn as sns

# ==============================================================================
def main():
	# ========== Open the csv results file ==========
	fn = "./results/AttributionResults.csv"
	df = pd.read_csv(fn, index_col=0)
	breakpoint()

# ==============================================================================
if __name__ == '__main__':
	main()