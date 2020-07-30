# ==============================================================================

__title__   = "ensenble significance"
__author__  = "Arden Burrell"
__version__ = "v1.0(23.06.2020)"
__email__   = "aburrell@whrc.org"

# ==============================================================================
import os
import sys
# ===== CHange the dir to the script location =====


if not os.path.dirname(sys.argv[0]) == "":
    os.chdir(os.path.dirname(sys.argv[0]))
    if not os.getcwd().endswith("NCdemo"):
    if "NCdemo" in os.getcwd():
        p1, p2, _ =  os.getcwd().partition("NCdemo")
        os.chdir(p1+p2)
import numpy as np
import pandas as pd
import xarray as xr
import bottleneck as bn
from scipy import stats
import matplotlib.pyplot as plt
# import statsmodels.stats.multitest as smsM

# ==============================================================================

def main()
	# The C3 run
	ds_c3 = xr.open_dataset("./results/TSSRattribution_Results_C3_125km.nc")
	# The C4 run
	ds_c4 = xr.open_dataset("./results/TSSRattribution_Results_C3_125km.nc")
	# The fraction of C4 plans 
	C4frac = xr.open_dataset("./data/AUSdemo_SYNMAP_C4Fraction_xrcoarsen_5win.nc")

	# ===== print one of the datasets =====
	print(ds_c3)

	# Set the name of the var and the pvar
	var  = "ClimateChange"
	pvar = "ClimateChange.Pvalue"
	# Calculate tthe adjusted 
	VARadj = Weighted_mean(ds_c3, ds_c4, C4frac, var)
	VARadj_pvalues = Combinee_Pvalues(ds_c3, ds_c4, C4frac, pvar)
	breakpoint()
# ==============================================================================

def Weighted_mean(ds_c3, ds_c4, C4frac, var):
    # pull out the dataarrays
    da_c3 = ds_c3[var]
    da_c4 = ds_c4[var]
    C4f   = C4frac["C4frac"].values

    #Calculate the weighted mean
    c4_adj = (da_c4 * C4f) + (da_c3 * (1-C4f))

    #return the results
    return c4_adj


def Combinee_Pvalues(ds_c3, ds_c4, C4frac, pvar):
	# Pull out the p values

	# Stack long and land to make 1d array
	da_c3p = ds_c3[pvar].stack(cord=('longitude', 'latitude'))
	da_c4p = ds_c4[pvar].stack(cord=('longitude', 'latitude'))
	C4f_p  = C4frac["C4frac"].stack(cord=('longitude', 'latitude'))
	C3f_p  = 1-C4f_p

	# stack into a 2d array so i can use  np.apply_along_axis 
	stacked = np.vstack([da_c3p.values, da_c4p.values, C4f_p.values, C3f_p.values])

	def _combine_2pvalue(array):
		pvals  = array[:2]
		wei    = array[2:]
		if bn.allnan(pvals):
			return np.NaN
		else:
			# Deal with pvalues that are= 0. It happens when the values are smaller than the is possible with float 32.  If they are pased to 
			# stats.combine_pvalues it will return nan
			pvals[pvals==0] = 0.000001
			#Deals with places where thing failed to get a p values
			pvals[np.isnan(pvals)] = 1 

			# Adjust the pvalues
			sta, pv = stats.combine_pvalues(pvals, method="stouffer", weights=wei)
			return np.array(pv)

	# ===== Apply the 2d pvalue combination function to calculate the new pvalues =====
	res = np.apply_along_axis(_combine_2pvalue, 0, stacked).reshape([1,-1])

	# Convert back intoa a datarray 
	pv_adj = da_c3p.copy(data=res)
	# Unstack the lats and lons then return the result
	return pv_adj.unstack().transpose('time', 'latitude', 'longitude').sortby("latitude", ascending=False)


# ==============================================================================

if __name__ == '__main__':
	main()