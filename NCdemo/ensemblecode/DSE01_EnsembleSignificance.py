# ==============================================================================

__title__   = "csv to netcdf and maps"
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
# ===== append that to the system path =====
sys.path.append(os.getcwd())

import numpy as np
import pandas as pd
import xarray as xr
import bottleneck as bn
from scipy import stats
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as smsM
import cartopy.crs as ccrs


# ==============================================================================
def main():

    # The ensemble result for one TSSRattribution category
    ds_ens = xr.open_dataset("./data/AUS_TSSRattribution_ensemble_example.nc")
    # ds_ens = xr.open_dataset("./data/workingset/AUS_ExportedTSSResults_ClimateChange_ensemble.nc")

    # ===== print out the datasets =====
    print(ds_ens)

    # Set the name of the var and the pvar
    var  = "ClimateChange"
    pvar = "ClimateChange.Pvalue"

    # ========== Pull out the var and the pvalues ==========
    raw_cv = ds_ens[var]
    raw_pv = ds_ens[pvar]
    raw_pv = raw_pv.where(raw_pv != 0, 0.000000001) 

    # ========== IPCC signifcance ==========
    fdr_ls = [FDRSignificanceCorrection(raw_pv.isel(run=rn).values) for rn in raw_pv.run]
    FDRbool = raw_pv.copy(data=np.stack(fdr_ls, axis=-1))
    IPCCsig  = xr.apply_ufunc(IPCC,raw_cv, FDRbool, input_core_dims=[["run"], ["run"]], vectorize=True)

    # ========== Combined pvalue significance ==========
    pvals    = xr.apply_ufunc(_combine_pvalue2d, raw_pv, input_core_dims=[["run"]], vectorize=True) #"time", "latitude", "longitude"
    FDRpvsig =  pvals.copy(data = FDRSignificanceCorrection(pvals.values))
    
    # ========== Calculate the ensembe mean ==========
    chmean = raw_cv.mean(dim="run")
    p = pvals.plot()
    plt.show()
    breakpoint()
# ==============================================================================

def _combine_pvalue2d(pvals):
    """
    takes an array and removes nans then combines
    args:
        pvalues:    np.array
            the pvalues to be combined 
    return:
        
    """

    if bn.allnan(pvals):
        return np.NAN
    else:
        ___, pv = stats.combine_pvalues(pvals[~np.isnan(pvals)])
        return pv

#============================================================================== 

def FDRSignificanceCorrection(da, FDRmethod = "fdr_bh", alpha = 0.10):
    """
    Takes the results of an existing trend detection aproach and modifies them to
    account for multiple comparisons.  
    args
        da: 3d Numpy array. 
            The time dime should have a size of 1
        FDRmethod: str
            The method to use for FDR
        alpha: float (9-1)
            Significance level
    returns:
        res: n-darray
            Same shape as da.  ! for significant, 0 for not significant. NaN for bad/no data
    """
    if FDRmethod == "fdr_by":
        t_method = "Adjusting for multiple comparisons using Benjamini/Yekutieli"
    elif FDRmethod == "fdr_bh":
        t_method = "Adjusting for multiple comparisons using Benjamini/Hochberg"
    else:
        raise ValueError("unknown MultipleComparisons method: %s, must be either fdr_by or fdr_bh " % FDRmethod)

    # ========== Work out the pvalue name =========== 
    shape = da.shape
    # ========== Locate the p values and reshape them into a 1d array ==========
    # ++++++++++ Find the pvalues ++++++++++
    ODimDa     = da.ravel()
    isnan      = np.isnan(ODimDa)
    
    # ++++++++++ pull out the non nan pvalus ++++++++++
    pvalue1d   = ODimDa[~isnan]
    
    # =========== Perform the MC correction ===========
    pvalue_adj = smsM.multipletests(pvalue1d, method=FDRmethod, alpha=alpha)
    # =========== Get the significance in bool ===========
    sigfloat   = pvalue_adj[0].astype(float)

    # make an empty dataarray
    re         = np.zeros_like(ODimDa)
    re[~isnan] = sigfloat
    re[isnan]  = np.NaN

    res = re.reshape(shape)
    # print(pvalue_adj[0], pvalue_adj[0].shape)
    return res

def IPCC(direction, significance):
    """
    Function to caluculate the significance using ipcc criteria
    args:
        direction: npA
            the sign of the change
        significance: npA
            FDR adjusted significance
    returns:
        ipcc_sig:   npA
            0 for not sig, 1 for sig
    """
    if bn.allnan(direction):
        return np.NaN
    else:
        totnum = bn.nansum((~np.isnan(direction)).astype(float))
        direc  = direction[~np.isnan(direction)]
        signi  = significance[~np.isnan(direction)]


        # Find the number of significant runs
        signum = np.sum(signi==1).astype(float)
        if signum == 0:
            # Remove the badly behaving zeros
            return 0

        # ========== IPCC significance testing ==========
        # +++++ 50% significant +++++
        sig50  = signum>=0.50
    
        # # +++++ 80% agreement +++++
        # agreement
        ag = [np.sum(np.logical_and((direc > 0), (signi == 1))), np.sum(np.logical_and((direc < 0), (signi == 1)))]
        agree = float(bn.nanmax(ag))

        agr80 = (agree/signum)>= 0.80
        return int(agr80)

#============================================================================== 

if __name__ == '__main__':
    main()