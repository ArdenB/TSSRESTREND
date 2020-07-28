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
# ===== append that to the system path =====
sys.path.append(os.getcwd())


# ========== Import packages ==========
import numpy as np
import pandas as pd
import argparse
# import datetime as dt
import warnings as warn
import xarray as xr
import dask
import bottleneck as bn
from numba import jit
from collections import OrderedDict
import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import palettable
import statsmodels.stats.multitest as smsM
import scipy.stats as sps

# ========== Load custom functions ==========
import CustomFunctions as cf 

# ==============================================================================
def main(args):
    # =========== Read the metadata file in ==========
    if args.use_archived is None:
        infofile = './results/infomation.json'
    else:
        './results/archive/infomation_%02d.json' % args.use_archived

    with open(infofile, "r+", errors='ignore') as f:
        info = json.load(f)
        # ========== fix the timedelta ==========
        if type(info["ComputeTime"]) == float:
            info["ComputeTime"] = pd.Timedelta(info["ComputeTime"], unit="sec")

    # ========== Open the csv results file ==========
    fn = "./results/AttributionResults.csv"
    df = pd.read_csv(fn, index_col=0)
    
    # ========== Fix the indexing and the error messaging ==========
    lonlat = np.array([_unpackLatLon(index) for index in df.index])
    df["longitude"] = lonlat[:, 0]
    df["latitude" ] = lonlat[:, 1]
    df = df.reset_index(drop=True).set_index(["longitude","latitude" ])
    df["errors"] = df["errors"].apply(_df_error)

    df["AnthropogenicClimateChange"] = df.ClimateChange + df.CO2
    # ========== function to calculate adjusted p values for ACC ==========
    pstack = np.stack((df["ClimateChange.Pvalue"].values, df["CO2.Pvalue"].values))
    df["AnthropogenicClimateChange.Pvalue"]= np.apply_along_axis(_combine_pvalue2d, 0, pstack)

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
    ds = ds.transpose('time', 'latitude', 'longitude').sortby("latitude", ascending=False)
    

    
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
    # ========== Write the dataset to a netcdf file ==========
    print("Starting write of data at:", pd.Timestamp.now())
    ds.to_netcdf("./results/TSSRattribution_Results_%s_%dkm.nc" % (info["photo"], (25*(info["coarsen"]-1)+25)), 
        format         = 'NETCDF4', 
        encoding       = encoding,
        unlimited_dims = ["time"])

    # ========== Make a series of maps ==========
    if args.plots:
        # ========== set the plot ==========
        plotpath = "./results/plots/"

        for va in ds.variables:
            # +++++ check its a variable to be skipped or mapped +++++
            if (va in ['longitude', 'latitude', "OtherFactorsValid", "errors", "time"] or va.endswith(".Pvalue")):
                continue

            # ========== Perform significance correction ==========
            if args.sig:
                signif, s_method = FDRSignificanceCorrection(ds, va, args.method)
            else:
                s_method = ""
                signif   = None
            # ========== Create the metadata for the plots ==========
            maininfo = "Plot from %s (%s):%s by %s, %s" % (__title__, __file__, 
                __version__, __author__, str(pd.Timestamp.now()))
            maininfo += s_method
            gitinfo   = cf.gitmetadata()
            
            print(va)
            # ========== Build, save and show the maps ==========
            MapSetup(ds, va, signif, maininfo, gitinfo, info, plotpath)

# ==============================================================================
# BUILD THE MAPS
# ==============================================================================
def MapSetup(ds, va, signif, maininfo, gitinfo, info, plotpath,  mask=None):
    """
    Builds the maps of the change attribution
    args:
        ensinfo
        ds          xrds: TSSRESTREND RESULTS
        va          str : the name of the change variable
        maininfo    str : plot string
        gitinfo     str : summary of the git header
        mask        xrda: Any additional boolean masks 
    """ 
    # ========== Create the mapdet object ==========
    mapdet = cf.mapclass("Australia", plotpath)
    mapdet.dpi  = 130
    mapdet.var  =  va 
    mapdet.desc =  "%s_%s_%dkm" % (va, info["photo"], (25*(info["coarsen"]-1)+25)) 
    # ========== Create the data array ==========
    da = ds[va].copy()
    ad = da.squeeze()
    # ========== Remove the non dryland areas ==========
    if not (mask is None):
        da *= mask
    
    # ========== Mask for FDR significnace ==========
    if not signif is None:
        da *= signif
        mapdet.desc += "_withNonSigMasked"
    # ========== calculate the tick position  ==========
    if info['annual']:
        mapdet.cblabel = r"$\times$10$^{-2}$ $\Delta$ NDVI$_{max}$ yr$^{-1}$" 
        # scale the NDVI to make it fit 
        da *= 1e2 
    else:
        mapdet.cblabel = r'$\Delta$ NDVI$_{max}$'


    if info['Nyears'] != 34.:
        warn.warn("The plot tick values were chosen for a time series of 34 years.  Other lengths may require differen colorbars")


    tks = np.array([-0.30, -0.24, -0.12, -0.06, -0.02, 0, 0.02, 0.06, 0.12, 0.24, 0.30])
    mapdet.cmin  = -0.3 #-0.1   # the min of the colormap
    mapdet.cmax =  0.3#0.1  # the max of the colormap

    # ========== Set the colormap ==========
    if va == "ObservedChange":
        cmapHex = palettable.colorbrewer.diverging.PRGn_10.hex_colors
    elif va in ["ObservedClimate", "ClimateChange"]:
        cmapHex     = palettable.colorbrewer.diverging.BrBG_10.hex_colors
    elif va == "ClimateVariability":
        cmapHex     = palettable.colorbrewer.diverging.RdBu_10.hex_colors
    elif va == "CO2":
        tks = np.array([0, 0.001, 0.02, 0.03, 0.04, 0.06, 0.24, 0.30])
        cmapHex = palettable.colorbrewer.sequential.Purples_7.hex_colors
        mapdet.cmin =  0.0#-0.1     # the min of the colormap
    elif va == "LandUse":
        cmapHex     = palettable.colorbrewer.diverging.PiYG_10.hex_colors
    elif va == "OtherFactors":
        # cmapHex = palettable.colorbrewer.diverging.RdBu_10.hex_colors
        cmapHex = palettable.cmocean.diverging.Curl_10_r.hex_colors
    elif va == "AnthropogenicClimateChange":
        cmapHex = palettable.colorbrewer.diverging.PuOr_10.hex_colors

        
    
    # ========== Set the variables ==========
    mapdet.cZero =  0.5     # the zero point of the colormap
    
    # ========== Add the Font info ========== 
    mapdet.gridalp  = 0.5

    # ========== Setup the cmap ==========
    if va == "CO2":
        spacing = 'uniform'
        cbounds = tks
        ticks   = tks
        # ========== create the colormap ==========
        cmap   = mpl.colors.ListedColormap(cmapHex)

        cmap.set_over(cmapHex[-1])
        
        norm   = mpl.colors.BoundaryNorm(cbounds, cmap.N)
        mapdet.extend  = "max"
    else:
        cmap, norm, ticks, cbounds, spacing = cf.ReplaceHexColor(
            cmapHex, mapdet.cmin, mapdet.cmax, ticks=tks, zeroR = 0.001)
        # print(cbounds)

    # cmap, norm, ticks, cbounds, spacing = pf.ReplaceHexColor(
    #   cmapHex, mapdet.cmin, mapdet.cmax, ticks=tks, zeroR = 0.001)

    cmap.set_bad('dimgrey',1.)
    
    mapdet.cmap    = cmap 
    mapdet.norm    = norm
    mapdet.ticks   = ticks  
    mapdet.cbounds = cbounds
    mapdet.spacing = spacing

    # ========== Make a map ========== 
    plotinfo, fname = cf.mapmaker(da, mapdet) 

    # ========== Save the metadata ========== 
    if fname:
        infomation = [maininfo, plotinfo, fname, gitinfo]
        cf.writemetadata(fname, infomation)
    # ipdb.set_trace()

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
        pnm = var + ".Pvalue"


    # ========== Locate the p values and reshape them into a 1d array ==========
    # ++++++++++ Find the pvalues ++++++++++
    try:
        ODimDa     = ds[pnm].stack(loc = ('time', 'latitude', 'longitude')).copy()
    except:
        breakpoint()
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
        attributes  Ordered Dictionary cantaining the attribute infomation
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
        index:  str
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

def _combine_pvalue2d(pvals):
    """
    takes an array and removes nans then combines
    args:
        pvalues:    np.array
    """
    if bn.allnan(pvals):
        return np.NAN
    elif np.sum(pvals) == 0:
        return 0
    else:
        ___, pv = sps.combine_pvalues(pvals[~np.isnan(pvals)])
        return pv
#===================
# ==============================================================================
if __name__ == '__main__':
    description='Arguments that can be passed to the saving and plotting script'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-p", "--plots", action="store_false", 
        help="Quiet plots: if selected plots will not be displayed ")
    parser.add_argument(
        "-s", "--sig", action="store_true", 
        help="Significance: Apply a zero mask using FDR adjustment and the Benjamini/Hochberg method")
    parser.add_argument(
        "--method", type=str,choices=["fdr_bh", 'fdr_by'], default="fdr_bh", help="The method used to adjust for False Discovery Rate. must be fdr_bh or fdr_by")
    parser.add_argument(
        "--use_archived", type=int, default=None, help="Use this argument to redo archived infomation.json files")
    args = parser.parse_args() 
    main(args)

