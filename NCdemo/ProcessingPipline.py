import xarray as xr 
import dask 
import bottleneck as bn
import pandas as pd
import os


fnV    = "./vod_AUS_GIMMS_monthly_1982_2015.nc"
fnoutV = "./AUSdemo_GIMMS_ndvi.nc"
fnP    = "./TerraClimate_ppt_AUS_remapped.nc"
fnoutP = "./AUSdemo_TERRACLIMATE_ppt.nc"
fnT    = "./TerraClimate_tmean_AUS_remapped.nc"
fnoutT = "./AUSdemo_TERRACLIMATE_tmean.nc"
coarsen = 10

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
breakpoint()