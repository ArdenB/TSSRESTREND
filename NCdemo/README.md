# Applying the TSS.RESTREND package to Netcdf File
This folder gives and example of how TSS-RESTREND R package can be used to do analysis on netcdf format spatial data.  It uses a mix of python and R because that is the way I process data. The implementation here is not computationally efficent and massive performace gains can be using pythons dask+xarray+scikit-learn to compute statistis on the entire dataset to calculate the observed change, the CO2, the climate change and climate varibility components.  

There are three netcdf files included in this repo, GIMMS derived NDVI (1982-2015), TERRACLIMATE derived precipitation (1960-2015) and TERRACLIMATE derived temperature (1960-2015). To make the datasets small enough that they can be stored on git, they have all been downsampled to approximatly 25km and their spatial coverage has been limited to Australia. Example datasets are located in [data](data). They are included for demonstaration purposes and should not be used for scientific analysis.  

## 1. Setup a programming environment
**Building a python environment using conda**

The following are instuctions on how to setup a conda python programming enviroment in any Ubuntu based linux distro including the Windows Subsystem for Linux (WSL).  Anaconda is availabe [here](https://www.anaconda.com/products/individual)

```bash

#  etting up Conda environment for the python parts of this code
conda config --add channels conda-forge
conda update -n base -c defaults conda

# Setup the new environment. These steps are optional but provent conflicts with existing python enviroments
conda create --name tssr
conda activate tssr

# install the relevant packages
conda install dask xarray pandas matplotlib palettable cartopy seaborn ipdb numba bottleneck netCDF4 webcolors gitpython geopandas shapely rasterio cdo python-cdo ipython statsmodels

# Install other usefull commmand line utils
sudo apt-get intall ncview

```

**Building an R programming environment**

The latest version of R is available from [CRAN](https://cran.r-project.org/mirrors.html). The TSSRESTREN package was built with [Rstudio](https://rstudio.com/products/rstudio/download/).

```R
# ========== Install dependincies ==========
install.packages(c("bfast", "broom", "RcppRoll", "forecast", "car", "rgl", "ggplot2", "mblm", "curl", "libcurl", "rjson", "optparse"))

# ========== Install TSS.RESTREND from cran ==========
install.packages("TSS.RESTREND")
```
The other option is to use this repo to build TSS-RESTREND. Rstudio project file TSS.RESTREND.Rproj is located [here](../TSS.RESTREND/)
```R
# The following packages are needed to build TSS.RESTREND in Rstudio
install.packages(c("devtools", "roxygen2"))

```

## 2. Example Ananlysis

To demonstrate how TSS-RESTREND can be applied to spatial data, this repo this repo contains four script that, if run in order will setup the analysis, process the data, perform TSS-RESTREND on the data, reassemble the results into a netcdf file and save maps of the results.  The results will be equivilant to a single run from Burrell et al., (2020).  

There are two scripts used to process files. The primary one is [S01_processingnetcdf.py](S01_processingnetcdf.py) which opens the vegetation, precipitation and temperature netcdf files then creates 3 .csv files that can be read in my R.  The second is the [ProcessingPipline.py](ProcessingPipline.py) script which is only usefull for downscaling 

Script opens the ndvi, precip and temp netcdfs provided, reshapes them and converts them to a dataframe which is saved as a csv.

