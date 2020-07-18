# Applying the TSS.RESTREND package to Netcdf File
This folder gives and example of how TSS-RESTREND R package can be used to do analysis on netcdf format spatial data.  It uses a mix of python and R because that is the way I process data. The implementation here is not computationally efficent and massive performace gains can be using pythons dask+xarray+scikit-learn to compute statistis on the entire dataset to calculate the observed change, the CO2, the climate change and climate varibility components.  Example datasets are located in [data](data).  

## Building a python environment using conda 

The following are instuctions on how to setup a conda python programming enviroment in any Ubuntu based linux distro including the Windows Subsystem for Linux (WSL).  Anaconda is availabe [here](https://www.anaconda.com/products/individual)

```bash

#  etting up Conda environment for the python parts of this code
conda config --add channels conda-forge
conda update -n base -c defaults conda

# Setup the new environment. These steps are optional but provent conflicts with existing python enviroments
conda create --name tssr
conda activate tssr

# install the relevant packages
conda install dask xarray pandas matplotlib palettable cartopy seaborn ipdb numba bottleneck netCDF4 webcolors gitpython geopandas shapely rasterio cdo python-cdo ipython

# Install other usefull commmand line utils
sudo apt-get intall ncview

```

## Building an R programming environment 

The latest version of R is available from [CRAN](https://cran.r-project.org/mirrors.html). I would also recomment [Rstudio](https://rstudio.com/products/rstudio/download/) as an IDE and for modifing the R packages that is included in this repo.

```R
# ========== Install dependincies ==========
install.packages(c("bfast", "broom", "RcppRoll", "forecast", "car", "rgl", "ggplot2", "mblm", "curl", "libcurl"))

# ========== Install TSS.RESTREND from cran ==========
install.packages("TSS.RESTREND")

# TSS.RESTREND can also be installed by building the  Rstudio project file TSS.RESTREND.Rproj in ./TSS.RESTREND
```

## Processing the netcdf files using python

There are two scripts used to process files. The primary one is [S01_processingnetcdf.py](S01_processingnetcdf.py) which opens the vegetation, precipitation and temperature netcdf files then creates 3 .csv files that can be read in my R.  The second is the [ProcessingPipline.py](ProcessingPipline.py) script which is only usefull for downscaling 

Script opens the ndvi, precip and temp netcdfs provided, reshapes them and converts them to a dataframe which is saved as a csv.

References for the original data: