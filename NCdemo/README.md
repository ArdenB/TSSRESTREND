# Applying the TSS.RESTREND package to Netcdf spatial data
This folder gives and example of how TSS-RESTREND R package can be used to do analysis on netcdf format spatial data.  It uses a mix of python and R because that is the way I process data. The implementation here is not computationally efficent and massive performace gains can be using pythons dask+xarray+scikit-learn to compute statistis on the entire dataset to calculate the observed change, the CO2, the climate change and climate varibility components.  

There are four netcdf files included in this repo, GIMMS derived NDVI (1982-2015), TERRACLIMATE derived precipitation (1960-2015), TERRACLIMATE derived temperature (1960-2015) and SYNMAP C4 vegetation fraction. To make the datasets small enough that they can be stored on git, they have all been downsampled to approximatly 25km and their spatial coverage has been limited to Australia. Example datasets are located in [data](data). They are included for demonstaration purposes and should not be used for scientific analysis.  

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

To demonstrate how TSS-RESTREND can be applied to spatial data, this repo this repo contains four script that, if run in order will setup the analysis, process the data, perform TSS-RESTREND on the data, reassemble the results into a netcdf file and save maps of the results.  The results will be equivilant to a single run from Burrell et al., (2020).  This example is split across four scripts that need to be run in order from S00 to S03.   

#### Setting up the metadata ####

The first script that needs to run is [S00_SetupMetadata.py](./S00_SetupMetadata.py).  This script saves a json file that is passed between all the following scripts with metadata about the run. It has six possible command line arguments:

 -  -h, --help      show help message and exit

 -  -c, --coarsen 	The size of the box used to downscale data. Defualt = zero. Must be an int. Passing a larger coarsen values will speed up analysis at the cost of resolution.  The defualt of 0 means the analysis will occur at 25km of the demonstration datasets and will take many hours.  Passing value of 10 should downscales the data to 250km pixels which allow all four demo scripts to be run in less than 30 minuts.  

 -  -y, --yearly	When calculating TSS-RESTREND, report values in change per year not Total Change. Defualt is Total Change

 -  --maxacp    	The maximum accumulation period in months. Must be an int. defulat = 12. See TSSRESTREND R package documentation. 

 -  --maxosp       	The maximim ofset period in months. Must be an int. defulat = 4. See TSSRESTREND R package documentation. 

 -  --photo 		The photosyenthetic pathyway to fit for calculating the CO2 effect size. Possilbe arguments: {"C3andC4","C3","C4"}. Defulat is "C3andC4" which uses the demo SYNMAP C4 vegetation fraction.  

 -  -a, --archive	Archive existing infomation.json file rather than overwriting it. The existing infomation.json file is moved to ./results/archive/.  With some tweaks this could be used to allow for multiple runs withour overwriting results.  


For this example we will use a coarsen value of 3 and the defualts for all other arguments.  

```sh
# Run in console
python S00_SetupMetadata.py -c 3
```

#### Reshape the demo data into two dimensional csv files ####


The [S01_processingnetcdf.py](S01_processingnetcdf.py) opens the demo vegetation, precipitation, temperature and C4 fraction netcdf files then creates 2 dimensional .csv files that can be read in by R. R can be very slow when worting with higher dimensional rasters. Using python to convert the 3D netcdf files to a 2D .csv file with the lat and long as the rowname and each column being a timestep reduces overall computation time. Which files are processed is determined using the infomation.json file produced in the previous step.  

- --use_archived		Use this argument to redo archived infomation.json files.  Must be an int which correpondes to the number of the desired ./data/archive/infomation_{int}.json file. 

```sh
# Run in console
python S01_processingnetcdf.py
```


