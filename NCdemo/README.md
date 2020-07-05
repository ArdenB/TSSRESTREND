# NC demonstration 
This folder gives and example of how TSS-RESTREND R package can be used to do analysis on netcdf format spatial data.  It uses a mix of python and R because that is the way I process data. The implementation here is not computationally efficent and massive performace gains can be using pythons dask+xarray+scikit-learn to compute statistis on the entire dataset to calculate the observed change, the CO2, the climate change and climate varibility components.  All scripts are to be run from the same dir that the provided datasets are located in

## process the netcdf files

S01_processingnetcdf.py

Script opens the ndvi, precip and temp netcdfs provided, reshapes them and converts them to a dataframe which is saved as a csv.


