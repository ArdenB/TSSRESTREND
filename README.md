# **Time Series Segmented Residual Trends (TSS-RESTREND)**

Time Series Segmented Residual Trends (TSS-RESTREND) is a method for the automated detection of land degradation from remotely sensed vegetation and climate datasets. TSS-RESTREND incorporates aspects of two existing degradation detection methods: RESTREND which is used to control for climate variability, and BFAST which is used to look for structural changes in the ecosystem.  

The full details of the testing and justification of the TSS-RESTREND method (version 0.1.02) are published in:      
* Burrell, Arden L., Jason P. Evans, and Yi Liu. 2017. “Detecting Dryland Degradation Using Time Series Segmentation and Residual Trend Analysis (TSS-RESTREND).” Remote Sensing of Environment 197 (August):43–57. <https://doi.org/10.1016/j.rse.2017.05.018> 
* Burrell, Arden L., Jason P. Evans, and Yi Liu. ‘The Impact of Dataset Selection on Land Degradation Assessment’. ISPRS Journal of Photogrammetry and Remote Sensing 146 (1 December 2018): 22–37. <https://doi.org/10.1016/j.isprsjprs.2018.08.017.>

The changes to the method included in version 0.2.15 focus on the inclusion of temperature as an additional climate variable. This allows for land degradation assessment in temperature limited drylands. The full details of this method were published in:
* Burrell, A. L., J. P. Evans, and Y. Liu. ‘The Addition of Temperature to the TSS-RESTREND Methodology Significantly Improves the Detection of Dryland Degradation’. IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 2019, 1–7. https://doi.org/10.1109/JSTARS.2019.2906466.

Version 0.3.0 is an expansion of the TSS-RESTREND package to perform additional vegetation change detection and attribution analysis.  It can be used to quantify vegetation change and then attribute it to some combination of CO2 fertilisation, climate change and climate variability.  For details on methods used to perform this attribution see:
* Burrell, A. L., J. P. Evans, and M. De Kauwe. ‘Anthropogenic Climate change has driven over 5 million km2 of drylands towards desertification’. 2020, in Press.  

## Logic and structure of the package 
This package was designed to work with a single pixel at a time, not 3D data structures like Rasters. This decision was made for a number of reasons the most significant of which is that processing spatial data, especially netcdf files, is much slower in R than using Python (using Xarray+Dask) or using the command line tool Climate Data Operators (CDO). It is computationally faster to turn spatial data into a 2D matrix/Dataframes then use the R package foreach and %do% to iterate through each row of the dataframe.  This also has the advantage of allowing the user to switch to a parallel workflow easily using backends like doParallel, doMC, doMPI or doSNOW.  

An example showing one way to process spatial data, then apply TSS-RESTREND to it, convert back to a spatial format and produce maps can be found [here](NCdemo).  

There are two main ways to use this package.  The easiest way is the new [TSSRattribution](TSS.RESTREND/man/TSSRattribution.Rd) function.  It takes a compete time seres ov monthly vegetation (CTSR,VI), precipitation (CTSR.RF), temperature (CTSR.TM) as well as the maximum accumulation and offset periods to consider.  

```R
library(TSS.RESTREND)
max.acp <- 12
max.osp <- 4
# Load in the data here 
# CTSR.VI <- LOAD YOUR VEGETATION DATA HERE. Must me a ts object with a monthly frequency
# CTSR.RF <- LOAD YOUR PRECIPITATION DATA HERE. Must me a ts object with a monthly frequency
# CTSR.TM <- LOAD YOUR TEMPERATURE DATA HERE. Must me a ts object with a monthly frequency

results = TSSRattribution(CTSR.VI, CTSR.RF, CTSR.TM, max.acp, max.osp)

print(results)
```

The original way to use this package was through the TSSRESTREND function. 
```R
library(TSS.RESTREND)
#Define the max accumuulation period
max.acp <- 12
#Define the max offset period
max.osp <- 4
#Create a table of every possible precipitation value given the max.acp and max.osp
ACP.table <- climate.accumulator(CTSR.VI, CTSE.RF, max.acp, max.osp)
ACT.table <- climate.accumulator(CTSR.VI, CTSE.TM, max.acp, max.osp, temperature=TRUE)

# perform the TSSRESTREND, retnonsig=FALSE replicates defualt behaviour is version<0.3.0
results <- TSSRESTREND(CTSR.VI, ACP.table, ACT.table, retnonsig=FALSE)
print(results)
plot(results, verbose=TRUE)
``` 

NOTE: TSSRESTREND was designed to demonstrate the methodology and is not optimised for perforce. If analysing a large area, major performance gains can be had by doing specific steps of the method on the entire dataset at once rather than on a per pixel basis.  

## FAQ
1. How do I make TSS-RESTREND work on spatial data.  
An example pipeline showing how this method can be applied to spatial data can be found [here](NCdemo).  

## **Change log**
### **TSS-RESTREND v0.3.0 (Release date: 2020-07-16)**
**Major changes:**

* Added a new wrapper function [TSSRattribution](TSS.RESTREND/man/TSSRattribution.Rd) that will perform the full attribution described in Burrell et al., (2020).  
* Added a function for removing the effects of eCO2 from a VI time series.  
* Added a function for detrending climate data using a moving window smoothing,  

**Minor tweaks and bug fixes:**
* Modified the [climate.accumulator]( TSS.RESTREND/man/climate.accumulator.Rd) function to allow rolling means of temperature data not just rolling sums
* Added the option to turn off the sigmasking in parts of the analysis. This means an estimate for total change will be calculated even when the underlying vegetation climate relationship is not significant.  Useful when doing multirun ensembles where a pixels is not considered in isolation
* In some rare cases TSSRESTREND function would return Total Change estimates outside of the range of the VI data passed.  I Was unable to replicate this error, but a catch statement has been added to detect these errors and flag them for the user. The TSSRESTREND function will now return result$summary$Method = "InvalidValueError" if it detects a total change larger than the range on the input vegetation data.  


### **TSS-RESTREND v0.2.15 (Release date: 2018-01-23)**
**Major changes:**

* Added the ability to include temperature as an additional climate variable. This involves changes to almost all parts of the package and new plotting functionality.   
* METHOD TWEAK: If the complete time series VPR/VCR relationship is not significant. TSS-RESTREND no longer returns indeterminate. It instead tests BFAST on the complete time series of Vegetation and uses those breakpoints.  This will cause minor differences in the results between earlier versions.

**Minor tweaks and bug fixes:**
* Large improvements to the internal documentation 
* Speed improvements made throughout the functions  
* Add the ability to vary BFAST window (h) 
* ACPcalculator updated to deal with places where there are negative CTS VPR relationship. Small tweaks to TSS.RESTREND master script and the Annual.Precipitation functions
* BUG FIX: in the rainfall.accumulator now has a check that makes sure the rainfall data has some variance (SD >0). No variance (SD == 0) breaks lm() calculations and causes failures
* BUG FIX: Fixed a bug in the indexing of AnnualRF.Cal (now AnnualClim.Cal) and the rainfall.accumulator (now climate.accumulator) that was causing downstream problems with the segVPR function. This bug may have slightly impacted some results in segVPR 



### Dependencies 
TSS.RESTREND depends on:
```R
stats, graphics, utils, bfast, broom, strucchange, ggplot2, RcppRoll
```
It also suggests two packages for 3D plots: 
```R
rgl, car
```
Because of the difficult of installing these packages on some unix based systems, these packages were not made dependencies
