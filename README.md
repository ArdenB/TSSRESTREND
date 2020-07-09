# **Time Series Segmented Residual Trends (TSS-RESTREND)**


Time Series Segmented Residual Trends (TSS-RESTREND) is a method for the automated detection of land degradation from remotely sensed vegetation and climate datasets. TSS-RESTREND incorporates aspects of two existing degradation detection methods: RESTREND which is used to control for climate variability, and BFAST which is used to look for structural changes in the ecosystem.  

The full details of the testing and justification of the TSS-RESTREND method (version 0.1.02) are published in:
      
* Burrell, Arden L., Jason P. Evans, and Yi Liu. 2017. “Detecting Dryland Degradation Using Time Series Segmentation and Residual Trend Analysis (TSS-RESTREND).” Remote Sensing of Environment 197 (August):43–57. <https://doi.org/10.1016/j.rse.2017.05.018> 
* Burrell, Arden L., Jason P. Evans, and Yi Liu. ‘The Impact of Dataset Selection on Land Degradation Assessment’. ISPRS Journal of Photogrammetry and Remote Sensing 146 (1 December 2018): 22–37. <https://doi.org/10.1016/j.isprsjprs.2018.08.017.>

The changes to the method included in version 0.2.15 focus on the inclusion of temperature as an additional climate variable. This allows for land degradation assessment in temperature limited drylands. The full details of this method were published in:

* Burrell, A. L., J. P. Evans, and Y. Liu. ‘The Addition of Temperature to the TSS-RESTREND Methodology Significantly Improves the Detection of Dryland Degradation’. IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 2019, 1–7. https://doi.org/10.1109/JSTARS.2019.2906466.


Version 0.3.0 is an expansion of the TSS-RESTREND package to perform additional vegetation change detection and attribution analysis.  It can be used to quantify vegetation change and then attribute it to some combination of CO2 fertilisation, climate change and climate varibility.  For details on methods used to perform this attribution see:

   
Basic overview of the attribution:
* Observed changes is ____

Conceptual requirements.
* For this approach to be valid _____

## Logic and structure of the package 

________________
The _____ function is the main way to this package. It calls ___ in order and returns ____

Versions 0.3.0 now has a total change detection function. it calls _____. 

It is designed to demonstrate the methodology and is not optimised for perforce. If analysiing a large area, large performace gains can be had by doing specific steps of the method on the entire dataset at once rather than on a per pixel basis.  



## FAQ

1. How do i make it work on spatial data.  ______________________

An example pipeline showing how this method can be applied to 

## Examples

A full list of examples and extra documentation are currently under development 

## **Change log**
* Added a function for removing the errects of eCO2 from a VI time series.  
* Modified the climate.accumulator function to allow rolling means of temperature data not just rolling sums
* Added the option to turn off the sigmasking in parts of the analysis. Usefull when doing multirun ensembles. 

###**TSS-RESTREND v0.2.15 (Release date: 2018-01-23)**


### Major changes: 
* Added the ability to include temperature as an additional climate variable. This involves changes to almost all parts of the package and new plotting functionality.   
* METHOD TWEAK: If the complete time series VPR/VCR relationship is not significant. TSS-RESTREND no longer returns indeterminate. It instead tests BFAST on the complete time series of Vegetation and uses those breakpoints.  This will cause minor differences in the results between earlier versions.

### Minor tweaks and bug fixes:
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