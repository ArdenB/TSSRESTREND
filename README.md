# Time Series Segmented Residual Trends (TSS-RESTREND)

Time Series Segmented Residual Trends (TSS-RESTREND) is a method for the automated detection of land degradation from remotely sensed vegetation and climate datasets. TSS-RESTREND incorporates aspects of two existing degradation detection methods: RESTREND which is used to control for climate variability, and BFAST which is used to look for structural changes in the ecosystem.  

The full details of the testing and justification of the TSS-RESTREND method (version 0.1.02) are published in:
      
  Burrell, Arden L., Jason P. Evans, and Yi Liu. 2017. “Detecting Dryland Degradation Using Time Series Segmentation and Residual Trend Analysis (TSS-RESTREND).” Remote Sensing of Environment 197 (August):43–57. <https://doi.org/10.1016/j.rse.2017.05.018> 
  
The changes to the method included in version 0.2.13 focus on the inclusion of temperature as an additional climate variable. This allows for land degradation assessment in temperature limited drylands. A paper that details this work is currently under review. There are also a number of bug fixes and speed improvements.  

   

## Change log

**TSS-RESTREND v0.2.03 (Release date: 2018-01-23)**


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


### Examples

A full list of examples and extra documentation are currently under development 


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