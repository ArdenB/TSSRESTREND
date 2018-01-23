TSS-RESTREND v0.2.03 (Release date: 2018-01-23)
==============

# Changes:
## Major: 
* Added the ability to include temperature as an additional climate variable. This involves changes to almost all parts of the package and new plotting functionality.   
* METHOD TWEAK: If the complete time series VPR/VCR relationship is not significant. TSS-RESTREND no longer returns indeterminate. It instead tests BFAST on the complete time series of Vegetation and uses those breakpoints.  This will cause minor differences in the results between earlier versions.

## Minor tweaks and bug fixes:
* Large improvements to the internal documentation 
* Speed improvements made throughout the functions  
* Add the ability to vary BFAST window (h) 
* ACPcalculator updated to deal with places where there are negative CTS VPR relationship. Small tweaks to TSS.RESTREND master script and the Annual.Precipitation functions
* BUG FIX: in the rainfall.accumulator now has a check that makes sure the rainfall data has some variance (SD >0). No variance (SD == 0) breaks lm() calculations and causes failures
* BUG FIX: Fixed a bug in the indexing of AnnualRF.Cal (now AnnualClim.Cal) and the rainfall.accumulator (now climate.accumulator) that was causing downstream problems with the segVPR function. This bug may have slightly impacted some results in segVPR 

