# Time Series Segmented Residual Trends (TSS-RESTREND)

Time Series Segmented Residual Trends (TSS-RESTREND) is a method for the automated detection of land degradation from remotely sensed vegetation and climate datasets. TSS-RESTREND incorporates aspects of two existing degradation detection methods: RESTREND which is used to control for climate variability, and BFAST which is used to look for structural changes in the ecosystem.  

The full details of the testing and justification of the TSS-RESTREND method (version 0.1.02) are published in:
      
  Burrell, Arden L., Jason P. Evans, and Yi Liu. 2017. “Detecting Dryland Degradation Using Time Series Segmentation and Residual Trend Analysis (TSS-RESTREND).” Remote Sensing of Environment 197 (August):43–57. \url{https://doi.org/10.1016/j.rse.2017.05.018.} (\url{http://www.sciencedirect.com/science/article/pii/S0034425717302171})

The changes to the method included in version 0.2.03 focus on the inclusion of temperature as an additional climate variable to allow for land degradation assessment in temperature limited drylands as well and a number of bug fixes and speed improvements. This work is currently under review. 

## To install from Github
```R
install.packages("devtools")
library(devtools)
install_github("ArdenB/TSS.RESTREND")

```
## to install from local package
...
R CMD INSTALL TSS.RESTREND_0.1.031.tar.gz 
...

## Dependencies 

This code was written and tested on both Linux (ubuntu 14.04) and Windows (Windows 10) in R 3.2.2. It depends on:
```R
library("bfast") # (1.5.7)
#library("RcppCNPy") #(0.2.5)
library("strucchange")#(1.5-1)
library("broom")#(0.4.1)
```
## Data Variables

All functions in this code use the same named data variables, 

*CTSR.VI - Complete Time Series Vegetation Index*     
An object of class "ts" that contains the monthly data set covering the entire time-series being studied.  

*CTSR.RF - Complete Time Series Accumulated Rainfall*   
An object of class "ts" that contains the optimally accumulated precipitation for every value in CTSR.VI. CTSR.RF must be the same length, cover the same time period and have the same temporal frequency as CTSR.VI.   
Note. The CTSR.RF must be calculated by the user.  

*anu.VI - Annual(Growing Season) Max Vegetation Index value*    
An object of class "ts" that contains the Annual Max VI value from the CTSR.VI time series.  

*acu.RF - Accumulated Rainfall*     
An object of class "ts" that contains the accumulated precipitation for the anu.VI. This should be calculated separately from the CTSR.RF variable.  this can be calculated using the _______ function.  acu.RF must be the same length, cover the same time period and have the same temporal frequency as anu.VI.   

*VI.index - the index in CTSR.VI where the anu.VI values occur*   
An object of class "ts". Note. Indexing in R starts from 1 rather than 0 like other languages.  

*rf.b4 - Accumulated Rainfall that is optimised for the period before a significant VPR breakpoint*     
An object of class "ts". rf.b4 must be the same length, cover the same time period and have the same temporal frequency as anu.VI.   

*rf.af - Accumulated Rainfall that is optimised for the period after a significant VPR breakpoint*    
An object of class "ts". rf.af must be the same length, cover the same time period and have the same temporal frequency as anu.VI.   

*ACCUM.TABLE - In Progress*     
to be added


## Functions

The functions in the code can be broken up into three types, the main function, the step functions and the demonstration functions.

###Main Function

*TSS.RESTREND*
```R
TSS.RESTREND(CTSR.VI, CTSR.RF, anu.VI, acu.RF, VI.index, rf.b4=FALSE, rf.af=FALSE, 
             sig=0.05, print=FALSE, plot=FALSE, details=FALSE)
```
######Parameters: 

######Returns:  

######Notes 



_Note. This is not the script used to produce the results in Burrell et. al., (2016). The code for that paper used both python and R, was designed specifically for the data sets used (can't be used with any dataset), and it was focused on parallelised computing. This code is much simplified and has not been speed-tested._   
