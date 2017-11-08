# Segmented RESTREND

This script is a a demonstration of the Segmented RESTREND method proposed in Burrell et. al., (2017). 

The current build (0.1.4) is under active development and has not been fully tested.  Version 0.1.32 is stable but doesn't handle some errors well.  

*Disclaimer. This code has been tested with the data provided int the `./demo_data` folder. Use caution when testing other dataset.*
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
