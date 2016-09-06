# SegmentedRESTREND

This script is a a demsonstration of the Segmented RESTREND method proposed in Burrell et. al., (2016?????). 

*Disclaimer. This code has been tested with the data provided int the `./demo_data` folder. Use caution when testing other dataset.*

## Dependencies 

This code was writen and tested on both linux (ubuntu 14.04) and Windows (Windows 10) in R 3.2.2. It depends on:
```R
library("bfast") # (1.5.7)
#library("RcppCNPy") #(0.2.5)
library("strucchange")#(1.5-1)
library("broom")#(0.4.1)
```
## Data Variables

All functions in this code use the same named data variables, 

*CTSR.VI - Complete Time Series Vegetation Index*     
An object of class "ts" that contains the monthly data set covering the entire timeseries being studied.  

*CTSR.RF - Complete Time Series Accumulated Rainfall*   
An object of class "ts" that contains the optimially accumulated precipitation for every value in CTSR.VI. CTSR.RF must be the same lenght, cover the same time period and have the same temporal frequency as CTSR.VI.   
Note. The CTSR.RF must be calculated by the user.  

*anu.VI - Annual(Growing Season) Max Vegetation Index value*  
An object of class "ts" that contains the Annual Max VI value from the CTSR.VI time series.  

*acu.RF - Accumulated Rainfall* 
An object of class "ts" that contains the accumulated precipitation for the anu.VI. This should be calculated seperatly from the CTSR.RF variable.  this can be calculated using the _______ function.  acu.RF must be the same lenght, cover the same time period and have the same temporal frequency as anu.VI.   

*VI.index - the index in CTSR.VI where the anu.VI values occur*
An object of class "ts". Note. Indexing in R starts from 1 rather than 0 like other languages.  

*rf.b4 - Accumulated Rainfall that is optimised for the period before a significant VPR breakpoint* 
An object of class "ts". rf.b4 must be the same lenght, cover the same time period and have the same temporal frequency as anu.VI.   

*rf.af - Accumulated Rainfall that is optimised for the period after a significant VPR breakpoint* 
An object of class "ts". rf.af must be the same lenght, cover the same time period and have the same temporal frequency as anu.VI.   

*ACCUM.TABLE - In Progress* 
to be added



## Functions




_Note. This is not the script used to produce the results in Burrell et. al., (2016). The code for that paper used both python and R, was designed specifically for the data sets used (can't be used with any dataset), and it was focused on parallelised computing. This code is much simplified and has not been speedtested._   
