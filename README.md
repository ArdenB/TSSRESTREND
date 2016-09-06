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
An object of class "ts". Monthly data set covering the entire timeseries being studied.  

*CTSR.RF - Complete Time Series Accumulated Rainfall*
An object of class "ts" that contains the optimially accumulated precipitation for every value in CTSR.VI




## Functions




_Note. This is not the script used to produce the results in Burrell et. al., (2016). The code for that paper used both python and R, was designed specifically for the data sets used (can't be used with any dataset), and it was focused on parallelised computing. This code is much simplified and has not been speedtested._   
