# Time Series Segmented Residual Trends (TSS-RESTREND)

Time Series Segmented Residual Trends (TSS-RESTREND) is a method for the automated detection of land degradation from remotely sensed vegetation and climate datasets. TSS-RESTREND incorporates aspects of two existing degradation detection methods: RESTREND which is used to control for climate variability, and BFAST which is used to look for structural changes in the ecosystem.  

The full details of the testing and justification of the TSS-RESTREND method (version 0.1.02) are published in:
      
  Burrell, Arden L., Jason P. Evans, and Yi Liu. 2017. “Detecting Dryland Degradation Using Time Series Segmentation and Residual Trend Analysis (TSS-RESTREND).” Remote Sensing of Environment 197 (August):43–57. <https://doi.org/10.1016/j.rse.2017.05.018> 
  
The changes to the method included in version 0.2.03 focus on the inclusion of temperature as an additional climate variable. This allows for land degradation assessment in temperature limited drylands. A paper that details this work is currently under review. There are also a number of bug fixes and speed improvements.  

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

### Installing the package
#### Installing in R

```R
# Install from CRAN
install.packages("TSS.RESTREND")

# Install the version on Github
install.packages("devtools")
library(devtools)
install_github("ArdenB/TSS.RESTREND")

```
#### Installing from local package
``` bash
R CMD INSTALL TSS.RESTREND_0.2.03.tar.gz 
```
### Change log
[NEWS](NEWS.md)
### Examples

A full list of examples and extra documentation are currently under development 
