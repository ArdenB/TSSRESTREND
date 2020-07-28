# Combining paired C3 and C4 runs 

This section describes a methods to combine TSS-RESTREND attribution runs done assuming a C3 and a C4 photosyenthetic pathway. The TSSRattribution function in Version 0.3.1 of the TSS.RESTREND has the optional argument C4frac which removes the need for this code. It is kept here for reference purposes only.  


### Build two runs

To build the two matched runs us the following:

```bash
# C3 run
python S00_SetupMetadata.py --coarsen 5  --photo "C4"
python S01_processingnetcdf.py
Rscript S02_TSSRESTRENDattribution.R  --ncores -1
python S03_MappingResults.py

# C4 run
python S00_SetupMetadata.py --coarsen 5  --photo "C3"
python S01_processingnetcdf.py
Rscript S02_TSSRESTRENDattribution.R  --ncores -1
python S03_MappingResults.py

```

After this code is complete, there will be two files in the results folder called "TSSRattribution_Results_C3_125km.nc" and "TSSRattribution_Results_C4_125km.nc".  