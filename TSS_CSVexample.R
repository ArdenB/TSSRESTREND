# Set the woking diroctory to the location
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# install the devtools package (contains the install_github)

# install.packages("devtools")

# library("devtools")

# install_github("ArdenB/SegmentedRESTREND_pub", subdir="TSS.RESTREND")
# library(TSS.RESTREND)

#Import the .csv that contains the monthly Vegetation data
# in.VI = read.csv("./demo_data/rabitVI.csv")
# in.VI = read.csv("./demo_data/mungoVI.csv")
# in.VI = read.csv("./demo_data/darlingVI.csv")
# in.VI = read.csv("./demo_data/desertVI.csv")
in.VI = read.csv("./demo_data/southwaVI.csv")

# turn that data into a time series object.
  #Starts january 1982 and ends december 2013 with a
  #monthly frequency (frequency = 12)
CTSR.VI = ts(in.VI, start=c(1982, 1), end=c(2013,12), frequency = 12)

# Import the associated rainfall series
  # Note that the rainfall sequences is longer to allow for the
  # accumulation period and the offset period. The start dates vary but the
  # end dates are the same.
# in.RF = read.csv("./demo_data/rabitRF.csv")
# in.RF = read.csv("./demo_data/mungoRF.csv")
# in.RF = read.csv("./demo_data/darlingRF.csv")
# in.RF = read.csv("./demo_data/desertRF.csv")
in.RF = read.csv("./demo_data/southwaRF.csv")

#turn the data into a time series object
rf.data = ts(in.RF, end=c(2013,12), frequency = 12)

# Define the max accumuulation period
max.acp <- 12

#Define the max offset period
max.osp <- 4

#Create a table of every possible precipitation value given the max.acp and max.osp
ACP.table <- rainfall.accumulator(CTSR.VI, rf.data, max.acp, max.osp)

# check the ACP.table is the right shape
print(dim(ACP.table))

#Pass the ACP.table and the CTSR.VI to the TSS.RESTREND
results <- TSS.RESTREND(CTSR.VI, ACP.table, print=TRUE, plot=TRUE)
print(results$summary)

