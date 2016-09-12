
CHOW <- function(anu.VI, acu.RF, VI.index, breakpoints, sig=0.05, print=FALSE){
  #test the data to make sure its valid
  if (class(anu.VI) != "ts") 
    stop("anu.VI Not a time series object")
  if (class(acu.RF) != "ts") 
    stop("acu.VI Not a time series object")
  ti <- time(anu.VI)
  f <- frequency(anu.VI)
  #check the two ts object cover the same time period
  ti2 <- time(acu.RF)
  f2 <- frequency(acu.RF)
  if (!identical(ti, ti2))
    stop("ts object do not have the same time")
  if (!identical(f, f2))
    stop("ts object do not have the same frequency")
  # need to test the breakpoints to make sure they are numeric
  if (class(breakpoints) != "numeric") 
    stop("Breakpoints are not class numeric")
  
  #count of all the breakpoints
  len <- length(breakpoints)
  # convert the breakpoints into year posistion
  #the breakpoint loc will be the last anaual max before the breakpoint
  
  #empty variables to add to a datafram
  empty.1 <- NaN
  empty.2 <- NaN 
  empty.3 <- NaN
  ind.df <- data.frame(abs.index=breakpoints, yr.index = empty.1, reg.sig=empty.2, VPR.bpsig = empty.3)
  for (bp in 1:length(breakpoints)){
    bpv = ind.df$abs.index[bp]
    for (n in 1:length(VI.index)){
      if (bpv>=VI.index[n] & bpv<=VI.index[n+1]){
        #print(n)}
        ind.df$yr.index[bp] = n
      }
    }
  }
  #create the lm for the VPR and test is significance, 
  #split VPR sig from VPR insignificant 
  VPR.fit <- lm(anu.VI ~ acu.RF)
  if (summary(VPR.fit)$coefficients[,4][2] > sig){
    if (print){print("VPR significance below critical threshold, Testing breakpoints in the VPR")}
    ind <- acu.RF #independent variable
    dep <- anu.VI #dependent variable
    Method = "seg.VPR"
  } else {
    ind <- ti
    dep <- VPR.fit$residuals
    Method = "seg.RESTREND"
  } 
  #Iterate over each of the breakpoints
  while (TRUE){
    for (bp.num in nrow(ind.df)){ #the breakpoints number, first bp is 1,  etc
      bp = ind.df$yr.index[bp.num]
      #start and ends
      if (identical(ind.df$yr.index[bp.num-1], numeric(0))){
        bp.start = 1
        #print("here0")
      }else{
        bp.start = ind.df$yr.index[bp.num-1]
        #print("here1")
      }
      if (is.na(ind.df$yr.index[bp.num+1])){
        bp.end = length(dep)
      }else{
        bp.end = ind.df$yr.index[bp.num+1]}
      #perform the chow test
      bkp = bp - (bp.start-1)
      print(bkp)
      chow <- sctest(dep[bp.start:bp.end] ~ ind[bp.start:bp.end], type = "Chow", point = bkp)
      
      
      ind.df$reg.sig[bp.num] = chow$p.value
      
    }
    if (nrow(ind.df)>1){ #*****This needs to be tested with mutiple breakpoints*****
      #delete breakpoint with the largest p values (lowest significance)
      ind.df <- ind.df[!(1:nrow(ind.df) %in% (which.max(ind.df$reg.sig))),]
    }else if (nrow(ind.df)==1){
      if (print){
        print(chow)
      }
      VPR.chow<- sctest(anu.VI ~ acu.RF, type = "Chow", point = ind.df$yr.index[1])
      ind.df$VPR.bpsig[1] = VPR.chow$p.value
      
      if (Method == "seg.VPR" & ind.df$reg.sig[1] > sig){ # cant chow non-sig residulas (bpRESID.chow = FALSE)
        ind.df$reg.sig[1] = NaN
        return(structure(list(n.Method = FALSE, bp.summary = ind.df,
                              bpRESID.chow = FALSE, bpVPR.chow=VPR.chow), class = "CHOW.Object"))
      }else if (Method == "seg.VPR" & ind.df$reg.sig[1] <= sig){
        ind.df$reg.sig[1] = NaN
        return(structure(list(n.Method = "seg.VPR", bp.summary = ind.df,
                              bpRESID.chow = FALSE, bpVPR.chow=VPR.chow), class = "CHOW.Object"))
      }else if (Method == "seg.RESTREND" & ind.df$reg.sig[1] > sig){
        return(structure(list(n.Method = "RESTREND", bp.summary = ind.df,
                              bpRESID.chow = chow, bpVPR.chow=VPR.chow), class = "CHOW.Object"))
      }else if (Method == "seg.RESTREND" & ind.df$reg.sig[1] <= sig){
        if (VPR.chow$p.value>sig){
          return(structure(list(n.Method = "seg.RESTREND", bp.summary = ind.df,
                                bpRESID.chow = chow, bpVPR.chow=VPR.chow), class = "CHOW.Object"))
        }else if(VPR.chow$p.value<=sig){
          return(structure(list(n.Method = "seg.VPR", bp.summary = ind.df,
                                bpRESID.chow = chow, bpVPR.chow=VPR.chow), class = "CHOW.Object"))
        }
      }else{
        print("Error in Method and df, exit point 1")
        return(FALSE)
      }
    }else{
      print("ind.df shape is wrong, exited to avoid infinit loop, exit point 2")
      return(FALSE)
    }
  }
}