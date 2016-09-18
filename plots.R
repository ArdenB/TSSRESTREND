  # Plots 
  #Print and plot
  if (print) { #Need better explination and a title but it will do for the moment
    print(summary(CTS.fit))
  }
  if (plot){
    #This plot need serius work but will do for the moment
    plot(resid.ts)
  }


  if (plot){
    #if plot is requested
    plot(bf.fit)
  }


  if (!bp){# no breakpoints detected by the BFAST
    test.Method = "RESTREND"
    if (plot){
      VImax.plot(anu.VI)
    }


        # if (plot){
    #   VImax.plot(anu.VI, brkp=brkp)
    # }