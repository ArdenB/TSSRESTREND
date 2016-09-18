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



  if (print){
    print(summary(bpanalysis))
  }

  if (plot){
    len <- length(VPR.resid)
    R.fit <- lm(VPR.residuals ~ year, segRES.df)
    R.fit0 <- lm(VPR.residuals[1:breakpoint] ~ year[1:breakpoint], segRES.df)
    R.fit1 <- lm(VPR.residuals[breakpoint+1:len] ~ year[breakpoint+1:len], segRES.df)

    start = as.integer(start(ti)[1])
    end = as.integer(end(ti)[1])
    xlim = c(start, end)
    m.range = 2*max(abs(bpanalysis$fitted.values))

    RESchow <- sctest(segRES.df$VPR.residuals ~ segRES.df$year, type = "Chow", point = breakpoint)
    plot(segRES.df$year[1:breakpoint], segRES.df$VPR.residuals[1:breakpoint], pch=16,
         xlab="time", ylab="Residuals", col="orange", xlim=xlim, ylim = c(-m.range, m.range))
    par(new=T)
    plot(segRES.df$year[breakpoint+1:len], segRES.df$VPR.residuals[breakpoint+1:len], pch=16,
         xlab="", ylab="", col="purple", main="",
         xlim=c(start, end), ylim = c(-m.range, m.range))
    abline(R.fit, col = "red",lwd = 2, lty = "dashed")


    par(new=T)
    bpa.fitts <- ts(bpanalysis$fitted.values, start=ti[1], end=tail(ti, 1), frequency = f)
    b4.bp = bpanalysis$coefficients[[1]]
    af.bp = bpanalysis$coefficients[[1]] + bpanalysis$coefficients[[3]]
    bpats2 <- append(bpa.fitts, c(b4.bp, af.bp), after = breakpoint)
    t2 <- append(ti, c(start+breakpoint-0.50001, start+breakpoint-0.49999), after=breakpoint)
    plot(t2, bpats2, pch=16, type = "l", lwd = 2,
         xlab="", ylab="", col="red", main="",
         xlim=c(start, end), ylim = c(-m.range, m.range))
    #add a breakpoint band
    abline(v=(breakpoint-0.5+start), col="white", lwd = 3, lty = "dotted")
    #Need to change the stastics that is shows here

    R.Fval = summary(bpanalysis)$f[[1]]
    R.pval = glance(bpanalysis)$p.value
    R.Rval = summary(bpanalysis)$r.squared
    rp = vector('expression',3)
    rp[1] = substitute(expression(italic(R^2) == R.Rval),
                       list(R.Rval = format(R.Rval,dig=3)))[2]
    rp[2] = substitute(expression(italic(F) == R.Fval),
                       list(R.Fval  = format(R.Fval,dig=3)))[2]

    rp[3] = substitute(expression(italic(p) == R.pval),
                       list(R.pval  = format(R.pval, digits = 3)))[2]
    legend('topleft', legend = rp, bty = 'n')
  }
