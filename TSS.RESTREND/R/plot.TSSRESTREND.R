#' @title  Plot function for class TSSRESTREND class
#'
#' @description Plot *******
#'
#' @import  bfast
#'
#' @param plots = "all"
#' User can create a single plot,
#'
#' @export

plot.TSSRESTREND <- function(x, verbose=FALSE, plots="all", sig=0.05, ...){
  #May add annual max VI vs RF
  #Plotting params
  breakpoint <- x$ols.summary$chow.sum$yr.index
  anu.VI <- x$ts.data$anu.VI
  len <- length(anu.VI)
  ti <- time(anu.VI)
  #Add a plot for VERBOSE (ctsr.vi vs VTSR.RF)
  if (plots=="all"||plots=="bfast"){
    plot(x$TSSRmodels$BFAST)
  # browser()
  if (plots=="all"||plots=="CTS"){
    par(mar=c(5,4,4,4))
    plot(x$ts.data$CTSR.VI, col="olivedrab", lwd = 2, pch=16, xlab="Year", ylab="")
    grid()
    mtext("VI",side=2,line=2,col="olivedrab")
    title("Complete Time Series of the VI and Rainfall")
    par(new=T)
    plot(x$ts.data$CTSR.RF, axes=F, col="steelblue2", lwd = 1, pch=16, xlab="", ylab="")
    axis(side=4)
    mtext("Accumulated Rainfall", side=4, line=2, col="steelblue2")
  }
  }

  if (verbose==TRUE||plots=="chow"){
    # chowbp <- x$TSSRmodels$BFAST
    VI.index <- x$ts.data$VI.index

    if (class(x$ols.summary$chow.ind) == "logical"){
      print("No Breakpoints to plot")
    }else if (dim(x$ols.summary$chow.ind)[1]<2){
      print("Insufficent Breakpoints to produce a chow plot")
    }else{
      colour <- c("orange", "purple", "olivedrab", "steelblue2", "lightslategrey", "hotpink1",
                  "lightgoldenrod3", "seagreen3", "sienna2")
      fit <- lm(x$ts.data$anu.VI ~ x$ts.data$acu.RF)
      len <- length(x$ts.data$anu.VI)
      ind <- c(0, x$ols.summary$chow.ind[,"yr.index"], len)
      t = seq(start(ti)[1], end(ti)[1])
      xlim <- c(start(ti)[1], end(ti)[1])
      m.range <- 1.5*max(abs(fit$residuals))
      ylim <- c(-m.range, m.range)
      R.fit <- lm(fit$residuals ~ t)
      plot(t, R.fit$fitted.values, xlab="Year", ylab="Residuals", main = "Chow Test on all Breakpoints",
           col="red", type = "l", lty="dashed", lwd=2, xlim=xlim, ylim = ylim)
      grid()
      bpn <- length(x$ols.summary$chow.ind[,"yr.index"])
      rp = vector('expression',bpn)
      for (n in 2:length(ind)){
        col.ind <- n-1
        st <- ind[n-1]
        en <- ind[n]
        par(new=T)
        df <- data.frame(z=t[(st+1):en], resi = fit$residuals[(st+1):en])
        r.sfit <- lm(df$resi ~ df$z)
        plot(df$z, df$resi, pch=16,
             xlab="", ylab="", col=colour[col.ind], main="",
             xlim=xlim, ylim = ylim)
        # par(new=T)
        lines(df$z, predict(r.sfit, df), col = colour[col.ind])


        if (n<length(ind)){
          abline(v=(ind[n]-0.5+start(ti)[1]), col="red", lty="dotted")
          RESchow <- sctest(fit$residuals[(st+1):ind[(n+1)]] ~ t[(st+1):ind[(n+1)]],
                              type = "Chow", point = (ind[n]-ind[(n-1)]))
          R.pval = RESchow$p.value
          # browser()
          rp[col.ind] = substitute(expression(italic(p) == R.pval),
                                   list(R.pval = format(R.pval, digits = 5)))[2]
        }

      }
      legend('topleft', legend = rp, bty = 'n')
      bkp.z <- as.numeric(x$ols.summary$chow.sum["yr.index"])
      RESchow.f <- sctest(fit$residuals ~ t, type = "Chow", point = bkp.z)
      plot(t, R.fit$fitted.values, xlab="Year", ylab="Residuals", main = "Chow Test on most Sig. Breakpoints",
           col="red", type = "l", lty="dashed", lwd=2, xlim=xlim, ylim = ylim)
      grid()
      par(new=T)
      plot(t[1:bkp.z], fit$residuals[1:bkp.z], pch=16,
           xlab="", ylab="", main="", col="orange",
           xlim=xlim, ylim = ylim)#main="RESTREND",
      new <- data.frame(b = t[1:bkp.z], 1)
      R.fit0 <- lm(fit$residuals[1:bkp.z] ~ t[1:bkp.z])
      lines(new$b, predict(R.fit0, new), col = "orange")
      par(new=T)
      plot(t[bkp.z+1:len], fit$residuals[bkp.z+1:len], pch=16,
           xlab="", ylab="", col="purple", main="",
           xlim=xlim, ylim = ylim)
      R.fit1 <- lm(fit$residuals[(bkp.z+1):len] ~ t[(bkp.z+1):len])

      new2 <- data.frame(c = t[(bkp.z+1):len])
      # pre <-
      lines(new2$c, predict(R.fit1, new2), col = "purple")
      abline(v=(bkp.z-0.5+start(ti)[1]), col="red", lty="dotted")
      r.Fval = RESchow.f$statistic
      r.pval = RESchow.f$p.value
      rp = vector('expression',2)
      rp[1] = substitute(expression(italic(F) == r.Fval),
                         list(r.Fval = format(r.Fval,dig=5)))[2]
      rp[2] = substitute(expression(italic(p) == r.pval),
                         list(r.pval = format(r.pval, digits = 5)))[2]
      legend('topleft', legend = rp, bty = 'n')
      # browser()

    }
  }


  if (plots=="all"||plots=="VPR"){
    if (x$summary$Method=="segmented.VPR"){
      StdVar.RF<-x$ts.data$StdVar.RF
      fit0 <- lm(anu.VI[1:breakpoint] ~ StdVar.RF[1:breakpoint])
      fit1 <- lm(anu.VI[(breakpoint+1):len] ~ StdVar.RF[(breakpoint+1):len])
      fitRES <- lm(anu.VI ~ StdVar.RF)
      # chow1 <- sctest(bpanalysis$residuals ~ t, type = "Chow", point = bkp)
      plt.ymin <- min(anu.VI)
      plt.ymax <- max(anu.VI)
      plt.xmin <- min(StdVar.RF)
      plt.xmax <- max(StdVar.RF)
      plot(StdVar.RF[1:breakpoint], anu.VI[1:breakpoint], pch=16,
           xlab="Rainfall Standard Variance", ylab="Annual max VI", col="orange",
           xlim=c(plt.xmin, plt.xmax), ylim=c(plt.ymin, plt.ymax))
      title("Segmented VPR")
      grid()
      par(new=T)
      plot(StdVar.RF[(breakpoint+1):len], anu.VI[(breakpoint+1):len], pch=16,
           xlab="", ylab="", col="purple", main="",
           xlim=c(plt.xmin, plt.xmax), ylim=c(plt.ymin, plt.ymax))
      par(new=T)
      # abline(fit, col = "red",lwd = 2, lty = "dashed")
      abline(fit0, col = "orange", lwd = 2)
      abline(fit1, col = "purple", lwd = 2)
      abline(fitRES, col="darkgrey", lwd=2, lty="dashed")
      top <- x$TSSRmodels$segVPR.fit$coefficients[[1]]
      bh <-  x$TSSRmodels$segVPR.fit$coefficients[[3]]
      bot <- top+bh
      # lines(x=c(0, 0), y=c(top, bot), col="red", lwd=2, pch=0)
      arrows(0, bot, x1=0, y1=top, length =.075,  angle = 90, code=3, col="red", lwd=2)
      R.Fval = summary(x$TSSRmodels$segVPR.fit)$f[[1]]
      R.pval = glance(x$TSSRmodels$segVPR.fit)$p.value
      R.Rval = summary(x$TSSRmodels$segVPR.fit)$r.squared
      rp = vector('expression', 4)
      rp[1] = substitute(expression(italic(BH) == bh),
                         list(bh  = format(bh, digits = 3)))[2]
      rp[2] = substitute(expression(italic(F) == R.Fval),
                         list(R.Fval  = format(R.Fval,dig=3)))[2]
      rp[3] = substitute(expression(italic(p) == R.pval),
                         list(R.pval  = format(R.pval, digits = 3)))[2]
      rp[4] = substitute(expression(italic(R^2) == R.Rval),
                         list(R.Rval = format(R.Rval,dig=3)))[2]

      legend('topleft', legend = rp, bty = 'n')
    }else{
      plot(as.numeric(x$ts.data$anu.VI) ~ as.numeric(x$ts.data$acu.RF), pch=16, xlab="Accumulated Rainfall (mm)",
           ylab="Annual VImax",  col="orange")
      title("VPR fit")
      abline(x$TSSRmodels$VPR.fit, col = "red",lwd = 2, lty = "dashed")
      grid()
      #Maybe add p values or r2 or something
    }
  }


  if (plots=="all"||plots=="anu.VI") {
    if (!breakpoint){
      yst <- start(anu.VI)[1]
      ynd <- end(anu.VI)[1]
      t = c(yst:ynd)
      plot(t, anu.VI, pch=16, xlab="Year", ylab="Annual max VI", col="orange")
      title("Annual VI max")
      grid()

    }else{
      yst <- start(anu.VI)[1]
      ynd <- end(anu.VI)[1]

      t = c(yst:ynd)
      plot(t[1:breakpoint], anu.VI[1:(breakpoint)], pch=16,xlab="Year",
           ylab="Annual max VI", col="orange", xlim=c(yst, ynd),
           ylim=c(min(anu.VI), max(anu.VI)) )
      grid()
      par(new=T)
      plot(t[(breakpoint+1):len], anu.VI[(breakpoint+1):len], pch=16,
           xlab="", ylab="", col="purple",main="", xlim=c(yst, ynd),
           ylim=c(min(anu.VI), max(anu.VI)))
      title("Annual VI max")
      abline(v=(breakpoint-0.5+yst), col="red", lty = "dashed")
    }
  }
  if(plots=="all"||plots=="final") {
    if (x$summary$Method=="segmented.RESTREND"){
      VPR.residuals <- x$TSSRmodels$VPR.fit$residuals
      len <- length(VPR.residuals)
      start = as.integer(start(ti)[1])
      end = as.integer(end(ti)[1])
      year = c(start:end)
      R.fit <- lm(VPR.residuals ~ year)
      R.fit0 <- lm(VPR.residuals[1:breakpoint] ~ year[1:breakpoint])
      R.fit1 <- lm(VPR.residuals[breakpoint+1:len] ~ year[breakpoint+1:len])


      xlim = c(start, end)
      m.range = 2*max(abs(x$TSSRmodels$resid.fit$fitted.values))

      RESchow <- sctest(VPR.residuals ~ year, type = "Chow", point = breakpoint)
      plot(year[1:breakpoint], VPR.residuals[1:breakpoint], pch=16,
           xlab="time", ylab="Residuals", col="orange", xlim=xlim, ylim = c(-m.range, m.range))
      title("Segmented RESTREND")
      par(new=T)
      plot(year[breakpoint+1:len], VPR.residuals[breakpoint+1:len], pch=16,
           xlab="", ylab="", col="purple", main="",
           xlim=c(start, end), ylim = c(-m.range, m.range))
      abline(R.fit, col = "darkgrey",lwd = 2, lty = "dashed")


      par(new=T)
      bpa.fitts <- ts(x$TSSRmodels$resid.fit$fitted.values, start=ti[1], end=tail(ti, 1), frequency = 1)
      b4.bp = x$TSSRmodels$resid.fit$coefficients[[1]]
      af.bp = x$TSSRmodels$resid.fit$coefficients[[1]] + x$TSSRmodels$resid.fit$coefficients[[3]]
      bpats2 <- append(bpa.fitts, c(b4.bp, af.bp), after = breakpoint)
      t2 <- append(ti, c(start+breakpoint-0.50001, start+breakpoint-0.49999), after=breakpoint)
      plot(t2, bpats2, pch=16, type = "l", lty = "dashed", lwd = 2,
           xlab="", ylab="", col="red", main="",
           xlim=c(start, end), ylim = c(-m.range, m.range))
      grid()
      #add a breakpoint band
      abline(v=(breakpoint-0.5+start), col="white", lwd = 3)#, lty = "dotted")
      #Need to change the stastics that is shows here

      R.Fval = summary(x$TSSRmodels$resid.fit)$f[[1]]
      R.pval = glance(x$TSSRmodels$resid.fit)$p.value
      R.Rval = summary(x$TSSRmodels$resid.fit)$r.squared
      # rp = vector('expression',3)

      top <- x$TSSRmodels$resid.fit$fitted.values[len]
      bot <- x$TSSRmodels$resid.fit$fitted.values[1]
      r.c <- top-bot
      arrows((end+0.5), bot, x1=(end+0.5), y1=top, length =.075,  angle = 90, code=3, col="red", lwd=2)

      rp = vector('expression', 3)
      rp[1] = substitute(expression(italic(rc) == r.c),
                         list(r.c  = format(r.c, digits = 3)))[2]
      rp[2] = substitute(expression(italic(R^2) == R.Rval),
                         list(R.Rval = format(R.Rval,dig=3)))[2]
      # rp[2] = substitute(expression(italic(F) == R.Fval),
                         # list(R.Fval  = format(R.Fval,dig=3)))[2]

      rp[3] = substitute(expression(italic(p) == R.pval),
                         list(R.pval  = format(R.pval, digits = 3)))[2]
      legend('topleft', legend = rp, bty = 'n')
    }else if(x$summary$Method=="RESTREND"){
      start = as.integer(start(ti)[1])
      end = as.integer(end(ti)[1])
      RES <- x$TSSRmodels$resid.fit
      m.range = 2*max(abs(RES$fitted.values))
      plot(c(start(ti)[1]:end(ti)[1]), x$TSSRmodels$VPR.fit$residuals, pch=16,xlab="Accumulated Rainfall",
           ylab="Annual VI", col="orange",main="RESTREND", ylim = c(-m.range, m.range))
      par(new=T)
      plot(c(start(ti)[1]:end(ti)[1]), RES$fitted.values, type = "l", lwd = 2, lty = "dashed", pch=16,
           xlab="", ylab="", col="red", main="", ylim = c(-m.range, m.range))
      grid()
      R.Fval = summary(RES)$f[[1]]
      R.Rval = summary(RES)$r.squared
      R.pval = glance(RES)$p.value
      # # lines(x=c(0, 0), y=c(top, bot), col="red", lwd=2, pch=0)
      top <- x$TSSRmodels$resid.fit$fitted.values[len]
      bot <- x$TSSRmodels$resid.fit$fitted.values[1]
      r.c <- top-bot
      arrows((end+0.5), bot, x1=(end+0.5), y1=top, length =.075,  angle = 90, code=3, col="red", lwd=2)

      rp = vector('expression', 3)
      rp[1] = substitute(expression(italic(rc) == r.c),
                         list(r.c  = format(r.c, digits = 3)))[2]
      rp[2] = substitute(expression(italic(R^2) == R.Rval),
                         list(R.Rval = format(R.Rval,dig=3)))[2]
      rp[3] = substitute(expression(italic(p) == R.pval),
                         list(R.pval = format(R.pval, digits = 3)))[2]
      # rp[3] = substitute(expression(italic(F) == R.Fval),
      #                    list(R.Fval = format(R.Fval,dig=3)))[2]

      legend('topleft', legend = rp, bty = 'n')
    }else if(x$summary$Method=="segmented.VPR"){
      # browser()
      # resid.raw<- x$TSSRmodels$segVPR.fit$residuals
      R2.BH <- x$summary$VPR.HeightChange
      # c(resid.raw[1:breakpoint], resid.raw[(breakpoint+1):len] + R2.BH)
      VPR.residuals <- x$TSSRmodels$segVPR.fit$residuals
      len <- length(VPR.residuals)
      start = as.integer(start(ti)[1])
      end = as.integer(end(ti)[1])
      year = c(start:end)
      RESchow <- sctest(VPR.residuals ~ year, type = "Chow", point = breakpoint)




      R.fit <- lm(VPR.residuals ~ year)
      R.fit0 <- lm(VPR.residuals[1:breakpoint] ~ year[1:breakpoint])
      R.fit1 <- lm(VPR.residuals[breakpoint+1:len] ~ year[breakpoint+1:len])


      xlim = c(start, end)
      m.range = 2*max(abs(x$TSSRmodels$resid.fit$fitted.values))


      plot(year[1:breakpoint], VPR.residuals[1:breakpoint], pch=16,
           xlab="time", ylab="Residuals", col="orange", xlim=xlim, ylim = c(-m.range, m.range))
      grid()
      title("Segmented VPR RESTREND")
      par(new=T)
      plot(year[breakpoint+1:len], VPR.residuals[breakpoint+1:len], pch=16,
           xlab="", ylab="", col="purple", main="",
           xlim=c(start, end), ylim = c(-m.range, m.range))
      # abline(R.fit, col = "red",lwd = 2, lty = "dashed")


      par(new=T)
      bpa.fitts <- ts(x$TSSRmodels$resid.fit$fitted.values, start=ti[1], end=tail(ti, 1), frequency = 1)
      b4.bp = x$TSSRmodels$resid.fit$coefficients[[1]]
      af.bp = x$TSSRmodels$resid.fit$coefficients[[1]] + x$TSSRmodels$resid.fit$coefficients[[3]]
      bpats2 <- append(bpa.fitts, c(b4.bp, af.bp), after = breakpoint)
      t2 <- append(ti, c(start+breakpoint-0.50001, start+breakpoint-0.49999), after=breakpoint)
      plot(t2, bpats2, pch=16, type = "l", lwd = 2,
           xlab="", ylab="", col="red", main="", lty = "dashed",
           xlim=c(start, end), ylim = c(-m.range, m.range))
      #add a breakpoint band
      abline(v=(breakpoint-0.5+start), col="white", lwd = 3)#, lty = "dotted")

      top <- x$TSSRmodels$resid.fit$fitted.values[len]
      bot <- x$TSSRmodels$resid.fit$fitted.values[1]
      r.c <- top-bot
      # # lines(x=c(0, 0), y=c(top, bot), col="red", lwd=2, pch=0)
      arrows((end+0.5), bot, x1=(end+0.5), y1=top, length =.075,  angle = 90, code=3, col="red", lwd=2)
      #Need to change the stastics that is shows here

      R.Fval = summary(x$TSSRmodels$resid.fit)$f[[1]]
      R.pval = glance(x$TSSRmodels$resid.fit)$p.value
      R.Rval = summary(x$TSSRmodels$resid.fit)$r.squared
      rp = vector('expression',4)
      rp[1] = substitute(expression(italic(rc) == r.c),
                         list(r.c  = format(r.c, digits = 3)))[2]
      rp[2] = substitute(expression(italic(R^2) == R.Rval),
                         list(R.Rval = format(R.Rval,dig=3)))[2]
      rp[3] = substitute(expression(italic(F) == R.Fval),
                         list(R.Fval  = format(R.Fval,dig=3)))[2]

      rp[4] = substitute(expression(italic(p) == R.pval),
                         list(R.pval  = format(R.pval, digits = 3)))[2]
      legend('topleft', legend = rp, bty = 'n')


    }
  }
}






