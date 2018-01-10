#' @title  Plot Function for ojects of the TSSRESTREND class
#'
#' @description Produces plots for class TSSRESTREND
#'
#' @import bfast
#' @import ggplot2
#' @importFrom  graphics axis grid lines mtext title
#' @importFrom stats predict
#'
#'
#' @param x
#'        Object of class TSSRESTREND
#' @param plots
#'        Defualts to "all", will produce the standard plots, plots can be called individually,
#'        "bfast", "chow", "VPR", "anu.VI", "final".
#' @param sig
#'        Significance
#' @param ...
#'        further arguments passed to the function.
#'
#' @export

plot.TSSRESTREND <- function(x, plots="all", sig=0.05, ...){
  #May add annual max VI vs RF
  #Plotting params
  if (class(x$ols.summary$chow.sum)=="logical"){
    if (plots=="all"||plots=="final"){
      stop("Incomplete TSSRESTREND object passed, This is a MISSING feature to be implemented in the next version")
      # print("Incomplete TSSRESTREND object passed. Producing the plots possible with available data.")
      # plots="final"
    }else{
      stop("Incomplete TSSRESTREND object passed, Can not produce the requested plot.")
    }
  }else{breakpoint <- x$ols.summary$chow.sum$yr.index}
  #pull out the annual max VI data
  anu.VI <- x$ts.data$anu.VI
  # Get the time and length
  len <- length(anu.VI)
  ti <- time(anu.VI)
  # Start plot. If statements allow for the creation of a single plot
  #BFAST
  #=========================================================================================================
  if (plots=="all"||plots=="bfast"){ #Plot the Bfast object
    plot(x$TSSRmodels$BFAST) #Uses defualt bfast plot
  }
  # Complete Time Series plot
  #=========================================================================================================
  if (plots=="all"||plots=="CTS"){
    # set the shape of the plot
    if (!is.null(x$ts.data$CTSR.TMraw)){
      #+++++ ADD THe Temperature data
      par(mar=c(5,4,4,6))
      plot(x$ts.data$CTSR.TMraw, axes=F, col="red2", lwd = 2, pch=16, xlab="", ylab="", lty=3)
      axis(side=4, col="red1", col.axis = 'red2')
      mtext("Temperature", side=4, line=2, col="red2")
      par(new=T)
      ofset = 2
    }else{
      par(mar=c(5,4,4,4))
      ofset = 0
    }
    #+++++ plot the Vegetation data
    plot(x$ts.data$CTSR.VI, col="olivedrab3", lwd = 2, pch=16, xlab="Year", ylab="",
         ylim=c((floor(range(x$ts.data$CTSR.VI)[1]*10)/10), (ceiling(range(x$ts.data$CTSR.VI)[2]*10)/10)))
    # Add a grid
    grid()
    # add Title and Grids
    mtext("VI",side=2,line=2,col="olivedrab3")
    title("Complete Time Series of the VI and Rainfall")
    #++++++ Plot the precipitation
    par(new=T)
    plot(x$ts.data$CTSR.RF, axes=F, col="royalblue2", lwd = 1.5, pch=16, xlab="", ylab="")
    axis(side=4, col="steelblue3", col.axis="royalblue2", line=(ofset+1))
    # add the axis title
    mtext("Accumulated Rainfall", side=4, line=ofset+3, col="royalblue2")
  }
  # Mulitbreak Chow plot
  #=========================================================================================================
  # Description
  # For pixels with multiple breakpoints detected in the VCR-residuals
  # this creates a plot of the residuals with all of the breakpoints.
  # To do:
  # Add additional comments so this is easier to read
  # Change the fit, maybe use the one saved in the result
  if (plots=="all"||plots=="chow"){
    # chowbp <- x$TSSRmodels$BFAST
    VI.index <- x$ts.data$VI.index
    #Check and see if there is something that needs to be plotted
    if (class(x$ols.summary$chow.ind) == "logical"){
      print("No Breakpoints to plot")
    }else if (dim(x$ols.summary$chow.ind)[1]<2){
      print("Insufficent (<2) Breakpoints to produce a chow plot")
    }else{
      #Set the segment colours
      colour <- c("orange", "purple", "olivedrab", "steelblue2", "lightslategrey", "hotpink1",
                  "lightgoldenrod3", "seagreen3", "sienna2")
      # Get a fit for the entire dataset
      browser() #This fit will not work if there is Temperature data
      fit <- lm(x$ts.data$anu.VI ~ x$ts.data$acu.RF)
      # Get the length of the data and the breakpoint indexes
      len <- length(x$ts.data$anu.VI)
      ind <- c(0, x$ols.summary$chow.ind[,"yr.index"], len)
      #Setup the X axis
      #Create a time series of dates for the x axis
      t = seq(start(ti)[1], end(ti)[1])
      xlim <- c(start(ti)[1], end(ti)[1])
      #Setup the values for the y axis
      m.range <- 1.5*max(abs(fit$residuals))
      ylim <- c(-m.range, m.range)
      # Get the fit of the residuals as a function of time
      R.fit <- lm(fit$residuals ~ t)
      #Start the plot
      plot(t, R.fit$fitted.values, xlab="Year", ylab="Residuals", main = "Chow Test on all Breakpoints",
           col="red", type = "l", lty="dashed", lwd=2, xlim=xlim, ylim = ylim)
      grid() # add grid
      # Get the number of breakpoints
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
  # Annual Max VI plot
  #=========================================================================================================
  # Description
  # Plots the Annual Maximum VI values with the Breakpoints marked
  # and colorcoded

  if (plots=="all"||plots=="anu.VI") {
    if (!breakpoint){
      #++++ Plot with no breakpoints
      # get the dates and build the time component
      yst <- start(anu.VI)[1]
      ynd <- end(anu.VI)[1]
      t = c(yst:ynd)
      plot(t, anu.VI, pch=16, xlab="Year", ylab="Annual max VI", col="orange")
      title("Annual VI max")
      grid()

    }else{
      #++++ Plot with a breakpoint
      # get the dates and build the time component
      yst <- start(anu.VI)[1]
      ynd <- end(anu.VI)[1]
      t = c(yst:ynd)

      # plot before the breakpoint
      plot(t[1:breakpoint], anu.VI[1:(breakpoint)], pch=16,xlab="Year",
           ylab="Annual max VI", col="orange", xlim=c(yst, ynd),
           ylim=c(min(anu.VI), max(anu.VI)) )
      grid()
      # plot after the breakpoint
      par(new=T)
      plot(t[(breakpoint+1):len], anu.VI[(breakpoint+1):len], pch=16,
           xlab="", ylab="", col="purple",main="", xlim=c(yst, ynd),
           ylim=c(min(anu.VI), max(anu.VI)))
      title("Annual VI max")
      # Add a line at the breakpoint
      abline(v=(breakpoint-0.5+yst), col="red", lty = "dashed")
    }
  }

  # Plot the VCR
  #=========================================================================================================
  # TO DO:
    # Add a Description
    # Change titles to VCR
    # Add in temperpature plot
    # ADd a dual fit
    # Add a temperature plot as well
    # ?? Add a 3D plot
  if (plots=="all"||plots=="VPR"){
    if (x$summary$Method=="segmented.VPR"){
      #++++ Plot a VPR/VCRwith a breakpoint (SEMGEMTED VPR)
      browser()
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
      arrows(0, bot, x1=0, y1=top, length =.075,  angle = 90, code=3, col="black", lwd=2)
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
    }else{ # No breakpoin the the VPR or VCR
      # Sort and see if i need to plot VPR or VCR (Is temperature a variable)
      if (is.null(x$ts.data$acu.TM)) { # no accumulated temperature data (may not be sig)
        # Plot the VPR
        plot(as.numeric(x$ts.data$anu.VI) ~ as.numeric(x$ts.data$acu.RF), pch=16, xlab="Accumulated Rainfall (mm)",
             ylab="Annual VImax",  col="orange")
        # Add a fit line for the VPR
        abline(x$TSSRmodels$VPR.fit, col = "red",lwd = 2, lty = "dashed")
        #Add title and Grid
        title("VPR fit")
        grid()
        #Get the Rsquared and pvalues and add them to the plot
        R.pval = glance(x$TSSRmodels$VPR.fit)$p.value
        R.Rval = summary(x$TSSRmodels$VPR.fit)$r.squared
        rp = vector('expression', 2)
        rp[1] = substitute(expression(italic(p) == R.pval),
                           list(R.pval  = format(R.pval, digits = 3)))[2]
        rp[2] = substitute(expression(italic(R^2) == R.Rval),
                           list(R.Rval = format(R.Rval,dig=3)))[2]
        # Set the location of the Rsquared an p values
        legend('topleft', legend = rp, bty = 'n')

      }else{ #Has temperature
        # Scatter plot with precip and Temperature as the x an y variables.
        # The VI values are expressed by the color of the points

        # Color points by the mpg variable
        # Put the relevant variables in a dataframe
        #Try a 3d version

        dfX <- data.frame(Precip=as.numeric(x$ts.data$acu.RF),
                          Temp=as.numeric(x$ts.data$acu.TM)/x$acum.df$tacp,
                          VI=as.numeric(x$ts.data$anu.VI))
        # Check if rgl is installed and openable
        if ((requireNamespace("rgl", quietly = TRUE)) && requireNamespace("car", quietly = TRUE)) {
          browser()
          # Build a 3d plot
          car::scatter3d(
            x=dfX$Precip,  y=dfX$VI, z=dfX$Temp, xlab="Accumulated Precipitation",
            ylab="Annual Max VI" , zlab = "Mean Monthly Accumulated Temperature",
            point.col="orange", sphere.size=0.5, surface.col="orange", axis.ticks=TRUE,
            neg.res.col="grey", pos.res.col="grey",
            axis.col=c("royalblue2", "olivedrab3", "red2")
            )
          #The
          browser()
          # pkg::f()
        } else {
          browser()
        }

        # the 2d version is a bit terrible
        # sp1<-ggplot(mtcars, aes(x=wt, y=mpg, color=mpg)) + geom_point()
        sp1<-ggplot(dfX, aes(x=Precip, y=Temp, color=VI)) + geom_point()
        sp1
        # Gradient between n colors
        sp1+scale_color_gradientn(colours = rainbow(5))
      }
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
      plot(c(start(ti)[1]:end(ti)[1]), x$TSSRmodels$VPR.fit$residuals, pch=16,xlab="Year",
           ylab="Residuals", col="orange",main="RESTREND", ylim = c(-m.range, m.range))
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






