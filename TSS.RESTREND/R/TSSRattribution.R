#' @title Vegetation change attribution using the Time Series Segmentation of Residual Trends (MAIN FUNCTION)
#'
#' @importFrom stats coef end frequency lm sd start time ts cor.test predict pchisq p.adjust
#' @importFrom mblm mblm
#' @importFrom utils tail read.csv
#' @importFrom broom glance
#' @importFrom RcppRoll roll_mean
#'
#' @description
#' This is a wrapper function around the TSS-RESTREND main function that dows additional attribution.
#' It measues the Observed vegetation change, land use, climate change and climate varibility.
#' Unlike TSSRESTREND function, this requires both temperature and precitation data fo work.
#' @author Arden Burrell, aburrell@whrc.org
#'
#' @inheritParams TSSRESTREND
#' @inheritParams climate.accumulator
#' @inheritParams franksCO2
#' @param C4frac
#'        The fraction of vegetation that follows the C4 photosynthetic pathway, between 0 and 1
#' @param returnmodels
#'        Return all the created models as well as the original data
#' @param AnnualRes
#'        Report results in change per year. Defualt is False. Instead reports total change from the
#'        start to the end of the time series.
#' @param SkipError
#'        Bool, If TRUE will handle most errors and return a dataframe filled with NA's.
#'        Usefull when processing large datasets to stop analysis failing on a single pixel.
#'        Use with caution. Defualt=TRUE
#' @param splitclim
#'        Bool, If TRUE Climate will be split into climate change and climate varibility as
#'        per Burrell et al., (2020). If FALSE. will just return climate as per IPCC:
#'        Chapter 3: Desertification in the IPCC Special Report on Climate Change, Desertification,
#'        Land Degradation, Sustainable Land Management, Food Security, and Greenhouse gas fluxes in Terrestrial Ecosystems.
#'        Defualt=True.
#' @return \bold{tacp}
#'        The accumulation period for the annual max time series temperature


#' @export

TSSRattribution <- function(
	CTSR.VI, CTSR.RF, CTSR.TM, max.acp, max.osp, C4frac = 0,  sig = 0.05, season = "none", exclude = 0,
	allow.negative = FALSE, allowneg.retest = FALSE, h = 0.15, returnmodels = FALSE, AnnualRes=FALSE,
  SkipError=TRUE, retnonsig=TRUE, splitclim=TRUE, cliwindow=20, CO2=FALSE, refyear=1980){

  # ========== Work out the Change unit  ==========
  if (AnnualRes){
    CU          = "ChangePerYear"
    ScaleFactor = 1
  }else{
    CU          = "TotalChange"
    ScaleFactor = length(CTSR.VI)/12
  }
  # ========== define the container to hold the results ==========
  overview <- data.frame(
    ObservedChange = NA, CO2 = NA, LandUse = NA, ClimateTotal = NA, ClimateChange = NA,
    ClimateVariability = NA, OtherFactors = NA, OtherFactorsValid = FALSE,
    Obs.Pvalue = NA, CO2.Pvalue = NA, LandUse.Pvalue = NA, ClimateTotal.Pvalue = NA,
    ClimateChange.Pvalue = NA, ClimateVariability.Pvalue = NA)
  models <- list(
    OBS = FALSE, CO2 = FALSE, LU = FALSE, CT = FALSE, CC = FALSE, CV = FALSE
  )

  # ========== Function to deal with failure points  ==========
  ret.fail <- function(overview, models, reason, SkipError, errormessage=FALSE){
    # The areguments are the current state of the overview
    if (SkipError){
      return(structure(list(summary = overview, models = models, errors = reason)))
    }else{
      # ========== fail ==========
      if (!errormessage==FALSE){
        print(errormessage)
      }
      stop(paste("Failed with reason: ", reason, "Use SkipError=TRUE to return NaN results instead of rasing exception"))
    }
  }
  # ========== Perfrom the data check ==========
  if (any(is.na(CTSR.VI))){
    return(ret.fail(overview, models, "NANinCTSR.VI", SkipError))
  }else if (any(is.na(CTSR.RF))){
    return(ret.fail(overview, models, "NANinCTSR.RF", SkipError))
  }else if (any(is.na(CTSR.TM))){
    return(ret.fail(overview, models, "NANinCTSR.TM", SkipError))
  }
  allerror <- try({

    # ========== Get the annual Max VI values ==========
    rawmax.df  <- AnMaxVI(CTSR.VI)
    CTSR.VIadj <- franksCO2(CTSR.VI, C4frac)
    adjmax.df  <- AnMaxVI(CTSR.VIadj)

    # ========== Calculate Observed Change ==========
    ti       <- time(rawmax.df$Max)
    obsVI    <- rawmax.df$Max
    OBScor   <- cor.test(ti, obsVI, method="spearman")
    OBStheil <- mblm(obsVI~ti, repeated=FALSE)
    models$OBS = list(SpearmanRho=OBScor, TheilSen=OBStheil)
    overview$ObservedChange = as.numeric(OBStheil[[1]][2]) *ScaleFactor
    overview$Obs.Pvalue     = glance(OBScor)$p.value

    # ========== Calculate CO2 Change ==========
    if (C4frac < 1){
      CO2change <- rawmax.df$Max - adjmax.df$Max
      CO2cor    <- cor.test(ti, CO2change, method="spearman")
      CO2theil  <- mblm(CO2change~ti, repeated=FALSE)
      models$CO2 = list(SpearmanRho=CO2cor, TheilSen=CO2theil)
      overview$CO2        = as.numeric(CO2theil[[1]][2]) *ScaleFactor
      overview$CO2.Pvalue = glance(CO2cor)$p.value
    }else{
      # Avoid cases where the CO2 values will all be zero
      overview$CO2        = 0.0
      overview$CO2.Pvalue = 1.0
    }

    # ========== Perform the TSSRESTREND to get the landuse ==========
  	# +++++ Calculate the accumulated temperature and Precipitation +++++
    # Create the ACP and ACT table
    ACP.table <- climate.accumulator(CTSR.VIadj, CTSR.RF, max.acp, max.osp)
    ACT.table <- climate.accumulator(CTSR.VIadj, CTSR.TM, max.acp, max.osp, temperature = TRUE)

    # +++++ Calculate the landuse component with error handling+++++
    results <- TSSRESTREND(CTSR.VIadj, ACP.table, ACT.table = ACT.table, retnonsig=retnonsig)
    overview$LandUse = results$summary$Total.Change*ScaleFactor/(length(CTSR.VI)/12)
    if (results$summary$Method != "segmented.VPR"){
      overview$LandUse.Pvalue = results$summary$residual.p
    }else{
      p = c(results$summary$residual.p, results$summary$VPRbreak.p)
      p[p == 0] = 0.00000001
      # ========== define fisher method for combining p-values ==========
      fisher.method <- function(pvalues){
        pvalues <- p.adjust(pvalues, "bonferroni")
        df      <- 2*length(pvalues)
        pchisq( -2*sum(log(pvalues), na.rm=TRUE), df, lower.tail=FALSE )
      }
      overview$LandUse.Pvalue = fisher.method(p)
    }

    # ========== Calculate the climate part ==========
    if (splitclim){
      # +++++ work out the pead veg location +++++
      uniqv <- unique(rawmax.df$Max.month)
      pveg  <- uniqv[which.max(tabulate(match(as.numeric(rawmax.df$Max.month), uniqv)))]

      # +++++ Get the new temperature and rainfall data tables +++++
      ACP.lr   <- climate.accumulator(CTSR.VIadj, CTSR.RF, max.acp, max.osp, cliwindow=cliwindow)
      ACT.lr   <- climate.accumulator(CTSR.VIadj, CTSR.TM, max.acp, max.osp, cliwindow=cliwindow, temperature = TRUE)
      VI.index <- seq(0, dim(ACP.lr)[2]-1, 12) + c(rep(pveg, cliwindow), rawmax.df$Max.month)
      if (results$summary$Method != "segmented.VPR"){
        # +++++ get the longer annual values +++++
        an.CC  <- roll_mean(ACP.lr[paste(c(results$acum.df$osp, results$acum.df$acp), collapse ="-"), VI.index], cliwindow)
        # ===== Check if Temperater data is included
        if (!is.na(results$acum.df$tacp)){
          # temperature data
          an.tw  <- roll_mean(ACT.lr[paste(c(results$acum.df$tosp, results$acum.df$tacp), collapse ="-"), VI.index], cliwindow)
          df.CC <- data.frame("acu.RF"=an.CC, "acu.TM"=an.tw)
        }else{
          # +++++ predict new values +++++
          df.CC <- data.frame("acu.RF"=an.CC)

        }
        # ========== Fix a dim error ==========
        if (dim(df.CC)[1] == (dim(rawmax.df)[1])+1){
          df.CC = tail(df.CC, -1)
        }else if (dim(df.CC)[1] != dim(rawmax.df)){
          stop("Unknown dim error")
        }
        VCRmod <- results$TSSRmodels$VPR.fit

        # browser()

      }else{
        # ========== Breakpoint the the VCR ==========
        # +++++ get the longer annual values +++++
        an.CCb4  <- roll_mean(ACP.lr[paste(c(results$acum.df$osp.b4, results$acum.df$acp.b4), collapse ="-"), VI.index], cliwindow)
        an.CCaf  <- roll_mean(ACP.lr[paste(c(results$acum.df$osp.af, results$acum.df$acp.af), collapse ="-"), VI.index], cliwindow)
        len      <- dim(rawmax.df)[1]
        if (length(an.CCb4) > len){
          an.CCb4  <- tail(an.CCb4, -1)
          an.CCaf  <- tail(an.CCaf, -1)
        }
        # +++++ calculate the before and after means +++++
        b4mean <- mean(ACP.table[paste(c(results$acum.df$osp.b4, results$acum.df$acp.b4), collapse ="-"), as.numeric(rawmax.df$index)])
        afmean <- mean(ACP.table[paste(c(results$acum.df$osp.af, results$acum.df$acp.af), collapse ="-"), as.numeric(rawmax.df$index)])
        b4std  <- sd(ACP.table[paste(c(results$acum.df$osp.b4, results$acum.df$acp.b4), collapse ="-"), as.numeric(rawmax.df$index)])
        afstd  <- sd(ACP.table[paste(c(results$acum.df$osp.af, results$acum.df$acp.af), collapse ="-"), as.numeric(rawmax.df$index)])
        CCb4z <- (an.CCb4-b4mean)/b4std
        CCafz <- (an.CCaf-afmean)/afstd
        breakpoint <- results$ols.summary$chow.sum$yr.index

        # ===== Setup the regression variables =====
        # Build a single adjusted precip vector
        adj.RF <- c(CCb4z[1:breakpoint], CCafz[(breakpoint + 1):len])
        # Create the dummy variable
        dummy <- rep(0, len)
        dummy[(breakpoint + 1):len] = 1

        # ===== Pull out the model =====
        VCRmod <- results$TSSRmodels$segVPR.fit

        # ========== do the same for temperature ==========
        if (!is.na(results$ols.summary$OLS.table["segVPR.fit", "temp.coef"])){

          # +++++ get the longer annual values +++++
          an.tCCb4  <- roll_mean(ACT.lr[paste(c(results$acum.df$tosp.b4, results$acum.df$tacp.b4), collapse ="-"), VI.index], cliwindow)
          an.tCCaf  <- roll_mean(ACT.lr[paste(c(results$acum.df$tosp.af, results$acum.df$tacp.af), collapse ="-"), VI.index], cliwindow)

          if (length(an.tCCb4) > len){
            an.tCCb4  <- tail(an.tCCb4, -1)
            an.tCCaf  <- tail(an.tCCaf, -1)
          }
          # # +++++ calculate the before and after means +++++
          tb4mean <- mean(ACT.table[paste(c(results$acum.df$tosp.b4, results$acum.df$tacp.b4), collapse ="-"), as.numeric(rawmax.df$index)])
          tafmean <- mean(ACT.table[paste(c(results$acum.df$tosp.af, results$acum.df$tacp.af), collapse ="-"), as.numeric(rawmax.df$index)])
          tb4std  <- sd(ACT.table[paste(c(results$acum.df$tosp.b4, results$acum.df$tacp.b4), collapse ="-"), as.numeric(rawmax.df$index)])
          tafstd  <- sd(ACT.table[paste(c(results$acum.df$tosp.af, results$acum.df$tacp.af), collapse ="-"), as.numeric(rawmax.df$index)])
          tCCb4z <- (an.tCCb4-tb4mean)/tb4std
          tCCafz <- (an.tCCaf-tafmean)/tafstd
          adj.TM <- c(tCCb4z[1:breakpoint], tCCafz[(breakpoint + 1):len])
          # ========== Build a dataframe ==========
          df.CC <- data.frame("sv.RF"=adj.RF, "sv.TM"=adj.TM, "breakpoint.var" = dummy)

        }else{
          # ========== Build a dataframe ==========
          df.CC <- data.frame("sv.RF"=adj.RF, "breakpoint.var" = dummy)
        }
      }
      # ========== Get the climate var removed VI estimates ==========
      CC.vi    <- predict(VCRmod, df.CC)
      CCcor    <- cor.test(ti, CC.vi, method="spearman")
      CCtheil  <-  mblm(CC.vi~ti, repeated=FALSE)
      models$CC = list(SpearmanRho=CCcor, TheilSen=CCtheil)
      overview$ClimateChange        = as.numeric(CCtheil[[1]][2]) *ScaleFactor
      overview$ClimateChange.Pvalue = glance(CCcor)$p.value

      # ========== Get the climate Varibility ==========
      CV.vi <- VCRmod$fitted.value - CC.vi
      CVcor    <- cor.test(ti, CV.vi, method="spearman")
      CVtheil  <-  mblm(CV.vi~ti, repeated=FALSE)
      models$CV = list(SpearmanRho=CVcor, TheilSen=CVtheil)
      overview$ClimateVariability        = as.numeric(CVtheil[[1]][2]) *ScaleFactor
      overview$ClimateVariability.Pvalue = glance(CVcor)$p.value
      # ========== Calculate the OtherFactors ==========
      overview$OtherFactors = overview$ObservedChange - (overview$CO2+overview$LandUse+overview$ClimateChange+overview$ClimateVariability)
    }else{
      # ========== No climate seperation ==========
      if (results$summary$Method != "segmented.VPR"){
        # +++++ Calculate the values +++++
        CLIres   <- results$TSSRmodels$VPR.fit$fitted.values
      }else{
        CLIres   <- results$TSSRmodels$segVPR.fit$fitted.values
      }
      CLIcor   <- cor.test(ti, CLIres, method="spearman")
      CLItheil <-  mblm(CLIres~ti, repeated=FALSE)
      overview$ClimateTotal        = as.numeric(CLItheil[[1]][2]) *ScaleFactor
      overview$ClimateTotal.Pvalue = glance(CLIcor)$p.value

      overview$OtherFactors = overview$ObservedChange-(overview$CO2+overview$LandUse+overview$ClimateTotal)
    }

    # ========== Check if Other Factors is valid ==========

    if ((overview$LandUse != 0) && (results$summary$Method[1] %in% c("RESTREND", "segmented.RESTREND", "segmented.VPR"))){
      overview$OtherFactorsValid = TRUE
    }
    # ========== Return Errors ==========
    return(structure(list(summary = overview, models = models, errors = "")))
  })
  # ========== Handle any exceptions ==========

  if (class(ret.fail) == "try-error"){
    print("deal with all error here ")
    browser()
    return(ret.fail(overview, models, "AttributionFailed", SkipError, errormessage=ret.fail))
  }

}

