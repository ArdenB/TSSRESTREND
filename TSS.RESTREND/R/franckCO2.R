#' @title Franks CO2 Vegetation correction
#'
#' @importFrom stats coef end frequency lm sd start time ts cor.test window
#' @importFrom graphics abline arrows legend par plot
#' @importFrom utils tail read.csv
#' @importFrom broom glance
#'
#' @description
#' Adjusts a vegetation time series to account for CO2 fertilisation
#'
#' @author Arden Burrell, aburrell@whrc.org
#'
#' @inheritParams TSSRESTREND
#' @inheritParams TSSRattribution
#' @param CO2
#'        A timeseries containg the CO2 concentration. The defualt is CMIP5 RCP8.5
#' @param refyear
#'        The Year that acts as a baseline for CO2. All vegetation values will be scaled #
#'		  to the CO2 concentration of this year. Defulat is 1980. THis function will pick the
#'		  first value in the selected year,
#'
#' @return \bold{CTSR.VIadj}
#'        A version of the CTSR.VI data that has been adjusted to account for CO2


#' @export

franksCO2 <- function(CTSR.VI, C4frac, CO2=FALSE, refyear=1980){
	# ========== Get the CO2 data ==========
	if (CO2 == FALSE){
		# ===== make the C)2 data =====
		CO2con = c(316.272,317.075,317.795,318.397,318.925,319.647,320.647,
			321.605,322.635,323.902,324.985,325.855,327.14,328.678,329.743,
			330.585,331.748,333.272,334.848,336.525,338.36,339.728,340.793,
			342.197,343.783,345.283,346.798,348.645,350.738,352.488,353.855,
			355.017,355.885,356.777,358.127,359.837,361.462,363.155,365.322,
			367.348,368.865,370.467,372.522,374.76,376.812,378.812,380.827,
			382.777,384.8,387.012,389.324,391.638,394.009,396.464,399.004,
			401.628,404.328,407.096,409.927,412.822,415.78,418.796,421.864,
			424.995,428.197,431.475,434.826,438.245,441.721,445.251,448.835,
			452.474,456.177,459.964,463.852,467.85,471.96,476.182,480.508,
			484.927,489.435,494.032,498.73,503.53,508.433,513.456,518.611,
			523.9,529.324,534.875,540.543)
		CO2 = ts(CO2con, start=c(1960, 1), end=c(2050,1), frequency = 1)
	}else if (class(CO2) != "ts")
	  stop("CTSR.VI Not a time series object. Please check the data")


  franks_FvC <- function(Ca){
	    theta = 0.7         # shape of the light response curve
	    gamma_star = 40.0   # CO2 compensation point
	    return ((theta * Ca - gamma_star) / (theta * Ca + 2.0 * gamma_star))
	 }

	# ========== Calculate the baseline ==========
	baseline_Ca = window(CO2, 1980, 1981)[[1]]
	baseAnet    = franks_FvC(baseline_Ca)

	# +++++ Get the start year and start month of the data +++++
	yst <- start(CTSR.VI)[1]
	# mst <- start(CTSR.VI)[2]
	yfn <- start(tail(CTSR.VI, n=1))[1]

	# # ========== get the CO2 concentration of the dataset years ==========
	# +++++ build a matching CO2 series +++++
	Ca <- as.numeric(window(CO2, yst, yfn))
	# +++++ Check the length +++++
	if (length(CTSR.VI)/12 == length(Ca)){
	  Ca <- rep(Ca, each=12)
	}else{
	  stop("Lengh problem in the CO2 adjustment, Ts must be either annual or monthly")
	}

	Anet           = franks_FvC(Ca)
	model_response = Anet / baseAnet

	# ========== Adjust the CTSRVI ==========
	VIadjC3 <- CTSR.VI / model_response

	# ========== Account for C4  ==========
	CTSR.VIadj <- (VIadjC3 * (1-C4frac)) + (CTSR.VI *C4frac)

	return(CTSR.VIadj)

}
