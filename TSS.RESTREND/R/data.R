#' @title
#' Rabbit Impacted Vegetation Precipitation Accumulation Table
#' @docType data
#' @description
#' Rainfall Accumulation Table form a part of the Simpson-sStrzelecki Dunefields bioregion
#' impacted by reduced rabbit predation after the release of RHD.
#'
#' @format A matrix containing the complete time series of every acp and ofp. See \code{\link{climate.accumulator}} for details
#' @source Awap data from \url{http://www.csiro.au/awap/cgi/awap2.pl}
#' @examples
#' data(rabbitACPtable)
"rabbitACPtable"

#' @title
#' Data frame containing the Complete Times Series data for a standard Restrend analysis
#' @description
#' contains CTSR.VI and CTSR.RF. See \code{\link{TSSRESTREND}}. Starts 1/1982 and extend to 12/2013
#' @docType data
#'
#' @format R data frame
#' @source a single pixel from \code{\link[gimms]{monthlyComposite}}
#' @seealso \code{\link{stdRESTREND}} for the annual values for the same pixel
"stdRESTRENDCTSR"

#' @title
#' Data frame containing the raw rainfall data set ending dec 2013 with a frquency of 12
#' @description
#' contains raw monthly precipitation
#' @docType data
#'
#' @format R data frame
#' @seealso \code{\link{stdRESTREND}} for the annual values for the same pixel
"stdRESTRENDctRF"


#' @title
#' Data frame containing the annual data for a standard Restrend analysis
#' @description
#' contains anu.VI, acu.RF and VI.index. Range 1982 - 2013
#' @docType data
#' @format R data frame
#' @source \code{\link[gimms]{gimms-package}}
"stdRESTREND"


#' @title
#' Precipitation Accumulation Table for the standard RESTREND demonstration pixel
#' @docType data
#' @description
#' Rainfall Accumulation Table
#'
#' @format A matrix containing the complete time series of every acp and ofp. See \code{\link{climate.accumulator}} for details
#' @source Awap data from \url{http://www.csiro.au/awap/cgi/awap2.pl}
"stdRESTRENDrfTab"

#' @title
#' Data frame containing the Complete Times Series data for a segmented VPR (VEGETATION PRECIPTATION RELATIONSHIP) analysis
#' @description
#' contains CTSR.VI and CTSR.RF. See \code{\link{TSSRESTREND}}. Starts 1/1982 and extend to 12/2013
#' @docType data
#'
#' @format R data frame
#' @source a single pixel from \code{\link[gimms]{monthlyComposite}}
#' @seealso \code{\link{stdRESTREND}} for the annual values for the same pixel
"segVPRCTSR"

#' @title
#' Data frame containing the raw rainfall data set for the segVPR data,  ending dec 2013 with a frquency of 12
#' @description
#' contains raw monthly precipitation
#' @docType data
#'
#' @source Awap data from \url{http://www.csiro.au/awap/cgi/awap2.pl}
#' @format R data frame
#' @seealso \code{\link{stdRESTREND}} for the annual values for the same pixel
"segVPRctRF"


#' @title
#' Data frame containing the annual data for a segVPR analysis
#' @description
#' contains anu.VI, acu.RF,  VI.index, rf.b4, rf.af. Range 1982 - 2013. Breakpoint for this pixel is 24 (2005)
#' @docType data
#' @format R data frame
#' @source \code{\link[gimms]{gimms-package}}
"segVPR"


#' @title
#' Precipitation Accumulation Table for the segVPR demonstration pixel
#' @docType data
#' @description
#' Rainfall Accumulation Table
#'
#' @format A matrix containing the complete time series of every acp and ofp. See \code{\link{climate.accumulator}} for details
#' @source Awap data from \url{http://www.csiro.au/awap/cgi/awap2.pl}
"segVPRrfTab"


#' @title
#' Data frame containing the Complete Times Series data for a segmented RESTREND analysis
#' @description
#' contains CTSR.VI and CTSR.RF. See \code{\link{TSSRESTREND}}. Starts 1/1982 and extend to 12/2013
#' @docType data
#'
#' @format R data frame
#' @source a single pixel from \code{\link[gimms]{monthlyComposite}}
#' @seealso \code{\link{stdRESTREND}} for the annual values for the same pixel
"segRESTRENDCTSR"

#' @title
#' Data frame containing the raw rainfall data set for the segRESTREND data,  ending dec 2013 with a frquency of 12
#' @description
#' contains raw monthly precipitation
#' @docType data
#'
#' @source Awap data from \url{http://www.csiro.au/awap/cgi/awap2.pl}
#' @format R data frame
#' @seealso \code{\link{stdRESTREND}} for the annual values for the same pixel
"segRESTRENDctRF"


#' @title
#' Data frame containing the annual data for a segRESTREND analysis
#' @description
#' contains anu.VI, acu.RF,  VI.index, rf.b4, rf.af. Range 1982 - 2013. Breakpoint for this pixel is 24 (2005)
#' @docType data
#' @format R data frame
#' @source \code{\link[gimms]{gimms-package}}
"segRESTREND"


#' @title
#' Precipitation Accumulation Table for the segRESTREND demonstration pixel
#' @docType data
#' @description
#' Rainfall Accumulation Table
#'
#' @format A matrix containing the complete time series of every acp and ofp. See \code{\link{climate.accumulator}} for details
#' @source Awap data from \url{http://www.csiro.au/awap/cgi/awap2.pl}
"segRESTRENDrfTab"
