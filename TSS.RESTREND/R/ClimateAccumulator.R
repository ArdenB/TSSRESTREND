#' @title Climate Accumulator
#'
#' @description
#' Takes the time series of rainfall and returns a rainfall accumulation table of every possible combination
#' of the max accumulation period and the max offset period.
#' @importFrom RcppRoll roll_sum roll_mean
#' @param CTSR.VI
#' Complete Time Series of Vegetation Index. An object of class \code{'ts'}. Monthly time series of VI values
#' @param clim.data
#' Complete Time Series of monthly rainfall or temperature. An object of class \code{'ts'}. Must
#' have the same end date as CTSR.VI and be longer than the CTSR.VI by more than the max acsumuation
#' period (max.acp) plus the max offset period.(max.ops)
#' @param max.acp
#' The max accumuation period. Must be an integer > 1.
#' @param max.osp
#' The max offset period. Must be an integer >1
#' @param temperature
#' Bool. If the clim.data being accumulated is temperature, will take a mean not a sum.  This makes
#' it easier to comapare regions with different accumulation and offset periods. This is new in v0.3.0.
#' Defualt=FALSE which replicates the behaviour of TSSRESTREND versions <0.2.16.
#' @param cliwindow
#' The size of the window in years to be used for calculating climate change.
#' @return ACP.table
#' A matrix with ever possible accumuated climate combination
#'
#' @export
#'
#' @examples
#' # Define the max accumuulation period
#' acp <- 12
#' #Define the max offset period
#' osp <- 4
#' rftable <- climate.accumulator(segRESTRENDCTSR$cts.NDVI, segRESTRENDctRF$precip, acp, osp)

climate.accumulator <- function(CTSR.VI, clim.data, max.acp, max.osp, temperature=FALSE, cliwindow=0){
  # ==============================================================================================
  # ========== Sanity check the input data ==========
  if (class(CTSR.VI) != "ts")
    stop("CTSR.VI Not a time series object")
  if (class(clim.data) != "ts")
    stop("clim.data Not a time series object")
  if (sd(clim.data) == 0)
    stop("The precipitation data has identical values (SD=0)")
  # +++++  Get the start and end dates of the precip and the VI +++++
  yst <- start(CTSR.VI)[1] - cliwindow
  mst <-  start(CTSR.VI)[2]
  y.en <- end(CTSR.VI)[1]
  m.en <- end(CTSR.VI)[2]
  clim.yend <- end(clim.data)[1]
  clim.mend <-  end(clim.data)[2]
  # +++++ Check to make sure they have no issues +++++
  if ((y.en != clim.yend) || clim.mend != m.en) {
    stop("clim.data does not end at the same time as CTSR.VI")}
  if (length(clim.data) < (length(CTSR.VI) + max.acp + max.osp+(cliwindow*12))) {
    stop("clim.data is not long enough for the set max.acp, max.ops and cliwindow")}

  # ==============================================================================================
  # ========== Build the accumulation table ==========
  # ===== Set up the table =====

  # make the row names
  row.nm <- rep(0, max.acp)
  if (max.osp > 1) {
    for (n in 1:(max.osp - 1)) {
      row.nm <- c(row.nm, rep(n, max.acp))
    }
  }
  # get the length of the matrix
  # len <- length(CTSR.VI) + cliwindow*12
  len <- (1+y.en-yst) * 12

  # +++++ Set up a blank matrix to write into +++++
  m <- matrix(nrow = (max.acp * max.osp), ncol = len)
  # set the row and col names
  rownames(m) <- paste(row.nm, rep(1:max.acp, max.osp), sep = "-")
  colnames(m) <- c(1:len)

  m2 <- matrix(nrow = (max.acp), ncol = length(clim.data))
  # Reverse the climate data to make it easy to get the start date
  rev.rf = rev(clim.data)

  # ===== pupulate the table with accumulation values =====
  # get the values for just the accumulation period
  for (n in 1:max.acp) {
    # +++++ perform the roll +++++
    if (temperature){
      # rolling mean the reversed climate data
      roll = (roll_mean(rev.rf, n))
      }else{
      # rolling sum the reversed climate data
      roll = (roll_sum(rev.rf, n))
      }
    if (n > 1) {
      # Add nans to the the end of the reversed climate data to make them the same size
      roll = c(roll, rep(NaN, (n - 1)))
    }
    # append the rolled and reversed data
    m2[n, ] = roll
  }
  # +++++ turn the suplied ts in a table +++++
  for (osp in 0:(max.osp - 1)) {
    # Extract the osp for each of the acp
    #   defence against index 1 i missing
    m3 <- (m2[,(1 + osp):(len + osp)])
    # flip the results back around the right way
    m4 <- (m3[, ncol(m3):1])
    # Add them to the table
    ind <- 1 + (osp*max.acp)
    m[ind:(ind + max.acp - 1),] <- m4

  }
  # Return the table
  return(m)
}
