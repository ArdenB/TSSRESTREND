#
# demo.stdRESTEND <- function(sig=0.05, print=TRUE, plot=TRUE, details=FALSE, mode="TSS.RESTREND"){
#   #set the environment variables
#   CTSR.VI <- stdRESTREND.CTSR$cts.NDVI
#   CTSR.RF <- stdRESTREND.CTSR$cts.precip
#   anu.VI <- stdRESTREND$max.NDVI
#   acu.RF <- stdRESTREND$acc.precip
#   VI.index <- stdRESTREND$index
#   if (mode == "TSS.RESTREND"){
#     TSSR.result <-TSS.RESTREND(CTSR.VI, CTSR.RF, anu.VI, acu.RF, VI.index, sig=sig, print=print, plot=plot, details = details)
#   }else{
#     #drops into a browser so the user can call the functions individually
#     print("loading standard RESTREND environment variables")
#     browser()
#   }
# }
#
# demo.segRESTEND <- function(sig=0.05, print=TRUE, plot=TRUE, details=FALSE, mode="TSS.RESTREND"){
#   #set the environment variables
#   CTSR.VI <- segRESTREND.CTSR$cts.NDVI
#   CTSR.RF <- segRESTREND.CTSR$cts.precip
#   anu.VI <- segRESTREND$max.NDVI
#   acu.RF <- segRESTREND$acc.precip
#   VI.index <- segRESTREND$index
#   if (mode == "TSS.RESTREND"){
#     TSSR.result <-TSS.RESTREND(CTSR.VI, CTSR.RF, anu.VI, acu.RF, VI.index, sig=sig, print=print, plot=plot, details = details)
#   }else{
#     #drops into a browser so the user can call the functions individually
#     print("loading segmented RESTREND environment variables")
#     browser()
#     #need to add the other modes
#   }
# }
#
# demo.segVPRD <- function(sig=0.05, print=TRUE, plot=TRUE, details=FALSE, mode="TSS.RESTREND"){
#   #set the environment variables
#   CTSR.VI <- segVPRD.CTSR$cts.NDVI
#   CTSR.RF <- segVPRD.CTSR$cts.precip
#   anu.VI <- segVPRD$max.NDVI
#   acu.RF <- segVPRD$acc.precip
#   VI.index <- segVPRD$index
#   rf.b4 <- segVPRD$acp.b4
#   rf.af <- segVPRD$acp.af
#
#   if (mode == "TSS.RESTREND"){
#     TSSR.result <-TSS.RESTREND(CTSR.VI, CTSR.RF, anu.VI, acu.RF, VI.index, rf.b4, rf.af, sig=sig, print=print, plot=plot, details = details)
#   }else{
#     #drops into a browser so the user can call the functions individually
#     print("loading segmented RESTREND environment variables")
#     # browser()
#     #need to add the other modes
#   }
# }
#
# demo.segVPRI <- function(sig=0.05, print=TRUE, plot=TRUE, details=FALSE, mode="TSS.RESTREND"){
#   #set the environment variables
#   CTSR.VI <- segVPRI.CTSR$cts.NDVI
#   CTSR.RF <- segVPRI.CTSR$cts.precip
#   anu.VI <- segVPRI$max.NDVI
#   acu.RF <- segVPRI$acc.precip
#   VI.index <- segVPRI$index
#   rf.b4 <- segVPRI$acp.b4
#   rf.af <- segVPRI$acp.af
#
#   if (mode == "TSS.RESTREND"){
#     TSSR.result <-TSS.RESTREND(CTSR.VI, CTSR.RF, anu.VI, acu.RF, VI.index, rf.b4, rf.af, sig=sig, print=print, plot=plot, details = details)
#   }else{
#     #drops into a browser so the user can call the functions individually
#     print("loading segmented RESTREND environment variables")
#     browser()
#     #need to add the other modes
#   }
# }
#
# #need to figure add a return
# # res <- demo.stdRESTEND()
# # browser()
# # res <- demo.segRESTEND()
# # browser()
# # res <- demo.segVPRD()
# # browser()
# # res <- demo.segVPRI()
# # browser()
#
#
# # print("hello World")
#
#
#
# # a<- BFAST.RESID(stdRESTREND.CTSR$cts.NDVI, stdRESTREND.CTSR$cts.precip, print=TRUE, plot=TRUE)
# # print(a)
# #
# #
# # #fin functions use class(a) to determine if its numeric or logical
# #
# # se<- BFAST.RESID(segRESTREND.CTSR$cts.NDVI, segRESTREND.CTSR$cts.precip, print=TRUE, plot=TRUE)
# # print(se)
# # #CHOW <- function(anu.VI, acu.RF, VI.index, se, sig=0.05, print=FALSE)
# #
# # #for testing the chow test
# # anu.VI <- segRESTREND$max.NDVI
# # acu.RF <- segRESTREND$acc.precip
# # VI.index <- segRESTREND$index
# #
# # se.chow <- CHOW(anu.VI, acu.RF, VI.index, se, sig=0.05, print=TRUE)
# # breakpoint = as.integer(se.chow$bp.summary[2])
# # se.RES <- seg.RESTREND(anu.VI, acu.RF, VI.index, breakpoint,  sig=0.05, print=TRUE, plot=TRUE)
# #
# # anu.VI <- stdRESTREND$max.NDVI
# # acu.RF <- stdRESTREND$acc.precip
# # VI.index <- stdRESTREND$index
# # res <- RESTREND(anu.VI, acu.RF, VI.index, sig=0.05, print=TRUE, plot=TRUE)
# #
# # res <- RESTREND(stdRESTREND$max.NDVI, stdRESTREND$acc.precip, stdRESTREND$index, sig=0.05, print=TRUE, plot=TRUE)
# #
