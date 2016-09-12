

VImax.plot <- function(anu.VI, brkp=FALSE){
  if (!brkp){
#     z <- rep(0, len)
#     z[(brkp+1):len] = 1
    yst <- start(anu.VI)[1]
    ynd <- end(anu.VI)[1]
    t = c(yst:ynd)
    plot(t, anu.VI, pch=16, xlab="Year", ylab="Annual max VI", col="orange") 

  }else{
    yst <- start(anu.VI)[1]
    ynd <- end(anu.VI)[1]
    len <- length(anu.VI)
    t = c(yst:ynd)
    plot(t[1:brkp], anu.VI[1:(brkp)], pch=16,xlab="Year", 
         ylab="Annual max VI", col="orange", xlim=c(yst, ynd),  
         ylim=c(min(anu.VI), max(anu.VI)) )
    par(new=T)    
    plot(t[(brkp+1):len], anu.VI[(brkp+1):len], pch=16, 
         xlab="", ylab="", col="purple",main="", xlim=c(yst, ynd), 
         ylim=c(min(anu.VI), max(anu.VI)) 
         )
    abline(v=(brkp-0.5+yst), col="red", lty = "dashed")
  }
}

