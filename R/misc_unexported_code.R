###################################
# package handlers
###################################

#' @useDynLib alignment
#' @importFrom Rcpp sourceCpp

#' @importFrom stats approx cor
#' @importFrom utils capture.output methods modifyList

#undefined globals
utils::globalVariables(c("index", "scores", "warp.y", "iteration"))


##################################
# to think about
##################################

# could move all unexported code to here ???
#    would be more obvious...

# think about unexported function naming
#    currently align_ for any common _align() code handler
#    and alignment_ for anythin that is for use with a alignment object
#        this is from when package was called align...
#            (maybe leave but it does not really work post name change...)


##################################
#cow_convert
##################################

#kr v.0.0.1

#converts cow outputs nSeg and BT
#into something that can be used
#by warp_frame

#probably not staying or not staying as is...

cow_convert <- function(nSeg, bT){
  bT <- as.vector(bT)
  d <- diff(nSeg)

  ans <- c(nSeg[1])
  for(i in 1:length(d)){
    ans <- c(ans,
             seq(bT[i]+1,bT[i+1], length.out=d[i]))
  }

  data.frame(x=min(nSeg):max(nSeg),
             y=ans)
}
