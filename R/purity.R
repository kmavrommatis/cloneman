#' Set the predicted purity for the sample
#'
#' @param X object of type cloneobj,
#' @param v purity (0-1)
#' @return the updated object X with the purity set
#' @export setPurity
#' @export
setPurity=function(X,v){
  UseMethod("setPurity",X)
}

#' @export setPurity.cloneobj
#' @export
setPurity.cloneobj=function(X,v){
  v=as.numeric(v)
  if(v > 1){ stop("Please provide a purity value between 0 and 1")}
  X$purity=v
  return(X)
}

#' Get the predicted purity for the sample
#'
#' @param X object of type cloneobj,
#' @return the purity of the sample
#' @export getPurity
#' @export
getPurity=function(X){
  #message("redirecting to getPurity for object of class", class(X))
  UseMethod("getPurity")
}

#' @export getPurity.cloneobj
#' @export
getPurity.cloneobj=function(X){
  #message("inside getPurity")
  return(X$purity)
}

#' @export getPurity.default
#' @export
getPurity.default=function(X){
  #message("inside default getPurity")
}
