
#' Set the predicted ploidy for the sample
#'
#' Two parameters are set:
#'
#' **ploidy** is the rounded ploidy (e.g. 2.1 will become 2)
#'
#' **real_ploidy** is the actual value provided by the user
#'
#' @param X object of type cloneobj,
#' @param v ploidy (float)
#' @return the updated object X with the ploidy set
#' @export setPloidy
#' @export
setPloidy=function(X,v){
  UseMethod("setPloidy")
}
#' @export setPloidy.default
#' @export
setPloidy.default=function(X,v){
  message("This is the default for method setPloidy")
}
#' @export setPloidy.cloneobj
#' @export
setPloidy.cloneobj=function(X,v){
  X$ploidy=round(as.numeric(v))
  X$real.ploidy=as.numeric(v)
  return(X)
}

#' Get the predicted ploidy for the sample
#'
#' Two parameters are returned:
#'
#' **ploidy** is the rounded ploidy (e.g. 2.1 will become 2)
#'
#' **real_ploidy** is the actual value provided by the user
#'
#' @param X object of type cloneobj,
#' @return a list with the _real.ploidy_ and _ploidy_
#' @export getPloidy
#' @export
getPloidy=function(X){
  UseMethod("getPloidy")
}

#' @export getPloidy.cloneobj
#' @export
getPloidy.cloneobj=function(X){
  return(list("ploidy"=X$ploidy,
              "real.ploidy"=X$real.ploidy))
}
