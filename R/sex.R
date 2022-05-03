#' Set the sex for the sample
#'
#' @param X object of type cloneobj,
#' @param v purity (0-1)
#' @return the updated object X with the purity set
#' @export
setSex=function(X,v){
  UseMethod("setSex")
}

setSex.cloneobj=function(X,v){
  v=as.character(v) %>% toupper()
  if(! v  %in% c("XX","XY")){ stop("Please provide a valid string for sex (either XX for female or XY for male)")}
  X$sex=v
  return(X)
}

#' Get the sex for the sample
#'
#' @param X object of type cloneobj,
#' @return the purity of the sample
#' @export getSex
#' @export
getSex=function(X){
  UseMethod("getSex")
}

#' @export getSex.cloneobj
#' @export
getSex.cloneobj=function(X){
  return(X$sex)
}



#' Infer the sex for a sample
#'
#' If events occur on the Y chromosome the sample is assigned to M
#' otherwise F
#' @param X object of type GRanges or GRangeslist
#' @return the sex of the sample as XX for female and XY for male
inferSex=function(gr){
  sex=NULL

  if( grepl( "GRangesList", class(gr)) ){
    male=which( seqnames(gr) %>% unlist() %>% unique() %in% c("Y","chrY")) %>% unlist
  }else{
    male=which( seqnames(gr) %in% c("Y",'chrY') )
  }
  if(length(male)>0){ sex='M'}else{sex='F'}
  return(sex)
}
