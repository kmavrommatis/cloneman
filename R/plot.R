#' plot the clone information
#'
#' @param X cloneobj
#' @param overlay
#' @export plot.cloneobj
#' @export
plot.cloneobj=function(X, overlay=FALSE,chromosomes=NULL){
  genome=X$genome
  X=cloneman::.finalize( X )
  #X=getSegmentFraction(X)
  gtrack <- Gviz::GenomeAxisTrack()
  ## subclone
  sub_track=getSubcloneTrack( X )

  hist_track=getHistTrack( X )

  #ref=getReferenceTrack(X)

  #if(!is.null(X$problems)){
  #  prob=getProblemsTrack( X )
  #}else{prob=GenomeAxisTrack()}


  #sl=GenomeInfoDb::sortSeqlevels( as(X$segments_list,"GRangesList") ) %>% GenomeInfoDb::seqlevels()
  sl=seqnames(X$segments_list) %>% unlist %>% unique()
  if(!is.null(chromosomes)){
    sl=intersect(sl, chromosomes)
  }
  if(overlay==FALSE){
    chroms <- data.frame(chromosome = sl)
    lattice::xyplot(1 ~ chromosome | chromosome, data = chroms, panel = function(x) {
      itrack <- Gviz::IdeogramTrack(genome = genome, chromosome = x)
      Gviz::plotTracks(c(itrack,
                   gtrack,
                   sub_track,
                   OverlayTrack(hist_track)
      ),chromosome = x, title.width = 2, frame=TRUE, add=TRUE)

    }, layout=c(round(sqrt(length(sl))),round(sqrt(length(sl)))), aspect = 1,  as.table=TRUE )
  }else{
    chroms <- data.frame(chromosome = sl)
    lattice::xyplot(1 ~ chromosome | chromosome, data = chroms, panel = function(chr) {
      itrack <- Gviz::IdeogramTrack(genome = genome, chromosome = chr)

      Gviz::plotTracks(list( itrack,
                       gtrack,
                       OverlayTrack( sub_track),
                       OverlayTrack( hist_track)
      ),chromosome = chr,legend=TRUE)

    }, layout=c(round(sqrt(length(sl))),round(sqrt(length(sl)))), aspect = 1,  as.table=TRUE)
  }
}

#' helper function to find the maximum number of CNV copies
#'
#' @param X cloneobj
#' @return Returns the maximum copy number among all the clones.
.maxCNV=function(X){
  ll=lapply(X$segments_list, function(K){
    max( K$total.allele.copies, na.rm=TRUE)
  })
  return(max(unlist(ll),na.rm = TRUE))
}

#' get the track for a subclone
#'
#' The track containst the segments with the various copy numbers
getSubcloneTrack=function(X,name=NULL){

  mx=.maxCNV(X)
  col=1
  sub_track=list()

  for( s in X$clone_info$clone.name ){
    if(is.null(X$segments_list[[s]])){next}
    jitter=col/20
    #message("Name is ",name, " s is ",s)
    #jitter=overlay
    if(s=="normal" | s=="germline"){next}
    n=s
    if(!is.null(name) ){n=name}

    sub_track[[ s ]]=Gviz::DataTrack( X$segments_list[[s]],
                                data=X$segments_list[[s]]$total.allele.copies+jitter * runif( 1, min=-1 , max=1) - X$ploidy,
                                type="histogram",
                                col.histogram=col,
                                fill.histogram="white",
                                baseline=0,
                                name = n,
                                showTitle=TRUE,
                                rotation.title=0,
                                cex.title=0.5,
                                ylim=c(-2,mx-2))
    col=col+1
  }
  return(sub_track)
}

#' get the CCF of a subclone
#'
#' The track contains the segments and the CCF of each of the segments

getHistTrack=function(X){
  # histogram with ccf of segment
  col=1
  hist_track=list()
  combcnv=getCombinedCNV(X)
  #for(cn in  X$cns){
    #label=paste0("CN_",cn)
    #message("Working with ", label)
    hist_track[[ col ]]=DataTrack( combcnv,
                                   data=mcols(combcnv)[ ,"CCF" ] ,
                                   name="ccf",
                                   type="histogram",
                                   baseline=0,
                                   showTitle=TRUE,
                                   rotation.title=0,
                                   cex.title=0.5,
                                   col.histogram=col,
                                   fill.histogram=col,
                                   alpha=0.5,
                                   alpha.title=1,
                                   ylim=c(0,1))

    col=col+1
  #}
  hist_track
}
