#' plot an integrated plot with CNV information
#'
#' Each plot is per chromosome
#' This plot contains at least four tracks
#' top: an ideogram of the chromosome
#' second from top: the location on the chromosome with coordinates
#' the following plots have the CNV of each segment
#' and the CCF for each segment
#'
#' If this function is provided with a reference cnvobj
#' an additional plot with the reference CNVs is added below
#'
#' If this function is used to compare two cnvobj objects
#' each panel has overlapping histograms with different colors
#'
#' @param X cloneobj. If a list of cnvobj is provided it is assumed that they will be compared.
#' @param chromosomes list of chromosomes to plot
#' @param filter.neutral filter regions that are NEUT cnv
#' @param ncols the number of columns to arrange the chromosome plots in a table format.
#' @export plot.cloneobj
#' @export
#'
plot.cloneobj=function(X, chromosomes=NULL, filter.neutral=TRUE, ncols=5){
  if( is.list(X) && !"cloneobj" %in% class(X)){
    # check all genomes
    genomes=sapply( X , function(x){x$genome}) %>% unique %>% length
    if(genomes !=1){stop("These objects don't have the same genome")}
    names(X)=sapply(X, function(x){x$name})

  }else{
    X=list( X )
    names(X)=X[[1]]$name
  }
  genome=X[[1]]$genome

  #X=getSegmentFraction(X)

  segments_l=lapply(1:length(X) , function(i){
    segments=GenomeInfoDb::keepStandardChromosomes(X[[i]]$segments, pruning.mode='coarse')
    sl=seqnames(segments) %>% unlist %>% unique() %>% as.character
    segments=segments[ which(segments$cnv.ccf >0)]
    segments
  })

  maxcnv=sapply(1:length(X), function(x){
    .maxCNV( X[[x]]$segments)
  }) %>% max()

  gtrack <- Gviz::GenomeAxisTrack()
  pallette=RColorBrewer::brewer.pal( n=8, name='Dark2')
  ## subclone
  sub_track_l=lapply(1:length(X), function(i){
    gr=segments_l[[i]]
    Gviz::DataTrack(
       gr,
       data=gr$total.allele.copies - X[[i]]$ploidy,
       type="histogram",
       col.histogram=colors()[5*i+2],
       fill.histogram=colors()[5*i],
       baseline=0,
       name = 'CNV diff',
       showTitle=TRUE,
       rotation.title=90,
       cex.title=0.5,
       alpha=0.75,
       alpha.title=1,
       ylim=c(-2,maxcnv-2)
    )
  })
  hist_track_l=lapply(1:length(X), function(i){
    gr=segments_l[[i]]
    if(filter.neutral==TRUE){gr=gr[which(gr$cnv.flag!="NEUT")]}
    Gviz::DataTrack(
      gr,
      data=gr$cnv.ccf,
      name="CCF",
      type="histogram",
      baseline=0,
      showTitle=TRUE,
      rotation.title=90,
      cex.title=0.5,
      col.histogram=colors()[5*i+2],
      fill.histogram=colors()[5*i],
      alpha=0.75,
      alpha.title=1,
      ylim=c(0,1)
    )

  })
  #ref=getReferenceTrack(X)

  #.DEBUG::: Gviz::plotTracks( list( sub_track,hist_track), chromosome='chr1')
  #sl=GenomeInfoDb::sortSeqlevels( as(X$segments,"GRangesList") ) %>% GenomeInfoDb::seqlevels()


  if(!is.null(chromosomes)){
    sl=intersect(sl, chromosomes)
  }

  grid::grid.newpage()
  grid::pushViewport( grid::viewport( layout=grid::grid.layout( ncols,ceiling( length(sl)/ncols ) )))

  for( i in seq_along(sl)){
    grid::pushViewport( grid::viewport( layout.pos.col = ((i - 1) %% ncols) + 1,
                                        layout.pos.row = (((i) - 1) %/% ncols) + 1)  )
    x=as.character(sl[i])
    itrack=Gviz::IdeogramTrack(genome = genome, chromosome = x)
    Gviz::plotTracks(
               c(
                 itrack,
                 gtrack,
                 OverlayTrack(sub_track_l),
                 OverlayTrack(hist_track_l)
               ),
               chromosome = x,
               frame=TRUE,
               add=TRUE,
               sizes=c(2,2,5,5)
    )
    grid::popViewport()

  }
  vcd::grid_legend(x = 'bottomright',
                   just = c(0,0),
                   col = colors()[ 1:length(X) * 5 ],
                   pch=15,
                   labels = names(X),
                   title = "CNV & CCF",
                   gp=gpar(lwd=2, cex=1),
                   hgap = unit(.8, "lines"),
                   vgap = unit(.9, "lines"))
  grid::popViewport()
}


#' helper function to find the maximum number of CNV copies
#'
#' @param X cloneobj
#' @return Returns the maximum copy number among all the clones.
.maxCNV=function(X){
  ll=max( X$total.allele.copies, na.rm=TRUE)

  return(ll)
}



#' plot an integrated plot with CNV information from multiple samples
#'
#' Each plot is per chromosome
#' This plot contains at least four tracks
#' top: an ideogram of the chromosome
#' second from top: the location on the chromosome with coordinates
#' the following plots have the CNV of each segment
#' and the CCF for each segment
#'
#' If this function is provided with a reference cnvobj
#' an additional plot with the reference CNVs is added below
#'
#' If this function is used to compare two cnvobj objects
#' each panel has overlapping histograms with different colors
#'
#' @param X is a list of cnvobj is provided it is assumed that they will be compared.
#' @param chromosomes list of chromosomes to plot
#' @param filter.neutral filter regions that are NEUT cnv
#' @param ncols the number of columns to arrange the chromosome plots in a table format.
#' @export plotCompare
#' @export
#'
plotCompare=function(X, chromosomes=NULL, filter.neutral=TRUE, ncols=5){
    cloneman::plot.cloneobj( X, chromosomes, filter.neutral, ncols=5)
}
#' get the track for a subclone
#'
#' The track containst the segments with the various copy numbers
getSubcloneTrack=function(segments,ploidy=2){
  mx=.maxCNV(segments)
  col=1
  sub_track=list()
  if(is.null(segments)){next}
  jitter=col/20
  sub_track=Gviz::DataTrack( segments,
                              data=segments$total.allele.copies+jitter * runif( 1, min=-1 , max=1) - ploidy,
                              type="histogram",
                              col.histogram=col,
                              fill.histogram="white",
                              baseline=0,
                              name = 'copy number',
                              showTitle=TRUE,
                              rotation.title=0,
                              cex.title=0.5,
                              ylim=c(-2,mx-2))
  col=col+1
  return(sub_track)
}

#' get the CCF of a subclone
#'
#' The track contains the segments and the CCF of each of the segments




#' plot the density of the CCF for an object
#'
#'
#' Depending if SNVs or CNVs are present in the object one can get
#' a density plot with a single line or two lines.
#' @param X cloneobj ,
#' @return an plot object from ggplot
#' @export
plotDensity=function(X){

  if( "segments" %in% names(X)){
    data=X$segments %>%
        as.data.frame() %>%
        dplyr::select('cnv.ccf') %>%
        dplyr::rename('ccf'='cnv.ccf')%>%
        dplyr::mutate( type='cnv')

  }
  if( "snvs" %in% names(X)){
    data=X$segments %>%
      as.data.frame() %>%
      dplyr::select('snv.ccf') %>%
      dplyr::rename('ccf'='snv.ccf')%>%
      dplyr::mutate( type='snv')
  }

  data %>%
    dplyr::filter(ccf>0)%>%
    ggplot2::ggplot( ggplot2::aes(x=ccf,fill=type)) +
    ggplot2::geom_density( bw=0.02, alpha=0.30)  +
    ggplot2::scale_x_continuous(breaks=seq(0,1.05,0.2))
}
