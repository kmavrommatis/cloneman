

#' CloneCNA claims that can efficiently and accurately identify somatic copy number alterations from heterogeneous tumor samples. (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1174-7)
#' CloneCNA can detect teh cellularity of each abberation, and from that extract the subclonal abundance. This comes in the column Cellularity.
#' When Cellularity is 0 it indicates normal
#' Presumably, segments with different values in that column correspond to different clones. As a result we can tag different clones based on these frequencies.
#' @param clonecna.fn input filename with CNV segments,
#' @return A GRanges object with the segments predicted by CloneCNA
#' @export
parseCloneCNA_file=function(
  clonecna.fn
){
  skip=grep('Chr', readLines(clonecna.fn)) -1
  clonecna=read.table( clonecna.fn,header=TRUE , skip=skip ,sep="\t")
  # keep only the ones that are not normal

  return( clonecna )
}


#' Generate a cnv Granges from a dataframe created by CloneCNA
#'
#' CloneCNA claims that can efficiently and accurately identify somatic copy number alterations from
#' heterogeneous tumor samples. (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1174-7)
#' CloneCNA can detect teh cellularity of each abberation, and from that extract the subclonal abundance. T
#' his comes in the column Cellularity.
#' When Cellularity is 0 it indicates normal
#' Presumably, segments with different values in that column correspond to different clones. As a result we can tag different clones based on these frequencies.
#' @param clonecna input filename with CNV segments or an already parsed dataframe,
#' @return A GRanges object with the segments predicted by CloneCNA
#' @export
parseCloneCNA=function(
  clonecna
){

  if( file.exists(clonecna) ){
    clonecna.df=parseCloneCNA_file( clonecna)
  }else if( is.data.frame( clonecna)){
    clonecna.df=clonecna
  }else{
    stop("Cannot recognize the input as either a file or a dataframe")
  }
  clonecna.gr=GenomicRanges::makeGRangesFromDataFrame( clonecna.df,
                                                       keep.extra.columns = TRUE,
                                                       ignore.strand = TRUE,
                                                       start.field = 'StartPos',
                                                       end.field = 'EndPos')%>%
    GenomeInfoDb::sortSeqlevels( ) %>%
    BiocGenerics::sort()
  GenomeInfoDb::seqlevelsStyle(clonecna.gr)='UCSC'

  #get the info from the file
  mcols_new=S4Vectors::mcols( clonecna.gr ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate( total.allele.copies=CN,
                   major.allele.copies=mCN,
                   minor.allele.copies=CN-mCN)

  flags=apply(mcols_new, 1, function(K){ cloneman:::cnvFlag(total.cn = K['total.allele.copies'], minor.cn = K['minor.allele.copies'] )}) # pass the total.cn to make sure that cases with uncertainty=-1 are parsed properly
  mcols_new = mcols_new %>%
    dplyr::mutate( cnv.flag=flags)
  S4Vectors::mcols( clonecna.gr)=mcols_new
  return( clonecna.gr )
}

#' Generate a CNV object with the segments from facets and additional information
#'
#' We keep all the events the program reports
#'
#'
#' @import mclust
#' @param clonecna    FACETS object (as provided by facets) or a file in .rds format with a FACETS object. This object is typically saved as "fit" and contains "fit$cncf", "fit$purity" and "fit$ploidy"
#' @param sample_name name of the sample to be used in the object
#' @export
loadCloneCNA=function(
    clonecna, sample_name=NULL
){
  if( file.exists(clonecna) ){
    clonecna.df=parseCloneCNA_file( clonecna)
  }else{
    stop("Cannot recognize the input as  a file")
  }
  segments=parseCloneCNA( clonecna )
  S4Vectors::mcols(segments)$cnv.ccf=S4Vectors::mcols(segments)$Cellularity

  keep=grep('---', readLines(clonecna)) ; keep=keep[2]-length(keep)
  header=read.table( clonecna, nrows=keep ,header=FALSE,comment.char='-' , sep=":", fill=TRUE)

  purity=header[ grep( 'Cellularity', header[,1]), 2]
  ploidy=header[ grep( 'Average copy number', header[,1]),2]
  obj=cloneobj(purity = purity,
               ploidy = ploidy,
               cnvlist =  segments,
               method='CloneCNA',
               sample_name=sample_name)
  obj

}
