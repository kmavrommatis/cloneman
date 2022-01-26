

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

  return( parseCloneCNA( clonecna.gr ) )
}


#' Generate a cnv Granges from a dataframe created by ControlFreec (file ending in .CNVs or contains pvalues.)
#' CloneCNA claims that can efficiently and accurately identify somatic copy number alterations from
#' heterogeneous tumor samples. (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1174-7)
#' CloneCNA can detect teh cellularity of each abberation, and from that extract the subclonal abundance. T
#' his comes in the column Cellularity.
#' When Cellularity is 0 it indicates normal
#' Presumably, segments with different values in that column correspond to different clones. As a result we can tag different clones based on these frequencies.
#' @param clonecna.df input filename with CNV segments,
#' @return A GRanges object with the segments predicted by CloneCNA
#' @export
parseCloneCNA=function(
  clonecna.df
){
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

  flags=apply(mcols_new, 1, function(K){ cnvFlag(total.cn = K['total.allele.copies'], minor.cn = K['minor.allele.copies'] )}) # pass the total.cn to make sure that cases with uncertainty=-1 are parsed properly
  mcols_new = mcols_new %>%
    dplyr::mutate( cnv.flag=flags)
  S4Vectors::mcols( clonecna.gr)=mcols_new
  return( clonecna.gr )
}
