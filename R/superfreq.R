#' Main function to parse the SNV segments and SNV events from  SuperFreq.
#'
#' Superfreq produces a variety of files, for the normal (if present) and tumor samples.
#' file data/CNAsegments{samplename}.tsv has a list of the CNV segments with their specific information
#' file rivers/{samplename}-river.tsv contains the list of the mutations & CNV segments assigned to clones
#' file data/clones_{samplename}.tsv contains the list of the clones with their clonality
#'
#' Our main file where we start all analyiss from is the rivers/{samplename}-river.tsv
#' and we use it as the scaffold, to add more information from the CNA segments and variants files.
#' main function for Superfreq
#' Create necessary objects from the superfreq directory
#' @param method_dir directory created by superfreq. Under it are expected the directories data and rivers
#'
#'
parseSuperFreq=function(method_dir){
  ll=.SuperFreqFilenames(method_dir)
  clone_info=.parseSuperFreq_info( ll$clones)
  subclone_events.gr=.parseSuperFreq_River( ll$subclone_events)
  CCF= clone_info[ match( subclone_events.gr$clone, rownames(clone_info)),
                  intersect( grep("stories", colnames( clone_info)), grep("normal",colnames(clone_info),invert = TRUE) )]
  subclone_events.gr$CCF=CCF

  superfreq_cnv.gr=.parseSuperFreq_CNV( subclone_events.gr,  ll$cnasegments)
  superfreq_snv.gr=.parseSuperFreq_SNV( subclone_events.gr,  ll$snvevents)
  return(list(
    "CNV"=superfreq_cnv.gr,
    "SNV"=superfreq_snv.gr,
    "clone_info"=clone_info
  ))
}

#' Create a cloneobj from SuperFreq
#'
#' @param method_dir directory created by superfreq. Under it are expected the directories data and rivers
#' @param method_name name of the method, defaults to SuperFreq
#' @param purity purity of sample, overrides the purity detected by the method, but does not change the CCF of clones
#' @param ploidy ploidy of sample, default is 2
#' @export
loadSuperFreq=function(method_dir, method_name='SuperFreq', purity=NULL, ploidy=2){
  ll=parseSuperFreq( method_dir)

  if(is.null(purity)){

    lab=  colnames(ll$clone_info)[grep(pattern = "stories", colnames(ll$clone_info))]
    lab=lab[ grep(pattern="normal", lab, invert=TRUE)]

    purity= ll$clone_info[ which(rownames(ll$clone_info) != "germline"),lab] %>% max(na.rm=TRUE)
  }
  sf=cloneobj(purity=purity,
              ploidy=ploidy,
              cnvlist=as(split( ll$CNV, ll$CNV$clone), "GRangesList"),
              snvlist=as(split( ll$SNV, ll$SNV$clone), "GRangesList"),
              method=method_name)
  sf
}


#' Pasre the CNV segments from SuperFreq
#' @param method_dir directory created by superfreq. Under it are expected the directories data and rivers
#' @return GRanges object with all the CNV segments
#' @export
parseSuperFreq_CNV=function(method_dir){
  a=parseSuperFreq(method_dir)
  return(a$CNV)
}

#' Parse the SNV events from SuperFreq
#' @param method_dir directory created by superfreq. Under it are expected the directories data and rivers
#' @return GRanges object with all the CNV segments
#' @export
parseSuperFreq_SNV=function(method_dir){
  a=parseSuperFreq(method_dir)
  return(a$SNV)
}


#' Search the directory of SuperFreq and get the list of files that are needed for hte analysis
#'
#' @param method_dir directory created by superfreq. Under it are expected the directories data and rivers
#' @return a list with the filenames of the filetypes needed
.SuperFreqFilenames=function(method_dir){
  if(! dir.exists( file.path(method_dir,"data")) |
     ! dir.exists( file.path(method_dir,"rivers"))){
    stop("Please provide the directory created by 'SuperFreq' with subdirectories data and rivers.")
  }
  sample_name=.SuperFreqSampleName( method_dir )
  subclone_events_info=list.files(file.path(method_dir, "rivers"), pattern=paste0(sample_name,".?-river.tsv"),full.names = TRUE)
  superfreq_info= list.files( file.path(method_dir,"data"), pattern=paste0("clones_.?",sample_name,".?.tsv"), full.names = TRUE)

  cnasegments=list.files( file.path(method_dir,"data"), pattern=paste0("CNAsegments_TP.*.tsv"), full.names=TRUE)

  snvevents=file.path(method_dir,"somaticVariants.csv")
  if( is.null(subclone_events_info) |
      is.null(superfreq_info) |
      is.null(cnasegments) |
      !file.exists( snvevents)){
    stop("Could not identify all the necessary files for Superfreq, is the ", method_dir, " the correct location of a Superfreq run output?\n")

  }

  return(list(
    "subclone_events"=subclone_events_info[1],
    "clones"=superfreq_info,
    "cnasegments"=cnasegments,
    "snvevents"=snvevents
  ))
}

.SuperFreqSampleName=function(method_dir){
  ff=list.files(file.path(method_dir, "data"), full.names = FALSE, pattern="clones_.*.tsv")
  ff=gsub("clones_","", ff)
  ff=gsub(".tsv","", ff)
  return(ff)
}


#' Parse the rivers/{samplename}-river.tsv file
#'
#' which contains the summary of CNV and SNV and their assignment to clones
#' @param fn the {samplename}-river.tsv file
#' @return a GRanges object with the events.
.parseSuperFreq_River=function( subclone_events_info){
  superfreq=read.table( subclone_events_info, header=TRUE, sep="\t") %>% tidyr::drop_na( c("start","end"))
  superfreq.gr=GenomicRanges::makeGRangesFromDataFrame( superfreq,
                                                        keep.extra.columns = TRUE,
                                                        ignore.strand = TRUE
  ) %>%
    GenomeInfoDb::sortSeqlevels( ) %>%
    BiocGenerics::sort()
  GenomeInfoDb::seqlevelsStyle(superfreq.gr)='UCSC'

  return( superfreq.gr)
}





#' Parse the info file that Superfreq produces
#'
#' This a file following the pattern data/clones_{samplename}.tsv
#' @param clone_info the file with the clone info (data/clones_{samplename}.tsv)
#' @return a dataframe with the clones information.
.parseSuperFreq_info=function( clone_info){

  info=read.table(clone_info, header=TRUE, row.names = 1 )
  return(info)
}



#' Parse the CNV segments from SuperFreq
#'
#' Generate a cnv Granges from a dataframe created by SuperFreq
#' it matches the CNVs in the suprefreq.gr to the CNVs listed in the CNA segments file
#' and transfers the information from the latter to the former
#' it also adds the clonal information of the clones

#' @param superfreq.gr Grages with SuperFreq output (loaded by parseSuperFreq_River)
#' @param clone_info dataframe with the information about the clones
#' @param cnasegments.fn file with the details of the CNV prediction
#' @return A GRanges object with the segments predicted by SuperFreq, with clonal ifnormation added
.parseSuperFreq_CNV=function( superfreq.gr, cnasegments.fn){

  cnvs=grep( "[KkMG]bp",superfreq.gr$name)
  superfreq_cnv.gr=superfreq.gr[ cnvs ]
  # Add the information from the CNVsegments file
  ss=read.table( cnasegments.fn , header=TRUE, sep="\t") %>%
    tidyr::drop_na(c("start","end")) %>%
    GenomicRanges::makeGRangesFromDataFrame( keep.extra.columns = TRUE, ignore.strand = TRUE) %>%
    GenomeInfoDb::sortSeqlevels( ) %>%
    BiocGenerics::sort()
  GenomeInfoDb::seqlevelsStyle(ss)='UCSC'

  cov=GenomicRanges::findOverlaps( subject=superfreq_cnv.gr, query=ss)
  cbind( S4Vectors::mcols( superfreq_cnv.gr)[ S4Vectors::subjectHits(cov),],
         S4Vectors::mcols( ss )[ S4Vectors::queryHits(cov),])
  ss_a=superfreq_cnv.gr[ S4Vectors::subjectHits(cov),]
  mc=cbind( S4Vectors::mcols( superfreq_cnv.gr)[ S4Vectors::subjectHits(cov),],
                    S4Vectors::mcols( ss )[ S4Vectors::queryHits(cov),])

  S4Vectors::mcols(ss_a)=mc

  # ss_a now is a Granges object that contains teh SNV detailed information with the clonality
  # Add the flags
  mcols_new=S4Vectors::mcols( ss_a ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate( total.allele.copies=stringr::str_length( call),
                 major.allele.copies=stringr::str_count( call, 'A'),
                 minor.allele.copies=stringr::str_count( call, 'B')
    )

  flags=apply(mcols_new, 1, function(K){ cnvFlag(total.cn = K['total.allele.copies'], minor.cn = K['minor.allele.copies'] )}) # pass the total.cn to make sure that cases with uncertainty=-1 are parsed properly
  mcols_new = mcols_new %>%
    dplyr::mutate( cnv.flag=flags)

  S4Vectors::mcols( ss_a )=mcols_new

  return(ss_a)
}




#' Extract the SNV events from the clones from SuperFreq
#'
#' @param superfreq.gr Grages with SuperFreq output (loaded by parseSuperFreq_River)
#' @param clone_info dataframe with the information about the clones
#' @param variants.fn file with the details of the CNV prediction
#' @return A GRanges object with the SNVs predicted by SuperFreq, with clonal information added (somaticVariants.tsv)
.parseSuperFreq_SNV=function(  superfreq.gr,variants.fn){


  snvs=grep( "[KkMG]bp",superfreq.gr$name, invert=TRUE)
  superfreq_snv.gr=superfreq.gr[ snvs ]
  # Add the information from the SNV somaticVariants file
  ss=read.table( variants.fn , header=TRUE, sep=",") %>%
    tidyr::drop_na(c("start","end")) %>%
    GenomicRanges::makeGRangesFromDataFrame( keep.extra.columns = TRUE, ignore.strand = TRUE) %>%
    GenomeInfoDb::sortSeqlevels( ) %>%
    BiocGenerics::sort()
  GenomeInfoDb::seqlevelsStyle(ss)='UCSC'
  keep=grep("normal", ss$sample, invert=TRUE)
  ss=ss[keep]


  cov=GenomicRanges::findOverlaps( subject=superfreq_snv.gr, query=ss)

  ss_a=superfreq_snv.gr[ S4Vectors::subjectHits(cov),]
  mc=cbind( S4Vectors::mcols( superfreq_snv.gr)[ S4Vectors::subjectHits(cov),],
            S4Vectors::mcols( ss )[ S4Vectors::queryHits(cov),])
  S4Vectors::mcols(ss_a)=mc



  return(ss_a)
}


