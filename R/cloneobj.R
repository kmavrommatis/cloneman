#' create the object
#'
#' cloneobj is a structure that holds information
#' for the clones
#'
#' This information includes
#'
#' **Method name** the name(s) of the method(s) used to create the data for this object
#'
#' **Purity**: the predicted purity of the sample (0-1)
#'
#' **Ploidy**: the predicted ploidy of the sample (default 2)
#'
#' **Copy number events** GRanges with segments from one or multiple clones (depending on the tool used). Each segment with its own CCF (if available)
#'
#' **SNV clones** GRanges of SNV events, with the SNV VAF (or CCF).
#'
#' **Sample Sex** either provided by the user or inferred'
#'
#' @param purity (float) Purity of the sample ,
#' @param ploidy (float) Ploidy of the sample,
#' @param cnvlist GRanges with cnv events,
#' @param snvlist GRanges with snv events,
#' @param method (list of string) Name of method used to predict these values
#' @param genome the name of the genome (hg19, hg38)
#' @return an object of type cloneobj
#' @export
cloneobj=function( purity=1, ploidy=2, cnvlist=NULL, snvlist=NULL,method='Unspecified', sample_name='sample',genome="hg38"){
  me=list(  purity=purity,
            ploidy=ploidy,
            sex=NULL,
            clone_info=data.frame("CCF"=NA,
                                  "CNV.segments"=NA,
                                  "CNV.chromosomes"=NA,
                                  "SNV.events"=NA,
                                  "SNV.chromosomes"=NA),
            segments=NULL,
            snv=NULL,
            is_finalized=FALSE,
            genome=genome,
            name=sample_name,
            methodName=method)
  class(me)="cloneobj"
  #if(!is.na(purity)){me=setPurity(me,purity)}
  #if(!is.na(ploidy)){me=setPloidy(me,ploidy)}
  # add the snvlist and cnvlist
  if(!is.null( cnvlist)){
    #message("Adding the CNV clones from a ", class(cnvlist), " object")
    me=addEvents(me, cnvlist,  type='cnv' )
  }

  if(!is.null( snvlist)){
    #message("Adding the SNV clones")
    me=addEvents(me, snvlist, type='snv')
  }
  me
}



#' fix the names of the clones in the list.
#' if a name vector is provided use that as the names of the list
#' else use the names that are assigned to the objects already
#' and if no names are assigned use the format 'cloneX'
#' @param gr GRangesList or GRanges object that will be assigned the names
#' @param clone_names a vector with names to replace the names of the elements
#' @return a vector with the names of the elements of the list
#'
.fixNames=function( gr, clone_names=NULL){
  # the input vector should have the length of the list
  if(!is.null( clone_names) & length(clone_names) != length(gr)){
    stop("Please provide a vector with names equal in size to the GRangesList object. Currently the names has size of ",length(names), " and the Granges object ", length(gr))
  }

  if("SimpleGRangesList" %in% class(gr) | "CompressedGRangesList" %in% class(gr)){
    # make the names if they are missing
    for(i in 1:length(gr)){
      if(is.null(names(gr)[i])){
        clone_names[i]=paste0("clone",i)
      }else{
        clone_names[i]=names(gr)[i]
      }
    }
  }
  if("GRanges" %in% class(gr)){
    if( is.null(clone_names) ){
      clone_names=paste0("clone", nrow(X$clone_info)+1)
    }
  }

  return(clone_names)
}


#' perform final tasks with an object
#'
#' @param X cloneobj
# .finalize=function(X){
#
#   if( X$is_finalized==TRUE){return(X)}
#
#   # convert the lists to Granges Lists
#   # if( "list" %in% class(X$segments)){ X$segments= as(X$segments,"GRangesList")}
#   # if( "list" %in% class(X$snvs)){ X$snvs= as(X$snvs,"GRangesList")}
#   # infer sex if necessary
#   if(is.null(X$sex)){ X$sex=cloneman:::inferSex( c(X$segments,X$snvs) %>% as("GRangesList") )}
#   # add the seq info
#   if( !is.null(X$segments)){
#     GenomeInfoDb::seqlevels( X$segments)=GenomeInfoDb::seqlevels( GenomeInfoDb::Seqinfo(genome=X$genome))
#     GenomeInfoDb::seqinfo( X$segments)=GenomeInfoDb::Seqinfo( genome=X$genome)
#     }
#   if(!is.null(X$snvs)){
#     GenomeInfoDb::seqlevels( X$snvs)=GenomeInfoDb::seqlevels( GenomeInfoDb::Seqinfo(genome=X$genome))
#     GenomeInfoDb::seqinfo( X$snvs)=GenomeInfoDb::Seqinfo( genome=X$genome)
#     }
#
#
#   X$is_finalized=TRUE
#   X
# }

#' add events (CNV and SNV) to the object
#' it uses a GRangesList as input, which contains all the SNV events for as sample
#' it expects a column named 'snv.ccf' which contains the CCF of the clone,
#' or a column names 'snv.vaf' with the VAF of the clone (where CCF=VAF *2 for most non amplified regions of the genome )
#' The overall CCF of this clone is the mean of this value.
#' For CNVs it expecst a GRanges object with all the CNV events.
#' with a column cnv.ccf with the ccf of each segment
#'
#' @param X cloneobj
#' @param events a GRanges object. If a list is provided the names of the list should correspond to the names of the clones
#' if no names are provided then they are created followign the pattern "cloneX" where X is 1, 2, ..
#' @param type either 'snv' or 'cnv' depending on the type of events that are added
#' @export
addEvents=function(X, events , type=c('snv','cnv') ){
  #clone_names=.fixNames( events, name)
  #message("inside addClones with object events of class ", class(events))
  if( ! "GRanges" %in% class(events)){
    stop("SVNs and CNVs should be a GRanges object. Instead they are ", class( events ))
  }
  if("GRanges" %in% class(events)){
      X=addClone( X, events, type=type)
  }
 # .finalize(X)

}

#' add a single CNV or SNV set of segments in the object
#'
#' @param X cloneobj
#' @param gr GRanges withe the segments of a specific clone. The column 'snv.ccf' or 'cnv.ccf' of this Granges object is used to get the CCF of this clone.
#' Alternatively the column 'snv.vaf' can be used only if CCF is not present for mutation objects
#' The CCF of the CNV segments are adjusted for purity, i.e. we estimate how much it would be if the sample was 100% pure, or - in other words - what is the CCF for the _tumor_ part only.
#' When segments are added the gaps of the chromosomes are filled as well
#'
#' @param type either 'snv' or 'cnv' depending on the type of events that are added
#' @return an updated cloneobj
addClone=function(X, gr,type=c("snv","cnv")){
  if( ! "GRanges" %in% class(gr) ){stop("Each clone must be a GRanges object")}
  if(type=='cnv'){
    if(! "cnv.ccf" %in% colnames(S4Vectors::mcols( gr ))){ stop("For a CNV clone the column 'cnv.ccf' must be present")}
    if(intersect( c("total.allele.copies","major.allele.copies", "minor.allele.copies","cnv.flag"),
                  colnames(S4Vectors::mcols( gr )))%>%length() !=4){
      stop("For a CNV clone the columns 'total.allele.copies','major.allele.copies','minor.allele.copies','cnv.flag' must be present")}

    # add the adjusted cnv.ccf
    S4Vectors::mcols( gr )$adj.cnv.ccf=S4Vectors::mcols( gr )$cnv.ccf/cloneman::getPurity(X)

    # add gaps to the segments
    grg=GenomicRanges::gaps(gr)

    grg=grg[ which(BiocGenerics::strand(grg)=='*'),]
    S4Vectors::mcols(grg)['cnv.ccf']=0
    S4Vectors::mcols(grg)['total.allele.copies']=2
    S4Vectors::mcols(grg)['cnv.flag']="NEUT"
    X$segments=c(grg, gr) %>% GenomeInfoDb::sortSeqlevels() %>% sort()
  }
  if(type=='snv'){
    if(intersect( c("snv.ccf", "snv.vaf"),
                  colnames(S4Vectors::mcols( gr )))%>%length() !=1){
      stop("For a SNV clone the column 'snv.ccf' OR the columne 'snv.vaf' must be present")
      }
    if('snv.vaf' %in% colnames( S4Vectors::mcols( gr )) && !'snv.ccf' %in% colnames( S4Vectors::mcols( gr )) ){ gr$snv.ccf=gr$snv.vaf * 2 }
    if('snv.ccf' %in% colnames( S4Vectors::mcols( gr )) && !'snv.vaf' %in% colnames( S4Vectors::mcols( gr )) ){ gr$snv.vaf=gr$snv.ccf / 2 }
    X$snvs=gr
  }
  #X$is_finalized=FALSE
  X
}





#' Print information for the object
#'
#' Provide a table with the clonal information of the object
#'
#' The table contains:
#'
#' known clones
#'
#' purity of the sample
#'
#' ploidy of the sample
#'
#' @export print.cloneobj
#' @export
print.cloneobj=function(x, ...){
  pl=getPloidy(x)
  if(is.null(x$method)){method="Unspecified"}else{method=x$method}
  cat("Method ", method , "\n")
  cat("Predicted ploidy (rounded)   ", pl$ploidy, "\n")
  cat("Predicted ploidy", pl$real.ploidy,"\n")
  cat("Predicted purity  ", x$purity,"\n")
  if(!is.null(x$segments)){
    cat("Contains ", length(x$segments)," CNV segments \n")
  }
  if(!is.null(x$snvs)){
    cat("Contains ", length(x$snvs)," SNV events \n")
  }

}



#' Override the GenomeInfoDb seqlevelstyle for cloneobj
#' @param X the cloneobj
#' @param style specify the style
#' @export 'seqlevelsStyle<-.cloneobj'
#' @export
'seqlevelsStyle<-' = function(X){
  UseMethod( 'seqlevelsStyle<-')
}

#' @export
'seqlevelsStyle<-.cloneobj' = function(X, style){
  if( !is.null(X$segments)){  GenomeInfoDb::seqlevelsStyle(X$segments)=style }
  if( !is.null(X$segments)){  GenomeInfoDb::seqlevelsStyle(X$snvs)=style }
}




#' Override the GenomeInfoDb seqlevelstyle for cloneobj
#' @param X the cloneobj
#' @return current style in the object
#' @export seqlevelsStyle.cloneobj
#' @export
seqlevelsStyle=function(X){
  UseMethod("seqlevelsStyle")
}

#' @export
seqlevelsStyle.cloneobj = function(X){
  s1=NULL
  s2=NULL
  if( !is.null(X$segments)){  s1=GenomeInfoDb::seqlevelsStyle(X$segments) ; res=s1}
  if( !is.null(X$snvs)){  s2=GenomeInfoDb::seqlevelsStyle(X$snvs) ; res=s2}
  if( !is.null( s1 ) && !is.null(s2) && s1 != s2){
    stop("SNVs and CNV segments in this object don't have the same sequence style. Please fix it and continue.")
  }

  return(
    res
  )
}



#' Filter regions
#'
#' Provide a list of regions that will be used to fragment the CNV and SVN of the object
#' This is useful to compare objects between them
#' and restrict the results on specific regions of the genome (e.g. the coordinates of an exome kit)
#'
#' @param X cloneobj to filter
#' @param gr GRanges object with the coordinates to filter
#' @return a new object with the filtered coordinates.
#' @export filter.cloneobj
#' @export
filter=function(X){
  UseMethod( "filter")
}

#' @export
filter.cloneobj=function(X, gr){
  for( x in c('segments', 'snvs')){
      if(is.null(X[[x]])){ next }
      ov=IRanges::findOverlapPairs(query=X[[x]], subject=gr)
      newlabel=paste( x, 'filtered',sep='.')
      X[[newlabel]]=pintersect(ov)
  }
  X[['filter']]=gr
  X
}



#' Set the name of the sample in a cloneobj
#'
#' The name is a string that can be used to identify the sample
#'
#' @param X cloneobj
#' @param sample_name the name of the object
#' @return an updated cloneobj
setName=function( X, sample_name){
  X$name=sample_name
  X
}

#' Get the name of the sample in a cloneobj
#'
#' The name is a string that can be used to identify the sample
#'
#' @param X cloneobj
#' @return the name of the object
getName=function( X){
  X$name
}

