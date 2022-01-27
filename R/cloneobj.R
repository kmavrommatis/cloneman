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
#' **Copy number clones** GRangesList with segments from one or multiple clones (depending on the tool used)
#'
#' **Sample Sex** either provided by the user or inferred
#'
#' **SNV clones** GRangesList of clones defined by the SNV VAF (or CCF)
#'
#'
#' @param purity (float) Purity of the sample ,
#' @param ploidy (float) Ploidy of the sample,
#' @param cnvlist GRangesList with cnv clones,
#' @param snvlist GRangesList with snv clones,
#' @param method (list of string) Name of method used to predict these values
#' @param genome the name of the genome (hg19, hg38)
#' @return an object of type cloneobj
#' @export
cloneobj=function( purity=NA, ploidy=NA, cnvlist=NULL, snvlist=NULL,method=NULL, genome="hg38"){
  me=list(  purity=numeric(),
            ploidy=numeric(),
            sex=NULL,
            clone_info=data.frame("clone.name"="germline",
                                  "CCF"=NA,
                                  "CNV.segments"=NA, "CNV.chromosomes"=NA,
                                  "SNV.events"=NA),
            segments_list=NULL,
            snv_list=NULL,
            is_finalized=FALSE,
            genome=genome,
            methodName=method)
  class(me)="cloneobj"
  if(!is.na(purity)){me=setPurity(me,purity)}
  if(!is.na(ploidy)){me=setPloidy(me,ploidy)}
  me$clone_info[1, 'CCF']=1-me$purity
  # add the snvlist and cnvlist
  if(!is.null( cnvlist)){

    #message("Adding the CNV clones from a ", class(cnvlist), " object")
    me=addClones(me, cnvlist, names(cnvlist), type='cnv' )
  }

  if(!is.null( snvlist)){
    #message("Adding the SNV clones")
    me=addClones(me, snvlist, names(snvlist), type='snv')
  }

  me$reference=NULL
  me$problems=NULL
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
.finalize=function(X){

  if( X$is_finalized==TRUE){return(X)}

  # convert the lists to Granges Lists
  if( "list" %in% class(X$segments_list)){ X$segments_list= as(X$segments_list,"GRangesList")}
  if( "list" %in% class(X$snv_list)){ X$snv_list= as(X$snv_list,"GRangesList")}
  # infer sex if necessary
  if(is.null(X$sex)){ X$sex=inferSex( c(X$segments_list,X$snv_list) %>% as("GRangesList") )}
  # add the seq info
  if( !is.null(X$segments_list)){GenomeInfoDb::seqinfo( X$segments_list)=GenomeInfoDb::Seqinfo( genome=X$genome)}
  if(!is.null(X$snv_list)){GenomeInfoDb::seqinfo( X$snv_list)=GenomeInfoDb::Seqinfo( genome=X$genome)}
  # add gaps to the segments list
  for( i in X$clone_info$clone.name){
    if( is.null(X$segments_list[[i]])){next}
    g=GenomicRanges::gaps(X$segments_list[[i]])
    S4Vectors::mcols(g)['CCF']=0
    S4Vectors::mcols(g)['total.allele.copies']=2
    S4Vectors::mcols(g)['cnv.flag']="NEUT"
    X$segments_list[[i]]=c( X$segments_list[[i]] ,g) %>% GenomeInfoDb::sortSeqlevels() %>% sort()
  }


  # some checks


  X$is_finalized=TRUE
  X
}

#' add clones to the object
#' it uses a GRangesList as input, which contains all teh CNV segments for each clone
#' it expects a column named 'cnv.ccf' which contains the CCF of the clone.
#' The overall CCF of this clone is the mean of this value.
#'
#' @param X cloneobj
#' @param cnv GRangesList or GRanges object. If a list is provided the names of the list should correspond to the names of the clones
#' if no names are provided then they are created followign the pattern "cloneX" where X is 1, 2, ..
#' @param names a vector with names that will be used to name the GRanges objects. Should have the same length as the list of GRanges.
#' This value overrides the names provided as names of the list.
#' @export
addCNVClones=function(X, cnv , name=NULL){
  addClones( X, cnv, name, type='cnv')
}



#' add clones to the object
#' it uses a GRangesList as input, which contains all teh SNV events for each clone
#' it expects a column named 'snv.ccf' which contains the CCF of the clone.
#' or a column names 'snv.vaf' with the VAF of the clone (where CCF=VAF *2 for most non amplified regions of the genome )
#' The overall CCF of this clone is the mean of this value.
#'
#' @param X cloneobj
#' @param snv GRangesList or GRanges object. If a list is provided the names of the list should correspond to the names of the clones
#' if no names are provided then they are created followign the pattern "cloneX" where X is 1, 2, ..
#' @param names a vector with names that will be used to name the GRanges objects. Should have the same length as the list of GRanges.
#' This value overrides the names provided as names of the list.
#' @export
addSNVClones=function(X, snv , name=NULL){
  addClones( X, snv, name, type='snv')
}



#' add clones to the object
#' it uses a GRangesList as input, which contains all teh SNV events for each clone
#' it expects a column named 'snv.ccf' which contains the CCF of the clone.
#' or a column names 'snv.vaf' with the VAF of the clone (where CCF=VAF *2 for most non amplified regions of the genome )
#' The overall CCF of this clone is the mean of this value.
#'
#' @param X cloneobj
#' @param events GRangesList or GRanges object. If a list is provided the names of the list should correspond to the names of the clones
#' if no names are provided then they are created followign the pattern "cloneX" where X is 1, 2, ..
#' @param names a vector with names that will be used to name the GRanges objects. Should have the same length as the list of GRanges.
#' This value overrides the names provided as names of the list.
#' @param type either 'snv' or 'cnv' depending on the type of events that are added
#' @export
addClones=function(X, events , name=NULL , type=c('snv','cnv') ){
  clone_names=.fixNames( events, name)
  #message("inside addClones with object events of class ", class(events))
  if("SimpleGRangesList" %in% class(events) |  "CompressedGRangesList" %in% class(events)){
    names(events)=clone_names
    for(i in names(events)){
      if(type == 'snv'){
        X=addSNVClone( X, events[[i]], clonename=i)
      }
      if(type == 'cnv'){
        X=addCNVClone( X, events[[i]], clonename=i)
      }
    }
  }

  if("GRanges" %in% class(events)){
    if( is.null(name) ){
      name=paste0("clone", nrow(X$clone_info)+1)
    }
    if(type == 'snv'){
      X=addSNVClone( X, events, clonename=name)
    }
    if(type == 'cnv'){
      X=addCNVClone( X, events, clonename=name)
    }
  }
  .finalize(X)

}

#' add a single CNV clone in the object
#'
#' @param X cloneobj
#' @param gr GRanges withe the segments of a specific clone. The column 'CCF' of this Granges object is used to get the CCF of this clone.
#' @param clonename the name of the clone that will be added in teh list of clones.
addCNVClone=function(X, gr, clonename=NA){
  addClone( X, gr, clonename, type='cnv')
}
#' add a single SNV clone in the object
#'
#' @param X cloneobj
#' @param gr GRanges withe the segments of a specific clone. The column 'CCF' of this Granges object is used to get the CCF of this clone.
#' @param clonename the name of the clone that will be added in teh list of clones.
addSNVClone=function(X, gr, clonename=NA){
  addClone( X, gr, clonename, type='snv')
}

#' add a single CNV clone in the object
#'
#' @param X cloneobj
#' @param gr GRanges withe the segments of a specific clone. The column 'CCF' of this Granges object is used to get the CCF of this clone. Alternatively the column 'VAF' can be used only if CCF is not present
#' @param clonename the name of the clone that will be added in teh list of clones.
#' @param type either 'snv' or 'cnv' depending on the type of events that are added
#' @return an updated cloneobj
addClone=function(X, gr, clonename=NA,type=c("snv","cnv")){
  if( ! "GRanges" %in% class(gr) ){stop("Each clone must be a GRanges object")}
  if(is.na(clonename)){ stop("Please provide a name for this clone")}
  gr$clone=clonename
  if(type=='cnv'){
    if(! "CCF" %in% colnames(S4Vectors::mcols( gr ))){ stop("For a CNV clone the column 'CCF' must be present")}
    if(intersect( c("total.allele.copies","major.allele.copies", "minor.allele.copies","cnv.flag"), colnames(S4Vectors::mcols( gr )))%>%length() !=4){ stop("For a CNV clone the columns 'total.allele.copies','major.allele.copies','minor.allele.copies','cnv.flag' must be present")}
    ccf=gr$CCF %>%mean()
    X$segments_list[[clonename]]=gr
  }
  if(type=='snv'){
    if(intersect( c("CCF", "VAF"), colnames(S4Vectors::mcols( gr )))%>%length() !=1){ stop("For a SNV clone the column 'CCF' OR the columne 'VAF' must be present")}
    if('VAF' %in% colnames( S4Vectors::mcols( gr ))){vaf=gr$VAF %>%mean(); ccf=vaf * 2}
    if('CCF' %in% colnames( S4Vectors::mcols( gr ))){ccf=gr$CCF %>%mean()}
    X$snv_list[[clonename]]=gr
  }
  X=addCloneEntry(X, name = clonename, ccf)
  X$is_finalized=FALSE
  X
}

#' add an entry for a new clone in the clone_info data frame
#'
#' @param X cloneobj
#' @param name the name of the clone
#' @param fraction the ccf of the clone
#' @return an updated cloneobj
addCloneEntry=function( X , name, fraction){
  # add an entry with information in clone_info
  name=as.character(name)
  fraction=as.numeric(fraction)
  cnv.count=NA
  cnv.chr.count=NA
  snv.count=NA
  if(name %in% names(X$segments_list)){
    cnv.count=length( X$segments_list[[name]])
    cnv.chr.count=length(GenomeInfoDb::seqnames( X$segments_list[[name]]) %>% unique)
  }
  if(name %in% names(X$snv_list)){
    snv.count=length( X$snv_list[[name]])
  }
  # message("Received name ", name, " with fraction ", fraction,", cnv count ", cnv.count, " and snv count ", snv.count)
  if(fraction<0 || fraction >1){stop("The CCF for a clone has to be between 0 and 1")}


  if( name %in% X$clone_info$clone.name){
    #message("Updating data for known clone ", name, " which already exists in ", X$clone_info$clone.name)
    X$clone_info[ X$clone_info$clone.name == name , c('clone.name','CCF','CNV.segments','CNV.chromosomes','SNV.events')]=c(name, fraction, cnv.count, cnv.chr.count, snv.count)
  }else{
    #message("Introducing a new clone ", name)
    X$clone_info= rbind( X$clone_info, c( "clone.name"=name, "CCF"=fraction, "CNV.segments"=cnv.count, "CNV.chromosomes"=cnv.chr.count, "SNV.events"=snv.count) )
  }
  X$clone_info=X$clone_info[ order( X$clone_info$CCF,decreasing = TRUE ), ]
  X$clone_info=unique(X$clone_info)
  X$is_finalized=FALSE
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
  cat("Predicted ploidy (rounded)   ", pl$ploidy, "\n")
  cat("Predicted ploidy", pl$real.ploidy,"\n")
  cat("Predicted purity  ", x$purity,"\n")
  if(!is.null(x$segments_list)){
    cat("Contains ", length(x$segments_list)," CNV GRanges for clones ", names(x$segments_list),"\n")
  }
  if(!is.null(x$snv_list)){
    cat("Contains ", length(x$snv_list)," SNV GRanges for clones ", names(x$snv_list),"\n")
  }
  cat("Clone information \n");

  print(x$clone_info); cat("\n")

}





#' total SNV events
#'
#' Merges the snv list and sums the VAF (or CCF) of events found in multiple clones
#' In the end the result contains the GRanges of the events
#' a column called 'clones' with the number of clones that have this event
#' and the overall VAF of the event
#'
#' @param X cnvobj
#' @return a GRanges with all the events in the SNV list
#' @export
getCombinedSNV=function(X){
  if(length(X$snv_list)==0){return(NULL)}
  #get a list of all the svn locations
  all_snvs=GenomicRanges::GRangesList( X$snv_list) %>% unlist()  %>% GenomicRanges::granges()
  #mcols(all_snvs)=mcols(all_snvs)[,c("reference","alternate")]
  all_snvs=all_snvs %>% unique()

  #all_snvs$constitutional.clones=0
  # now check which ones map to a clone and get the VAF
  all_snvs$VAF=0
  all_snvs$CCF=0
  all_snvs$clone=NA
  for( i in names(X$snv_list) ){
    ov=findOverlaps( query=all_snvs, subject=X$snv_list[[i]])
    if("VAF" %in% colnames(S4Vectors::mcols(X$snv_list[[i]]))){
      all_snvs$VAF[ queryHits(ov) ]= all_snvs$VAF[ queryHits(ov) ] + X$snv_list[[i]]$VAF[ subjectHits(ov)]
    }
    if("CCF" %in% colnames(S4Vectors::mcols(X$snv_list[[i]]))){
      all_snvs$CCF[ queryHits(ov) ]= all_snvs$CCF[ queryHits(ov) ] + X$snv_list[[i]]$CCF[ subjectHits(ov)]
    }
    all_snvs$clone[ queryHits(ov) ]=i
  }
  if(all( all_snvs$CCF == 0)){ S4Vectors::mcols(all_snvs)= S4Vectors::mcols(all_snvs) %>% tibble::as_tibble() %>% dplyr::select( -CCF)}
  if(all( all_snvs$VAF == 0)){ S4Vectors::mcols(all_snvs)= S4Vectors::mcols(all_snvs) %>% tibble::as_tibble() %>% dplyr::select( -VAF)}
  all_snvs
}


#' total CNV segments
#'
#'create a list of all the segments and their fraction
#' e.g. if a segment is in two clones, its fraction is the sum of the clones
#' at the end the result contains the CNV segments,
#' the column named 'clones' with the number of clones containing the event
#' and the column CCF with the total CCF of the event
#'
#' @param X cloneobj
#' @return a Granges object with combined CNV segments
#' @export

getCombinedCNV=function(X){
  # we use the reference as  the basis of comparison with each sublone
  all_cns=c()
  if(length(X$segments_list)==0){return(NULL)}
  #first create a list of all the segments
  all_cns=X$segments_list[[1]][ which(X$segments_list[[1]]$cnv.flag != "NEUT")]
  for(i in 2:length(X$segments_list)){
    a1=GenomicRanges::setdiff( all_cns, X$segments_list[[i]][ which(X$segments_list[[i]]$cnv.flag != "NEUT")])
    a2=GenomicRanges::setdiff(  X$segments_list[[i]][ which(X$segments_list[[i]]$cnv.flag != "NEUT")] ,all_cns )
    a3=GenomicRanges::intersect( X$segments_list[[i]][ which(X$segments_list[[i]]$cnv.flag != "NEUT")],all_cns)
    all_cns=c( a1, a2, a3)
  }
  all_cns=all_cns %>% GenomeInfoDb::sortSeqlevels( ) %>% sort()
  all_cns$CCF=0
  all_cns$constitutional.clones=0
  for(i in names( X$segments_list )){
    ov=findOverlaps( query=all_cns, subject=X$segments_list[[i]])
    all_cns$CCF[ queryHits(ov) ]= all_cns$CCF[ queryHits(ov) ] + X$segments_list[[i]][ subjectHits(ov)]$CCF

    all_cns$clone=i
  }
  all_cns
}



'seqlevelsStyle<-.cloneobj' = function(X, style){

  GenomeInfoDb::seqlevelsStyle(X$segments_list)=style
  GenomeInfoDb::seqlevelsStyle(X$snv_list)=style
}

seqlevelsStyle.cloneobj = function(X){
  return(
    GenomeInfoDb::seqlevelsStyle(X$segments_list)
  )
}
