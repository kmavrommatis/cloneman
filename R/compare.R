#' Compare the CNVs between two or more cloneobj
#'
#'
#' Comparing the CNV segments between two methods is trick
#' Segments with slightly different coordinates may be overlapping by most
#' of their length, but they are not identical
#' In the case of WES data, we can compare the CNV status of each of the targets
#' This assumes that the dataset has been filtered to select only those regions
#'
#'
#' This function can compare two cloneobj and find how many segments have
#' 1. the same status (regardless of copy number), e.g. two segments with flag AMP
#' will be considered the same even if one is 4 copies and the other 8
#' 2. the same status and same copy number.
#' 3. the same status and same CCF. In this case the CCF is normalized for purity
#'
#' @param query the query cloneobj
#' @param subject the subject cloneobj
#' @param method The method of comparison to employ (based on the description above). Options are 'status','cnv','ccf'
#' @param wiggle The differnce between the CCF when compared that is considered acceptable
#' @param plot if set to TRUE draw an alluvial plot that shows how events are distributed between samples.
#' @return a table with the results
#' @export
compareCNV=function( query, subject, method='status', wiggle=0.1){

  if( !"cloneobj" %in% class( query) ||  !"cloneobj" %in% class( subject) ){
    stop("query and subject arguments have to be of class 'cloneobj'")
  }
  df=cloneman:::compareStatus( query, subject, wiggle=wiggle)

  stats=computeStatistics(df, method='status')

  return(stats)
}


#' Internal method to calculate statistic
#' @param the dataframe with the status information
#' @return a data frame with statistics
computeStatistics=function(df, method='status'){


  stats=list(
  'TP'=length( which( df$status == "TP")),
  'FP'=length( which( df$status == 'FP')),
  'FN'=length( which( df$status == 'FN')),
  'filtered TP'=length( which( df$status == 'TP' & df$query !='NEUT')),
  'total'=nrow( df )
  )

  stats[['Precision']]=stats[['TP']]/(stats[['TP']] + stats[['FP']])
  stats[['Recall']]   =stats[['TP']]/(stats[['TP']] + stats[['FN']])
  stats[['F-measure']]=2 * (( stats[['Precision']] * stats[['Recall']] )/( stats[['Precision']] + stats[['Recall']]) )

  return(stats)
}


#' Internal method to plot an alluvial with distribution of events
#'
#' @param query  query the query cloneobj
#' @param subject the subject cloneobj
#' @return a figure with distribution of events
#' @export
plotAlluvialCNV=function( query, subject,wiggle=0.1){
  if( !"cloneobj" %in% class( query) ||  !"cloneobj" %in% class( subject) ){
    stop("query and subject arguments have to be of class 'cloneobj'")
  }

  df=cloneman:::compareStatus( query, subject, wiggle=wiggle)
  dt=df[,c(1:2,7:9)]

  dt=dt %>%
    dplyr::mutate( cnv=ifelse( cnv == 'FN', 'underestimate',ifelse(cnv=='FP','overestimate','correct'))) %>%
    dplyr::mutate( ccf=ifelse( ccf == 'FN', 'underestimate',ifelse(ccf=='FP','overestimate','correct')))
  easyalluvial::alluvial_wide( data=dt, max_variables=5, fill_by='all_flows',colorful_fill_variable_stratum=FALSE)
}

#' Compare two cloneobj by status
#'
#' We use the cnv.flag to do the comparisons
#' Importatn assumptions: both methods have the same segments
#'
#' For the CNV comparison in mode status:
#'   TP: same in query and subject
#'   FP: NEUT in query and something else in subject
#'   FN: something else in query NEUT in subject
#'
#' For the CNV comparison in mode cnv (compare the actual segment copy number)
#' The comparison is restricted to the TP segments of the 'status' part
#'   TP: same in query and subject
#'   FP: query has less copies than subject (more copies in subject)
#'   FN: query has more copies than subject (less copies in subject)
#'
#' For the CNV comparison in mode ccf (compare the CCF of each segment)
#' The comparison is restricted to the TP segments of the 'status' part
#'   TP: |query CCF - subject CCF| < wiggle
#'   FP:
#'
#' @param query cloneobj
#' @param subject cloneobj
#' @return a dataframe with the status, cnv and ccf for each overlapping segment
#' @export
compareStatus=function(query,subject, wiggle=0.1){
  if( !"cloneobj" %in% class( query) ||  !"cloneobj" %in% class( subject) ){
    stop("query and subject arguments have to be of class 'cloneobj'")
  }
  query.gr=query$segments
  subject.gr=subject$segments
  query.name=ifelse(is.null(query$name),'query',getName(query))
  subject.name=ifelse( is.null(subject$name),'subject', getName(subject))


  queryOv=pintersect( findOverlapPairs( query.gr, subject.gr ))
  subjectOv=pintersect( findOverlapPairs( subject.gr, query.gr ))
  keep=intersect( which(subjectOv$hit==TRUE),which(queryOv$hit==TRUE) )
  queryOv=queryOv[keep]
  subjectOv=subjectOv[keep]
  # all the ranges
  all=length( subjectOv )
  all_bp=sum(width(subjectOv))

  # find the common ranges (same cnv.flag)
  tt=cbind(
    query=queryOv$cnv.flag,
    subject=subjectOv$cnv.flag,
    query_CNV=queryOv$total.allele.copies,
    subject_CNV=subjectOv$total.allele.copies,
    query_CCF=queryOv$cnv.ccf,
    subject_CCF=subjectOv$cnv.ccf
    ) %>% as.data.frame()
  tt=tt %>%
    dplyr::mutate(across( c(query_CNV,subject_CNV,query_CCF, subject_CCF) , as.numeric ) )  %>%
    dplyr::mutate( status= ifelse( query == subject, "TP", ifelse( query == "NEUT", "FP" ,"FN"))) %>%
    dplyr::mutate( cnv =   ifelse( status == "TP", ifelse( query_CNV == subject_CNV, 'TP', ifelse( query_CNV < subject_CNV, 'FP', 'FN'))  , NA) ) %>%
    dplyr::mutate( ccf =   ifelse( status == "TP", ifelse( abs(query_CCF - subject_CCF)<wiggle, 'TP',ifelse( query_CCF< subject_CCF, 'FP', 'FN'  )), NA)  ) %>%
    dplyr::rename_with( ~stringr::str_replace( .x,pattern='query', replacement = query.name), starts_with( 'query' ),
                        ~stringr::str_replace( .x,pattern='subject', replacement = subject.name), starts_with( 'subject' ))

  return(tt)
}
