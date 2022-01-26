
#' get the flag for a CNV,
#' according to the table
#' |---------|----------|--------------|---------------------------------------------------------------------|
#' |total cn | minor cn |	type         | meaning of type                                                     |
#' |---------|----------|--------------|---------------------------------------------------------------------|
#' |0        |	0 	    | HOMDEL       | Homozygous deletion                                                 |
#' |0        |	NA 	    | HOMDEL       | Homozygous deletion                                                 |
#' |1        |	0 	    | HETDEL       | Heterozygous deletion (1 out of 2 alleles)                          |
#' |1        |	NA 	    | HETDEL       | Heterozygous deletion (this is an edge case that should not occur)  |
#' |2        |	0 	    | NEUT-LOH     | Loss of heterozygosity without change in the number of alleles      |
#' |2        |	1 	    | NEUT         | Normal                                                              |
#' |2        |	NA 	    | NEUT/Unknown | Normal or undetermined by the algorithm                             |
#' |2+       |	0 	    | AMP-LOH      | Amplification with loss of heterozygosity                           |
#' |2+       |	1+ 	    | AMP          | Amplification                                                       |
#' |2+       |	NA 	    | AMP-[LOH?]   | Amplification with (possible) loss of heterozygosity                |
#' |---------|----------|--------------|---------------------------------------------------------------------|
cnvFlag=function(major.cn=NA, minor.cn=0, total.cn=NA){
  #message(major.cn, " ", class(major.cn)," " ,minor.cn," ",class(minor.cn), " " , total.cn)
  if(is.na( total.cn ) & is.na(major.cn)){
    stop("Please provide either the total number of alleles or the major allele copies")
  }
  major.cn=as.numeric(major.cn)
  minor.cn=as.numeric(minor.cn)
  total.cn=as.numeric(total.cn)
  if(is.na(total.cn)){
    total.cn=ifelse( is.na(minor.cn), major.cn, major.cn + minor.cn)
  }

  flag="UNKNOWN"
  # set the flag
  if( total.cn == 0 & (minor.cn ==0 | is.na(minor.cn))){ flag="HOMDEL"}
  if( total.cn == 1 & (minor.cn ==0 | is.na(minor.cn))){ flag="HETDEL"}
  if( total.cn == 2 & (!is.na(minor.cn) & minor.cn ==0 )){ flag="NEUT-LOH"}
  if( total.cn == 2 & (!is.na(minor.cn) & minor.cn ==1 )){ flag="NEUT"}
  if( total.cn == 2 & is.na(minor.cn) ){ flag="NEUT/UNKNOWN"}
  if( total.cn >  2 & (!is.na(minor.cn) & minor.cn ==0 )){ flag="AMPLOH"}
  if( total.cn >  2 & (!is.na(minor.cn) & minor.cn >= 1 )){ flag="AMP"}
  if( total.cn >  2 & is.na(minor.cn)  ){ flag="AMP/LOH"}

  return( flag )

}
