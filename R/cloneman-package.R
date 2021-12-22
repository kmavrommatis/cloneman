#' @keywords internal
#' @md
#' @title
#' cloneman: a package to parse and manage information from Copy Number calling and Clone Inference tools.
#'
#' @description
#' There exist many programs that are used to infer the copy number abberations in samples
#' They typically infer the copy number segments (either as integer numbers, e.g. 1, 3, 4 copies etc)
#' or as log ratios to the "normal" coverage.
#' Each program contains produces its own specific output format, in one or multiple files.
#'
#' This package was conceived in an attempt to normalize these outputs in a standard R format
#' with a set of predetermined columns and objects that can be used for downstream analysis.
#'
#' @section CNV segments:
#'
#' The simplest parsed output of these tools is a list of segments.
#' More specifically for CNVs a GRanges object is produced with all the segments and their
#' associated information
#'
#' The following columns are added to the metadata columns of the GRanges
#'
#' **total.allele.copies**: the total number of alleles of that segment,
#'
#' **major.allele.copies**: the number of copies of the _major_ allele,
#'
#' **minor.allele.copies**: the number of copies of the _minor_ allele,
#'
#' **normal.allele.copies**: the number of copies expected for a normal sample (typically coming from the normal ploidy of the sample)
#'
#' **cnv.flag**: For each segment depending on the combination of the above values a flag is added to indicate the type and
#' value of the aberration according to the schema:
#'
#' |**total.allele.copies** | **minor.allele.copies** |	**cnv.flag**         |
#' |---------|----------|--------------|
#' |0        |	0 	    | HOMDEL       |
#' |0        |	NA 	    | HOMDEL       |
#' |1        |	0 	    | HETDEL       |
#' |1        |	NA 	    | HETDEL       |
#' |2        |	0 	    | NEUT-LOH      |
#' |2        |	1 	    | NEUT         |
#' |2        |	NA 	    | NEUT/Unknown |
#' |2+       |	0 	    | AMP-LOH      |
#' |2+       |	1+ 	    | AMP          |
#' |2+       |	NA 	    | AMP-[LOH?]   |
#'
#' For methods where the ploidy is predicted it is parsed and added in the column _ploidy_.

#' Since this package was developed with human genome in mind a column caled _normal.allele.copies_ is added
#' with the default value of '2'
#' and for the Y (chrY) chromosome this is set to 1
#'
#' @section Clone Abundance:
#'
#' The clone abundance can be calculated by the CNV calling method, but can also be computed from SNVs, using tools such as Sciclone, Pyclone, Pyclone-vi, Fastclone etc.
#' There are two ways of expressing the clone abundance
#' one as _VAF_ (variant allele frequency) and as _CCF_ (cellular clone frequency).
#' Since we typically assume a diploid baseline, the VAF is typically the frequence of the mutation observed on a single allele and it is half of the CCF.
#'
#' To avoid this confusion we include both in the results
#'
#' **variant.allele.frequency**: represents the VAF of an event. Unless there is a way to compute it from other data, it is defaulted to NA.
#'
#' **cellular.abundance**: represents the CCF of an event.
#'
#' In all cases they are expressed in fractions that are ranging between 0 and 1.
#'
#' @section list of supported tools:
#'
#' |**tool**    |**type of data**      | **comments**|
#' |------------|----------------------|--------------|
#' |[FACETS](https://pubmed.ncbi.nlm.nih.gov/27270079/)      | CNV/CI/purity/ploidy |              |
#' |[CloneCNA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1174-7)    | CNV/CI/purity/ploidy |              |
#' |[ControlFreec](http://bioinformatics.oxfordjournals.org/content/28/3/423.long)| CNV/CI/purity/ploidy | CI not working well, ignored in this version|
#'
#' @section cloneman functions:
#' **Simple parsers**
#'
#' [parseFACETS()] : parse a dataframe produced by FACETS.
#'
#' [parseFACETS_rds()]: load an Rdata file produced by FACETS. This is the preferred method in order to retain additional information.
#'
#' [parseFACETS_file()]: load a tab delimited file with segments produced by FACETS. Additional information (e.g. purity) is not retained
#'
#' [parseCloneCNA()]: parse a dataframe produced by CloneCNA
#'
#' [parseCloneCNA_file()]: parse a file produced by CloneCNA. This is the preferred method in order to retain additional information.
#'
#' [parseCfreec()]: parse a dataframe with segments produced by ControlFreec.
#'
#' [parseCfreec_file()]: parse a _CNV or _CNV.pvalue.txt file produced by ControlFreec
#'
#' @docType package
#' @name cloneman
"_PACKAGE"
## usethis namespace: start
## usethis namespace: end
NULL
