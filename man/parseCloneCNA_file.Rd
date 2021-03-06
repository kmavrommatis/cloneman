% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CloneCNA.R
\name{parseCloneCNA_file}
\alias{parseCloneCNA_file}
\title{CloneCNA claims that can efficiently and accurately identify somatic copy number alterations from heterogeneous tumor samples. (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1174-7)
CloneCNA can detect teh cellularity of each abberation, and from that extract the subclonal abundance. This comes in the column Cellularity.
When Cellularity is 0 it indicates normal
Presumably, segments with different values in that column correspond to different clones. As a result we can tag different clones based on these frequencies.}
\usage{
parseCloneCNA_file(clonecna.fn)
}
\arguments{
\item{clonecna.fn}{input filename with CNV segments,}
}
\value{
A GRanges object with the segments predicted by CloneCNA
}
\description{
CloneCNA claims that can efficiently and accurately identify somatic copy number alterations from heterogeneous tumor samples. (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1174-7)
CloneCNA can detect teh cellularity of each abberation, and from that extract the subclonal abundance. This comes in the column Cellularity.
When Cellularity is 0 it indicates normal
Presumably, segments with different values in that column correspond to different clones. As a result we can tag different clones based on these frequencies.
}
