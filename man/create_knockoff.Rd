% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_knockoff.R
\name{create_knockoff}
\alias{create_knockoff}
\title{Create knockoff genotype data}
\usage{
create_knockoff(
  dat.hap,
  pos,
  M = 10,
  maxcor = 0.7,
  maxbp = 80000,
  phasing.dad = NULL,
  phasing.mom = NULL
)
}
\arguments{
\item{dat.hap}{A 6n*p matrix for the haplotype data, in which n is the number of trios and p is the number of variants. Each trio must consist of father, mother, and offspring (in this order). The haplotypes must be coded as 0 or 1. Missing haplotypes are not allowed.}

\item{pos}{A numeric vector of length p for the position of p variants.}

\item{M}{A positive integer for the number of knockoffs. The default is 10.}

\item{maxcor}{A real number in a range of [0,1] for the correlation threshold in hierarchical clustering, such that variants from two different clusters do not have a correlation greater than maxcor when constructing knockoff parents. The default is 0.7.}

\item{maxbp}{A positive integer for the size of neighboring base pairs used to generate knockoff parents. The default is 80000.}

\item{phasing.dad}{A numeric vector of length n that contains 1 or 2 to indicate which paternal haplotype was transmitted to offspring in each trio. If NA, the function will calculate the phasing information based on the input haplotype matrix.}

\item{phasing.mom}{A numeric vector of length n that contains 1 or 2 to indicate which maternal haplotype was transmitted to offspring in each trio. If NA, the function will calculate the phasing information based on the input haplotype matrix.}
}
\value{
A 3n*p*M array for the knockoff genotype data.
}
\description{
Create knockoff genotype data using phased haplotype data.
}
\examples{
data(KnockoffTrio.example)
dat.ko<-create_knockoff(KnockoffTrio.example$dat.hap,KnockoffTrio.example$pos,M=10)
}
