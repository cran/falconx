\name{getASCN.x}
\alias{getASCN.x}
\title{Getting Allele-specific DNA Copy Number}
\description{Given a set of breakpoints where parent-specific copy number changes, this function estimates the parent-specific copy number for each segment, and the haplotype for the major chromosome on segments where the two homologous chromosomes have different copy numbers.  You are recommended to specify the parameter "rdep", the case-control genome-wide average coverage ratio.  Usually, a good estimate of rdep is (total mapped reads in tumor)/(total mapped reads in normal).  
}
\usage{
getASCN.x(readMatrix, biasMatrix, tauhat=NULL, threshold=0.15, COri=c(0.95,1.05), 
error=1e-5, maxIter=1000, independence=TRUE, pos=NULL, readlength=NULL)
}
\arguments{
  \item{readMatrix}{A data frame with four columns and the column names are "AN", "BN", "AT"
and "BT". They are A-allele coverage in the tumor (case) sample, B-allele coverage
in the tumor (case) sample, A-allele coverage in the normal (control) sample,
and B-allele coverage in the normal (control) sample, respectively.}
  \item{biasMatrix}{A data frame with two columns and the column names are "sN", "sT".  They are the site-specific bias in total coverage for normal (control) sample and tumor (case) sample, respectively.}
  \item{tauhat}{The estimated break points.  If it is not specified (NULL), then this function will first estimate the break points by calling the function "getChangepoints.x", and then estimate the parent-specific DNA copy number for each segment.}
  \item{threshold}{The estimated copy number are set to be 1 if it differs from 1 by less than this threshold.}
 \item{COri, error, maxIter}{Parameters used in estimating the success probabilities of the mixed binomial distribution.  See the manuscript by Chen and Zhang for more details.  "pOri" provides the initial success probabilities.  The two values in pOri needs to be different.  "error" provides the stopping criterion.  "maxIter" is the maximum iterating steps if the stopping criterion is not achieved.}
   \item{independence}{The model assumes reads are conditionally independent.  If this argument is FALSE, the pruning approach will be performed.}
  \item{pos}{The locations (in base pair) of the heterozygous sites. This information is needed when "independence=FALSE".}
  \item{readlength}{The length of read if the data is from single-end sequencing, and the maximum span of read pairs if the data if from paired-end sequencing.  This information is needed when "independence=FALSE".}

}
\value{
  \item{tauhat}{A vector holding the estimated break points in terms of the index in the coverage vectors.}
  \item{ascn}{The estimated parent-specific copy numbers in the segments between the break points in tauhat.}
  \item{Haplotype}{The estimated haplotype for the major chromosome (the chromosome has a higher copy number compared to its homologous chromosome) on segments where the two homologous chromosomes have different copy numbers.}
}
\seealso{
  \code{\link{getChangepoints.x}}, \code{\link{view}}
}
\examples{
data(Example) 
cn = getASCN.x(readMatrix, biasMatrix, tauhat=tauhat)
 # cn$tauhat would give the indices of change-points.  
 # cn$ascn would give the estimated allele-specific copy numbers for each segment.
 # cn$Haplotype[[i]] would give the estimated haplotype for the major chromosome in segment i 
 # if this segment has different copy numbers on the two homologous chromosomes.
}

