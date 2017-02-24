\name{getChangepoints.x}
\alias{getChangepoints.x}
\title{Getting Change-points}
\description{This function estimates the change-points where one or both parent-specific copy numbers change.  It uses a circular binary segmentation approach to find change-points in a binomial mixture process.  The output of the function is the set of locations of the break points.  If the whole genome is analyzed, it is recommended to run this function chromosome by chromosome, and the runs on different chromosomes can be done in parallel to shorten the running time.
}
\usage{
getChangepoints.x(readMatrix, biasMatrix, verbose=TRUE, COri=c(0.95,1.05), error=1e-5, 
maxIter=1000, independence=TRUE, pos=NULL, readlength=NULL)
}
\arguments{
  \item{readMatrix}{A data frame with four columns and the column names are "AN", "BN", "AT"
and "BT". They are A-allele coverage in the tumor (case) sample, B-allele coverage
in the tumor (case) sample, A-allele coverage in the normal (control) sample,
and B-allele coverage in the normal (control) sample, respectively.}
  \item{biasMatrix}{A data frame with two columns and the column names are "sN", "sT".  They are the site-specific bias in total coverage for normal (control) sample and tumor (case) sample, respectively.}
  \item{verbose}{Provide progress messages if it is TRUE.  This argument is TRUE by default.  Set it to be FALSE if you want to turn off the progress messages.}
  \item{COri, error, maxIter}{Parameters used in estimating the success probabilities of the mixed binomial distribution.  See the manuscript by Chen and Zhang for more details.  "pOri" provides the initial success probabilities.  The two values in pOri needs to be different.  "error" provides the stopping criterion.  "maxIter" is the maximum iterating steps if the stopping criterion is not achieved.}
  \item{independence}{The model assumes reads are conditionally independent.  If this argument is FALSE, the pruning approach will be performed.}
  \item{pos}{The locations (in base pair) of the heterozygous sites. This information is needed when "independence=FALSE".}
  \item{readlength}{The length of read if the data is from single-end sequencing, and the maximum span of read pairs if the data if from paired-end sequencing.  This information is needed when "independence=FALSE".}
}
\seealso{
  \code{\link{getASCN.x}}
}
\examples{
data(Example) 
 # tauhat = getChangepoints.x(readMatrix, biasMatrix)  # uncomment this to run the function. 
}

