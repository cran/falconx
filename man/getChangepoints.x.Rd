\name{getChangepoints.x}
\alias{getChangepoints.x}
\title{Getting Change-points}
\description{This function estimates the change-points where one or both parent-specific copy numbers change.  It uses a circular binary segmentation approach to find change-points in a binomial mixture process.  The output of the function is the set of locations of the break points.  If the whole genome is analyzed, it is recommended to run this function chromosome by chromosome, and the runs on different chromosomes can be done in parallel to shorten the running time.
}
\usage{
getChangepoints.x(readMatrix, biasMatrix, verbose=TRUE, COri=c(0.95,1.05), error=1e-5, 
maxIter=1000)
}
\arguments{
  \item{readMatrix}{A data frame with four columns and the column names are "AN", "BN", "AT"
and "BT". They are A-allele coverage in the tumor (case) sample, B-allele coverage
in the tumor (case) sample, A-allele coverage in the normal (control) sample,
and B-allele coverage in the normal (control) sample, respectively.}
  \item{biasMatrix}{A data frame with two columns and the column names are "sN", "sT".  They are the site-specific bias in total coverage for normal (control) sample and tumor (case) sample, respectively.}
  \item{verbose}{Provide progress messages if it is TRUE.  This argument is TRUE by default.  Set it to be FALSE if you want to turn off the progress messages.}
  \item{COri, error, maxIter}{Parameters used in estimating the success probabilities of the mixed binomial distribution.  See the manuscript by Chen and Zhang for more details.  "pOri" provides the initial success probabilities.  The two values in pOri needs to be different.  "error" provides the stopping criterion.  "maxIter" is the maximum iterating steps if the stopping criterion is not achieved.}
}
\seealso{
  \code{\link{getASCN.x}}
}
\references{
Chen, H., Jiang, Y., Maxwell, K., Nathanson, K. and Zhang, N. (2016). Allele-specific copy number estimation by whole Exome sequencing.
}
\examples{
data(Example) 
 # tauhat = getChangepoints.x(readMatrix, biasMatrix)  # uncomment this to run the function. 
}

