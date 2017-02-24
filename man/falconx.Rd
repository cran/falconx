\name{falconx}
\alias{falconx}
\title{Finding Allele-Specific Copy Number in Whole-Exome Sequencing Data}
\description{This library contains a set of tools for allele-specific DNA copy number profiling using whole exome sequencing.  Given the allele-specific coverage and site biases at the variant loci, this program segments the genome into regions of homogeneous allele-specific copy number.  It requires, as input, the read counts for each variant allele in a pair of case and control samples, as well as the site biases. For detection of somatic mutations, the case and control samples can be the tumor and normal sample from the same individual.
}
\author{
Hao Chen and Nancy R. Zhang

Maintainer: Hao Chen (hxchen@ucdavis.edu)
}

\seealso{
  \code{\link{getChangepoints.x}}, \code{\link{getASCN.x}}, \code{\link{view}}
}
\examples{
data(Example) 
 # tauhat = getChangepoints.x(readMatrix, biasMatrix)  # uncomment this to run the function. 
cn = getASCN.x(readMatrix, biasMatrix, tauhat=tauhat)
 # cn$tauhat would give the indices of change-points.  
 # cn$ascn would give the estimated allele-specific copy numbers for each segment.
 # cn$Haplotype[[i]] would give the estimated haplotype for the major chromosome in segment i 
 # if this segment has different copy numbers on the two homologous chromosomes.
view(cn)
}

