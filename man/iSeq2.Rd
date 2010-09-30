\name{iSeq2}
\alias{iSeq2}
\title{Bayesian hierarchical modeling of ChIP-seq data through hidden Ising models}
\description{
  iSeq2 implements the method that models the bin-based tag counts using Poisson-Gamma distribution
  and the hidden states of the bins using a hidden high-order Ising model. 
}
\usage{
iSeq2(Y,gap=300,burnin=500,sampling=2000,winsize=2,ctcut=5,
      a0=1,b0=1,a1=5,b1=1,k=3,verbose=FALSE)
}
\arguments{
  \item{Y}{Y should be a data frame containing the first 4 columns of
    the data frame returned by function 'mergetag()'.  The columns 1-4 of Y
    are chromosome IDs, start positions of the bins, end positions of the
    bins, tag counts in the bins.  For one-sample analysis, the tag
    counts can be the number of forward and reverse tags falling in the
    bins. For two-sample analysis, tag counts are the adjusted counts of
    ChIP samples, which are obtained by subtracting the control tag counts from
    corresponding ChIP tag counts for each bin.  If the user provides
    his/her own Y, Y must be firstly sorted by the chromosome ID, then by
    the start position, and then by the end position.}
  \item{gap}{gap is the average length of the sequenced DNA fragments.  If the
    distance between two nearest bins is greater than 'gap', a bin with
    0 tag count is insterted into the two neighboring bins for modeling.}
  \item{burnin}{The number of MCMC burn-in iterations.}
  \item{sampling}{The number of MCMC sampling iterations.  The
    posterior probability of enriched and non-enriched state is calculated
    based on the samples generated in the sampling period. }
  \item{winsize}{The parameter to control the order of interactions
    between genomic regions.  For example, winsize = 2, means that
    genomic region i interacts with regions i-2,i-1,i+1 and i+2. A
    balance between high sensitivity and low FDR could be achieved by
    setting winsize = 2.}
  \item{ctcut}{A value used to set the initial state for each
    genomic bin. If tag count of a bin is greater than ctcut, its state
    will be set to 1, otherwise -1. }
  \item{a0}{The scale hyper-parameter of the Gamma prior, alpha0.}
  \item{b0}{The rate hyper-parameter of the Gamma prior, beta0.}
  \item{a1}{The scale hyper-parameter of the Gamma prior, alpha1.}
  \item{b1}{The rate hyper-parameter of the Gamma prior, beta1.}
  \item{k}{The parameter used to control the strength of interaction
    between neighboring bins, which must be a positive value (k>0).  The
    larger the value of k, the stronger iterations between neighboring
    bins. The value for k may not be too small (e.g. < 1.0).  Otherwise, the Ising
    system may not be able to reach a super-paramagnetic state.}
  \item{verbose}{A logical variable.  If TRUE, the number of completed MCMC
    iterations is reported.}
}
\seealso{
  \code{\link{iSeq1}}, \code{\link{peakreg}},\code{\link{mergetag}},\code{\link{plotreg}}
}
\value{
  A list with the following elements.
  \item{pp}{The posterior probabilities of the bins in the enriched state.}
  \item{lambda0}{The posterior samples of the model parameter lambda0}
  \item{lambda1}{The posterior samples of the model parameter lambda1.}
}

\examples{
data(nrsf)
chip = rbind(nrsf$chipFC1592,nrsf$chipFC1862,nrsf$chipFC2002)
mock = rbind(nrsf$mockFC1592,nrsf$mockFC1862,nrsf$mockFC2002)
tagct = mergetag(chip=chip,control=mock,winsize=50)
tagct22 = tagct[tagct[,1]=="chr22",]
res2 = iSeq2(Y=tagct22[,1:4],gap=300, burnin=100,sampling=500,winsize=2,ctcut=5,
  a0=1,b0=1,a1=5,b1=1,k=1.0,verbose=FALSE)

}

\author{Qianxing Mo \email{moq@mskcc.org}}

\references{
  Qianxing Mo, Faming Liang.  (2010). Bayesian modeling of ChIP-chip data through
  a high-order Ising model. \emph{Biometrics}, 2010 Jan 29 [Epub ahead
  of print]. DOI: 10.1111/j.1541-0420.2009.01379.x

  Qianxing Mo (2010). A fully Bayesian hidden Ising model for ChIP-seq
  data analysis. (submitted).
}
\keyword{models}

