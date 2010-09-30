\name{iSeq1}
\alias{iSeq1}
\title{Bayesian modeling of ChIP-seq data through hidden Ising models}
\description{
  iSeq1 implements the method that models the bin-based tag counts using Poisson-Gamma distribution
  and the hidden states of the bins using a standard 1D Ising model. 
}
\usage{
iSeq1(Y,gap=300,burnin=500,sampling=2000,ctcut=5,a0=1,b0=1,a1=5,b1=1,
  k0=3,mink=0,maxk=10,normsd=0.1,verbose=FALSE)
}
\arguments{
  \item{Y}{Y should be a data frame containing the first 4 columns of
    the data frame returned by function 'mergetag()'.  The columns 1-4 of Y
    are chromosome IDs, start position of the bin, end position of the
    bin, tag counts in the bins.  For one-sample analysis, the tag
    counts can be the number of forward and reverse tags falling in the
    bins. For two-sample analysis, tag counts are the adjusted counts of
    ChIP samples, which are obtained by subtracting the control tag counts from
    corresponding ChIP tag counts for each bin.  If the user provides
    his/her own Y, Y must be firstly sorted by the chromosome ID, then by
    the start position, and then by the end position.   
  }
  \item{gap}{gap is the average length of the sequenced DNA fragments.  If the
    distance between two nearest bins is greater than 'gap', a bin with
    0 tag count is insterted into the two neighboring bins for modeling.}
  \item{burnin}{The number of MCMC burn-in iterations.}
  \item{sampling}{The number of MCMC sampling iterations.  The
    posterior probability of enriched and non-enriched state is calculated
    based on the samples generated in the sampling period. }
  \item{ctcut}{A value used to set the initial state for each
    window/bin.  If tag count of a bin is greater than ctcut, its state
    will be set to 1, otherwise -1. }
  \item{a0}{The scale hyper-parameter of the Gamma prior, alpha0.}
  \item{b0}{The rate hyper-parameter of the Gamma prior, beta0.}
  \item{a1}{The scale hyper-parameter of the Gamma prior, alpha1.}
  \item{b1}{The rate hyper-parameter of the Gamma prior, beta1.}
  \item{k0}{The initial parameter used to control the strength of interaction
    between neighboring bins, which must be a positive value (k0>0).  A larger value of
    kappa represents a stronger interaction between neighboring bins.}
  \item{mink}{The minimum value of k(kappa) allowed.}
  \item{maxk}{The maximum value of k(kappa) allowed.}
  \item{normsd}{iSeq1 uses a Metropolis random walk proposal for
    sampling from the posterior distributions of the model
    parameter kappa. The proposal distribution is a normal distribution with
    mean 0 and standard deviation specified by normsd.}
  \item{verbose}{A logical variable.  If TRUE, the number of completed MCMC
    iterations is reported.}
}
\seealso{
  \code{\link{iSeq2}},\code{\link{peakreg}},\code{\link{mergetag}},\code{\link{plotreg}}
}
\value{
  A list with the following elements.
  \item{pp}{The posterior probabilities of bins in the enriched state.}
  \item{kappa}{The posterior samples of the interaction parameter of the
    Ising model.}
  \item{lambda0}{The posterior samples of the model parameter lambda0}
  \item{lambda1}{The posterior samples of the model parameter lambda1.}
}

\examples{
 data(nrsf)
 chip = rbind(nrsf$chipFC1592,nrsf$chipFC1862,nrsf$chipFC2002)
 mock = rbind(nrsf$mockFC1592,nrsf$mockFC1862,nrsf$mockFC2002)
 tagct = mergetag(chip=chip,control=mock,winsize=50)
 tagct22 = tagct[tagct[,1]=="chr22",]
 res1 = iSeq1(Y=tagct22[,1:4],gap=300,burnin=200,sampling=1000,
 ctcut=3,a0=1,b0=1,a1=5,b1=1,k0=3,mink=0,maxk=10,normsd=0.1,verbose=FALSE)

 reg1 = peakreg(tagct22[,1:3],tagct22[,5:6]-tagct22[,7:8],res1$pp,0.5,
        method="ppcut",maxgap=300)

reg2 = peakreg(tagct22[,1:3],tagct22[,5:6]-tagct22[,7:8],res1$pp,0.05,
         method="fdrcut",maxgap=300)

 ID = (reg1[1,4]):(reg1[1,5])
 plotreg(tagct22[ID,2:3],tagct22[ID,5:6],tagct22[ID,7:8],peak=reg1[1,6])
}

\author{Qianxing Mo \email{moq@mskcc.org}}

\references{
  Qianxing Mo (2010). A fully Bayesian hidden Ising model for ChIP-seq
  data analysis. (submitted).
}
\keyword{models}

