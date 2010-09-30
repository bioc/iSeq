\name{mergetag}
\alias{mergetag}
\title{Aggregate sequence tags into genomic windows/bins and count the number of
  tags in the windows/bins.}
\description{
  A function to aggregate sequence tags into genomic windows/bins with certain length
  specified by the user and count the number of tags falling in the windows/bins.
}
\usage{
mergetag(chip,control,winsize=50)
}

\arguments{
  \item{chip}{A n by 3 matrix or data frame. The Rows correspond to sequence
    tags. chip[,1] contains chromosome IDs; chip[,2] contains the genomic
    positions of sequence tags matched to the reference genome. For each
    tag, in order to accurately infer the true binding sites, we suggest
    using the middle positions of the tags as the tags' positions on the
    chromosomes.  Note a genomic position must be an integer.
    chip[,3] contains the direction indicators of the sequence tags.  The
    user can basically use any symbols to represent the forward or reverse
    chains. Function 'mergetag' use integer 1 and 2 to represent the
    directions of the chains by doing as.numeric(as.factor(chip[,3])).
    Therefore, the user should know the directions referred by integer 1
    and 2. For example, if the forward and reverse chains are
    represented by 'F' and 'R', respectively, then chains 1 and 2 will refer to
    the forward and reverse chain, respectively.  In the output, the tag
    counts are summarized for chains 1 and 2, respectively (see the below for details). }
  \item{control}{A n by 3 matrix or data frame. The column names of control
    must be the same as the column names of chip.}
  \item{winsize}{The window size for the genomic window/bin into which
    sequence tags are aggregated.}
}

\seealso{
  \code{\link{iSeq1}}, \code{\link{iSeq2}}, \code{\link{peakreg}},\code{\link{plotreg}}
}
\value{
  A data frame with rows corresponding to the bins and columns corresponding to the following: 
  \item{chr}{Chromosome IDs.}
  \item{gstart}{The start position of the bin.}
  \item{gend}{The end position of the bin.}
  \item{ct12}{For one-sample analysis, where only the ChIP data are
    available, ct12 = ipct1 + ipct2. For two-sample analysis, where both the ChIP
    and control data are available.  ct12 = maximum(ipct1+ipct2-conct1-conct2,0).}
  \item{ipct1}{The number of sequence tags for the chain 1 of the ChIP data.}
  \item{ipct2}{The number of sequence tags for the chain 2 of the ChIP data.}
  \item{conct1}{The number of sequence tags for the chain 1 of the control data.}
  \item{conct2}{The number of sequence tags for the chain 2 of the control data.}
}

\examples{
data(nrsf)
chip = rbind(nrsf$chipFC1592,nrsf$chipFC1862,nrsf$chipFC2002)
mock = rbind(nrsf$mockFC1592,nrsf$mockFC1862,nrsf$mockFC2002)

tagct = mergetag(chip=chip,control=mock,winsize=50)

}

\author{Qianxing Mo \email{moq@mskcc.org}}

\references{

  Qianxing Mo (2010). A fully Bayesian hidden Ising model for ChIP-seq
  data analysis. (submitted).
  
}
\keyword{models}

