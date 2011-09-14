\name{plotreg}
\alias{plotreg}
\title{A function used to plot enriched genomic regions}
\description{
  A function used to plot enriched genomic regions.
}
\usage{
plotreg(gpos,ipct,conct,peak,col=c("yellow","green","grey0","blue"))
}
\arguments{
  \item{gpos}{A n by 2 matrix or data frame. The rows correspond to
    genomic bins. The first and second columns contain the start and end positions
    of the genomic windows/bins, respectively.}
  \item{ipct}{A n by 2 matrix containing the ChIP tag counts
    corresponding to the bins in gpos. ipct[,1] contains the counts for
    the chain 1 (usually the forward chain); ipct[,2] contains the
    counts for the chain 2 (usually the reverse chain).}
  \item{conct}{A n by 2 matrix containing the control tag counts
    corresponding to the bins in gpos. ipct[,1] contains the counts for
    the chain 1 (usually the forward chain); ipct[,2] contains the
    counts for the chain 2 (usually the reverse chain).}
  \item{peak}{A vector containing the peak (center) positions of the
    genomic regions.}
  \item{col}{The colors used to fill the rectangles. col[1] is used for
    ipct[,1], col[2] for ipct[,2], col[3] for conct[,1] and col[4] for
    conct[,2], respectively.}
}
\seealso{
  \code{\link{iSeq1}}, \code{\link{iSeq2}}, \code{\link{peakreg}},\code{\link{mergetag}}
}
\value{
  No value returned. 
}

\examples{
#see the example in iSeq1 

}

\author{Qianxing Mo \email{moq@mskcc.org}}

\references{
   
  Qianxing Mo (2011). A fully Bayesian hidden Ising model for ChIP-seq
  data analysis. \emph{Biostatistics}, Advance Access published September 13,
  2011. doi:10.1093/biostatistics/kxr029
   
}
\keyword{models}

