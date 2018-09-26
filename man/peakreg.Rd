\name{peakreg}
\alias{peakreg}
\title{Call and merge enriched genomic windows/bins.}
\description{
  A function used to call and merge enriched bins using the posterior
  probability calculated by iSeq1 or iSeq2 functions at certain
  posterior probability and false discovery rate (FDR) cutoffs. 
}
\usage{
peakreg(chrpos,count,pp,cutoff,method=c("ppcut","fdrcut"),maxgap=300)
}

\arguments{
  \item{chrpos}{A n by 3 matrix or data frame. The rows correspond to 
    genomic bins. The first column contains chromosome
    IDs; the second and third columns contain the start and end positions
    of the bin, respectively.}
  \item{count}{A n by 2 matrix containing the number of sequence tags in
    the bins specified by chrpos.  The first column contains the tag
    counts for chain 1 (usually the forward chain), and the second
    column contains the tag counts for chain 2 (usually the reverse
    chain).  See the document of the function 'mergetag' for the
    definition of chain 1 and 2.  The function uses the information in
    'count' to find the center of the enriched regions, where the true
    binding sites are usually located.}
  \item{pp}{A vector containing the posterior probabilities of bins
    in the enriched state returned by functions iSeq1 or iSeq2. }
  \item{cutoff}{The cutoff value (a scalar) used to call enriched
    bins.  If use posterior probability as a criterion
    (method="ppcut"), a bin is said to be enriched if its pp is
    greater than the cutoff.  If use FDR as a criterion
    (method="fdrcut"), bins are said to be enriched
    if the bin-based FDR is less than the cutoff.  The FDR is
    calculated using a direct posterior probability approach (Newton et
    al., 2004).}
  \item{method}{'ppcut' or 'fdrcut'.}
  \item{maxgap}{The criterion used to merge enriched bins.  If the
    genomic distance of adjacent bins is less than maxgap, the bins will
    be merged into the same enriched region.}
}
\seealso{
  \code{\link{iSeq1}}, \code{\link{iSeq2}}, \code{\link{mergetag}},\code{\link{plotreg}}
}
\value{
  A data frame with rows corresponding to enriched regions and columns
  corresponding to the following: 
  \item{chr}{Chromosome IDs.}
  \item{gstart}{The start genomic position of the enriched region.}
  \item{gend}{The end genomic position of the enriched region.}
  \item{rstart}{The row number for gstart in chrpos.}
  \item{rend}{The row number for gend in chrpos.}
  \item{peakpos}{The inferred center (peak) of the enriched region.}
  \item{meanpp}{The mean posterior probability of the merged
    regions/bins.}
  \item{ct1}{total tag counts for the region from gstart to gend for the
    chain corresponding to count[,1]; ct1=sum(count[rstart:rend,1])}
  \item{ct2}{total tag counts for the region from gstart to gend for the
    chain corresponding to count[,2]; ct2=sum(count[rstart:rend,2])}
  \item{ct12}{ct12 = ct1 + ct2}
  \item{sym}{A parameter used to measure if the forward and reverse tag
    counts are symmetrical (or balanced) in enriched regions. The values
    range from 0.5 (perfect symmetry) to 0 (complete asymmetry).}
%  \item{cp}{The number of change points in the merged enriched
%    region. A change point is a genomic position where the signs of
%    the net tag counts of the the neighboring bins change.  For each
%    bin, the net tag counts = tag counts of chain 1 - tag counts of
%    chain 2. For example, if the signs of the net tag counts for two
%    neighboring bins are +  - (or, -  +), then we say there is a
%    change point. A typical enriched region has one change point.}
}

\examples{

data(nrsf)
chip = rbind(nrsf$chipFC1592,nrsf$chipFC1862,nrsf$chipFC2002)
mock = rbind(nrsf$mockFC1592,nrsf$mockFC1862,nrsf$mockFC2002)
tagct = mergetag(chip=chip,control=mock,maxlen=80,minlen=10,ntagcut=20)
tagct22 = tagct[tagct[,1]=="chr22",]
res1 = iSeq1(Y=tagct22[,1:4],gap=200,burnin=200,sampling=500,ctcut=0.95,a0=1,b0=1,
 a1=5,b1=1, k0=3,mink=0,maxk=10,normsd=0.1,verbose=FALSE)

reg1 = peakreg(tagct22[,1:3],tagct22[,5:6]-tagct22[,7:8],res1$pp,0.5,
        method="ppcut",maxgap=200)

reg2 = peakreg(tagct22[,1:3],tagct22[,5:6]-tagct22[,7:8],res1$pp,0.05,
         method="fdrcut",maxgap=200)
}

\author{Qianxing Mo \email{qianxing.mo@moffitt.org}}

\references{

  Qianxing Mo. (2012). A fully Bayesian hidden Ising model for ChIP-seq
  data analysis. \emph{Biostatistics} 13(1), 113-28.
  
  Newton, M., Noueiry, A., Sarkar, D., Ahlquist, P. (2004). Detecting
  differential gene expression with a semiparametric hierarchical mixture method.
  \emph{Biostatistics}  5 , 155-176.
}
\keyword{models}

