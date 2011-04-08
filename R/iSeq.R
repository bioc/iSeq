## The latest R functions for iChip software
# Quincy Mo, moq@mskcc.org
# Department of Epidemiology and Biostatistics,Memorial Sloan-Kettering Cancer center
# New York, NY 10065

# high-order Ising model for ChIP-chip data
## Y[,1] = chromosome,Y[,2]=start position of the bin, Y[,3]=end position of the bin,
#Y[,4] = counts or adjusted counts
iSeq2 = function(Y,gap=300,burnin=500,sampling=2000,winsize=2,ctcut=0.95,a0=1,b0=1,a1=5,b1=1,k=3,verbose=FALSE){

  if(missing(Y)){
    stop("Argument Y is missing!\n")
  }

  if((burnin < 1) || (sampling < 1)){
    stop("burin or sampling is too small!\n")
  }
  if(winsize < 1){
    stop("winsize must be >= 1\n")
  }

  if((ctcut < 0.5) || (ctcut > 1)){
    stop("Error: ctcut outside (0.5, 1)\n")
  }
  
  if(ctcut < 0.9){
    warning("Warning: ctcut outside [0.9, 0.99)\n")
  }
  if(k < 0){
    stop("kappa must be greater than 0! \n")
  }
  
  chrn = as.integer(as.factor(Y[,1]))
  rown = nrow(Y)
  ID = .C("fillgap",as.integer(gap),as.integer(rown),chrn,as.integer(Y[,2]),as.integer(Y[,3]),
    ipos=as.integer(rep(0, rown)),fpos=as.integer(rep(0, rown)),fn = as.integer(0))
  ipos = ID$ipos+1
#  fpos = ID$fpos[1:ID$fn]+1
  newn = rown + ID$fn
  count = rep(0,newn)
  count[ipos] = Y[,4]
  
  burn = as.integer(burnin)
  Size = as.integer(sampling)
  rowdim = as.integer(length(count))
  idata = as.integer(count) # exact genomic position is not needed
  post = as.double(rep(0,rowdim))
  halfw = as.integer(winsize)
  a0 = as.double(a0)
  b0 = as.double(b0)
  a1 = as.double(a1)
  b1 = as.double(b1) 
  X = as.integer(rep(0,rowdim))
  ctCut = as.integer(quantile(Y[,4],probs=ctcut))
  kp = as.double(k)
  lambda0 = as.double(rep(0,(burnin+sampling)))
  lambda1 = lambda0
  verb = as.integer(0)
  if(verbose == TRUE){verb = as.integer(1)}
  res = .C("iSeq2",burn,Size,rowdim,idata,halfw,ctCut,kp,a0,b0,a1,b1,postX=post,X=X,
    lambda0=lambda0,lambda1=lambda1,verbose=verb,PACKAGE="iSeq")
#  list(pp=res$postX,lambda0=res$lambda0,lambda1=res$lambda1)
  list(pp=res$postX[ipos],kappa=res$pKappa,lambda0=res$lambda0,lambda1=res$lambda1)
}


## standard one-dimensional Ising model for ChIP-chip data 
iSeq1 = function(Y,gap=300,burnin=500,sampling=2000,ctcut=0.95,a0=1,b0=1,a1=5,b1=1,
  k0=3,mink=0,maxk=10,normsd=0.1,verbose=FALSE){

  if(missing(Y)){
    stop("Argument Y is missing!\n")
  }
  if((burnin < 1) || (sampling < 1)){
    stop("burin or sampling is too small!\n")
  }
  if(k0 < 1){
    warning("k0 may be too small. You may set k0 between 1 and 4.\n")
  }
  
  if(mink < 0){
    stop("mink must be greater than 0!\n")
  }
  if(mink >= maxk){
    stop("mink must be less than maxk! \n")
  }
  if((ctcut < 0.5) || (ctcut > 1)){
    stop("Error: ctcut outside (0.5, 1)\n")
  }
  
  if(ctcut < 0.9){
    warning("Warning: ctcut outside [0.9, 0.99)\n")
  }
  
  if(normsd < 0){
    stop("normsd must be greater than 0. \n")
  }

  chrn = as.integer(as.factor(Y[,1]))
  rown = nrow(Y)
  ID = .C("fillgap",as.integer(gap),as.integer(rown),chrn,as.integer(Y[,2]),as.integer(Y[,3]),
    ipos=as.integer(rep(0, rown)),fpos=as.integer(rep(0, rown)),fn = as.integer(0))
  ipos = ID$ipos+1
#  fpos = ID$fpos[1:ID$fn]+1
  newn = rown + ID$fn
  count = rep(0,newn)
  count[ipos] = Y[,4]
  
  burn = as.integer(burnin)
  Size = as.integer(sampling)
  rowdim = as.integer(length(count))
  idata = as.integer(count)
  ctCut = as.integer(quantile(Y[,4],probs=ctcut))
  a0 = as.double(a0)
  b0 = as.double(b0)
  a1 = as.double(a1)
  b1 = as.double(b1) 
  kStart = as.double(k0)
  minK = as.double(mink)
  maxK = as.double(maxk)
  ransd = as.double(normsd)
  postx = as.double(rep(0,rowdim))
  X = as.integer(rep(0,rowdim))
  lambda0 = as.double(rep(0,(burnin+sampling)))
  lambda1 = lambda0
  pK = lambda0
  verb = as.integer(0)
  if(verbose == TRUE){verb = as.integer(1)}
  res = .C("iSeq1",burn,Size,rowdim,idata,ctCut,kStart,minK,maxK,ransd,postX=postx,
    X=X,pKappa=pK,a0,b0,a1,b1,lambda0=lambda0,lambda1=lambda1,verb,PACKAGE="iSeq")
  list(pp=res$postX[ipos],kappa=res$pKappa,lambda0=res$lambda0,lambda1=res$lambda1)
}

#function to find the binding regions based on posterior probability and selected cutoff
#pos[,1]- chromosome; pos[,2] - genomic start position; pos[,3] - genomic end position
#pos[,4] - countF; pos[,5] - countR
#pp is the posterior probability, cutoff is selection criteria for pp
#pos and pp are the results of the whole genome. 
#probes are merged if the genomic distance between neighboring probes are less than the maxgap
peakreg = function(chrpos,count,pp,cutoff,method=c("ppcut","fdrcut"),maxgap=300){

  if(nrow(chrpos) != length(pp)){
    stop("nrow(chrpos) must be equal to length(pp) \n")
  }
  if(cutoff < 0 || cutoff >1){
    stop("Error: cutoff should be in region [0, 1].\n")
  }
  if(missing(cutoff)){
    stop("Error: cutoff must be specified!")
  }
  if(maxgap < 1){
    stop("maxgap must be greater than 1! \n")
  }
  
  type = match.arg(method)
  ppcut = FALSE
  fdrcut = FALSE
  if(type == "ppcut"){ppcut = TRUE}
  else if(type == "fdrcut"){fdrcut = TRUE}
  
  chrf = as.factor(chrpos[,1])
  chrn = as.integer(chrf)
  rowID=1:nrow(chrpos)
  id = NULL
  if(ppcut){
    id = (pp > cutoff)
  }else if(fdrcut){
    sig = dirFDR(pp,cutoff)
    id = (sig$id > 0)
  }else{
    stop("Error: arg for method must be one of 'ppcut','fdrcut'.")
  }
  nsig = sum(id)
  if(nsig == 0){
    message("- No enriched region found - \n")
    br=NA
    return(br)
  }
  x = NULL
  if(sum(id) != 1){
    x = cbind(chrn[id],chrpos[id,2:3],count[id,1:2],rowID[id])
  }else{
    x = matrix(c(chrn[id],chrpos[id,2:3],count[id,1:2],rowID[id]),nrow=1)
  }
  xlen = nrow(x) #length
  ######### chrom start end countF coutR row.s row.e #############
  if(xlen == 1){
    br = cbind(x[,1:3],x[,6],x[,6],floor((x[,2]+x[,3])/2),round(pp[id],2),count[id,],sum(count[id,]),round(min(count[id,])/sum(count[id,]),2))
    colnames(br) = c("chr","gstart","gend","rstart","rend","peakpos","meanpp","ct1","ct2","ct12","sym")
    br = as.data.frame(br)
    br[1,1] = chrpos[br[1,4],1] #use the original chromosome label
    return(br)
  }
  zero = as.integer(rep(0,xlen))
  res = .C("mergeReg",as.integer(x[,1]),as.integer(x[,2]),as.integer(x[,3]),as.integer(x[,4]),as.integer(x[,5]),
    as.integer(x[,6]),as.integer(xlen),as.integer(maxgap),ochr=zero,ogstart=zero,ogend=zero,orstart=zero,
    orend=zero,opeak=zero,obimode=zero,nregion=as.integer(0),PACKAGE="iSeq")

  y = cbind(chr=res$ochr,gstart=res$ogstart,gend=res$ogend,rstart=res$orstart,rend=res$orend,peakpos=res$opeak)
                                        #cp=res$obimode)
  y = y[1:res$nregion,]
  
  if(res$nregion==1){
    y = matrix(y,ncol=6)
  }
  
  meanFun = function(y,pp){
    mean(pp[y[4]:y[5]])
  }

  countFun = function(y,ct){
    sum(ct[y[4]:y[5]])
  }
  
  mpp = round(apply(y,1,meanFun,pp=pp),2)
  countF = apply(y,1,countFun,ct=count[,1])
  countR = apply(y,1,countFun,ct=count[,2])
  countFR = countF + countR
  symtry = round(pmin(countF,countR)/countFR,2)
  symtry[countFR <= 0] = 0
  br = data.frame(y,meanpp=mpp,ct1=countF,ct2=countR,ct12=countFR,sym=symtry)
  br[,1] = chrpos[br[,4],1] #use the original chromosome label
  return(br)
}

dirFDR <- function(pp, fdrcut){
  beta.g = 1-pp
  k = sort(unique(beta.g))
  res = .C("fdr",as.integer(length(k)),as.double(k),as.integer(length(beta.g)),
    as.double(beta.g),efdr=as.double(rep(0,length(k))),PACKAGE="iSeq")
  cutoff = k[sum(res$efdr<=fdrcut)] #find the biggest k
  pcut = 1 - cutoff  ## this pcut should be >= pcut
  id = rep(0,length(pp))
  id[beta.g <= cutoff] = 1
  return(list(id=id,pcut=pcut))
}


### function used for iSeq package ####
### IP[,1] - chromosome, IP[,2]: genomic position, IP[,3], chain indicator 
mergetag = function(chip,control,maxlen=80,minlen=10,ntagcut=10){
  if(missing(chip)){stop("Argument chip is missing")}
  
  Y = NULL
  message("- Sort the data -")
  IP = chip[order(chip[,1],as.numeric(chip[,2])),]

  message("- Merge the ChIP sequence tags - ")
  xrow = nrow(IP)
  chr = as.integer(rep(0,xrow))  #a vector of 0 for initialization
  chrinput = as.numeric(as.factor(IP[,1])) #actual chromosome
  chain = as.numeric(as.factor(IP[,3])) 
  res = .C("binning",as.integer(chrinput),as.integer(IP[,2]),as.integer(chain),as.integer(xrow),
    as.integer(maxlen),as.integer(minlen),as.integer(ntagcut),chro=chr,gstart=chr,gend=chr,gend2=chr,
    count1=chr,count2=chr,chroID=chr,nregion=as.integer(0),PACKAGE="iSeq")
  Y = data.frame(gstart=res$gstart,gend=res$gend,ct12=(res$count1+res$count2),
    ct1=res$count1,ct2=res$count2,gend2=res$gend2)[1:res$nregion,]
  Y = data.frame(chr=IP[res$chroID[1:res$nregion],1],Y)

  if(!(missing(control))){
    message("- Background substraction - ")
    CON = control[order(control[,1],as.numeric(control[,2])),]
    IPchr = as.factor(Y[,1])
    CONchr = as.factor(CON[,1])
    CONchain = as.numeric(as.factor(CON[,3]))
    if(length(attributes(IPchr)$levels) != length(attributes(CONchr)$levels)){
      message("Error: chip[,1] and control[,1] have different numbers/types of chromosomes!\n")
      stop("The numbers/types of chromosomes in chip and control must be the same.\n")
    }
    IPchr = as.integer(as.numeric(IPchr))
    CONchr = as.integer(as.numeric(CONchr))
    IProw = nrow(Y)
    CONrow = nrow(CON)
    res2 = .C("subBkg",IPchr,as.integer(Y[,2]),IPend=as.integer(Y[,3]),as.integer(Y$gend2),as.integer(Y[,5]),
      as.integer(IProw),as.integer(maxlen),as.integer(ntagcut),CONchr,as.integer(CON[,2]),as.integer(CONchain),
      as.integer(CONrow),countIP=as.integer(Y[,4]),countFCON=as.integer(rep(0,IProw)),
      countRCON=as.integer(rep(0,IProw)),PACKAGE="iSeq")
    Y = data.frame(chr=Y[,1],gstart=Y[,2],gend=res2$IPend,adjct=res2$countIP,ipct1=Y[,5],ipct2=Y[,6],
               conct1=res2$countFCON,conct2=res2$countRCON)
  }
  return(Y)
}

# a function to plot genomic regions, gpos[,1]:start position, gpos[,2]: end position 
plotreg = function(gpos,ipct,conct,peak,col=c("yellow","green","grey0","blue")){
  if(is.vector(gpos)){
    gpos = matrix(c(gpos[,1],gpos[,2]),ncol=2)
    ipct = matrix(ipct,ncol=2)
    if(nrow(gpos) != nrow(ipct)){
      stop("Error: the dimensions for gpos and ipct are not the same!\n")
    }
    if(!missing(conct)){
      conct = matrix(conct,ncol=2)
    }
  }
  xmin = min(gpos[,1])
  xmax = max(gpos[,2])
  ymin = NA
  ymax = NA
  if(!missing(conct)){
    ymin = -max(c(ipct[,2],conct[,2]))
    ymax = max(c(ipct[,1],conct[,1]))
  }else{
    ymin = -max(ipct[,2])
    ymax = max(ipct[,1])
  }
  nreg = nrow(gpos)
  plot(c(xmin,xmax),c(ymin,ymax),type="n",xlab="Genomic Gposition",ylab="Counts")
  rect(gpos[,1],rep(0,nreg),gpos[,2],ipct[,1],col=col[1],border=TRUE)
  rect(gpos[,1],rep(0,nreg),gpos[,2],-ipct[,2],col=col[2],border=TRUE)
  if(!missing(conct)){
    rect(gpos[,1],rep(0,nreg),gpos[,2],conct[,1],col=col[3],border=TRUE)
    rect(gpos[,1],rep(0,nreg),gpos[,2],-conct[,2],col=col[4],border=TRUE)
  }
  if(!missing(peak)){
    abline(v=peak,lty=2,lwd=2)
  }
  abline(h=0,lwd=2)
}

## binobj is the object returned by binning function
# control[,1] -- chromosome, control[,2] - genomic position, control[,3] chain indicator
# control must be sorted, firstly by chromosome and then by genomic position
subcon = function(binobj,control){
  IPchr = as.factor(binobj$reg[,1])
  CONchr = as.factor(control[,1])
  CONchain = as.numeric(as.factor(control[,3]))
  if(length(attributes(IPchr)$levels) != length(attributes(CONchr)$levels)){
    message("Error: binobj$reg[,1] and control[,1] have different numbers/types of chromosomes!\n")
    stop("The numbers/types of chromosomes in binobj and control must be the same.\n")
  }
  IPchr = as.integer(as.numeric(IPchr))
  CONchr = as.integer(as.numeric(CONchr))
  IProw = nrow(binobj$reg)
  CONrow = nrow(control)
  res = .C("subBkg",IPchr,as.integer(binobj$reg[,2]),IPend=as.integer(binobj$reg[,3]),
    as.integer(IProw),as.integer(binobj$maxlen),CONchr,as.integer(control[,2]),as.integer(CONchain),
    as.integer(CONrow),countIP=as.integer(binobj$reg[,4]+binobj$reg[,5]),countFCON=as.integer(rep(0,IProw)),
    countRCON=as.integer(rep(0,IProw)),PACKAGE="iSeq")
  data.frame(chr=binobj$reg[,1],gstart=binobj$reg[,2],gend=res$IPend,adjct=res$countIP,
             confct=res$countFCON,conrct=res$countRCON)
}
