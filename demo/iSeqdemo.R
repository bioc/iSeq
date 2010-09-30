
library(iSeq)

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

### plot an enriched region ############
ID = (reg1[1,4]):(reg1[1,5])
plotreg(tagct22[ID,2:3],tagct22[ID,5:6],tagct22[ID,7:8],peak=reg1[1,6])


#res2 = iSeq2(Y=tagct22[,1:4],gap=300, burnin=200,sampling=1000,winsize=2,ctcut=5,
#  a0=1,b0=1,a1=5,b1=1,k=1.0,verbose=FALSE)
