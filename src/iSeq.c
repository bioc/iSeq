/* Codes used for Hidden Ising modeling of ChIP-seq data */
/* Coding start date: 2/18/2010.  Complete the package on Sept. 14, 2010 */
/* Qianxing Mo (qmo@bcm.edu),
   Division of Biostatistics, Dan L. Duncan Cancer Center, Baylor College of Medicine */
/*
#include <stdlib.h>
#include <stdio.h>
*/

#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void iSeq2(int *burning,int *size,int *nrow,int *mydt,int *halfwin,int *countCutoff,
	   double *kappa,double *alpha0,double *beta0,double *alpha1,double *beta1,
	   double *postX,int *X,double *lambda0,double *lambda1,int *verbose){

  int sampleSize,i,j,r,n0,n1,state,wtsum;
  double p,s0,s1;

  sampleSize = (*burning) + (*size);

  s0 = 0;
  s1 = 0;
  n1 = 0;
  for(i=0; i<(*nrow); i++){
    if(mydt[i] > (*countCutoff)){
      X[i] = 1;
      s1 += mydt[i];
      n1 += 1;
    }else {
      X[i] = -1;
      s0 += mydt[i];
    }
  }
  n0 = (*nrow) - n1;
  /*  Rprintf("lambda0 = %f, lambda1 = %f \n",s0/n0,s1/n1); */

  /* Gibbs sampling */
  GetRNGstate();
  for(r=0; r<sampleSize;r++){
    lambda0[r] = rgamma((s0+(*alpha0)),1.0/(n0+(*beta0)));
    lambda1[r] = rgamma((s1+(*alpha1)),1.0/(n1+(*beta1)));

    /* for the first halfwin probes */
    for(i=0; i<(*halfwin); i++){
      state = X[i];
      wtsum = 0 - X[i];
      /* count the number of 1 around probe[i] including probe[i] itself,so initialize wtsum = -X[i] */
      for(j=0; j<=(i+(*halfwin)); j++){
	wtsum = wtsum + X[j];
      }
      /* p = 1.0/(1.0+R_pow_di(lambda0[r]/lambda1[r],mydt[i])*exp(lambda1[r]-lambda0[r]-(*kappa)*wtsum*2.0));*/
      p = 1.0/(2+expm1(mydt[i]*log(lambda0[r]/lambda1[r])+lambda1[r]-lambda0[r]-(*kappa)*wtsum*2.0));
      if(runif(0,1) < p){
	X[i] = 1;
      }else {
	X[i] = -1;
      }
      if(state != X[i]){
	/*if state change from -1 to 1 */
	if(state == -1){
	  s0 = s0 - mydt[i];
	  n0 = n0 - 1;
	  s1 = s1 + mydt[i];
	  n1 = n1 + 1;
	}else{
	  s0 = s0 + mydt[i];
	  n0 = n0 + 1;
	  s1 = s1 - mydt[i];
	  n1 = n1 - 1;
	}
      }
      /* identifiability constraint: let lambda0 < lambda1 for; if not, -1 is binding site*/
      if(r >= (*burning)){
	if(lambda0[r] < lambda1[r]){
	  if(X[i] == 1){
	    postX[i] = postX[i] + 1;
	  }
	}else {
	  if(X[i] == -1){
	    postX[i] = postX[i] + 1;
	  }
	}
      }
    }

    /*for the probes between halfwin and n-halfwin */
    for(i=(*halfwin); i<((*nrow)-(*halfwin)); i++){
      state = X[i];
      wtsum = 0 - X[i];
      /* count the number of 1 around probe[i] including probe[i] itself,so initialize wtsum = -X[i] */
      for(j=(i-(*halfwin)); j<=(i+(*halfwin)); j++){
	wtsum = wtsum + X[j];
      }
      /* p = 1.0/(1.0+R_pow_di(lambda0[r]/lambda1[r],mydt[i])*exp(lambda1[r]-lambda0[r]-(*kappa)*wtsum*2.0));*/
      p = 1.0/(2+expm1(mydt[i]*log(lambda0[r]/lambda1[r])+lambda1[r]-lambda0[r]-(*kappa)*wtsum*2.0));
      if(runif(0,1) < p){
	X[i] = 1;
      }else {
	X[i] = -1;
      }
      if(state != X[i]){
	/*if state change from -1 to 1 */
	if(state == -1){
	  s0 = s0 - mydt[i];
	  n0 = n0 - 1;
	  s1 = s1 + mydt[i];
	  n1 = n1 + 1;
	}else{
	  s0 = s0 + mydt[i];
	  n0 = n0 + 1;
	  s1 = s1 - mydt[i];
	  n1 = n1 - 1;
	}
      }
      /* identifiability constraint: let lambda0 < lambda1; if not, -1 is binding site*/
      if(r >= (*burning)){
	if(lambda0[r] < lambda1[r]){
	  if(X[i] == 1){
	    postX[i] = postX[i] + 1;
	  }
	}else {
	  if(X[i] == -1){
	    postX[i] = postX[i] + 1;
	  }
	}
      }
    }
    
    /* for the last halfwin probes */
    for(i=((*nrow)-(*halfwin)); i<(*nrow); i++){
      state = X[i];
      wtsum = 0 - X[i];
      /* count the number of 1 around probe[i] including probe[i] itself,so initialize wtsum = -X[i] */
      for(j=(i-(*halfwin)); j<(*nrow); j++){
	wtsum = wtsum + X[j];
      }
      /*p = 1.0/(1.0+R_pow_di(lambda0[r]/lambda1[r],mydt[i])*exp(lambda1[r]-lambda0[r]-(*kappa)*wtsum*2.0));*/
      p = 1.0/(2+expm1(mydt[i]*log(lambda0[r]/lambda1[r])+lambda1[r]-lambda0[r]-(*kappa)*wtsum*2.0));
      if(runif(0,1) < p){
	X[i] = 1;
      }else {
	X[i] = -1;
      }
      if(state != X[i]){
	/*if state change from -1 to 1 */
	if(state == -1){
	  s0 = s0 - mydt[i];
	  n0 = n0 - 1;
	  s1 = s1 + mydt[i];
	  n1 = n1 + 1;
	}else{
	  s0 = s0 + mydt[i];
	  n0 = n0 + 1;
	  s1 = s1 - mydt[i];
	  n1 = n1 - 1;
	}
      }
      /* identifiability constraint: let lambda0 < lambda11; if not, -1 is binding site*/
      if(r >= (*burning)){
	if(lambda0[r] < lambda1[r]){
	  if(X[i] == 1){
	    postX[i] = postX[i] + 1;
	  }
	}else {
	  if(X[i] == -1){
	    postX[i] = postX[i] + 1;
	  }
	}
      }
    }

    if(r%2000 == 0){
      R_CheckUserInterrupt();
      if(*verbose == 1){
	Rprintf("%d  ",r);
      }
    }
  }
  PutRNGstate();
  
  if(*verbose == 1){
    Rprintf("\n");
  }

  /*  Rprintf("n0 = %d, n1 = %d \n",n0,n1); */

  if((n0 < 1) || (n1 < 1)){
    Rprintf("Warning: all bins are in the same state at the last MCMC iteration.\n NO enriched region is found!\n");
  }

  for(i=0; i<(*nrow); i++){
    postX[i] = (double)postX[i]/(*size);
  }

}


/* one parameter standard ising model */
/* treat kappa as random variable, use MH algorithm to estimate it */
void iSeq1(int *burning,int *Size,int *nrow,int *dt,int *countCutoff,double *kappaStart,double *minKappa,
	   double *maxKappa,double *ransd,double *postX,int *X,double *postKappa,double *alpha0,double *beta0,
	   double *alpha1,double *beta1,double *lambda0,double *lambda1,int *verbose){

  int sampleSize,i,j,r,n0,n1,state,nrowm1,sumxx;
  double p,s0,s1,kappa,newKappa,ratio;

  sampleSize = (*burning) + (*Size);
  nrowm1 = *nrow - 1;
 
  s1 = 0.0;
  s0 = 0.0;
  n1 = 0;
  for(i=0; i<(*nrow); i++){
    if(dt[i] > (*countCutoff)){
      X[i] = 1;
      s1 += dt[i];
      n1 += 1;
    }else {
      X[i] = -1;
      s0 += dt[i];
    }
  }

  n0 = (*nrow) - n1;
  /*  Rprintf("lambda0 = %f, lambda1 = %f \n",s0/n0,s1/n1); */
  kappa = (*kappaStart);

  /* Gibbs sampling */
  GetRNGstate();
  for(r=0; r<sampleSize;r++){
    lambda0[r] = rgamma((s0+(*alpha0)),1.0/(n0+(*beta0)));
    lambda1[r] = rgamma((s1+(*alpha1)),1.0/(n1+(*beta1)));

    /* for the first probe */
    i = 0;
    state = X[i];
    /*    p = 1.0/(1+R_pow_di(lambda0[r]/lambda1[r],dt[i])*exp(lambda1[r]-lambda0[r]-2.0*X[i+1]*kappa)); */
    p = 1.0/(2+expm1(dt[i]*log(lambda0[r]/lambda1[r])+lambda1[r]-lambda0[r]-2.0*X[i+1]*kappa));
    if(runif(0,1)<p){
      X[i] = 1;
    }else {
      X[i] = -1;
    }
  
    if(state != X[i]){
      /*if state change from -1 to 1 */
      if(state == -1){
	s0 = s0 - dt[i];
	n0 = n0 - 1;
	s1 = s1 + dt[i];
	n1 = n1 + 1;
      }else {
	s0 = s0 + dt[i];
	n0 = n0 + 1;
	s1 = s1 - dt[i];
	n1 = n1 - 1;
      }
    }
    /*
      if(X[i] == 1 && r >= (*burning)){
      postX[i] = postX[i] + X[i];
      }
    */
    /* identifiability constraint: let lambda0 < lambda1; if not, -1 is binding site*/
    if(r >= (*burning)){
      if(lambda0[r] < lambda1[r]){
	if(X[i] == 1){
	  postX[i] = postX[i] + 1;
	}
      }else {
	if(X[i] == -1){
	  postX[i] = postX[i] + 1;
	}
      }
    }
    
    for(i=1; i<nrowm1; i++){
      state = X[i];
      /*      p = 1.0/(1+R_pow_di(lambda0[r]/lambda1[r],dt[i])*exp(lambda1[r]-lambda0[r]-2.0*(X[i-1]+X[i+1])*kappa)); */
      p = 1.0/(2+expm1(dt[i]*log(lambda0[r]/lambda1[r])+lambda1[r]-lambda0[r]-2.0*(X[i-1]+X[i+1])*kappa));
      if(runif(0,1)<p){
	X[i] = 1;
      }else {
	X[i] = -1;
      }
      
      if(state != X[i]){
	/*if state change from -1 to 1 */
	if(state == -1){
	  s0 = s0 - dt[i];
	  n0 = n0 - 1;
	  s1 = s1 + dt[i];
	  n1 = n1 + 1;
	}else {
	  s0 = s0 + dt[i];
	  n0 = n0 + 1;
	  s1 = s1 - dt[i];
	  n1 = n1 - 1;
	}
      }
      /*
	if(X[i] == 1 && r >= (*burning)){
	postX[i] = postX[i] + 1;
	}
      */
      /* identifiability constraint: let lambda0 < lambda1; if not, -1 is binding site*/
      if(r >= (*burning)){
	if(lambda0[r] < lambda1[r]){
	  if(X[i] == 1){
	    postX[i] = postX[i] + 1;
	  }
	}else {
	  if(X[i] == -1){
	    postX[i] = postX[i] + 1;
	  }
	}
      }
    }
    /* For the last probe */
    i = nrowm1;
    state = X[i];
    /*    p = 1.0/(1+R_pow_di(lambda0[r]/lambda1[r],dt[i])*exp(lambda1[r]-lambda0[r]-2.0*X[i-1]*kappa)); */
    p = 1.0/(2+expm1(dt[i]*log(lambda0[r]/lambda1[r])+lambda1[r]-lambda0[r]-2.0*X[i-1]*kappa));
    if(runif(0,1)<p){
      X[i] = 1;
    }else {
      X[i] = -1;
    }
    
    if(state != X[i]){
      /*if state change from -1 to 1 */
      if(state == -1){
	s0 = s0 - dt[i];
	n0 = n0 - 1;
	s1 = s1 + dt[i];
	n1 = n1 + 1;
      }else {
	s0 = s0 + dt[i];
	n0 = n0 + 1;
	s1 = s1 - dt[i];
	n1 = n1 - 1;
      }
    }
    /*
      if(X[i]==1 && r >= (*burning)){
      postX[i] = postX[i] + X[i];
      }   
    */
    /* identifiability constraint: let lambda0 < lambda1; if not, -1 is binding site*/
    if(r >= (*burning)){
      if(lambda0[r] < lambda1[r]){
	if(X[i] == 1){
	  postX[i] = postX[i] + 1;
	}
      }else {
	if(X[i] == -1){
	  postX[i] = postX[i] + 1;
	}
      }
    }

    /* Metropolis random walk update kappa */
    newKappa = kappa + rnorm(0, *ransd);
    if(newKappa > (*minKappa) && newKappa < (*maxKappa)){
      sumxx = 0;
      for(i=0; i<nrowm1; i++){
	sumxx = sumxx + X[i]*X[i+1];
      }
      ratio = nrowm1*(log(cosh(kappa))-log(cosh(newKappa)))+(newKappa-kappa)*sumxx;
      if(ratio >= 0){
	kappa = newKappa;
      }else {
	if(runif(0,1) < exp(ratio)){
	  kappa = newKappa;	
	}
      }
    }
    postKappa[r] = kappa;

    if(r%2000 == 0){
      R_CheckUserInterrupt();
      if(*verbose == 1){
	Rprintf("%d  ",r);
      }
    }
  }
  PutRNGstate();

  /*  Rprintf("n0 = %d, n1 = %d \n",n0,n1); */

  if(*verbose == 1){
    Rprintf("\n");
  }

  if((n0 < 1) || (n1 < 1)){
    Rprintf("Warning: all bins are in the same state at the last MCMC iteration. \n NO enriched region is found!\n");
  }

  for(i=0; i<(*nrow); i++){
    postX[i] = (double)postX[i]/(*Size);
  }
}

/* input: 
   chr: chromosomes; pos: genomic positions; len: length of the vector
   chain: 0 for positive chain, 1 for negative chain. regsize: size of binning region;
   output: all vector should be initialized with 0s
   chro: chromosome; gstart: genomic start position; gend: genomic end position; 
   count1: count for + chain; count2: count for negative chain.
   chroID: the index used for replace the numberic chromosome with the original chromosome ID
   because character chromosomes are changed to nummeric chromosome for compuation efficiency
*/
void binning2(int *chr,int *pos,int *chain,int *len,int *regsize,int *maxtag,int *chro,int *gstart,int *gend,
	     int *count1,int *count2,int *chroID,int *nregion){
  int i, j;  /* column start from 0 */

  /* copy the first row of x and initiate the number of reads to 1 at that point */
  j = 0;
  chro[j] = chr[j];
  chroID[j] = j+1;     /*because R index is from 1:n */
  gstart[j] = pos[j];      /* region start position */ 
  gend[j] = pos[j];  /* region end position */
  if(chain[j]==1){       /* count for positive strand in the region */
    count1[j] = 1; 
  }else{
    count2[j] = 1;
  }

  for(i=1; i < (*len); i++){
    if(chr[i] == chro[j]){ /* if the same chromosome */
      /*      if((pos[i]-gstart[j]) < (*regsize)){ */
      if(((pos[i]-gstart[j]) < (*regsize)) && (count1[j]<(*maxtag)) && (count2[j]<(*maxtag))){ 
	gend[j] = pos[i];  /* increase the end position of the region */
	if(chain[i]==1){
	  count1[j] += 1;       /* count for positive strand in the region */
	}else{
	  count2[j] += 1;   /* count for positive strand in the region */
	}
      }else{
	j = j+1;
	chro[j] = chr[i];
	chroID[j] = i+1;
	gstart[j] = pos[i];
	gend[j] = pos[i];  /* keep the end position of the region */
	if(chain[i]==1){
	  count1[j] = 1;   
	}else{
	  count2[j] = 1;
	}
      }
    }else{
      j = j+1;    
      chro[j] = chr[i];
      chroID[j] = i+1;
      gstart[j] = pos[i];
      gend[j] = pos[i];  /* keep the end position of the region */
      if(chain[i]==1){
	count1[j] = 1;   
      }else{
	count2[j] = 1;
      }
    }    
  }
  *nregion = j + 1; /*because R index is from 1:n */
}


/* f: foreground signal == IP enriched 
   b: background signal == control samples 
   input: chrf, posfStart, posfEnd, lenf(length of chrf), count1f (tag counts for chain 1),
          regsize - enriched region size, must equal to the size used in binning function
	  chrb,posb,chainb,lenb(length of chrb)
   output: countf - IP total count (for chain 1 and 2); coutf is also input
           count1b - count for backgroud chain 1 
           count2b - count for backgroud chain 2 
   the region can be changed. that is, posfEnd can be increased as long as posfEnd[i] <= posfStart[i]+regsize 
*/
void subBkg2(int *chrf,int *posfStart,int *posfEnd,int *count1f,int *lenf,int *regsize, int *maxtag,int *chrb,int *posb,
	    int *chainb,int *lenb,int *countf,int *count1b,int *count2b){
  int i, j,lenfm1,newEnd;
  i = 0;    /* ID for foreground */
  j = 0;    /* ID for background */
  lenfm1 = (*lenf) - 1;

  /* decide the end position for bin i */
  /* Whenever i is increased, it needs to calculate newEnd */
  if((count1f[i] == (*maxtag)) || ((countf[i] - count1f[i])==(*maxtag))){
    newEnd = posfEnd[i];
  }else{
    newEnd = posfStart[i]+(*regsize)-1; /* the bin size == regsize, note <= is changed to < in binning function */
  }
  while(j < (*lenb) && i < (*lenf)){
    /* if posb fall in the IP enriched region */
    if(chrb[j]==chrf[i]){
      if(posb[j]>=posfStart[i] && posb[j]<=newEnd){
	if(posb[j] > posfEnd[i]){
	  posfEnd[i] = posb[j];
	}
	if(countf[i] > 0){ /*total count truncate at 0 */
	  countf[i] -= 1;
	}
	if(chainb[j] == 1){
	  count1b[i] += 1;
	}else{
	  count2b[i] += 1;
	}
      }else if(posb[j]>newEnd && i<lenfm1){
	/* that is, if in the same chr,posb[j]<posfStart[i],except j++, do not do anything else */
	do{
	  i++;
	  if((count1f[i] == (*maxtag)) || ((countf[i] - count1f[i])==(*maxtag))){
	    newEnd = posfEnd[i];
	  }else{
	    newEnd = posfStart[i]+(*regsize)-1; /* the bin size == regsize, note <= is changed to < in binning function */
	  }
	  /*	  newEnd = posfStart[i] + (*regsize); */
	  if(posb[j]>=posfStart[i] && posb[j]<= newEnd && chrb[j]==chrf[i]){
	    if(posb[j] > posfEnd[i]){
	      posfEnd[i] = posb[j];
	    }
	    if(countf[i] > 0){
	      countf[i] -= 1;
	    }
	    if(chainb[j] == 1){
	      count1b[i] += 1;
	    }else{
	      count2b[i] += 1;
	    }
	  }
	}while(posb[j]>newEnd && chrb[j]==chrf[i] && i<(*lenf)); 
      } /*After the do loop chrf[i] = chrb[j] + 1 or posb[j] < newEnd, and posb has been count, need to j++ */
      j++; 
    }else if(chrb[j-1] != chrb[j]){ /*means chrb[j] > chrf[i] */
      do{
	i++;
	if((count1f[i] == (*maxtag)) || ((countf[i] - count1f[i])==(*maxtag))){
	  newEnd = posfEnd[i];
	}else{
	  newEnd = posfStart[i]+(*regsize)-1; /* the bin size == regsize, note <= is changed to < in binning function */
	}
      }while(chrb[j] != chrf[i] &&  i<(*lenf));
    }else if(chrf[i-1] != chrf[i]){ /* mean chrf[i] > chrb[j]*/
      do{j++;}while(chrb[j]!=chrf[i] && j<(*lenb));
    }
  }
  /*  Rprintf("\n IProw=%d,CONrow=%d, \n",i,j); */
}


void binning(int *chr,int *pos,int *chain,int *len,int *regsize,int *minsize, int *maxtag,int *chro,
	     int *gstart,int *gend,int *gend2,int *count1,int *count2,int *chroID,int *nregion){
  int i, j, jSize;  /* column start from 0 */

  /* copy the first row of x and initiate the number of reads to 1 at that point */
  j = 0;
  chro[j] = chr[j];
  chroID[j] = j+1;     /*because R index is from 1:n */
  gstart[j] = pos[j];      /* region start position */
  gend[j] = pos[j];  /* actual region end position */
  gend2[j] = pos[j] + (*regsize)-1;  /* predefined region end position */
  if(chain[j]==1){       /* count for positive strand in the region */
    count1[j] = 1;
  }else{
    count2[j] = 1;
  }
  
  jSize = *regsize;
  for(i=1; i < (*len); i++){
    if(chr[i] == chro[j]){ /* if the same chromosome */
      if((pos[i]-gstart[j]) < jSize){
        gend[j] = pos[i];  /* increase the end position of the region */
        if(chain[i]==1){
          count1[j] += 1;       /* count for positive strand in the region */
        }else{
          count2[j] += 1;   /* count for positive strand in the region */
        }
      }else{
        j = j+1;
        chro[j] = chr[i];
        chroID[j] = i+1;
        gstart[j] = pos[i];
        gend[j] = pos[i];  /* keep the end position of the region */
        if(chain[i]==1){
          count1[j] = 1;
        }else{
          count2[j] = 1;
        }
	if((count1[j-1]+count2[j-1]) >= (*maxtag)){
	  if(jSize >= (*minsize)*2){
	    jSize = jSize/2;
	  }
	}else{
	  if(jSize <= (*regsize/2)){
	    jSize = 2*jSize;
	  }
	}
	gend2[j] = pos[i] + jSize-1;
      }
    }else{
      j = j+1;
      chro[j] = chr[i];
      chroID[j] = i+1;
      gstart[j] = pos[i];
      gend[j] = pos[i];  /* keep the end position of the region */
      if(chain[i]==1){
        count1[j] = 1;
      }else{
        count2[j] = 1;
      }
      if((count1[j-1]+count2[j-1]) >= (*maxtag)){
	if(jSize >= (*minsize)*2){
	  jSize = (jSize/2);
	}
      }else{
	if(jSize <= (*regsize/2)){
	  jSize = 2*jSize;
	}
      }
      gend2[j] = pos[i] + jSize-1;
    }
  }
  *nregion = j + 1; /*because R index is from 1:n */
}

/* This function is used in mergetag function. Note added on 09/25/2012 */
void subBkg(int *chrf,int *posfStart,int *posfEnd, int *posfEnd2, int *count1f,int *lenf,int *regsize, int *maxtag,
	    int *chrb,int *posb,int *chainb,int *lenb,int *countf,int *count1b,int *count2b){
  int i, j,lenfm1;
  i = 0;    /* ID for foreground */
  j = 0;    /* ID for background */
  lenfm1 = (*lenf) - 1;

  while(j < (*lenb) && i < (*lenf)){
    /* if posb fall in the IP enriched region */
    if(chrb[j]==chrf[i]){
      if(posb[j]>=posfStart[i] && posb[j]<=posfEnd2[i]){
        if(posb[j] > posfEnd[i]){
          posfEnd[i] = posb[j];
        }
        if(countf[i] > 0){ /*total count truncate at 0 */
          countf[i] -= 1;
        }
        if(chainb[j] == 1){
          count1b[i] += 1;
        }else{
          count2b[i] += 1;
        }
      }else if(posb[j]>posfEnd2[i] && i<lenfm1){
        /* that is, if in the same chr,posb[j]<posfStart[i],except j++, do not do anything else */
        do{
          i++;
          if(posb[j]>=posfStart[i] && posb[j]<= posfEnd2[i] && chrb[j]==chrf[i]){
            if(posb[j] > posfEnd[i]){
              posfEnd[i] = posb[j];
            }
            if(countf[i] > 0){
              countf[i] -= 1;
            }
            if(chainb[j] == 1){
              count1b[i] += 1;
            }else{
              count2b[i] += 1;
            }
          }
        }while(posb[j]>posfEnd2[i] && chrb[j]==chrf[i] && i<(*lenf));
      } /*After the do loop chrf[i] = chrb[j] + 1 or posb[j] < newEnd, and posb has been count, need to j++ */
      j++;
    }else if(chrb[j-1] != chrb[j]){ /*means chrb[j] > chrf[i] */
      do{
        i++;
      }while(chrb[j] != chrf[i] &&  i<(*lenf));
    }else if(chrf[i-1] != chrf[i]){ /* mean chrf[i] > chrb[j]*/
      do{j++;}while(chrb[j]!=chrf[i] && j<(*lenb));
    }
  }
  /*  Rprintf("\n IProw=%d,CONrow=%d, \n",i,j); */
}


/* Input: chr, gstart, gend, count1 and count2 could be the counts for IP samples or subtracted IP samples
   rowID: the ID for the merged region, xrow=length(chr),maxgap is the maximum gap between the end of a region
   and the start of next region. 
   output: ochr,ogstart,ogend,orstart,orend,opeak(peak position - the center of between the peaks of + and - chain),
   changepoint:number of change points, nregion - number of merged region
*/
void mergeReg(int *chr,int *gstart,int *gend,int *count1,int *count2,int *rowID,int *xrow,int *maxgap,
	      int *ochr,int *ogstart,int *ogend,int *orstart,int *orend,int *opeak,int *changepoint,int *nregion){
 
  int i,j,cdiff0,cdiff1,maxpos1,maxpos2,idist,idistP,idistC,idistN;
  float maxden1,maxden2,iden1,iden2;
  int *istart,*iend; /*for keep track of the position of input vectors (eg. gstart, gend) for output region*/

  istart = (int *)R_alloc(*xrow,sizeof(int));
  iend = (int *)R_alloc(*xrow,sizeof(int));
  if((istart==NULL) || (iend==NULL)){
    error("Error: Fail to allocate memory for istart or iend!\n");
  }
  /* Just copy the data of the first region */
  ochr[0] = chr[0];
  ogstart[0] = gstart[0];
  ogend[0] = gend[0];
  orstart[0] = rowID[0];
  orend[0] = rowID[0];
  opeak[0] = (gstart[0]+gend[0])/2;
  cdiff0 = count1[0] - count2[0];
  maxden1 = 1.0*count1[0]/(gend[0]-gstart[0]+1);
  maxden2 = 1.0*count2[0]/(gend[0]-gstart[0]+1);
  maxpos1 = opeak[0];
  maxpos2 = opeak[0];
  changepoint[0] = 0;
  istart[0] = 0;
  iend[0] = 0;

  j = 0;
  for(i=1; i < *xrow; i++){
    idist = gend[i]-gstart[i]+1;
    iden1 = 1.0*count1[i]/idist;
    iden2 = 1.0*count2[i]/idist;
    if(chr[i] == ochr[j]){ /* if the same chromosome */
      if((gstart[i]-ogend[j]) < (*maxgap)){ /*beginning of the next region - the end of current region */
        ogend[j] = gend[i];
        orend[j] = rowID[i];
	iend[j] = i;
	cdiff1 = count1[i] - count2[i];
	/* count the number of change points */
	if(cdiff0 * cdiff1 <= 0){
	  changepoint[j] += 1;
	}
	cdiff0 = cdiff1;

	if(iden1 >= maxden1){
	  maxden1 = iden1;
	  maxpos1 = (gstart[i] + gend[i])/2;
	  opeak[j] = (maxpos1 + maxpos2)/2; 
	}

	if(iden2 > maxden2){
	  maxden2 = iden2;
	  maxpos2 = (gstart[i] + gend[i])/2;
	  opeak[j] = (maxpos1 + maxpos2)/2; 
	}
      }else {
	j = j+1;
	ochr[j] = chr[i];
	ogstart[j] = gstart[i];
	ogend[j] = gend[i];
	orstart[j] = rowID[i];
	orend[j] = rowID[i];
	opeak[j] = (gstart[i]+gend[i])/2;
	maxpos1 = opeak[j];
	maxpos2 = opeak[j];
	maxden1 = iden1;
	maxden2 = iden2;
	cdiff0 = count1[i] - count2[i];
	changepoint[j] = 0;
	istart[j] = i;
	iend[j] = i;
      }
    }else{
      j = j+1;    
      ochr[j] = chr[i];
      ogstart[j] = gstart[i];
      ogend[j] = gend[i];
      orstart[j] = rowID[i];
      orend[j] = rowID[i];
      opeak[j] = (gstart[i]+gend[i])/2;
      maxpos1 = opeak[j];
      maxpos2 = opeak[j];
      maxden1 = iden1;
      maxden2 = iden2;
      cdiff0 = count1[i] - count2[i];
      changepoint[j] = 0;
      istart[j] = i;
      iend[j] = i;
    }    
  }

  *nregion = j;
  /* recaculate the peak position using moving average approach */ 
  for(j=0; j<(*nregion+1);j++){
    if((iend[j]-istart[j])>=3){ /*redefine the peak for the regions with >=4 bins */
      idistC = gend[istart[j]]-gstart[istart[j]]+1;
      idistN = gend[istart[j]+1]-gstart[istart[j]+1]+1;
      maxden1 = (1.0*count1[istart[j]]/idistC + 1.0*count1[istart[j]+1]/idistN)/2.0;
      maxden2 = (1.0*count2[istart[j]]/idistC + 1.0*count2[istart[j]+1]/idistN)/2.0;
      maxpos1 = (gstart[istart[j]] + gend[istart[j]])/2;
      maxpos2 = maxpos1;
      idistP = idistC;
      idistC = idistN;
      for(i=(istart[j]+1); i<iend[j];i++){
	idistN = gend[i+1]-gstart[i+1]+1;
	iden1 = (1.0*count1[i-1]/idistP + 1.0*count1[i]/idistC + 1.0*count1[i+1]/idistN)/3.0;
	iden2 = (1.0*count2[i-1]/idistP + 1.0*count2[i]/idistC + 1.0*count2[i+1]/idistN)/3.0;
	if(iden1 >= maxden1){
	  maxden1 = iden1;
	  maxpos1 = (gstart[i] + gend[i])/2;
	}

	if(iden2 > maxden2){
	  maxden2 = iden2;
	  maxpos2 = (gstart[i] + gend[i])/2;
	}
	idistP = idistC;
	idistC = idistN;
      }
      /* now i = iend[j] */
      iden1 = (1.0*count1[i-1]/idistP + 1.0*count1[i]/idistC)/2.0;
      iden2 = (1.0*count2[i-1]/idistP + 1.0*count2[i]/idistC)/2.0;
      if(iden1 >= maxden1){
	maxden1 = iden1;
	maxpos1 = (gstart[i] + gend[i])/2;
      }
      
      if(iden2 > maxden2){
	maxden2 = iden2;
	maxpos2 = (gstart[i] + gend[i])/2;
      }
      opeak[j] = (maxpos1 + maxpos2)/2;
    }
  }

  *nregion = *nregion + 1;
  /* *nregion = j + 1; because R index is from 1:n */
}


void fdr(int *klen, double *kp, int *blen, double *beta, double *efdr){
  int k,b;
  int *J;
  J = (int *)R_alloc(*klen,sizeof(int));
  if(J==NULL){
    error("Error: Fail to allocate memory!\n");
  }

  for(k=0; k<(*klen); k++){
    J[k] = 0;
  }

  for(k=0; k<(*klen); k++){
    for(b=0; b<(*blen);b++){
      if(beta[b]<=kp[k]){
	J[k] = J[k] + 1;
	efdr[k] = efdr[k] + beta[b];
      }      
    }
    efdr[k] = efdr[k]/J[k];
  }
}
/* gap is the length of genomic gap between two regions 
   n is the length of ichr,istart,iend,ipos,fpos
   ichr, istart,iend is the input chromose, start and end positions
   ipos is the position of ichr, istart, iend in the augamented matrix, where the gap is filled with zero ct
   fpos is the position of the gaps in the augamented matrix
   fn is the total number of gaps needed to be filled

 */
void fillgap(int *gap, int *n, int *ichr, int *istart, int *iend, int *ipos, int *fpos, int *fn){
  int i,fid,newn;
  ipos[0] = 0;
  fid = 0;
  newn = 0;
  for(i=1; i<(*n); i++){
    newn += 1;
    if(ichr[i]==ichr[i-1]){
      if((istart[i]-iend[i-1])> (*gap)){
        fpos[fid] = newn;
        fid = fid + 1;
        newn += 1;
      }
    }
    ipos[i] = newn;
  }
  *fn = fid;
}



R_NativePrimitiveArgType iSeq2Args[16] = {INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,REALSXP,REALSXP,REALSXP,
					    REALSXP,REALSXP,REALSXP,INTSXP,REALSXP,REALSXP,INTSXP};
R_NativePrimitiveArgType iSeq1Args[19] = {INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,REALSXP,REALSXP,REALSXP,
					  REALSXP,REALSXP,INTSXP,REALSXP,REALSXP,REALSXP,REALSXP,
					  REALSXP,REALSXP,REALSXP,INTSXP};
R_NativePrimitiveArgType mergeRegArgs[16] =  {INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,
					      INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP};
R_NativePrimitiveArgType subBkgArgs[15] =  {INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,
					INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP};
R_NativePrimitiveArgType binningArgs[15] =  {INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,
					    INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP};
R_NativePrimitiveArgType fdrArgs[5] =  {INTSXP,REALSXP,INTSXP,REALSXP,REALSXP};
R_NativePrimitiveArgType fillgapArgs[8] =  {INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP};

static const R_CMethodDef CEntries[] = {
  {"iSeq2", (DL_FUNC)&iSeq2, 16, iSeq2Args},
  {"iSeq1", (DL_FUNC)&iSeq1, 19, iSeq1Args},
  {"mergeReg", (DL_FUNC)&mergeReg,16,mergeRegArgs},
  {"subBkg", (DL_FUNC)&subBkg,15,subBkgArgs},
  {"binning", (DL_FUNC)&binning,15,binningArgs},
  {"fdr", (DL_FUNC)&fdr, 5, fdrArgs},
  {"fillgap",(DL_FUNC)&fillgap,8,fillgapArgs},
  {NULL, NULL, 0}
}; 

void R_init_iSeq(DllInfo *dll){
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
}
