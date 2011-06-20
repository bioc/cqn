#############################################################
##                                                         ##                      
##  R package mclust and nor1mix are required.             ##
##                                                         ##                      
#############################################################
##
##  THIS IS SQN allowing both tails to be extended. 
#############################################################
## Usage: SQN(y,ctrl.id,N.mix=5,model.weight=.9) 
##            y : matrix of data to be normalized
##      ctrl.id : index of control probes that the normalization is based on
##        N.mix : the number of normal mixture distribution used to estimate the control CDF
## model.weight : weight given to parametric normal mixture
##
##
## Example:
##
##
##  
##  mu0a=mu0b=rnorm(1000,5,1.8) ## control probes 
##  mu1a=rnorm(5000,5,0.8);mu1b=mu1a+2  ## probes for signal
##  ya=c(rnorm(1000,mu0a,.3),rnorm(5000,mu1a,.3))
##  yb=c(rnorm(1000,mu0b,.3),rnorm(5000,mu1b,.3))+1 ## a systematic bias introduced
##  Y=cbind(ya,yb)   #before normalization
##  Ynorm=SQN(Y,ctrl.id=1:1000)  #after normalization
##
##  par(mfrow=c(1,2))
##  boxplot(Y,main="before normalization")
##  boxplot(Y[1:1000,],add=T,col=3,boxwex=.4)

##  boxplot(Ynorm,main="after normalization")
##  boxplot(Ynorm[1:1000,],add=T,col=3,boxwex=.4)
##  legend(1,10,legend=c("probes for signal","negative control probes"),text.col=c(1,3),bg="white")
## LTrule=1: leave it unnormalized.
## 

SQN2 <-function(y,N.mix=5,ctrl.id,model.weight=.9,max.q=.95,min.q=.01,LTrule=1){
  QE=apply(y[ctrl.id,],2,sort)
  QN=apply(QE,1,median)
  mix.param=suppressWarnings(Mclust(QN,G=N.mix)$parameters)
  mix.param=norMix(mu=mix.param$mean,sig2=mix.param$variance$sigmasq,w=mix.param$pro)

  qq=seq(1/(2*length(QN)),1-1/(2*length(QN)),1/length(QN))
  qq=qnorMix(qq,mix.param)
  QN1=QN*(1-model.weight)+qq*model.weight
  
  ynorm=apply(y,2,mix.qn,ctrl.id,NQ=QN1,mix.param=mix.param,max.q=max.q,min.q=min.q,
        LTrule)
  model.par=c("N.mix"=N.mix,"model.weight"=model.weight,"max.q"=max.q,"min.q"=min.q,"LTrule"=LTrule)
  
  list("yout"=ynorm,"QN"=QN1,par=list("mix.param"=mix.param,"model.par"=model.par))
  
}  

mix.qn=function(y0,ctrl.id,NQ,mix.param,max.q=0.95,min.q=.01,LTrule=1){ 
  ##
  ##NQ:normalized quantiles
  ##
  ECDF=ecdf(y0[ctrl.id])
  Q=ECDF(y0) ## if outside the range, 0 or 1
  id0=which(Q<min.q)

  id2=which(Q>max.q)
  B=length(id2)
  Q2=max.q+(1-max.q)*(1:B)/(B+1)    
  y2o=order(y0[id2])
  Q[id2][y2o]=Q2
  
  ynorm=vector(mode="numeric",length=length(y0))
  #if (LTrule==1) ## keep original
  #  ynorm[id0]= y0[id0]
  #if (LTrule==2) ## use model to predict
  #  ynorm[id0]= qnorMix(Q[id0],mix.param)
  #if (LTrule==3) ## trucate
  #  ynorm[id0]= quantile(QN1,min.q)

  ynorm[id0]= y0[id0]-quantile(y0[ctrl.id],min.q)+quantile(NQ,min.q)
  ynorm[-c(id0,id2)]=quantile(NQ,Q[-c(id0,id2)])
  ynorm[id2]= qnorMix(Q[id2],mix.param)
  ynorm
}


