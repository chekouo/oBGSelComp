oBGSelCompAlgo=function(outcome=outcome,X=X,Pathway=Pathway,sample=10000, burnin=1000,
                        hyperprob=c(1,1),hypersig=c(5,0.01),hyperlamb=c(5,2), seed=1){
q=ncol(outcome);
K=nrow(Pathway);
n=nrow(outcome);
p=ncol(X);
outc1=as.vector(t(as.matrix(outcome)));
X=as.matrix(X);
outcome=as.matrix(outcome)
yinit1=as.vector(t(as.matrix(OutcomeInit(outcome=outcome,X=X))));
X1=as.vector(t(as.matrix(X)));
path1=as.vector(t(as.matrix(Pathway)));

a=hyperprob[1];b=hyperprob[2];
alpha0=hypersig[1];beta0=hypersig[2];
a.lambda=hyperlamb[1];b.lambda=hyperlamb[2]
result <- .C("mainMCMCFunction",q1=as.integer(q),K1=as.integer(K),n1=as.integer(n),p1=as.integer(p),
             outc1=outc1,yinit1=yinit1,X1=X1,path1=as.integer(path1),burninsample1=as.integer(burnin),nbrsample1=as.integer(sample),
             a=as.double(a),b=as.double(b),alpha01=as.double(alpha0),beta01=as.double(beta0),al1=as.double(a.lambda), 
             bll=as.double(b.lambda),seed1=as.integer(seed),gamMean1=as.vector(rep(0,q*K)),
             BetaSample1=as.vector(rep(0,q*sample*p)),PostPredSample1=as.double(rep(0,n*q*sample)),
             logpost=as.double(rep(0,sample+burnin)));

#double* a, double* b, double* alpha01, double* beta01
gamMean=matrix(result$gamMean1,q,K,byrow=T);
BetaSample=list();PostPredSample=array(0,c(n,q,sample))
l1=0;
for (l in 1:q){
BetaSample[[l]]=matrix(result$BetaSample1[(1+l1):(sample*p+l1)],sample,p,byrow=T)
l1=l1+sample*p;
}
l2=0;
for (i in 1:n){
PostPredSample[i,,]=matrix(result$PostPredSample1[(1+l2):(sample*q+l2)],q,sample,byrow=T)
l2=l2+sample*q;
}
outc=matrix(outc1,n,q,byrow = T)
return(list(PostProbGrp=gamMean,BetaSample=BetaSample,PostPredSample=PostPredSample,logpost=result$logpost,outcome=outc))
}

