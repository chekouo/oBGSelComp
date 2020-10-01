OutcomeInit=function(outcome=outcome,X=X){
set.seed(1)
#library(sirt);
#library(glmnet);
y=outcome
nr=nrow(y);
q=ncol(y);
y0=apply(y,1,function(x) sum(x==0))
ind=seq(1,length(y0))[y0>0]
yno=y
if (length(ind)>0){
  yno[ind,]= (y[ind,]*(nr-1)+1/q)/nr;
}
ynew=(y*(nr-1)+1/q)/nr;
alphaS=matrix(0,nr,q)
for (i in 1:nr){
  y1=ynew[i,];
### We duplicate some data
modd=c(0,0.001,0.01,0.05);vv=matrix(runif(q*length(modd)),q)*modd;
 yy=rbind(yno[i,],(y1+vv)/(1+apply(vv,2,sum)))
  yi=matrix(y[i,],30,4,byrow = T)
  yy=rbind(yi,yy)
  alphaS[i,]=log(dirichlet.mle(yy,maxit=10000)$alpha);
}
alphaSS=alphaS;

#print(head(alphaS))
alphaSS
}
