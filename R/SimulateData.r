SimulateData<-function(n=100,betaAbs=1,sig=.1,r2=0.5,propOverl=0,seed=1,K=20,p=200,q=4,nbrpath=4){
  ##See Section 5.2 manuscript
#library(gtools)
if (!propOverl%in%c(0,10,50)){
    stop("We should enter for propOverl 0, 10 or 50.")
}
  if (p>200){
    if ( propOverl!=0){
      print("Warning: When p>200, we only provided the option to generate non-overlap groups")
    }
    propOverl=0
  }
#K=20;p=200;q=4;nbrpath=4;
p1=p;
set.seed(1)
PathKnow=matrix(0,q,K)
PathKnow[,1:nbrpath]=1
mRNA=matrix(0,n,p1);
Path0=BuildPath(p1,K);
Path=BuildPath(p1, K, propOverl)
uk=sqrt(r2)*rnorm(K)
for (j in 1:p1){
  mRNA[,j]=c(t(Path0[,j])%*%(Path0[,j]*uk))+rnorm(n)
}
meaM=apply(mRNA,2,mean)
sdM=apply(mRNA,2,sd)
mRNA=t((t(mRNA)-meaM)/sdM)
Beta=matrix(0,ncol(mRNA),q)
Rgene=rep(0,ncol(mRNA))
for (l in 1:q){
  selpath=seq(1,K)[PathKnow[l,]==1]
for (k in 1: length(selpath)){
  selGene=seq(1,p1)[Path[selpath[k],]==1]
  unif=runif(1)
  if ((propOverl==10)||(propOverl==50)){
    if (k==1){
      sel=rep(0,2)
      sel[1]=selGene[1];sel[2]=selGene[length(selGene)];
    }
    if (k==2){
      sel=rep(0,2)
      selGene1=which(Path[selpath[k-1],]*Path[selpath[k],]==1)
      sel[1]=selGene1[length(selGene1)];sel[2]=selGene1[length(selGene1)]+1;
    }
    if (k==3){
      sel=rep(0,2)
      selGene1=which(Path[selpath[k-1],]*Path[selpath[k],]==1)
      sel[1]=selGene1[length(selGene1)]+1;sel[2]=selGene[length(selGene)];
    }
    if (k==4){
      sel=rep(0,3)
      selGene1=which(Path[selpath[k-1],]*Path[selpath[k],]==1)
      sel[1]=selGene1[length(selGene1)];sel[2]=selGene1[length(selGene1)]+1;sel[3]=selGene[length(selGene)];
    }
  } else if (propOverl==0){
    sel=c(selGene[1],selGene[length(selGene)])
  }
  if (unif<0.5){
    Beta[sel,l]=betaAbs
  } else {
    Beta[sel,l]=-betaAbs
  }
  Rgene[sel]=1
}
}
alpha=matrix(0,nrow(mRNA),q)
for (l in 1:q){
  alpha[,l]=as.vector(mRNA%*%t(t(Beta[,l])))+sig*rnorm(nrow(mRNA))
}
set.seed(seed)
  y=matrix(0,nrow(mRNA),q)
  for (i in 1:nrow(mRNA)){
    y[i,]=rdirichlet(1,exp(alpha)[i,])
  }
list(y=y,X=mRNA,Pathway=Path,GamKnown=PathKnow,Beta=Beta)
}
