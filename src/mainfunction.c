#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include "header.h"

/* Compute logposterior*/

void logposterior(int n,int p,int K, int NbrOut,int ** njg,double ** alphaS, double ** outc,double * s2, double ** beta,
             double *beta0l,double ** X, double *** Tau,double *lambda2,_Bool** gampath,double * q,
             double ab, double bb,double al, double *bl,double alpha0, double beta0,double sigma20,double *logpost){
  double logp=0;
  for (int i=0;i<n;i++){
    double AlphaI[NbrOut]; 
    for (int l=0;l<NbrOut;l++){
      AlphaI[l]=exp(alphaS[l][i]);
    }
    logp+=gsl_ran_dirichlet_lnpdf(NbrOut, AlphaI,outc[i]);
    for (int l=0;l<NbrOut;l++){
      double mu=0;
     for (int j=0;j<p;j++){
        if (njg[l][j]>0)
        mu+=X[i][j]*beta[l][j];
      }
      mu+=beta0l[l];
      logp+=log(gsl_ran_gaussian_pdf (alphaS[l][i]-mu,sqrt(s2[l])));
    }
  }
  /*
  for (int l=0;l<NbrOut;l++){
    for (int j=0;j<p;j++){
      if  (njg[l][j]>0){
        //Prior for beta[l,j]
        logp+=log(gsl_ran_gaussian_pdf (beta[l][j],sqrt(s2[l]*Tau[l][j][1]))); 
        //Prior for Tau2
        logp+=log(lambda2[l]/2)-Tau[l][j][1]*lambda2[l]/2;
      }
    }
    //Prior for gam_l
    for (int k=0;k<K;k++){
    logp+=gampath[l][k]*log(q[l])+(1-gampath[l][k])*log(1-q[l]);
    }
    
    //Prior for p_l (or q_l)
    logp+=log(gsl_ran_beta_pdf(q[l], ab, bb));
    //Prior for lambda2_l
    logp+=log(gsl_ran_gamma_pdf(lambda2[l], al, 1/bl[l]));
    //Prior for sigma2_l
    logp+=-log(pow(s2[l],2))+log(gsl_ran_gamma_pdf(1/s2[l], alpha0, 1/beta0));
    //Prior for beta0_l
    logp+=log(gsl_ran_gaussian_pdf (beta0l[l],sqrt(sigma20)));
  }
   */
    *logpost=logp;
}

/*Generate q*/
double geneBeta(gsl_rng * r,int K,_Bool * gampath,double ab, double bb){
double ab1=0;
int k;
for (k=0;k<K;k++){
ab1+=gampath[k];
}
double bb1=K-ab1;
return gsl_ran_beta (r, ab1+ab,  bb1+bb);
}
//***** Generate Tau

void geneTau(gsl_rng * r, int p, int K, int * njg, double ** Tau, double * beta, double s2,double lambda2){
double lamb=lambda2; double mu=0;
int j;
//printf("\nLambda=%lf\n ",lambda2);
//printf("\nTau=\n");
for (j=0;j<p;j++){
if ((njg[j]>0)&&(fabs(beta[j])>pow(10,-10))){
//if ((njg[j]>0)){
mu=sqrt(lambda2*s2)/fabs(beta[j]); 
Tau[j][1]=1.0/inverseGaussian(r,  mu, lamb);

if (Tau[j][1]<0){
//printf("Tau=%lf\n",Tau[j][1]); 
}
if (((Tau[j][1]-Tau[j][1])!=0)||(Tau[j][1]<0)){
Tau[j][1]=1/mu+1/lamb;
//printf("J=%d, MU=%lf, lamb=%lf, S2=%lf, Beta=%lf", j,mu,lamb,s2,beta[j]);
}
} else {
Tau[j][1]=gsl_ran_exponential (r, 2.0/lambda2);
}
//printf("%lf ",Tau[j][1]);
}


}

double geneLambda(gsl_rng * r, int p, int * njg, double ** Tau, double al, double bl){
int j;
double bl1=bl;
double al1=al;
int p1=0;
for (j=0;j<p;j++){
if (njg[j]>0){
bl1+=Tau[j][1]/2;
p1+=1;
}
}
al1+=p1;
return gsl_ran_gamma (r, al1, 1/bl1)+0.0001;
}

double * Sigmag(int k,int n,int * njg,double ** tau,double **Ps,int *IDX){
double *Sig=malloc(n*n*sizeof(double));
int i,j,l;

for (i=0;i<n;i++){
for (j=0;j<=i;j++){
double a=0;
for (l=0;l<k;l++){
int l1=0;
if (njg[IDX[l]]>0)
l1=1;
else 
l1=0;
a+=tau[IDX[l]][l1]*Ps[i][l]*Ps[j][l];
}
if (i==j){
a+=1.0;
}
Sig[i*n+j]=Sig[j*n+i]=a;
//printf("%f %d %d \n",Ainv[i*k+j],i,j);
}
}
return Sig;
}



//double logmarggam1(int n, p,int *nx,double * Yk,double ** X,double * Xy,int * IDX, double ** Tau, gsl_matrix * A,_Bool* path,double alpha0, double beta0){

gsl_vector *Muk(int n, int np,const gsl_matrix * L,double ** X,double * Yk,int *IDX,double *Xy){
int i,j;
int np1=np;
if (np==0) np1=1;
gsl_vector * muk=gsl_vector_alloc (np1);
if (np>0){ 
for (j=0; j<np;j++){
Xy[j]=0;
for (i=0; i<n;i++){
Xy[j]+=X[i][IDX[j]]*Yk[i];
}
}
gsl_vector_view b= gsl_vector_view_array (Xy, np);
gsl_linalg_cholesky_solve (L, &b.vector, muk);
} else {
gsl_vector_set_zero (muk);
}
 return muk;
}

GSLMatVect logmarggam1(int n, int p, int *np,double * Yk,double ** X,int * IDX,_Bool * path,int * njg, double ** Tau,double alpha0, double beta0, double *logmarg){
int i,j;
GSLMatVect LMu;
double logdet=0;
double logRR=0;
double XyMu=0;
double ** PG=PGG(n,p,path,X,njg, IDX,np);
double * Ainv= Ainbeta(*np,n,njg,Tau,PG,IDX);
//gsl_matrix * L;gsl_vector *Mu;
if (*np>0){
/*
 * printf("Tau=\n");
for (j=0;j<*np;j++){
printf("%d \n",IDX[j]);
printf("%lf ",Tau[IDX[j]][1]);
}
printf("\nAinv=\n");
int j1;
for (j=0; j<*np;j++){
for (j1=0; j1<=j;j1++){
printf("A=%lf J=%d, J1=%d", Ainv[j+j1* (*np)],j,j1);
}
printf("\n");
}
*/

gsl_matrix_view m = gsl_matrix_view_array (Ainv, *np,*np);
gsl_linalg_cholesky_decomp (&m.matrix);
for (j=0; j<*np;j++){
int jj=IDX[j];
logdet+=2*log(gsl_matrix_get(&m.matrix,j,j))+log(Tau[jj][1]);
}
double *Xy=malloc(*np*sizeof(double));
//gsl_vector *muk=Muk(n, *np,&m.matrix,X,Yk,IDX,Xy);
gsl_vector * muk=gsl_vector_alloc (*np);
for (j=0; j<*np;j++){
Xy[j]=0;
for (i=0; i<n;i++){
Xy[j]+=X[i][IDX[j]]*Yk[i];
}
}
gsl_vector_view b= gsl_vector_view_array (Xy, *np);
gsl_linalg_cholesky_solve (&m.matrix, &b.vector, muk);

for (j=0;j<*np;j++){
XyMu+=Xy[j]*gsl_vector_get(muk,j);
}
free(Xy);
LMu.Mat=gsl_matrix_alloc (*np, *np);
gsl_matrix_memcpy (LMu.Mat, &m.matrix);
LMu.Vec=gsl_vector_alloc (*np);
gsl_vector_memcpy (LMu.Vec, muk);
gsl_vector_free(muk);
} else {
LMu.Mat=gsl_matrix_alloc (1, 1);LMu.Vec=gsl_vector_alloc (1);
gsl_matrix_set_zero (LMu.Mat);gsl_vector_set_zero (LMu.Vec);
}

for (i=0; i<n;i++){
logRR+=Yk[i]*Yk[i];
}
logRR=-0.5*logdet-((n+2*alpha0)/2.0)*log(0.5*(logRR-XyMu)+beta0);;
*logmarg=logRR;
free(Ainv);free_dmatrix(PG,0,n-1,0,*np-1);
//gsl_matrix_memcpy (LMu.Mat, L);gsl_vector_memcpy (LMu.Vec, Mu);
return LMu;
}

void gampat1(int n, int p, int K, double * y,double ** X, gsl_rng * r,_Bool ** path,int * njg,_Bool * gampath, double ** Tau, double * beta, double alpha0, double beta0,double q,double s2){
int i,j,k;
double Yk[n];
_Bool gampath1[K];
for (k=0;k<K;k++){
gampath1[k]=gampath[k];
}
int ** IDX=imatrix(0,K-1,0,p-1);
int * np=ivector(0,K-1);
int ** IDX2=imatrix(0,K-1,0,p-1);
int * np2=ivector(0,K-1);
double log1,log2;
for (k=0;k<K;k++){
for (i=0; i<n;i++){
double xb=0;
for (j=0; j<p;j++){
if ((path[k][j]==0)&&(njg[j]>0))
xb+=X[i][j]*beta[j];
}
Yk[i]=y[i]-xb;
}
//log1=logmarggam(n,p,&np[k],Yk,X,IDX[k],path[k],njg,Tau,alpha0,beta0);
/* Newwwwwww ****/
//if (np[k]>0) 
//printf("NPP=%d\n",np[k]);

//double logmarggam1(int n, int nx,double * Yk,double * Xy,int * IDX, double ** Tau, const gsl_matrix * A,double *muk,double alpha0, double beta0)
GSLMatVect LMu= logmarggam1(n, p, &np[k],Yk,X,IDX[k],path[k],njg,Tau, alpha0,  beta0, &log1);
gampath1[k]=1-gampath[k];
NJG(p,K,path,gampath1,njg);
//log2=logmarggam(n,p,&np2[k],Yk,X,IDX2[k],path[k],njg,Tau,alpha0,beta0);
GSLMatVect LMunew= logmarggam1(n, p, &np2[k],Yk,X,IDX2[k],path[k],njg,Tau, alpha0,  beta0, &log2);
double rat=log1-log2;
double ratT=0;
_Bool a=0;_Bool b=0;
if (gampath[k]==1){
a=1;
ratT=log(q)-log(1-q)+rat;
} else {
ratT=log(q)-log(1-q)-rat;
}
double  un=gsl_ran_flat (r, 0, 1);
if (log(un)-log(1-un)<ratT){
b=1;
gampath[k]=1;
} else gampath[k]=0;
if (a==b) NJG(p,K,path,gampath,njg);
gampath1[k]=gampath[k];
/*** We update beta coeffocient******///
int npk=np[k];
if (a!=b) npk=np2[k];
if (npk>=1){
gsl_vector *result=gsl_vector_alloc (npk);
for (j=0; j<npk; j++) gsl_vector_set(result, j, sqrt(s2)*gsl_ran_ugaussian(r) );
if (a==b){
gsl_blas_dtrsv (CblasLower, CblasTrans, CblasNonUnit,LMu.Mat, result);
for (j=0; j<npk; j++){
int id=IDX[k][j];
beta[id]= gsl_vector_get(LMu.Vec,j)+gsl_vector_get(result,j);
}
} else{
gsl_blas_dtrsv (CblasLower, CblasTrans, CblasNonUnit,LMunew.Mat, result);
for (j=0; j<npk; j++){
int id=IDX2[k][j]; 
beta[id]= gsl_vector_get(LMunew.Vec,j)+gsl_vector_get(result,j);
}
}
gsl_vector_free(result);
} 
gsl_vector_free(LMunew.Vec);gsl_matrix_free(LMunew.Mat);
gsl_vector_free(LMu.Vec);gsl_matrix_free(LMu.Mat);
int npp;
findc(p,path[k],0,IDX[k], &npp);
for (j=0; j<npp; j++){ 
if (njg[IDX[k][j]]==0) beta[IDX[k][j]]=0; 
}

}

/*for (j=0; j<p;j++){
if (njg[j]==0)
beta[j]=0;
}
*/
free(np);free(np2);
free_imatrix(IDX,0,K-1,0,p-1);free_imatrix(IDX2,0,K-1,0,p-1);
}











double logmarggam(int n, int p, int *np,double * Yk,double ** X,int * IDX,_Bool * path,int * njg, double ** Tau,double alpha0, double beta0){
double logRR=0;
int i;
double ** PG=PGG(n, p,path,X,njg, IDX,np);
double *SigmaGam=Sigmag(*np,n,njg,Tau,PG,IDX);
gsl_matrix_view m
    = gsl_matrix_view_array (SigmaGam, n,n);
if (*np>0){
gsl_linalg_cholesky_decomp (&m.matrix);
gsl_vector_view b= gsl_vector_view_array (Yk, n);
gsl_vector *xr=gsl_vector_alloc (n);
gsl_linalg_cholesky_solve (&m.matrix, &b.vector, xr);

double logdet=0;
for (i=0; i<n;i++){
logRR+=0.5*Yk[i]*gsl_vector_get(xr,i);
logdet+=2*log(gsl_matrix_get(&m.matrix,i,i));
}
logRR=-0.5*logdet-((n+2*alpha0)/2)*log(logRR+beta0);
gsl_vector_free(xr);
} else {
for (i=0; i<n;i++){
logRR+=0.5*Yk[i]*Yk[i];
}
logRR=-((n+2*alpha0)/2)*log(logRR+beta0);
}
free_dmatrix(PG,0,n-1,0,*np-1);
free(SigmaGam);
return logRR;
}
double sigma2(gsl_rng * r,int n,int p,int *njg,double ** tau, double *beta,double **X,double *y, double alpha0, double beta0,int sample){
double alphag=alpha0+(double)n/2.0;
double bet=0;
double tauSum=0;
int p1=0;int i,j;
for (j=0; j<p;j++){
if ((njg[j]>0)&& fabs(beta[j])>0){
p1+=1;
bet+=pow(beta[j],2)/tau[j][1];
tauSum+=tau[j][1];
}
}
bet=bet/2;
alphag=alphag+(double)p1/2.0;
double yx=0;
for (i=0;i<n;i++){
double xb=0;
for (j=0; j<p;j++){
xb+=X[i][j]*beta[j];
}
yx+=pow(y[i]-xb,2);
}
yx=yx/2;
yx=yx+bet+beta0;

double inv=1/yx;
//if (inv<pow(10,-10))
//inv=pow(10,-10);
double gam=gsl_ran_gamma (r, alphag, inv);
//if (sample>18300){
//printf("Alphag=%lf , YX=%lf ,InvSIGMA=%lf\n",alphag,yx,gam);
//}
return  1/gam;
}

void gampat(int n, int p, int K, double * y,double ** X, gsl_rng * r,_Bool ** path,int * njg,_Bool * gampath, double ** Tau, double * beta, double alpha0, double beta0,double q,double s2){
int i,j,k;
double Yk[n];
_Bool gampath1[K];
for (k=0;k<K;k++){
gampath1[k]=gampath[k];
}
int ** IDX=imatrix(0,K-1,0,p-1);
int * np=ivector(0,K-1);
int ** IDX2=imatrix(0,K-1,0,p-1);
int * np2=ivector(0,K-1);
double log1,log2;
for (k=0;k<K;k++){
for (i=0; i<n;i++){
double xb=0;
for (j=0; j<p;j++){
if ((path[k][j]==0)&&(njg[j]>0))
xb+=X[i][j]*beta[j];
}
Yk[i]=y[i]-xb;
}
log1=logmarggam(n,p,&np[k],Yk,X,IDX[k],path[k],njg,Tau,alpha0,beta0);
gampath1[k]=1-gampath[k];
NJG(p,K,path,gampath1,njg);
log2=logmarggam(n,p,&np2[k],Yk,X,IDX2[k],path[k],njg,Tau,alpha0,beta0);
double rat=log1-log2;
double ratT=0;
if (gampath[k]==1){
ratT=log(q)-log(1-q)+rat;
} else {
ratT=log(q)-log(1-q)-rat;
}
double  un=gsl_ran_flat (r, 0, 1);
if (log(un)-log(1-un)<ratT){
gampath[k]=1;
} else gampath[k]=0;
NJG(p,K,path,gampath,njg);
gampath1[k]=gampath[k];
/*** We update beta coeffocient******///
double ** PG=PGG(n,p,path[k],X,njg, IDX[k],&np[k]);
if (np[k]>=1){
double * Ainv= Ainbeta(np[k],n,njg,Tau,PG,IDX[k]);
gsl_matrix_view m = gsl_matrix_view_array (Ainv, np[k],np[k]);
    gsl_linalg_cholesky_decomp (&m.matrix);
double xy[np[k]];
for (j=0; j<np[k];j++){
xy[j]=0;
for (i=0; i<n;i++){
xy[j]+=X[i][IDX[k][j]]*Yk[i];
}
}
gsl_vector_view b= gsl_vector_view_array (xy, np[k]);
gsl_vector *xr=gsl_vector_alloc (np[k]);
gsl_linalg_cholesky_solve (&m.matrix, &b.vector, xr);
gsl_vector *result=gsl_vector_alloc (np[k]);
for (j=0; j<np[k]; j++)
        gsl_vector_set(result, j, gsl_ran_ugaussian(r) );
gsl_blas_dtrsv (CblasLower, CblasTrans, CblasNonUnit,&m.matrix, result);
for (j=0; j<np[k]; j++){
beta[IDX[k][j]]= gsl_vector_get(xr,j)/s2+gsl_vector_get(result,j);
}

free(Ainv);
gsl_vector_free(result);gsl_vector_free(xr);
} 
int npp;
findc(p,path[k],0,IDX[k], &npp);
for (j=0; j<npp; j++){ 
if (njg[IDX[k][j]]==0) beta[IDX[k][j]]=0; 
}

free_dmatrix(PG,0,n-1,0,np[k]-1);
}

/*for (j=0; j<p;j++){
if (njg[j]==0)
beta[j]=0;
}
*/
free(np);free(np2);
free_imatrix(IDX,0,K-1,0,p-1);free_imatrix(IDX2,0,K-1,0,p-1);
}

double * Ainbeta(int k,int n,int * njg,double ** tau,double **Ps,int * IDX){
double *Ainv=malloc(k*k*sizeof(double));
int i,j,l;

for (i=0;i<k;i++){
for (j=0;j<=i;j++){
double a=0;
for (l=0;l<n;l++)
a+=Ps[l][i]*Ps[l][j];
if (i==j){
a+=(1.0/tau[IDX[i]][1]);
//printf("A===%lf ",a);
}

Ainv[i*k+j]=Ainv[j*k+i]=a;
//printf("%f %d %d \n",Ainv[i*k+j],i,j);
}
}
return Ainv;
}

double **PGG(int n, int p,_Bool * path,double **X,int *njg, int *IDX, int *np){
int i,j;
_Bool *path1=malloc(p*sizeof(_Bool));
for (j=0;j<p;j++){
path1[j]=path[j];
if (njg[j]==0){
path1[j]=0;
}
}
findc(p,path1,0,IDX, np);
int np1=*np;
//printf(" X=%d",np); 
double **PG=dmatrix(0,n-1,0,np1-1);
for (i=0;i<n;i++){
for (j=0;j<np1;j++){
PG[i][j]=X[i][IDX[j]];
//printf("%lf ",PG[i][j]);
}
}
free(path1);
return PG;
}
/* Alpha parameters for dirichlet*/

 static double lanczos_7_c[9] = {
     0.99999999999980993227684700473478,
     676.520368121885098567009190444019,
    -1259.13921672240287047156078755283,
     771.3234287776530788486528258894,
    -176.61502916214059906584551354,
     12.507343278686904814458936853,
    -0.13857109526572011689554707,
     9.984369578019570859563e-6,
     1.50563273514931155834e-7
  };

double logGamma (int N,double* alpha){
int g=7;
int i;
double maxalpha=MAX(max(N,alpha),log(g-1+0.5));
double logsum=0;
double sumealpha=0;
for (i=0;i<N;i++){
logsum+=exp(alpha[i]-maxalpha);
sumealpha+=exp(alpha[i]);
}
logsum+=exp(log(g-1+0.5)-maxalpha);
logsum=(sumealpha-1+0.5)*(maxalpha+log(logsum));
int k;
double Ag = lanczos_7_c[0];
 for(k=1; k<=8; k++) { Ag += lanczos_7_c[k]/(sumealpha-1+k); }
return logsum+log(Ag)-(sumealpha-1+g+0.5)+0.5*log(2*PI);;
}

void alphaStar1(gsl_rng * r,int n,int p, int NbrOut,double ** alphaS, double ** outc,double * s2, double ** beta,double *beta0l,double ** X,int *acceptRate,int sample){
int i,l,j;

double sig0=0.5;
for (i=0;i<n;i++){
double sumalp=0;
for (l=0;l<NbrOut;l++){
alphaS[l][i]=alphaS[l][i]+beta0l[l];
//AlphaI[l]=alphaS[l][i];
//AlphaINew[l]=alphaS[l][i];
sumalp+=exp(alphaS[l][i]);
}
//printf("%lf ",sumalp);
for (l=0;l<NbrOut;l++){
double mu=0;
for (j=0;j<p;j++){
mu+=X[i][j]*beta[l][j];
}
mu+=beta0l[l];
double alphaold=exp(alphaS[l][i]);
//printf("alphaOld=%lf\n", alphaold);
double logalphaold=alphaS[l][i];
//double logpostold= gsl_sf_lngamma (sumalp)-gsl_sf_lngamma (exp(logalphaold))+(exp(logalphaold)-1)*log(outc[i][l])-0.5*log(s2[l])-0.5*pow(logalphaold-mu,2)/s2[l];
double meannew=s2[l]*alphaold*(gsl_sf_psi(sumalp)-gsl_sf_psi(alphaold)+log(outc[i][l]))+mu;
//printf("II=%d,LL=%d,logOut=%lf\n",i,l,log(outc[i][l]));
double logalphanew=meannew+sig0*gsl_ran_ugaussian(r);

double alphanew=exp(logalphanew);
//printf("logalphaNew=%lf\n", logalphanew);
//AlphaINew[l]=logalphanew;
double sumalpn=sumalp-alphaold+alphanew;
double meanold=s2[l]*alphanew*(gsl_sf_psi(sumalpn)-gsl_sf_psi(alphanew)+log(outc[i][l]))+mu;
double logpostnew= gsl_sf_lngamma (sumalpn)-gsl_sf_lngamma (alphanew)+(alphanew-1)*log(outc[i][l])-0.5*log(s2[l])-0.5*pow(logalphanew-mu,2)/s2[l];
double lognum=log( gsl_ran_gaussian_pdf (logalphaold-meanold,sig0))+logpostnew;
double logpostold= gsl_sf_lngamma (sumalp)-gsl_sf_lngamma (alphaold)+(alphaold-1)*log(outc[i][l])-0.5*log(s2[l])-0.5*pow(logalphaold-mu,2)/s2[l];
double logdeno=log( gsl_ran_gaussian_pdf (logalphanew-meannew,sig0))+logpostold;
double accept=lognum-logdeno;
double uni=gsl_ran_flat (r, 0, 1);
if (log(uni)<accept){
alphaS[l][i]=logalphanew;
sumalp=sumalpn;
*acceptRate=*acceptRate+1;
//printf("%lf",alphaS[l][i]);
} 
}
for (l=0;l<NbrOut;l++){
//if (l==1)
//printf("all= %lf\n ",alphaS[l][i]);
alphaS[l][i]=alphaS[l][i]-beta0l[l];
}
}

}



void alphaStar2(gsl_rng * r,int n,int p, int NbrOut,double ** alphaS, double ** outc,double * s2, double ** beta,double *beta0l,double ** X,int *acceptRate,int sample){
int i,l,j;

int niter=20;
double sig20=0.5;
for (i=0;i<n;i++){
double sumalp=0;
for (l=0;l<NbrOut;l++){
alphaS[l][i]=alphaS[l][i]+beta0l[l];
sumalp+=exp(alphaS[l][i]);
}
//printf("%lf ",sumalp);
for (l=0;l<NbrOut;l++){
double mu=0;
for (j=0;j<p;j++){
mu+=X[i][j]*beta[l][j];
}
mu+=beta0l[l];
double alphaold=exp(alphaS[l][i]);
//printf("alphaOld=%lf\n", alphaold);
//printf("SUMALPHA=%lf\n", sumalp);
//printf("II=%d,LL=%d,logOut=%lf\n",i,l,log(outc[i][l]));
double meannew=alphaold;
int it=0;
while (it<niter){
//printf("Iter=%d,MeanNew=%lf\n", it,meannew);
//printf("Mu=%lf, Rat=%lf, AA=%lf\n",mu,-(log(meannew)-mu)/(s2[l]*meannew),gsl_sf_psi(sumalp)+log(outc[i][l])-(log(meannew)-mu)/(s2[l]*meannew));
meannew=InvDigamma(gsl_sf_psi(sumalp)+log(outc[i][l])-1/meannew-(log(meannew)-mu)/(s2[l]*meannew));
it+=1;
}
//printf("II=%d,LL=%d,logOut=%lf\n",i,l,log(outc[i][l]));
//printf("MeanNew=%lf\n", meannew);
sig20=-1/(gsl_sf_psi_1(sumalp)-gsl_sf_psi_1(meannew)+1/pow(meannew,2)+(log(meannew)-mu-1)/(s2[l]*pow(meannew,2)));
//printf("Sigma2=%lf\n", sig20);
alphaS[l][i]=log(Truncate(meannew,sqrt(sig20), 0.0, r));
sumalp=sumalp-alphaold+exp(alphaS[l][i]);
}
for (l=0;l<NbrOut;l++){
alphaS[l][i]=alphaS[l][i]-beta0l[l];
}
}

}


void alphaStar(gsl_rng * r,int n,int p, int NbrOut,double ** alphaS, double ** outc,double * s2, double ** beta,double *beta0l,double ** X,int *acceptRate,int sample){
int i,l,j;

//double sig0=0.1;
double sig0=0.4;
for (i=0;i<n;i++){
//double sumalp=0;
double AlphaI[NbrOut],AlphaINew[NbrOut];
for (l=0;l<NbrOut;l++){
alphaS[l][i]=alphaS[l][i]+beta0l[l];
AlphaI[l]=alphaS[l][i];
AlphaINew[l]=alphaS[l][i];
//sumalp+=exp(alphaS[l][i]);
}
//printf("%lf ",sumalp);
for (l=0;l<NbrOut;l++){
double mu=0;
for (j=0;j<p;j++){
mu+=X[i][j]*beta[l][j];
}
mu+=beta0l[l];
double logalphaold=alphaS[l][i];
//double logpostold= gsl_sf_lngamma (sumalp)-gsl_sf_lngamma (exp(logalphaold))+(exp(logalphaold)-1)*log(outc[i][l])-0.5*log(s2[l])-0.5*pow(logalphaold-mu,2)/s2[l];
double logalphanew=logalphaold+sig0*gsl_ran_ugaussian(r);


double alphan=exp(logalphanew);
AlphaINew[l]=logalphanew;
//if (l==0){
//if (sample>18200){
//printf("YX1=%lf ,Exp1=%lf",logalphanew, alphan);
//printf("YX2=%lf ,Exp2=%lf",alphaS[l][i], exp(logalphaold));
//}
//}
//double sumalpn=sumalp-exp(logalphaold)+alphan;
//double logpostnew= gsl_sf_lngamma (sumalpn)-gsl_sf_lngamma (alphan)+(alphan-1)*log(outc[i][l])-0.5*log(s2[l])-0.5*pow(logalphanew-mu,2)/s2[l];
//double accept=gsl_sf_lngamma (sumalpn)-gsl_sf_lngamma (sumalp)-gsl_sf_lngamma (alphan)+gsl_sf_lngamma (exp(logalphaold))+(alphan-exp(logalphaold))*log(outc[i][l])-0.5*pow(logalphanew-mu,2)/s2[l]+0.5*pow(logalphaold-mu,2)/s2[l];
//printf("Accept1= %lf\n ",accept);
double accept=logGamma(NbrOut,AlphaINew)-logGamma(NbrOut,AlphaI)-logGamma(1,&logalphanew)+logGamma(1,&logalphaold)+(alphan-exp(logalphaold))*log(outc[i][l])-0.5*pow(logalphanew-mu,2)/s2[l]+0.5*pow(logalphaold-mu,2)/s2[l];
//printf("Accept2= %lf\n ",accept);

//printf("%lf ",accept);
//if (accept>0)
//accept=0;
double uni=gsl_ran_flat (r, 0, 1);
if (log(uni)<accept){
alphaS[l][i]=logalphanew;
AlphaI[l]=alphaS[l][i];
//sumalp=sumalpn;
*acceptRate=*acceptRate+1;
//printf("%lf",alphaS[l][i]);
} else {
AlphaINew[l]=alphaS[l][i];
}
}
for (l=0;l<NbrOut;l++){
//if (l==1)
//printf("all= %lf\n ",alphaS[l][i]);
alphaS[l][i]=alphaS[l][i]-beta0l[l];
}
}

}

 
void Beta0(gsl_rng *r,int n,int p,double *betao, double s2,double sigma20,double * beta,double *y,double ** X,int t){
int i,j;
double sigma2beta=1/(n/s2+1/sigma20);
//if (t>18200){
//printf("SigmaBeta0=%lf ,S2=%lf",sigma2beta, s2);
//}
double meany=0;
double newy[n];
for (i=0;i<n;i++){
y[i]=y[i]+*betao;
double xb=0;
for (j=0;j<p;j++){
xb+=X[i][j]*beta[j];
}
newy[i]=y[i]-xb;
meany+=newy[i];
}
meany=meany/n;
//if (t>18200){
//printf("Meany=%lf ",meany);
//}
//printf("Meany==%lf",meany);

*betao=sigma2beta*n*meany/s2+sqrt(sigma2beta)*gsl_ran_ugaussian(r);
for (i=0;i<n;i++){
y[i]=y[i]-*betao;
}
}

