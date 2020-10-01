#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_randist.h>
 #include <gsl/gsl_statistics.h>
 #include <gsl/gsl_math.h>
#include <time.h>
#include "header.h"

//R CMD SHLIB mainAlgo.c  utils.c mainfuntion2.c -lgsl -lgslcblas -o PathSel.so

void mainMCMC(_Bool update_lambda,double ** y,double ** yinit,int burninsample, int nbrsample,int n, int p,int NbrOutc, 
              int K,_Bool ** path,_Bool** gampath,double ** gamMean, double **outc, double ** X,double al, 
              double *bl,double *lambda2S,double ab, double bb, double alpha0, double beta0, 
              double *** BetaSample){
int l,k,t,j;
long seed1=1;

//_Bool** gampath=bmatrix(0,NbrOutc-1,0,K-1);
double *s2=malloc(NbrOutc*sizeof(double));
double * q=malloc(NbrOutc*sizeof(double));
double sigma20=1e10;
double *lambda2=malloc(NbrOutc*sizeof(double));
//double alpha0=0.001;double beta0=0.001;
//alpha0=4.0;
//beta0=0.01;

gsl_rng * r = gsl_rng_alloc (gsl_rng_rand48);
gsl_rng_set (r, seed1);

double q0=0.2;
//double ab=2;double bb=20;
double *** Tau=malloc(NbrOutc*sizeof(double **));
for (l=0;l<NbrOutc;l++){
Tau[l]=dmatrix(0,p-1,0,2-1);
lambda2S[l]=0;
lambda2[l]=al/bl[l];
s2[l]=1;
q[l]=q0;
}
//double ** gamMean=dmatrix(0,NbrOutc-1,0,K-1);
double **beta =dmatrix(0,NbrOutc-1,0,p-1);
//double *** BetaSample=malloc(NbrOutc*sizeof(double **));
double ** njgMean =dmatrix(0,NbrOutc-1,0,p-1);
int ** njg=imatrix(0,NbrOutc-1,0,p-1);
double *s2mean=malloc(NbrOutc*sizeof(double));
double *beta0l=malloc(NbrOutc*sizeof(double));
double *beta0lMean=malloc(NbrOutc*sizeof(double));
double *RGmean=malloc(p*sizeof(double));

for (l=0;l<NbrOutc;l++){
s2mean[l]=0;
beta0l[l]=0;
beta0lMean[l]=0;
//BetaSample[l]=dmatrix(0,nbrsample-1,0,p-1);
for (j=0;j<p;j++){
RGmean[j]=0;
beta[l][j]=0;
for (k=0;k<2;k++){
Tau[l][j][k]=0.1;
if (k==0)
Tau[l][j][k]=0;
}
}
for (j=0;j<K;j++){
//double  uni=gsl_ran_flat (r, 0, 1);
//  if (uni<0.1)  gampath[l][j]=1; else gampath[l][j]=0;
gamMean[l][j]=0;
}
NJG(p,K,path,gampath[l],njg[l]);
for (j=0;j<p;j++){
if (njg[l][j]!=0)
beta[l][j]=0.1;
njgMean[l][j]=0;
}
}
int i;
for (l=0;l<NbrOutc;l++){
for (i=0;i<n;i++){
y[l][i]=yinit[l][i];
}
}


int N=burninsample+nbrsample;
int acceptRate=0;
for (t=0;t<N;t++){
for (l=0;l<NbrOutc;l++){
q[l]=geneBeta(r,K,gampath[l],ab,bb);
Beta0(r,n,p,&beta0l[l], s2[l],sigma20,beta[l],y[l],X,t);
//printf("Beta000=%.2lf ",beta0l[l]);
geneTau(r,p,K, njg[l], Tau[l], beta[l],s2[l],lambda2[l]);
gampat1(n,p, K, y[l],X,r,path,njg[l], gampath[l], Tau[l],  beta[l], alpha0,  beta0,q[l],s2[l]);
//s2[l]=0.01;
s2[l]=sigma2(r,n,p,njg[l],Tau[l], beta[l],X,y[l], alpha0, beta0,t);
lambda2[l]=geneLambda(r, p, njg[l], Tau[l], al, bl[l]);
}
alphaStar(r,n,p, NbrOutc,y, outc,s2, beta,beta0l,X,&acceptRate,t);


if (t>=burninsample){
for (l=0;l<NbrOutc;l++){
lambda2S[l]+=lambda2[l]/nbrsample;
}
}
if (update_lambda==0){
if (t%(N/5)==1){
printf("\n");
printf("The number of mcmc  iterations is %d\n\n",t);

}
//printf("Tau=\n");
//for (j=0;j<p;j++){
//for (i=1;i<2;i++){
//printf("%.2lf ",Tau[j][i]);
//}
//printf("\n");
//}
if (t>=burninsample){
for (l=0;l<NbrOutc;l++){
for (j=0;j<K;j++){
gamMean[l][j]+=(double) gampath[l][j]/(nbrsample);
}
s2mean[l]+=s2[l]/nbrsample;
beta0lMean[l]+=beta0l[l]/nbrsample;
//printf("NJG=\n");
}
for (j=0;j<p;j++){
double njx=0;
for (l=0;l<NbrOutc;l++){
njgMean[l][j]+=(double) njg[l][j]/(nbrsample);
BetaSample[l][t-burninsample][j]=beta[l][j];
if (njg[l][j]>0)
njx+=1;
}
if (njx>0)
RGmean[j]+=1.0/nbrsample;
}


}
}
}// End of MCMC 

if (update_lambda==0){

printf("AcceptRateAlphaStar=%lf\n",1.0*acceptRate/((burninsample+nbrsample)*n*NbrOutc));

printf("\n");
printf("Sigma2Mean=");
for (l=0;l<NbrOutc;l++){
printf("%.2lf ",s2mean[l]);
}
printf("\n");
printf("Beta0Mean=");
for (l=0;l<NbrOutc;l++){
printf("%.2lf ",beta0lMean[l]);
}
printf("\n");
}
//char *RepSav=malloc(100*sizeof(char));

//free(RepSav);

//free_dmatrix(gamMean,0,NbrOutc-1,0,K-1);
free_imatrix(njg,0,NbrOutc-1,0,p-1);
free_dmatrix(njgMean,0,NbrOutc-1,0,p-1);
gsl_rng_free (r);
//free_bmatrix(gampath,0,NbrOutc-1,0,K-1);
for (l=0;l<NbrOutc;l++){
free_dmatrix(Tau[l],0,p-1,0,2-1);
//free_dmatrix(BetaSample[l],0,nbrsample-1,0,p-1);
}
free_dmatrix(beta,0,NbrOutc-1,0,p-1);
//free(BetaSample);
free(Tau);
//free_bmatrix(pathNoOv,0,K-1,0,p-1);
free(s2mean);
free(s2);free(lambda2);
free(beta0l);
free(beta0lMean);
free(RGmean);
free(q);
}
