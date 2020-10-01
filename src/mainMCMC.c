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

//R  CMD SHLIB  mainAlgo.c MCMC.c  utils.c mainfuntion2.c -lgsl -lgslcblas -o PathSel.so

//R CMD SHLIB mainAlgo.c  utils.c mainfuntion2.c -lgsl -lgslcblas -o PathSel.so
void mainMCMCFunction(int* q1, int* K1,int *n1,int* p1,double * outc1,double * yinit1,double * X1,int * path1,int* burninsample1,int* nbrsample1,double* a, double* b, double* alpha01, double* beta01,double * al1, double* bll,int * seed1, double * gamMean1,double * BetaSample1){
//printf("Hello World!\n");


setvbuf(stdout, NULL, _IONBF, 0);
clock_t t = clock();
int burninsample=burninsample1[0];
int nbrsample=nbrsample1[0];
int p=p1[0];
printf("Number of markers is %d\n",p);
int n=n1[0];
printf("Number of samples is %d\n",n);
int K=K1[0];
printf("Number of pathways/groups is %d\n",K);
//int q=4; //nbr of outcomes
int q=q1[0];
double ** X=dmatrix(0,n-1,0,p-1);
double ** outc=dmatrix(0,n-1,0,q-1);
double ** yinit=dmatrix(0,q-1,0,n-1);
int i,l,j,k;
for (i=0;i<n;i++){
for (j=0;j<p;j++){
X[i][j]=X1[i*p+j];
}
for (l=0;l<q;l++){
outc[i][l]=outc1[i*q+l];
yinit[l][i]=yinit1[i*q+l];
}
}

_Bool ** path=bmatrix(0,K-1,0,p-1);
for (k=0;k<K;k++){
for (j=0;j<p;j++){
path[k][j]=(path1[k*p+j]==1);
}
}
double *** BetaSample=malloc(q*sizeof(double **));
for (l=0;l<q;l++){
BetaSample[l]=dmatrix(0,nbrsample-1,0,p-1);
}
_Bool** gampath=bmatrix(0,q-1,0,K-1);
double ** gamMean=dmatrix(0,q-1,0,K-1);

//double ** alphaSta=dmatrix(0,q-1,0,n-1);
double ** y=dmatrix(0,q-1,0,n-1);
long seed=(long) seed1[0];
gsl_rng * r = gsl_rng_alloc (gsl_rng_rand48);
gsl_rng_set (r, seed);




// We impute zero values
OutcomeInpute(outc,n,q);
/*
Initialization of gampath **/
for (l=0;l<q;l++){
for (k=0;k<K;k++){
double  uni=gsl_ran_flat (r, 0, 1);
  if (uni<0.1)  gampath[l][k]=1; else gampath[l][k]=0;
}
}

for (l=0;l<q;l++){
VectorCentered(n,yinit[l]);
}
/* Set hyperparameters*/
double ab=a[0]; double bb1=b[0];
double alpha0=alpha01[0]; double beta0=beta01[0];
//double bll=2;
double bl[q];
for (l=0;l<q;l++){
//bl[l]=1;
bl[l]=bll[0];
}
printf("\n");
double al=al1[0];
//double al=8;
double *lambda2S=malloc(q*sizeof(double));
//printf("ALLLL is %lf",al);
mainMCMC(0,y,yinit,burninsample, nbrsample,n, p,q,K,path,gampath,gamMean, outc,X,al, bl,lambda2S,ab, bb1, alpha0, beta0,BetaSample);
int kk=0;int ti;
int bb=0;
for (l=0;l<q;l++){
for (ti=0;ti<nbrsample;ti++)
for (j=0;j<p;j++)
BetaSample1[bb++]=BetaSample[l][ti][j];
for (k=0;k<K;k++){
gamMean1[kk++]=gamMean[l][k];
}
}

t = clock() - t;
double  time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("\n\nTime taken in seconds is %f\n",time_taken);
    printf("\nTime taken in minutes is %f\n",time_taken/60);
   // printf("\nTime taken in hours is %f\n",time_taken/3600);

free_bmatrix(path,0,K-1,0,p-1);
free_dmatrix(X,0,n-1,0,p-1);
free_dmatrix(outc,0,n-1,0,q-1);
free_dmatrix(gamMean,0,q-1,0,K-1);
free_bmatrix(gampath,0,q-1,0,K-1);
free_dmatrix(y,0,q-1,0,n-1);
free_dmatrix(yinit,0,q-1,0,n-1);
gsl_rng_free (r);
free(lambda2S);
for (l=0;l<q;l++){
free_dmatrix(BetaSample[l],0,nbrsample-1,0,p-1);
}
free(BetaSample);
}

