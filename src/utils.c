#include <string.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
 #include <gsl/gsl_statistics.h>
 #include <gsl/gsl_math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#include "header.h"
double Truncate(double mu,double sd, double lower, const gsl_rng * r){
double lowern=(lower-mu)/sd;
double alphaopt=(lowern+sqrt(pow(lowern,2)+4))/2;
double z=lowern+gsl_ran_exponential (r, 1/alphaopt);
double qz=exp(-pow(z-alphaopt,2)/2);
double u=gsl_ran_flat (r, 0,1);
while (u>qz){
z=lowern+gsl_ran_exponential (r, 1/alphaopt);
qz=exp(-pow(z-alphaopt,2)/2);
u=gsl_ran_flat (r, 0,1);
}
return z*sd+mu;
}

double InvDigamma(double y){
int niter=5;
double gam=-gsl_sf_psi(1);
double x=exp(y)+0.5;
if (y<-2.22)
x=-1.0/(y+gam);
int i=0;
//printf("XX=%lf\n",x);
//printf("YY=%lf\n",y);
for (i=0;i<niter;i++){
x=x-(gsl_sf_psi(x)-y)/gsl_sf_psi_1(x);
}
 return x;
}


void OutcomeInpute(double ** outc,int n,int NbrOutc){
int i,l;
double ZeroOu[n];
for (i=0;i<n;i++){
ZeroOu[i]=0;
for (l=0;l<NbrOutc;l++){
if (outc[i][l]==0) {
ZeroOu[i]+=1;
}
}
}
for (l=0;l<NbrOutc;l++){
for (i=0;i<n;i++){
if (ZeroOu[i]>0)
outc[i][l]=(outc[i][l]*(n-1)+1.0/NbrOutc)/n;
}
}
}
double max(int n, double * x){
double xmax=x[0];
int i;
for (i=0;i<n;i++){
if (x[i]>=xmax)
xmax=x[i];
}
return xmax;
}


void NJG(int p, int K, _Bool ** path, _Bool * gampath, int* njg){
int k,j;
for (j=0;j<p;j++){
njg[j]=0;
for (k=0;k<K;k++){
if ((gampath[k]==1)&&(path[k][j]==1)){
njg[j]+=1;
}
}
}
}
void findc(int n,_Bool R[n],int a,int * IDX, int *nx){
int ii_data[n];
int idx = 0;
int  ii = 0;
_Bool  exitg2 = 0;
_Bool guard2=0;
  while ((exitg2 == 0) && (ii < n)) {
    guard2 = 0;
    if (R[ii]!= a) {
      idx++;
      ii_data[idx - 1] = ii;
      if (idx >= n) {
        exitg2 = 1;
      } else {
        guard2 = 1;
      }
    } else {
      guard2 = 1;
    }

    if (guard2 == 1) {
      ii++;
    }
  }

int loop_ub=idx;
 for (idx = 0; idx < loop_ub; idx++) {
    IDX[idx] = ii_data[idx];
  }
*nx=loop_ub;

}

void sort(int n,double *x,int *idx)
{
int i,j;
double a;
int id;
for (i = 0; i < n; i++)
idx[i]=i;
for (i = 0; i < n; ++i)
    {
        for (j = i + 1; j < n; ++j)
        {
            if (x[i] <= x[j])
            {
                a =  x[i];
             id=idx[i];
                idx[i]=idx[j];
                x[i] = x[j];
                idx[j]=id;
                x[j] = a;
            }
        }
    }

}


double auc(int n, double * esti,_Bool class[n]){
double fpr[n+2],tpr[n+2];
double auc1=0;
int P=0;//P=positive instances
int i,j;
double esti1[n];
for (i=0;i<n;i++){
esti1[i]=esti[i];
if (class[i]==1)
P+=1;
}
int idx[n];
sort(n,esti1,idx);
fpr[n+1]=1;tpr[n+1]=1;
fpr[0]=0;tpr[0]=0;
for (i=n;i>=1;--i){
double af=0;double at=0;
for (j=0;j<n;j++){
if (esti[j]>esti1[i-1]){
if (class[j]==0){
af+=1;
}
else {
at+=1;
} } }
tpr[i]=at/P;
fpr[i]=af/(n-P);
auc1+=(fpr[i+1]-fpr[i])*(tpr[i+1]+tpr[i]);
}
auc1+=(fpr[1]-fpr[0])*(tpr[1]+tpr[0]);
auc1=0.5*(auc1);
return auc1;
}





void VectorCentered(int nR,double * x){
double mean1 =mean(nR,x);
//double var1 =var(nR,x);
int i;
for (i=0;i<nR;i++){
x[i]=(x[i]-mean1);
}
}

void NormalizeVector(int nR,double * x){
double mean1 =mean(nR,x);
double var1=1;
//double var1 =var(nR,x);
int i;
for (i=0;i<nR;i++){
x[i]=(x[i]-mean1)/sqrt(var1);
}
}




double inverseGaussian(gsl_rng * r, double mu, double lambda) {
double v=gsl_ran_gaussian (r, 1);  // sample from a normal distribution with a mean of 0 and 1 standard deviation
        double y = v*v;
        double x = mu + (mu*mu*y)/(2*lambda) - (mu/(2*lambda)) * sqrt(4*mu*lambda*y + mu*mu*y*y);
double test=gsl_ran_flat (r, 0, 1);    // sample from a uniform distribution between 0 and 1
        if (test <= (mu)/(mu + x))
               return x;
        else
               return (mu*mu)/x;
}



double mean(int n,double * x){
int i;
double me=0;
for (i=0;i<n;i++)
me+=x[i];
return me/n;
}

double var(int n,double * x){
int i;
double me=mean(n,x);
double va=0;
for (i=0;i<n;i++)
va+=(x[i]-me)*(x[i]-me);
return va/(n-1);
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;
    
	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;
    
	for(i=nrl;i<=nrh;i++) 
   {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

_Bool **bmatrix(int nrl, int nrh, int ncl, int nch)
{
        int i;
        _Bool **m;

        m=(_Bool **) malloc( (nrh-nrl+1)*sizeof(_Bool*));
        if (!m) nrerror("allocation failure 1 in dmatrix()");
        m -= nrl;

        for(i=nrl;i<=nrh;i++)
   {
                m[i]=(_Bool *) malloc((nch-ncl+1)*sizeof(_Bool));
                if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
                m[i] -= ncl;
        }
        return m;
}
void free_bmatrix(_Bool **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((_Bool*) (m[i]+ncl));

	free((_Bool*) (m+nrl));
}


int **imatrix(int nrl, int nrh, int ncl, int nch)
{
        int i;
        int **m;

        m=(int **) malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
        if (!m) nrerror("allocation failure 1 in dmatrix()");
        m -= nrl;

        for(i=nrl;i<=nrh;i++)
   {
                m[i]=(int *) malloc((unsigned) (nch-ncl+1)*sizeof(int));
                if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
                m[i] -= ncl;
        }
        return m;
}
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
        int i;

        for(i=nrh;i>=nrl;i--) free((int*) (m[i]+ncl));

        free((int*) (m+nrl));
}

int *ivector(int nl, int nh)
{
        int *v;

        v=(int *)malloc( (nh-nl+1)*sizeof(int));
        if (!v) nrerror("allocation failure in dvector()");
        return v-nl;
}
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));

	free((char*) (m+nrl));
}

void nrerror(char error_text[])
{
	printf("Utils run-time error...\n");
	printf("%s\n",error_text);
	printf("...now exiting to system...\n");
	exit(1);
}

 void save2d(char *filename,int n,int p,double ** data)
         {
         FILE *qfile;
         qfile=fopen(filename,"wb");
          if(qfile==NULL){
                printf("\nUnable to open file.");
                      exit(1);
                        } else {
        int i, j;
        for (i = 0; i < n; i++) {
              for (j = 0; j < p; j++) {
        fprintf(qfile, "%lf ", data[i][j]);
        }
        fprintf(qfile,"\n");
                }
                fclose(qfile);
                 }
}

void save2db(char *filename,int n,int p,_Bool ** data)
         {
         FILE *qfile;
         qfile=fopen(filename,"wb");
          if(qfile==NULL){
                printf("\nUnable to open file.");
                      exit(1);
                        } else {
        int i, j;
        for (i = 0; i < n; i++) {
              for (j = 0; j < p; j++) {
        fprintf(qfile, "%d ", data[i][j]);
        }
        fprintf(qfile,"\n");
                }
                fclose(qfile);
                 }
}
void readDoubleArray(char *filename, int nRows, int nCols, double ** data )
{

   FILE *fp = fopen (filename, "r");
   if (fp==NULL)
   {
      printf("We can't open the file (%s).\n", filename);
      exit(1);
   }
   else
   { int iR,iC;
      for (iR = 0; iR < nRows; ++iR )
      {  
         for (iC = 0; iC < nCols; ++iC )
         {  
            fscanf(fp, "%lf" , &data[iR][iC] );
                     }
                           }
                              fclose(fp);
                             }
}

void readBoolArray(char *filename, int nRows, int nCols, _Bool ** data )
{

   FILE *fp = fopen (filename, "r");
   if (fp==NULL)
   {
      printf("We can't open the file (%s).\n", filename);
      exit(1);
   }
   else
   { int iR,iC;
      for (iR = 0; iR < nRows; ++iR )
      {
         for (iC = 0; iC < nCols; ++iC )
         {
int x;
          fscanf(fp, "%d" , &x );
_Bool bb=(x!=0);
data[iR][iC]=bb;
                     }
                      }

                              fclose(fp);
}
}

void save1d(char *filename,int n,double *data)
         {
         FILE *qfile;
         qfile=fopen(filename,"wb");
          if(qfile==NULL){
                printf("\nUnable to open file.");
                      exit(1);
                        } else {
        int i;
        for (i = 0; i < n; i++) {
        fprintf(qfile, "%f ", data[i]);
                }
                fclose(qfile);
                 }
}
