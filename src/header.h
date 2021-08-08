#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define PI 3.14159265358979323846
typedef struct
{
  gsl_matrix * Mat;
  gsl_vector * Vec;
} GSLMatVect;


void mainMCMCFunction(int* q1, int* K1,int *n1,int* p1,double * outc1,double * yinit1,double * X1,
                      int * path1,int* burninsample1,
                      int* nbrsample1,double* a, double* b, double* alpha01, 
                      double* beta01,double * al1, double* bll,int * seed1, 
                      double * gamMean1,double * BetaSample1, double *PostPredSample1,double * logpost);
void logposterior(int n,int p,int K, int NbrOut,int ** njg,double ** alphaS, double ** outc,double * s2, double ** beta,
             double *beta0l,double ** X, double *** Tau,double *lambda2,_Bool** gampath,double * q,
             double ab, double bb,double al, double *bl,double alpha0, double beta0,double sigma20,double *logpost);
double geneBeta(gsl_rng * r,int K,_Bool * gampath,double ab, double bb);
double gigrnd(gsl_rng * rr,double P, double a, double b);
double InvDigamma(double y);
double Truncate(double mu,double sd, double lower, const gsl_rng * r);
void VectorCentered(int nR,double * x);
void OutcomeInpute(double ** outc,int n,int NbrOutc);
double max(int n, double * x);
void alphaStar(gsl_rng * r,int n,int p, int NbrOut,double ** alphaS, double ** outc,double * s2, double ** beta,double *beta0l,double **X,int *acceptRate,int sample);
void alphaStar2(gsl_rng * r,int n,int p, int NbrOut,double ** alphaS, double ** outc,double * s2, double ** beta,double *beta0l,double ** X,int *acceptRate,int sample);
void alphaStar1(gsl_rng * r,int n,int p, int NbrOut,double ** alphaS, double ** outc,double * s2, double ** beta,double *beta0l,double ** X,int *acceptRate,int sample);
void save2db(char *filename,int n,int p,_Bool ** data);
void save2d(char *filename,int n,int p,double ** data);
double auc(int n, double * esti,_Bool class[n]);
void sort(int n,double *x,int *idx);
double geneLambda(gsl_rng * r, int p, int * njg, double ** Tau, double al, double bl);
void NormalizeVector(int nR,double * x);
double sigma2(gsl_rng * r,int n,int p,int *njg,double ** tau, double *beta,double **X,double *y, double alpha0, double beta0,int sample);
void geneTau(gsl_rng * r, int p, int K, int * njg, double ** Tau, double * beta, double s2,double lambda2);
double inverseGaussian(gsl_rng * r, double mu, double lambda);
double logmarggam(int n, int p, int *np,double * Yk,double ** X,int * IDX,_Bool * path,int * njg, double ** Tau,double alpha0, double beta0);
void gampat(int n, int p, int K, double * y,double ** X, gsl_rng * r,_Bool ** path,int * njg,_Bool * gampath, double ** Tau, double * beta, double alpha0, double beta0,double q,double s2);
//void pathNoOverl(int p, int K,_Bool ** path,_Bool ** path1);
void gampat1(int n, int p, int K, double * y,double ** X, gsl_rng * r,_Bool ** path,int * njg,_Bool * gampath, double ** Tau, double * beta, double alpha0, double beta0,double q,double s2);


GSLMatVect logmarggam1(int n, int p, int *np,double * Yk,double ** X,int * IDX,_Bool * path,int * njg, double ** Tau,double alpha0, double beta0, double *logmarg);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void nrerror(char error_text[]);
void free_bmatrix(_Bool **m, int nrl, int nrh, int ncl, int nch);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
double var(int n,double * x);
double mean(int n,double * x);
_Bool **bmatrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
void nrerror(char error_text[]);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
int *ivector(int nl, int nh);
void findc(int n,_Bool R[n],int a,int * IDX, int *nx);
void NJG(int p, int K, _Bool ** path, _Bool * gampath, int* njg);
double * Ainbeta(int k,int n,int * njg,double ** tau,double **Ps,int * IDX);
double **PGG(int n, int p,_Bool * path,double **X,int *njg, int *IDX, int *np);
void Effects(int n, int p, int K, double * y,double ** X, gsl_rng * r,_Bool ** path,int * njg, double ** Tau, double * beta, double s2);
void generatedata(gsl_rng * r,int n, int p,int K,_Bool ** path,_Bool * gampath, double *y, double ** X, int propOverl);
void mainMCMC(_Bool update_lambda,double ** y,double ** alphaInit,int burninsample, int nbrsample,int n, int p,int NbrOutc, int K,_Bool ** path,
              _Bool** gampath,double ** gamMean, double **outc, double ** X,double al, double *bl,double *lambda2S,double ab, double bb, 
              double alpha0, double beta0,double *** BetaSample,double ***PostPredSample,double * logpost);
void GenerateDirichlet(double **alphaSta,double ** outcom, int q, gsl_rng * r,int n, int p,int K,_Bool ** path,_Bool ** gampath, double ** X, int propOverl,double rr);
void readBoolArray(char *filename, int nRows, int nCols, _Bool ** data );
void readDoubleArray(char *filename, int nRows, int nCols, double ** data );
void save1d(char *filename,int n,double *data);
void Beta0(gsl_rng *r,int n,int p,double *betao, double s2,double sigma20,double * beta,double *y,double ** X,int t);
gsl_vector *dirichlet_mle(gsl_matrix *D);
double _ipsi(double y);
