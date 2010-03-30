#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "strsubs.h" 

void vsp(double *a, double *b, double c, int n);
void vst(double *a, double *b, double c, int n);
void vvt(double *a, double *b, double *c, int n);
void vvp(double *a, double *b, double *c, int n);
void vvm(double *a, double *b, double *c, int n);
void vvd(double *a, double *b, double *c, int n);
void vsqrt(double *a, double *b,  int n) ;
void vinvert(double *a, double *b,  int n) ;
void vabs(double *a, double *b,  int n)  ;
void vlog(double *a, double *b,  int n)  ;
void vlog2(double *a, double *b,  int n)  ;
void vexp(double *a, double *b,  int n)  ;
void vclear(double *a,  double c, int n) ;
void vzero(double *a, int n) ;
void ivvp(int *a, int *b, int *c, int n);
void ivvm(int *a, int *b, int *c, int n);
void ivsp(int *a, int *b, int c, int n);
void ivst(int *a, int *b, int c, int n);
void ivclear(int *a,  int c, int n) ;
void ivzero(int *a, int n) ;
void cclear(unsigned char *a,  unsigned char c, int n) ;
double clip(double x, double lo, double hi)  ;
void vclip(double *a, double *b,double loval, double hival,int n)   ;
void vmaxmin(double *a, int n, double *max, double *min)  ;
void vlmaxmin(double *a, int n, int *max, int *min)  ;
void ivmaxmin(int *a, int n, int *max, int *min)  ;
void ivlmaxmin(int *a, int n, int *max, int *min)  ;
void getdiag(double *a, double *b, int n)  ;
void flipiarr(int *a, int *b, int n) ;
void fliparr(double *a, double *b, int n)  ;

void copyarr(double *a,double *b,int n) ;
void copyiarr(int *a,int *b,int n) ;
void copyiparr(int **a,int **b,int n) ;

void dpermute(double *a, int *ind, int len)  ;
void ipermute(int *a, int *ind, int len)  ;
void dppermute(double **a, int *ind, int len)  ;
void ippermute(int **a, int *ind, int len)  ;

double  asum(double *a, int n) ;
double  asum2(double *a, int n) ;
int     intsum(int *a, int n) ;
int     idot(int *a, int *b, int n)  ;
double  aprod(double *a, int n) ;
double  vdot(double *a, double *b, int n)  ;
double  corr(double *a, double *b, int n)  ;
double  variance(double *a, int n)  ;
double trace(double *a, int n) ;
int  nnint(double a) ;
void countcat(int *tags, int n,int *ncat,int nclass)  ;
void rowsum(double *a, double *rr, int n) ;
void colsum(double *a, double *cc, int n) ;
void rrsum(double *a, double *cc, int m, int n)  ;
void ccsum(double *a, double *cc, int m, int n)  ;
void printmat(double *a, int m, int n) ;
void printmatw(double *a, int m, int n, int w) ;
void printmatl(double *a, int m, int n) ;
void printmatwl(double *a, int m, int n, int w) ;
void floatit(double *a, int *b, int n) ;
void fixit(int  *a, double *b, int n) ;
void printimat(int *a, int m, int n) ;
int  findfirst(int *a, int n, int val) ;
int  findlast(int *a, int n, int val) ;
int  binsearch(int *a, int n, int val) ;
void idperm(int *a, int n)  ;
double NPlog2(double y) ;
double log2fac(int  n) ;
double logfac(int n)  ;
double logbino(int n, int k)  ;
double logmultinom(int *cc, int n) ;
double addlog(double a, double b) ;
double vldot(double *x, double *y, int n) ;
double pow10 (double x) ;
double vpow10 (double *a, double *b, int n) ;
double vlog10 (double *a, double *b, int n) ;
/* matrix transpose */
void transpose(double *aout, double *ain, int m, int n)  ;
void addouter(double *out, double *a, int n) ;
void subouter(double *out, double *a, int n) ;

/* storage allocation */
int **initarray_2Dint(int numrows, int numcolumns, int initval);
void free2Dint(int ***xx, int numrows) ;
double **initarray_2Ddouble(int numrows, int numcolumns, double initval);
void clear2D(double ***xx, int numrows, int numcols, double val)   ;
void iclear2D(int ***xx, int numrows, int numcols, int val)  ;
void free2D(double ***xx, int numrows) ;

double bal1 (double *a, int n)  ;
void vcompl(double *a, double *b, int n) ;
void setidmat(double *a, int n) ; 

int stripit(double *a, double *b, int *x, int len) ;
int istripit(int *a, int *b, int *x, int len) ;
int cstripit(char **a, char **b, int *x, int len) ;

void mapit(int *a, int *b, int n, int inval, int outval) ;
int ifall(int n, int k)  ;  // falling factorial = n (n-1) (n-2) ... (n-k+1)


