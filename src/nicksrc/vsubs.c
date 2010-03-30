#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "strsubs.h" 
#include "vsubs.h" 

/** 
 tiny routines BLAS? 
 a small library to do simple arithmetic
 on 1D vectors with no skips 
*/
void 
vsp(double *a, double *b, double c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] + c ;
}
void 
vst(double *a, double *b, double c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] * c ;
}
void 
vvt(double *a, double *b, double *c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] * c[i] ;
}
void 
vvp(double *a, double *b, double *c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] + c[i] ;
}
void 
vvm(double *a, double *b, double *c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] - c[i] ;
}
void 
vvd(double *a, double *b, double *c, int n) 
{
   int i ;
   for (i=0; i<n; i++)  {
    if (c[i] == 0.0) 
      fatalx("(vvd): zero value in denominator\n") ;
    a[i] = b[i] / c[i] ;
   }
}
void 
vsqrt(double *a, double *b,  int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    if (b[i]<0.0) 
     fatalx("(vsqrt): negative value %g\n",b[i]) ;
    if (b[i] == 0.0) {
     a[i] = 0.0 ;
     continue ;
    }
    a[i] = sqrt(b[i]) ;
   }
}
void 
vinvert(double *a, double *b,  int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    if (b[i] == 0.0) 
     fatalx("(vinvert): zero value\n") ;
    a[i] = 1.0 / b[i] ;
   }
}
void 
vabs(double *a, double *b,  int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    a[i] = fabs(b[i]) ;
   }
}
void 
vlog(double *a, double *b,  int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    if (b[i]<=0.0) 
     fatalx("(vlog): negative or zero value %g\n",b[i]) ;
    a[i] = log(b[i]) ;
   }
}
void 
vlog2(double *a, double *b,  int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    if (b[i]<=0.0) 
     fatalx("(vlog2): negative or zero value %g\n",b[i]) ;
    a[i] = NPlog2(b[i]) ;
   }
}

void 
vexp(double *a, double *b,  int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    a[i] = exp(b[i]) ;
   }
}
void 
vclear(double *a,  double c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] =  c ;
}
void 
vzero(double *a,  int n) 
{
  vclear(a, 0.0, n) ;
}
void 
cclear(unsigned char *a,  unsigned char c, int n) 
/** 
 be careful nothing done about NULL at end
*/
{
   int i ;
   for (i=0; i<n; i++)  {
    a[i] =  c ;
   }
}
void 
ivvp(int *a, int *b, int *c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] + c[i] ;
}
void 
ivvm(int *a, int *b, int *c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] - c[i] ;
}
void 
ivsp(int *a, int *b, int c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] + c ;
}
void 
ivst(int *a, int *b, int c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] * c ;
}
void 
ivclear(int *a,  int c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] =  c ;
}
void
ivzero(int *a,  int n) 
{
  ivclear(a, 0, n) ;
}
double clip(double x, double lo, double hi) 
/* clip off values to range [lo,hi] */
{
  if (x<lo) return lo ;
  if (x>hi) return hi ;
  return x ;
}

void 
vclip(double *a, double *b,double loval, double hival,int n)  
{
/* clip off values to range [loval,hival] */
   int i ;
   double t ;

   for (i=0; i<n; i++) {
    t = MAX(b[i],loval) ;
    a[i] = MIN(t,hival) ;
   }
}
void vmaxmin(double *a, int n, double *max, double *min)  
{

    int i ;
    double tmax, tmin ;

    tmax = tmin = a[0] ;
    for (i=1; i<n; i++) {
      tmax = MAX(tmax, a[i]) ;
      tmin = MIN(tmin, a[i]) ;
    }
    if (max != NULL) *max = tmax ;
    if (min != NULL) *min = tmin ;
}
void vlmaxmin(double *a, int n, int *pmax, int *pmin)  
/** 
 return location 
*/
{

    int i ;
    double tmax, tmin ;
    double lmax, lmin ;

    tmax = tmin = a[0] ;
    lmax = lmin = 0 ;
    for (i=1; i<n; i++) {
      if (a[i]>tmax) {
       tmax = a[i] ;
       lmax=i ;
      }
      if (a[i]<tmin) {
       tmin = a[i] ;
       lmin=i ;
      }
    }
    if (pmax != NULL) *pmax = lmax ;
    if (pmin != NULL) *pmin = lmin ;
}
void ivmaxmin(int *a, int n, int *max, int *min)  
{

    int i ;
    int tmax, tmin ;

    tmax = tmin = a[0] ;
    for (i=1; i<n; i++) {
      tmax = MAX(tmax, a[i]) ;
      tmin = MIN(tmin, a[i]) ;
    }
    if (max != NULL) *max = tmax ;
    if (min != NULL) *min = tmin ;
}
void ivlmaxmin(int *a, int n, int *pmax, int *pmin)  
/** 
 return location 
*/
{

    int i ;
    int tmax, tmin ;
    int lmax, lmin ;

    tmax = tmin = a[0] ;
    lmax = lmin = 0 ;
    for (i=1; i<n; i++) {
      if (a[i]>tmax) {
       tmax = a[i] ;
       lmax=i ;
      }
      if (a[i]<tmin) {
       tmin = a[i] ;
       lmin=i ;
      }
    }
    if (pmax != NULL) *pmax = lmax ;
    if (pmin != NULL) *pmin = lmin ;
}
double 
vdot(double *a, double *b, int n) 
{
   int i; 
   double ans=0.0 ;
   for (i=0; i<n; i++) 
    ans += a[i]*b[i] ;

   return ans ;

}
double corr(double *a, double *b, int n) 
{
  double v12, v11, v22, y1, y2, y ;
  double *aa, *bb ;
  ZALLOC(aa, n, double) ;
  ZALLOC(bb, n, double) ;
  y1 = asum(a,n)/ (double) n ;
  y2 = asum(b,n)/ (double) n ;

  vsp(aa, a, -y1, n) ;
  vsp(bb, b, -y2, n) ;

  v12 = vdot(aa, bb, n) ;
  v11 = asum2(aa, n) ;
  v22 = asum2(bb, n) ;

  y = v11*v22 ;
  if (y==0.0) fatalx("(corr) constant vector\n") ;


  free(aa) ;  free(bb) ;
  return (v12/sqrt(y)) ;

}

double variance(double *a, int n) 
{

  double *aa ;
  double y1, y2 ;

  ZALLOC(aa, n, double) ;
  y1 = asum(a,n)/ (double) n ;
  vsp(aa, a, -y1, n) ;

  y2 = asum(aa,n)/ (double) n ;

  free(aa) ;
  return y2 ;

}

void 
getdiag(double *a, double *b, int n)  
/* extract diagonal */
{
   int i, k ;

   for (i=0; i<n; i++) {
    k = i*n+i ;
    a[i] = b[k] ;
   }
}

void copyarr(double *a,double *b,int n) 
{
  int i ;
  for (i=0; i<n; i++) {
   b[i] = a[i] ;
  }
}

void copyiarr(int *a,int *b,int n) 
{
  int i ;
  for (i=0; i<n; i++) {
   b[i] = a[i] ;
  }
}
void copyiparr(int **a,int **b,int n) 
{
  int i ;
  for (i=0; i<n; i++) {
   b[i] = a[i] ;
  }
}
void dpermute(double *a, int *ind, int len) 
{

  int i , k ;
  double *rrr ;

  ZALLOC(rrr, len, double) ;

  for (i=0; i<len; i++) {
   rrr[i] = a[i] ;
  }

  for (i=0; i<len; i++) {
   k = ind[i] ;
   a[i] = rrr[k] ;
  }

  free (rrr) ;
}
void ipermute(int *a, int *ind, int len) 
{

  int i , k ;
  int *rrr ;

  ZALLOC(rrr, len, int) ;

  copyiarr(a, rrr, len) ;

  for (i=0; i<len; i++) {
   k = ind[i] ;
   a[i] = rrr[k] ;
  }

  free (rrr) ;
}
void dppermute(double **a, int *ind, int len) 
{

  int i , k ;
  double **rrr ;

  ZALLOC(rrr, len, double *) ;

  for (i=0; i<len; i++) {
   rrr[i] = a[i] ;
  }

  for (i=0; i<len; i++) {
   k = ind[i] ;
   a[i] = rrr[k] ;
  }

  free (rrr) ;
}
void ippermute(int **a, int *ind, int len) 
{

  int i , k ;
  int **rrr ;

  ZALLOC(rrr, len, int *) ;

  for (i=0; i<len; i++) {
   rrr[i] = a[i] ;
  }

  for (i=0; i<len; i++) {
   k = ind[i] ;
   a[i] = rrr[k] ;
  }

  free (rrr) ;
}

double  asum(double *a, int n) 
{
   int i; 
   double ans=0.0 ;
   for (i=0; i<n; i++) 
    ans += a[i] ;

   return ans ;
}
int  intsum(int *a, int n) 
{
   int i; 
   int ans=0 ;
   for (i=0; i<n; i++) 
    ans += a[i] ;

   return ans ;
}
int idot(int *a, int *b, int n) 
{
   int i; 
   int ans=0.0 ;
   for (i=0; i<n; i++) 
    ans += a[i]*b[i] ;

   return ans ;

}


double  aprod(double *a, int n) 
/* overflow not checked */
{
   int i; 
   double ans=1.0 ;
   for (i=0; i<n; i++) 
    ans *= a[i] ;

   return ans ;
}
double  asum2(double *a, int n) 
{
   int i; 
   double ans=0.0 ;
   for (i=0; i<n; i++) 
    ans += a[i]*a[i] ;

   return ans ;
}

double trace(double *a, int n) 
{
   double *diags, t ;
   ZALLOC(diags,n,double) ;
   getdiag(diags,a,n) ; /* extract diagonal */
   t = asum(diags,n) ;
   free(diags) ;
   return t ;
}

int nnint(double x)
{
   long int lrint(double x) ;  
// double round(double x) ;
   return (int) lrint(x) ;
}
void 
countcat(int *tags, int n,int *ncat,int nclass)  
/* simple frequency count of integer array */

{
   int i, k; 
   ivzero(ncat, nclass) ;
   for (i=0 ; i<n ; i++)  {
    k = tags[i] ;
    if ( (k<0) || (k >= nclass)) 
     fatalx("(countcat) bounds error\n") ;
    ++ncat[k] ;
   }
}
void rowsum(double *a, double *rr, int n) 
{
   int i,j ;
   vclear(rr,0.0,n) ;
   for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
     rr[j] += a[i+j*n] ;
    }
   }
}
void colsum(double *a, double *cc, int n) 
{
   int i,j ;
   vclear(cc,0.0,n) ;
   for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
     cc[i] += a[i+j*n] ;
    }
   }
}
void rrsum(double *a, double *rr, int m, int n) 
{
   int i,j ;
   vclear(rr,0.0,n) ;
   for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
     rr[j] += a[i+j*m] ;
    }
   }
}
void ccsum(double *a, double *cc, int m, int n) 
{
   int i,j ;
   vclear(cc,0.0,m) ;
   for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
     cc[i] += a[i+j*m] ;
    }
   }
}
void printmat(double *a, int m, int n) 
/** 
 print a matrix n wide m rows  
*/
{
      printmatw(a, m, n, 5) ;
}
void printmatw(double *a, int m, int n, int w) 
/** 
 print a matrix n wide m rows  w to a row
*/
{
      int i,j, jmod ;
      for (i=0; i<m; i++)  {
        for (j=0; j<n; j++)  {
         printf("%9.3f ", a[i*n+j]) ;
         jmod = (j+1) % w ;
         if ((jmod == 0)  && (j<(n-1))) {
          printf("  ...\n") ;
         }
        }
        printf("\n") ;
      }
}
void printmatl(double *a, int m, int n) 
/** 
 print a matrix n wide m rows  
*/
{
      printmatwl(a, m, n, 5) ;
}
void printmatwl(double *a, int m, int n, int w) 
/** 
 print a matrix n wide m rows  w to a row
 15.9f format
*/
{
      int i,j, jmod ;
      for (i=0; i<m; i++)  {
        for (j=0; j<n; j++)  {
         printf("%15.9f ", a[i*n+j]) ;
         jmod = (j+1) % w ;
         if ((jmod == 0)  && (j<(n-1))) {
          printf("  ...\n") ;
         }
        }
        printf("\n") ;
      }
}
void floatit(double *a, int *b, int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    a[i] = (double) b[i] ;
   } 
}
void printimat(int  *a, int m, int n) 
/** 
 print a matrix n wide m rows  
*/
{
      int i,j, jmod ;
      for (i=0; i<m; i++)  {
        for (j=0; j<n; j++)  {
         printf("%5d ", a[i*n+j]) ;
         jmod = (j+1) % 10 ;
         if ((jmod == 0)  && (j<(n-1))) {
          printf("  ...\n") ;
         }
        }
        printf("\n") ;
      }
}

void fixit(int  *a, double *b, int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    a[i] = nnint(b[i]) ;
   } 
}
int findfirst(int *a, int n, int val) 
{
   int i ;
   for (i=0 ; i<n; i++)  {
    if (a[i] == val) return i ;
   }
   return -1 ;  
}
int findlast(int *a, int n, int val) 
{
   int i ;
   for (i=n-1 ; i>=0; i--)  {
    if (a[i] == val) return i ;
   }
   return -1 ;  
}
int binsearch (int *a, int n, int val)
// binary search.  a sorted in ascending order
{
#define TINYS 12
  int  x,  m, h, v ;

  if (n<=TINYS) return findfirst(a, n, val) ;
  if (val<a[0]) return -1 ;
  if (val>a[n-1]) return -1 ;
  h = n/2 ;
  v = a[h] ;
  if (val<v) return binsearch (a, h, val) ;
  if (val==v) return h ;
  m = (n-1) - (h+1) + 1 ;
  x = binsearch (a+h+1, m, val) ;
  if (x<0) return -1 ;
  return x + h + 1 ;
}

void idperm(int *a, int n) 
{
     int i ;
     for (i=0; i<n; i++) 
      a[i] = i ;
}
double NPlog2(double y) 
{
  if (y<=0.0) fatalx("(NPlog2) negative argument\n") ;
  return (log(y)/log(2.0)) ;
}

double logfac(int n) 
/** 
 log (factorial n))
*/
{
  double y, x ;
  x = (double) (n+1) ;
  y = lgamma(x) ;
  return (y) ;
}

double logbino(int n, int k) 
/* log n choose k */
{
  double top, bot ;
   
  top = logfac(n) ;
  bot = logfac(n-k) + logfac(k) ;

  return top-bot ;
}

double log2fac(int n) 
/** 
 log base2 (factorial n))
*/
{
  double y, x ;
  x = (double) (n+1) ;
  y = lgamma(x) ;
  return (y/log(2.0)) ;
}

double addlog(double a, double b) 
{
 /* given a = log(A)
          b = log(B)
    returns log(A+B) 
    with precautions for overflow etc
*/
    double x, y, z ;

    x = MIN(a,b) ;
    y = MAX(a,b) ;

/** 
 answer is log(1 + A/B) + log (B)  
*/
    z = x-y ;  
    if (z < -50.0)  return y ;
    z = 1.0+exp(z) ;
    z = log(z) + y ;
    return (z) ;

}



double vldot(double *x, double *y, int n) 
/** 
 x. log(y) 
*/
{
  double *z, ans ;
  double tiny = 1.0e-19 ;
  int i ;

  ZALLOC(z, n, double) ;
  vsp(z, y, 1.0e-20, n) ;
  vlog(z, z, n) ;

  ans = 0.0 ;
  for (i=0; i< n ; i++) {  
   if (x[i]>tiny) ans += x[i]*z[i] ;
  }
  free (z) ;
  return ans ;
}

double pow10 (double x) 
{
     return exp(x*log(10.0)) ;
}


double vpow10 (double *a, double *b, int n) 
{
     int i ;
     for (i=0; i< n; i++)  
      a[i] = exp(b[i] * log(10.0)) ;
}

double vlog10 (double *a, double *b, int n) 
{
     int i ;
     for (i=0; i< n; i++)  
      a[i] = log10(b[i]) ;
}

void transpose(double *aout, double *ain, int m, int n) 
/** 
 aout and ain must be identical or not overlap 
 does matrix transpose 

 input  m vectors of length n  (m x n) 
 output n vectors of length m 
*/
{
    double *ttt  ;
    int i, j, k1, k2 ;
    if (aout == ain) {  
     ZALLOC(ttt, m*n, double) ;
    }
    else ttt = aout ;

    for (i=0; i<m; i++) 
     for (j=0; j<n; j++)  {
      k1 = i*n+j ;
      k2 = j*m+i ;
      ttt[k2] = ain[k1]  ;
    }
    if (aout == ain) {  
     copyarr(ttt, aout, m*n) ;
     free(ttt) ;
    }
}
int **initarray_2Dint(int numrows, int numcolumns, int initval)
{
  int i,j;
  int **array;

   		
  ZALLOC(array, numrows, int *) ;
  for (i=0; i<numrows; i++) {
    ZALLOC(array[i],numcolumns,int);
    if (initval != 0) 
     ivclear(array[i], initval, numcolumns) ;
  }
  return array;
}

void free2Dint (int ***xx, int numrows) 
{
   int **array ;
   int i ;
   array = *xx ;

  for (i=numrows-1; i>=0; i--) {
    free(array[i]) ;
  }
  free(array) ;
}




double **initarray_2Ddouble(int numrows, int numcolumns, double initval)
{
  int i,j;
  double **array;

   		
  ZALLOC(array, numrows, double *) ;
  for (i=0; i<numrows; i++) {
    ZALLOC(array[i],numcolumns,double);
    if (initval != 0.0) 
     vclear(array[i], initval, numcolumns) ;
  }
  return array;
}

void clear2D(double ***xx, int numrows, int numcols, double val)  
{
   double **array ;
   int i ;
   array = *xx ;

  for (i=numrows-1; i>=0; i--) {
    vclear(array[i], val, numcols) ;
  }

}

void iclear2D(int ***xx, int numrows, int numcols, int val)  
{
   int **array ;
   int i ;

   array = *xx ;

  for (i=numrows-1; i>=0; i--) {
    ivclear(array[i], val, numcols) ;
  }

}


void free2D (double ***xx, int numrows) 
{
   double **array ;
   int i ;
   array = *xx ;

  for (i=numrows-1; i>=0; i--) {
    free(array[i]) ;
  }
  free(array) ;
}

void addouter(double *out, double *a, int n) 
/* 
 add outerprod(a)  to out
 trival to recode to make ~ 2 * faster
*/
{
  int i,j ;
  for (i=0; i<n; i++) {
   for (j=0; j<n; j++) {
    out[i*n+j] += a[i]*a[j] ;
   }
  }
}

void subouter(double *out, double *a, int n) 
/* 
 subtract outerprod(a)  to out
 trival to recode to make ~ 2 * faster
*/
{
  int i,j ;
  for (i=0; i<n; i++) {
   for (j=0; j<n; j++) {
    out[i*n+j] -= a[i]*a[j] ;
   }
  }
}

double bal1 (double *a, int n) 
// WARNING a is input and output
{
   double y ;

   y = asum(a, n) ; 
   if (y<=0.0) fatalx("bad bal1\n") ;
   vst(a, a, 1.0/y, n) ;
   return y ; 

}

double logmultinom(int *cc, int n) 
/* log multinomial */
{
   int t, k, i ;
   double y, ytot ;

   if (n<=1) return 0.0 ;  
   t = intsum(cc, n) ;  
   if (t==0) return 0.0 ;
   ytot = 0 ;
   for (i=0; i<n-1; i++) {  
    k = cc[i] ;
    y = logbino(t,k) ;
    ytot += y ;
    t -= k ;
   }
   return ytot ;
}
void flipiarr(int *a, int *b, int n) 
// reverse array 
{
   int *x, k ;
   ZALLOC(x, n, int) ;  

   for (k=0; k<n; ++k)  {  
    x[n-1-k] = b[k] ;
   }

   copyiarr(x, a, n) ;  

   free(x) ;


}

void fliparr(double *a, double *b, int n) 
{
   double *x ;
   int  k ;

   ZALLOC(x, n, double) ;  

   for (k=0; k<n; ++k)  {  
    x[n-1-k] = b[k] ;
   }

   copyarr(x, a, n) ;  

   free(x) ;

}
void vcompl(double *a, double *b, int n) 
// a <- 1 - b 
{
   double *x ;
   ZALLOC(x, n, double) ;  

   vvm(x, x, b, n) ;
   vsp(x, x, 1.0, n) ; 
 
   copyarr(x, a, n) ;

   free(x) ;
}
void setidmat(double *a, int n) 
// a <- identity matrix
{
  int i ;  
  vzero(a, n*n) ;
  for (i=0; i<n; i++) { 
   a[i*n+i] = 1.0 ;
  }
}
   

int stripit(double *a, double *b, int *x, int len)
// copy b to a leave out elems where x < 0
{
   int k, n ;

   n = 0 ;
   for (k=0; k<len; ++k) {
    if (x[k] >= 0)  {
     a[n] = b[k] ;
     ++n ;
    }
   }
   return n ;
}

int istripit(int *a, int *b, int *x, int len)
// copy b to a leave out elems where x < 0
{
   int k, n ;

   n = 0 ;
   for (k=0; k<len; ++k) {
    if (x[k] >= 0)  {
     a[n] = b[k] ;
     ++n ;
    }
   }
   return n ;

}
int cstripit(char **a, char **b, int *x, int len)
// copy b to a leave out elems where x < 0
{
   int k, n ;

   n = 0 ;
   for (k=0; k<len; ++k) {
    if (x[k] >= 0)  {
     a[n] = b[k] ;
     ++n ;
    }
   }
   return n ;
}

void mapit(int *a, int *b, int n, int inval, int outval) 
{
  int k ;

  copyiarr(b, a, n) ;
  for (k=0; k<n; ++k) {  
   if (a[k]==inval)  a[k] = outval ;
  }
}

int ifall(int n, int k)  
// falling factorial
{

 int prod = 1, t=n, j ;

 for (j=0; j<k; ++j) { 
  prod *= t ;
  --t ;
 }
 return prod ;
}



