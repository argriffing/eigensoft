#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "vsubs.h"  
#include "strsubs.h" 
#include "linsubs.h" 

void bal(double *a, double *b, int n) 
/** 
 normalize mean 0 s.d 1 
*/
{
    double t ;
    t = asum(b,n)/ (double) n ;
    vsp (a, b, -t, n) ;

    t = asum2(a,n)/ (double) n ;
    vst (a, a, 1.0/sqrt(t), n) ;
}

void mulmat(double *a, double *b, double *c, int a1, int a2, int a3) 
/* b is a1 x a2 , c a2 x a3 so a is a1 x a3  */  
{
    double *t ;
    int i,j,k ;
    ZALLOC(t, a1*a3, double) ;

    for (i=0; i<a1; i++)  
     for (j=0; j<a3; j++)  
      for (k=0; k<a2; k++)  
       t[i*a3+j] += b[i*a2+k]*c[k*a3+j] ;

    copyarr(t, a, a1*a3) ;

    free (t) ;

}

void pdinv(double *cinv, double *coeff, int n) 
{
   double *tt;
   double *p ;
   double t, sum ;
   int i,j, k ;

/**
   pmat(coeff, n) ;
*/
   ZALLOC (tt, n*n, double);
   ZALLOC (p, n, double );
   

  copyarr(coeff,tt,n*n); 
  
  choldc (tt, n, p) ;

 
  for (i=0; i<n; i++) {
    tt[i*n+i] = 1.0/p[i] ;
    for (j=i+1; j<n; j++) {
      sum=0.0 ;
      for (k=i; k<j; k++) {
        sum -= tt[j*n+k]*tt[k*n+i] ;
      }
      tt[j*n+i] = sum/p[j] ;

    }
  }

   for (i=0; i<n; i++) 
    for (j=i; j<n; j++) {
     sum=0.0 ;
     for (k=j; k<n; k++) {
      sum += tt[k*n+j]*tt[k*n+i] ;
     
     }
     cinv[i*n+j] = cinv[j*n+i] = sum ;

    }
   

   free(tt) ;
   free(p) ;

}

int            
solvit (double *prod, double *rhs, int n, double *ans) 
{ 
  //The coefficient matrix should be positive definite

  /*AT : changed this code to take in matrix in a linear array form*/
  double *ttt;                   
  double *b; 
  double *p; 
  int i ; 
  int ret ;
 
 
  ZALLOC (ttt, n*n, double); 
  ZALLOC (p, n, double); 
  ZALLOC(b,n,double); 
 
  copyarr(prod,ttt,n*n); 
  copyarr(rhs,b,n); 
 
  ret = choldc (ttt, n, p); 
  if (ret<0) return -1 ;  // not pos def
  cholsl (ttt, n, p, b, ans); 
   
  free (ttt) ;  
  free(b);           
  free (p) ; 
 
  return 1 ;
} 

void 
cholsl (double *a, int n, double p[], double b[], double x[])

/** 
 Numerical Recipes.  Must change 
*/

{
  /*AT: Changing the code*/

 int i, k; 
  double sum; 
     
 
  for (i = 0; i < n; i++) 
    { 
      sum = b[i]; 
      for (k = i - 1; k >= 0; k--) 
        sum -= a[i*n+k] * x[k]; 
      x[i] = sum / p[i]; 
    } 
 
  for (i = (n-1); i >= 0; i--) 
    { 
      sum = x[i]; 
      for (k = i + 1; k < n; k++) 
        sum -= a[k*n+i]* x[k]; 
      x[i] = sum / p[i]; 
    }        


}

int 
choldc (double *a, int n, double p[])
{
   int i, j,k;         
  double sum; 
     
  
  for (i = 0; i < n; i++) 
    {           
      for (j = i; j < n; j++) 
        {                 
          sum = a[i*n+j]; 
          for  (k = i - 1; k >= 0; k--) 
            sum -= a[i*n+k] * a[j*n+k]; 
          if (i == j) 
            {          
	      /**                                       
              printf("zzchol %d %20.10f %9.3f\n",i, sum, a[i][i]) ; 
	      */                                  
              if (sum <= 0.0) {                  
                return -1 ; // not pos def
              }  
              p[i] = sqrt (sum); 
            }   
          else  
            { 
              a[j*n+i] = sum / p[i]; 
               
            } 
        }                                    
    } 
  
    return 1 ;
                
}

void pmat(double *mat, int n)  

/** 
 print square matrix 
*/

{
 int i,j ;
 double *diag ;

 ZALLOC(diag,n, double) ;
 getdiag(diag, mat, n) ;

 printf("pmat:\n") ;
 
 for (i=0; i<n; i++) {
  printf("diag %5d %9.3f\n",i, diag[i]) ;
  for (j=0; j<n; j++) {
   if ((j%10) == 9)  printf("\n") ;
   if ((n%10) != 0) printf("%9.3f ",mat[i*n+j]) ;
  }
  printf("\n") ;
 }
 printf("\n") ;
 printf("\n") ;

 free(diag) ;
}

void cholesky(double *cf, double *a, int n) 
{  
  int i, j, k ;
  double *tt ;
  double *p ;
  
  ZALLOC(tt, n*n, double) ;
  ZALLOC(p, n, double) ; 
  copyarr(a,tt,n*n);


   choldc(tt, n, p ) ;

   vzero(cf, n*n) ;

   for (i = 0; i < n; i++) {
    tt[i*n+i] = p[i] ;
    
    for (j=0; j <= i ; j++) {  
     k = (i)*n+(j) ;
     cf[k] = tt[i*n+j] ;
    }
   }
  
   free(tt) ; 
   free(p) ;
}
