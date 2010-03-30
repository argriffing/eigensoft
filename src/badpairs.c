#include <stdio.h>
#include <string.h>

#include "admutils.h" 
#include "mcio.h"  
#include "badpairs.h"




void
dobadpairs(char *badpairsname, SNP **snpm, int numsnps)
{

   FILE *fff ;
   int i, n  ;
   SNP *cupt1, *cupt2 ;
   char line[MAXSTR] ;
   char *spt[MAXFF], *sx ;
   int nsplit, num=0, k1, k2, *g1, *g2  ;
   double y1, y2 ;
   int tot = 0, t, thash = 0 ;


   if (badpairsname == NULL)  { 
      printf("*** warning no badpairsname\n") ;
      return ;
   }
   else openit(badpairsname, &fff, "r") ;

   while (fgets(line, MAXSTR, fff) != NULL)  {
    nsplit = splitup(line, spt, MAXFF) ; 
    if (nsplit ==0 ) continue ;
    sx = spt[0] ;
    if (sx[0] == '#')  { 
     freeup(spt, nsplit) ;
     continue ;
    }
    k1 = snpindex(snpm, numsnps, spt[0]) ;
    k2 = snpindex(snpm, numsnps, spt[1]) ;
    if (k1<0) {  
     printf("*** warning.  snp %s in badpairsname but not in main snp file\n", spt[0]) ;
     freeup(spt, nsplit) ;
     continue ;
    }
    if (k2<0) {  
     printf("*** warning.  snp %s in badpairsname but not in main snp file\n", spt[1]) ;
     freeup(spt, nsplit) ;
     continue ;
    }
    cupt1 = snpm[k1] ;
    cupt2 = snpm[k2] ;
    if ( (cupt1 -> ignore) || (cupt2 -> ignore)) {  
     freeup(spt, nsplit) ;
     continue ;
    }
    if ( (cupt1 -> isfake) || (cupt2 -> isfake)) {  
     freeup(spt, nsplit) ;
     continue ;
    }

    if ( (cupt1 -> isrfake) || (cupt2 -> isrfake)) {  
     freeup(spt, nsplit) ;
     continue ;
    }

    n = cupt1 -> ngtypes ;
    if (cupt1 -> score > cupt2 -> score)  { 
     cupt1 = snpm[k2] ;
     cupt2 = snpm[k1] ;
    }
    for (i=0; i<n; i++) { 
     if (getgtypes(cupt1, i) < 0) continue ;
     if (getgtypes(cupt2, i) < 0) continue ;
     putgtypes(cupt1, i, -1) ; 
     ++tot ;
    }
     freeup(spt, nsplit) ;
   }

  fclose(fff)  ;  
  printf("dobadpairs: removed %d genotypes\n", tot) ;
}



void loadbadpsc(SNP **snpm, int numsnps, int randombadpairs, char *goodsnpname) 
{

  SNP *cupt ;
  double dd[4] ;
  double theta = 0.2 ;
  double tt[4], pa, pe, qa, qe, pt, qt  ;
  double xa, xe, xx ;
  double score ;
  int num, t ;

  pt = theta ; qt = 1.0 - pt ;
  for (num= 0; num < numsnps; ++num) {

   cupt = snpm[num] ;
   if (cupt -> ignore) continue  ;

   tt[0] = cupt -> af_nn[0] ;
   tt[1] = cupt -> af_nn[1] ;
   tt[2] = cupt -> cauc_nn[0] ;
   tt[3] = cupt -> cauc_nn[1] ;

   xa = asum(tt,2) ;
   xe = asum(tt+2,2) ;
   pa = tt[0]/xa ;
   pe = tt[2]/xe ;
   qa = 1.0 - pa ;
   qe = 1.0 - pe ;

   dd[0] = qt*pa ;
   dd[1] = qt*qa ;
   dd[2] = pt*pe ;
   dd[3] = pt*qe ;

   if (randombadpairs) {  
    score = 1000000.0*DRAND()  ;
   }
   else {
    t  = nnint(1000000.0*mutx(dd)) ;
    score = (double) t ;
    score +=  .1 * gauss() ; //  dither 
   }
   cupt -> score = score ; 
  }
  dogood(goodsnpname, snpm, numsnps) ;
}
void 
dogood(char *goodsnpname, SNP **snpm, int numsnps) 
{
   FILE *fff ;
   int i, n  ;
   SNP *cupt ;
   char line[MAXSTR] ;
   char *spt[MAXFF], *sx ;
   int nsplit, num=0, k1, k2, *g1, *g2  ;
   double y1, y2 ;
   int tot = 0 ;


   if (goodsnpname == NULL)  { 
      return ;
   }
   else openit(goodsnpname, &fff, "r") ;

   while (fgets(line, MAXSTR, fff) != NULL)  {
    nsplit = splitup(line, spt, MAXFF) ; 
    if (nsplit ==0 ) continue ;
    sx = spt[0] ;
    if (sx[0] == '#')  { 
     freeup(spt, nsplit) ;
     continue ;
    }
    k1 = snpindex(snpm, numsnps, spt[0]) ;
    if (k1<0) {  
     printf("*** warning.  snp %s in goodpairsname but not in main snp file\n", spt[0]) ;
     freeup(spt, nsplit) ;
     continue ;
    }
    cupt = snpm[k1] ;
    cupt -> score += 1.0e7 ;
    ++tot ;
   }

  fclose(fff) ;
  printf("number of goodsnps:  %d\n", tot) ;

}

double mutxx(double *dd, int m, int n) 
{
   double *aa, *bb ; 
   double y ;

   ZALLOC(aa, m, double) ;
   ZALLOC(bb, n, double) ;

   rrsum(dd, bb, m, n) ;
   ccsum(dd, aa, m, n) ;
/**
   if ( ((m*n)<=4)  && verbose) {  
    printf("m: %d n: %d\n",m,n) ;
    printmat(dd,m,n) ;
    printf("\n") ;
    printmat(aa,1,m) ;
    printmat(bb,1,n) ;
   }
*/
   y = entrop(aa,m) + entrop(bb,n) - entrop(dd,m*n) ;

  free(aa) ;  free(bb) ;
  return y ;

}
double mutx(double *dd) 
/** mutual inf from 2 x 2 */
{
    double y ;
    y = mutxx(dd, 2,2 ) ;
    return y  ;
}

void
printsnpsc(char *snpscname, SNP **snpm, int numsnps, double *sc) 
{

   FILE *fff ;
   int i, score  ;
   SNP *cupt ;


   if (snpscname == NULL)  fff = stdout ;  
   else openit(snpscname, &fff, "w") ;

   for (i=0; i<numsnps; i++)  {  
    cupt = snpm[i] ;
    if (cupt -> ignore) continue ;
    if (sc[i] < 0.0) continue ;
    score  = nnint(1.0e6*sc[i]) ;
    fprintf(fff, "%20s  %8d\n", cupt -> ID, score) ;
   }
   if (snpscname != NULL)  fclose(fff)  ;  
}


int 
killsnps(Indiv **indivmarkers, SNP **snpmarkers, int numsnps, int mincasenum) 
{
  int i, t, tot=0 ; 
  SNP *cupt ;

  for (i=0; i< numsnps; i++) {  
   cupt = snpmarkers[i] ;
   if (cupt -> isfake) continue ; 
   if (cupt -> ignore) continue ;
   t =  numvalidgtx(indivmarkers, cupt, 1)  ;          
   if (t < mincasenum) { 
     printf(" removing SNP %s\n", cupt -> ID) ;
     cupt-> ignore = YES ;
     ++tot ;
   }
  }
  return tot ;
}
