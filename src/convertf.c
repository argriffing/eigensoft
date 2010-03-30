#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>

#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  

#define WVERSION   "2401" 
/** 
 reformats files.             
 pedfile junk (6, 7 cols, ACGT added)
 packped (.bed) added  input and output 
 pedfile can now be in random order
 tersemode supported 
 ped -> ancestrymap outputs alleles
 various fixups for ped files in computing variance, reference alleles 
 alleles can (say) be A X if second allele unknown
 allele flipped for bed format
 huge ascii file trapped now
 fastdup added
 r2thresh added
 killr2 added to convertf (default is NO)  
 bug in processing C0 SNPS in ped files fixed

 outputall added  (all indivs, snps output)
 genolistname added (list of packed ancestrymap files one per line)
 support for enhanced eigenstrat format and genolist can also be enhanced eigenstrat
 Chromosome # 1-89 now supported.  Improved treatment of funny chromosomes. 
 MT -> 90 
 XY -> 91

 Various bugfixes and checks for files size  
 inpack2 no longer needs big memory 
 hash table lookup for snpindex
 check size more carefully
 update to mcio.c  (snprawindex)
 maxmissfracsnp and maxmissfracind added
*/


#define MAXFL  50   
#define MAXSTR  512
#define NO_PENDING 3

char *trashdir = "/var/tmp" ;
int verbose = NO ;
int qtmode = NO ;
Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 

char  *genotypename = NULL ;
char  *genotypelist = NULL ;

char  *snpname = NULL ;
char *indoutfilename = NULL ;
char *snpoutfilename = NULL ;
char  *genooutfilename = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;

double r2thresh = -1.0 ;  
double r2genlim = 0.01 ; // Morgans 
double r2physlim = 5.0e6 ;     
double maxmissfracsnp = 1000.0 ; // no thresh
double maxmissfracind = 1000.0 ; // no thresh
int killr2 = NO ;  

int packout = -1 ;
int tersem  = YES ;
extern enum outputmodetype outputmode  ;
extern int checksizemode ;
char *omode = "packedancestrymap" ;
extern int packmode ;
int ogmode = NO ;
int fastdup = NO ;
int fastdupnum = 0  ;

int xchrom = -1 ;
int lopos = -999999999 ; 
int hipos = 999999999 ;
int minchrom = 1 ; 
int maxchrom = 97 ;


char  unknowngender = 'U' ;

void setomode(enum outputmodetype *outmode, char *omode)  ;
void readcommands(int argc, char **argv) ;
void outfiles(char *snpname, char *indname, char *gname, SNP **snpm, 
  Indiv **indiv, int numsnps, int numind, int packem, int ogmode) ;


int main(int argc, char **argv)
{

  int **snppos ;
  int *snpindx ;
  char **snpnamelist, **indnamelist ;
  char **eglist ;
  int  lsnplist, lindlist, numeg ;
  int i,j; 
  SNP *cupt, *cupt1, *cupt2, *cupt3 ;
  Indiv *indx ;
  double gpos1,gpos2,cpos1,cpos2,gd, cd, gd100 ;
  double rthresh, zt ;
  int mpflag, ret, nvalid;

  int ch1, ch2 ;
  int fmnum , lmnum ;
  int num, n1, n2 ;
  int nkill = 0 ;
  int t, k ;

  int nindiv = 0, e, f, lag=1  ;
  double xc[9], xd[4], xc2[9] ;
  double ychi, zscore, zthresh = 20.0 ;
  double y1, y2 ; 
  int nignore, numrisks = 1 ;

  char **genolist ;
  int numgenolist ;
  int maxmiss ; 

  malexhet = YES ;    // convertf default is don't change the data
  tersem = YES ;     // no snp counts

  readcommands(argc, argv) ;

  setomode(&outputmode, omode) ;
  packmode = YES ;
  settersemode(tersem) ;

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;


  for (i=0; i<numsnps; i++)  {  
   if (xchrom == -1) break ;  
   cupt = snpmarkers[i] ; 
   if (cupt -> chrom != xchrom) cupt -> ignore = YES ; 
   if (cupt -> ignore) continue ; 
   t = nnint(cupt -> physpos) ; 
   if ( (t< lopos) || (t >hipos)) cupt -> ignore = YES ;
  }

  nignore = 0 ;
  for (i=0; i<numsnps; i++)  {  
   cupt = snpmarkers[i] ; 
   if (cupt -> chrom > maxchrom) cupt -> ignore = YES ;  
   if (cupt -> chrom < minchrom) cupt -> ignore = YES ;  
   if (cupt -> ignore) ++nignore ;
  }

  if (numsnps == nignore) fatalx("no valid snps\n") ;

/**
  cupt = snpmarkers[0] ;
  printf("zz2: %d %d %d %20s: %d\n", numsnps, nignore, cupt -> chrom, cupt -> ID, cupt -> ignore) ;
*/

  numindivs = getindivs(indivname, &indivmarkers) ;

  if (genotypelist!= NULL) {  
    getgenos_list(genotypelist, snpmarkers, indivmarkers, 
     numsnps, numindivs, nignore) ;
  }

  else {
   setgenotypename(&genotypename, indivname) ;
   getgenos(genotypename, snpmarkers, indivmarkers, 
     numsnps, numindivs, nignore) ;
  }

  if (outputall) {
   outfiles(snpoutfilename, indoutfilename, genooutfilename, 
    snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode) ;

   printf("##end of convertf run (outputall mode)\n") ;
   return 0 ;
  }

  if (killr2) {
   nkill = killhir2(snpmarkers, numsnps, numindivs, r2physlim, r2genlim, r2thresh) ;
   if (nkill>0) printf("killhir2.  number of snps killed: %d\n", nkill) ;
  }

  setstatus(indivmarkers, numindivs, "Case") ;

  /******************************************************************/
  /* removesubthreshold(indivmarkers, snpmarkers, numindiv, numsnps, 
    maxmissfracind, maxmissfracsnp); */
  /******************************************************************/

  if (fastdup)  {  
   if (fastdupnum > 0) setfastdupnum(fastdupnum) ;
   fastdupcheck(snpmarkers, indivmarkers, numsnps, numindivs) ;  
  }

  if (decim>0) {  
   snpdecimate(snpmarkers, numsnps, decim, dmindis, dmaxdis) ;
  }

  outfiles(snpoutfilename, indoutfilename, genooutfilename, 
   snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode) ;

  printf("##end of convertf run\n") ;
  return 0 ;
}
void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  char *parname = NULL ;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n ;

  while ((i = getopt (argc, argv, "p:vV")) != -1) {

    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'v':
	printf("version: %s\n", WVERSION) ; 
	break; 

      case 'V':
	verbose = YES ;
	break; 

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
   pcheck(parname,'p') ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

/**
DIR2:  /fg/nfiles/admixdata/ms2
SSSS:  DIR2/outfiles 
genotypename: DIR2/autos_ccshad_fakes
eglistname:    DIR2/eurlist  
output:        eurout
*/
   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "genotypelist:", &genotypelist) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "badsnpname:", &badsnpname) ;
   getstring(ph, "indoutfilename:", &indoutfilename) ;
   getstring(ph, "indivoutname:", &indoutfilename) ; /* changed 11/02/06 */
   getstring(ph, "snpoutfilename:", &snpoutfilename) ;
   getstring(ph, "snpoutname:", &snpoutfilename) ; /* changed 11/02/06 */
   getstring(ph, "genooutfilename:", &genooutfilename) ; 
   getstring(ph, "genotypeoutname:", &genooutfilename) ; /* changed 11/02/06 */
   getstring(ph, "outputformat:", &omode) ;
   getstring(ph, "outputmode:", &omode) ;
   getint(ph, "checksizemode:", &checksizemode) ;

   getint(ph, "outputgroup:", &ogmode) ;
   getint(ph, "malexhet:", &malexhet) ;
   getint(ph, "nomalexhet:", &malexhet) ; /* changed 11/02/06 */
   getint(ph, "tersemode:", &tersem) ;
   getint(ph, "familynames:", &familynames) ;
   getint(ph, "packout:", &packout) ; /* now obsolete 11/02/06 */
   getint(ph, "decimate:", &decim) ; 
   getint(ph, "dmindis:", &dmindis) ; 
   getint(ph, "dmaxdis:", &dmaxdis) ; 
   getint(ph, "fastdup:", &fastdup) ;
   getint(ph, "fastdupnum:", &fastdupnum) ;
   getint(ph, "killr2:",  &killr2) ;
   getint(ph, "hashcheck:", &hashcheck) ;
   getint(ph, "outputall:", &outputall) ;
   getint(ph, "sevencolumnped:", &sevencolumnped) ;

   getdbl(ph, "r2thresh:", &r2thresh) ;
   getdbl(ph, "r2genlim:", &r2genlim) ;
   getdbl(ph, "r2physlim:", &r2physlim) ;

   getint(ph, "chrom:", &xchrom) ;
   getint(ph, "lopos:", &lopos) ;
   getint(ph, "hipos:", &hipos) ;

   getint(ph, "qtmode:", &qtmode) ;

   getint(ph, "minchrom:", &minchrom) ;
   getint(ph, "maxchrom:", &maxchrom) ;
   getdbl(ph, "maxmissfracsnp:", &maxmissfracsnp) ;
   getdbl(ph, "maxmissfracind:", &maxmissfracind) ;
   

   writepars(ph) ;
   closepars(ph) ;

}
