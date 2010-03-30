HOMEL=$(PWD)
TOP=../bin
## binary for install  
# DEBUG_OPTIONS= -g -pg
DEBUG_OPTIONS= 
BIN=$(HOMEL)/../bin
PERLDIR=$(HOMEL)/perlsrc

# "make eigenstrat" to make eigenstrat program (ditto for other programs)
# "mv eigenstrat ../bin" (or "make install" to copy all C programs to ../bin)
# "make clean" to clean up extra files in this directory
# "make clobber" to clobber all files and subdirectories except source files
#    so as to enable recompiling from scratch.

NLIB=$(HOMEL)/smartlib/nicklib.a
IDIR=$(HOMEL)/smartinclude
IDIR3=$(HOMEL)/eigensrc
VPATH=.:nicksrc:perlsrc:eigensrc
BLAS = blas
# may need to change to BLAS = blas (depends on blas/lapack installation)

# CC=/util/bin/gcc 
CC=gcc
FF=g77     

CFLAGS= -c -pthread -I$(IDIR) -I$(IDIR3) -Wimplicit
OBJ=strsubs.o sortit.o vsubs.o statsubs.o linsubs.o getpars.o xsearch.o gauss.o	gds.o
TWTAB=\"$(HOMEL)/smarttables/twtable\"


eigx.o:		eigensrc/eigx.f
	$(FF)  -c -o eigx.o  eigensrc/eigx.f


statsubs.o:     nicksrc/statsubs.c
	$(CC)  $(CFLAGS) -DTWTAB=$(TWTAB) -o statsubs.o nicksrc/statsubs.c

M1=smartpca
M1O=smartpca.o  twsubs.o mcio.o qpsubs.o admutils.o egsubs.o regsubs.o eigsubs.o  eigx.o 

M1Z=smartrel
M1ZO=smartrel.o  twsubs.o mcio.o qpsubs.o admutils.o egsubs.o regsubs.o eigsubs.o  eigx.o 

M2=convertf
M2O=convertf.o  mcio.o admutils.o

M3=twstats
M3O=twstats.o  

M4=eigenstrat
M4O=eigenstrat.o 

M5=eigenstratQTL
M5O=eigenstratQTL.o

M6=pca
M6O=pca.o eigsubs.o eigx.o

M7=pcatoy
M7O=pcatoy.o eigsubs.o eigx.o

M8=baseprog
M8O=baseprog.o  mcio.o admutils.o  badpairs.o

P4=mergeit
P4O=mergeit.o mcio.o admutils.o

M9=smarteigenstrat
M9O=smarteigenstrat.o mcio.o admutils.o

PROGS= smartpca smartrel convertf twstats eigenstrat eigenstratQTL smarteigenstrat \
   pca baseprog mergeit pcatoy

PERLPROGS= fixgreen ploteig 

all: nicklib $(PROGS) $(PERLPROGS)

install:	all	tables
	mkdir -p $(TOP)
	cp $(PROGS) $(TOP)
	cp $(PERLDIR)/fixgreen $(TOP)
	cp $(PERLDIR)/ploteig  $(TOP)

uninstall:
	rm -f $(NLIB) tables
	cd  $(TOP) 
	rm -f $(PROGS) $(PERLPROGS)

$(M1): nicklib $(M1O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -pthread -o $(M1) $(M1O) $(NLIB) -L/broad/tools/lib64 -lg2c -lm -llapack -latlas -l$(BLAS) -Wimplicit 

$(M1Z): nicklib $(M1ZO)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -pthread -o $(M1Z) $(M1ZO) $(NLIB)  -L/broad/tools/lib64 -lm -lg2c -latlas -llapack -l$(BLAS) -Wimplicit

$(M2): nicklib $(M2O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M2) $(M2O) $(NLIB) -Wimplicit

$(M3): nicklib $(M3O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M3) $(M3O) $(NLIB) -Wimplicit

$(M4): nicklib $(M4O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M4) $(M4O) $(NLIB) -Wimplicit

$(M5): nicklib $(M5O)
	gcc  -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M5) $(M5O) $(NLIB) -Wimplicit

$(M6): nicklib $(M6O)
	gcc  -I$(IDIR) $(DEBUG_OPTIONS) -pthread -o $(M6) $(M6O) $(NLIB) -L/broad/tools/lib64 -lg2c -lm -latlas -llapack -l$(BLAS) -Wimplicit

$(M7): nicklib $(M7O)
	gcc  -I$(IDIR) $(DEBUG_OPTIONS) -pthread -o $(M7) $(M7O) $(NLIB) -L/broad/tools/lib64 -lg2c -lm -latlas -llapack -l$(BLAS)  -Wimplicit

$(M8): nicklib $(M8O)
	gcc  -I$(IDIR) $(DEBUG_OPTIONS) -o $(M8) $(M8O) $(NLIB) -lg2c -lm  -Wimplicit

$(M9): nicklib $(M9O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M9) $(M9O) $(NLIB) -Wimplicit

$(P4): nicklib $(P4O)
	gcc  -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(P4) $(P4O) $(NLIB) -Wimplicit

libnick.a:	dirs tables  $(OBJ)
	ar -r libnick.a $(OBJ)

nicklib:	dirs  tables  libnick.a 
	cp libnick.a  $(NLIB)

tables:    
	echo "tables made"  > tables
	cp twtable  $(HOMEL)/smarttables
##	cp hg17affy.map  $(HOMEL)/smarttables
	
dirs:	
	mkdir -p  $(HOMEL)/smartlib
	mkdir -p  $(HOMEL)/smarttables
	mkdir -p  $(HOMEL)/smartinclude
	mkdir -p  $(BIN)
	cp  *.h  $(IDIR)
	cp  eigensrc/*.h  $(IDIR)
	cp  nicksrc/*.h  $(IDIR)

clean: 
	rm -f *.o 
	rm -f *junk*
	rm -f core
	rm -f libnick.a
	rm -f $(PROGS)
	rm -f mcmcsrc/*.o
	rm -f eigensrc/*.o
	rm -f nicksrc/*.o

clobber: clean rmdirs uninstall


rmdirs: 
	rm -rf $(HOMEL)/smartlib 
	rm -rf $(HOMEL)/smarttables 
	rm -rf $(HOMEL)/smartinclude 

