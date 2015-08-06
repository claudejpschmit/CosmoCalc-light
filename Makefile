#################################################
# This file should just make a bunch of objects,#
# to simulate what class does.                  #
#################################################

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build
LIBRARIES = lib/
.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

vpath %.cpp src:main:$(LIBRARIES)ALGLIB_source:$(LIBRARIES)GLOBAL21CM_source
vpath %.o build
vpath .base build

# Compiler required for c++ code.
# including -ffast-math may not be as bad as anticipated.
CXX = g++ -Wall -std=c++11 -ffast-math -s -Wno-deprecated 
# Compiler required for fortran RECFAST
FF = gfortran

OPTFLAG = -O4
ARMAFLAGS = -larmadillo

# This decreases the precision of the bessel functions significantly if 
# added to the compilation of the files containing boost->sph_bess(l,x).
OPTFLAG_CLASS = -ffast-math
OMPFLAG = -fopenmp
CCFLAG = -g -fPIC
LDFLAG = -g -fPIC

# leave blank to compile without HyRec, or put path to HyRec directory
# (with no slash at the end: e.g. CLASS_hyrec or ../CLASS_hyrec)
# automatically add external programs if needed. First, initialize to blank.
EXTERNAL =

INCLUDES = -I../include
INCLUDES += -I../$(LIBRARIES)ALGLIB_include
INCLUDES += -I../$(LIBRARIES)GLOBAL21CM_include

# These lines seem to be unnecessary, but I leave them in anyways.
INCLUDES += -I/usr/include/boost
LINKER = -L/usr/include/boost #-lboost_filesystem

%.o: %.cpp .base
	cd $(WRKDIR);$(CXX) $(OPTFLAG) $(OMPFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o $(ARMAFLAGS)

# This line creates the CLASS objects.
%.o: %.c .base
	cd $(WRKDIR);$(CC) $(OPTFLAG) $(OPTFLAG_CLASS) $(OMPFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

ALGLIB = alglibinternal.o alglibmisc.o ap.o dataanalysis.o diffequations.o fasttransforms.o integration.o interpolation.o linalg.o optimization.o solvers.o specialfunctions.o statistics.o
GLOBAL21CM = dnumrecipes.o dcomplex.o dcosmology.o astrophysics.o twentyonecm.o spline.o spline2D.o
SRC = CosmoBasis.o CosmologyCalculatorClass.o FisherClass.o CAMB_interface.o Global21cmInterface.o LevinIntegrator.o
MAIN = Main.o

all: calc 

calc: $(SRC) $(SOURCE) $(TOOLS) $(OUTPUT) $(EXTERNAL) $(ALGLIB) $(GLOBAL21CM) $(ODE) $(MAIN) 
	cd $(MDIR);$(CXX) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(LINKER) -o calc $(addprefix build/, $(notdir $^)) -lm $(ARMAFLAGS)

install:
	$(FF) -o $(LIBRARIES)GLOBAL21CM_dependencies/RECFAST_CODE/recfast $(LIBRARIES)GLOBAL21CM_dependencies/RECFAST_CODE/recfast.for
	cd $(MDIR)
	make all

clean_integration: .base
	rm $(WRKDIR)/Integrator.o;
	rm $(WRKDIR)/Main.o;
	rm $(WRKDIR)/CosmologyCalculatorClass.o;
	rm calc

clean: .base
	rm -rf $(WRKDIR);
	rm calc;
