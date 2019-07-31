###################################################################################
# NSCOUETTE by Jose Manuel Lopez (PPE formulation), Liang Shi, et al.
#
# Please cite:  
#   A hybrid MPI-OpenMP parallel implementation for pseudospectral 
#   simulations with application to Taylor-Couette flow
#   Shi, L.; Rampp, M.; Hof, B.; Avila, M. 
#   Computers & Fluids, 106, 1-11 (2015)   
#   (arXiv:1311.2481) 
#
###################################################################################
#
# Makefile for a standard x86_64 Linux software stack with Intel/GNU/PGI compilers:
#  a number of platform-specific configurations are provided in ARCH/make.arch.*
#
# requirements:
#
#  1) Compiler: a modern Fortran compiler which is OpenMP 3 compliant
#               and additionally a corresponding C-compiler
#               (tested with ifort/icc 12-15 or later, GCC 4.7-4.9, PGI 14)
#     -> select with: make COMPILER=[Intel|GNU|PGI] 
#
#     for free software see: 
#        https://gcc.gnu.org/
#
#  2) MPI: an MPI library which supports thread-level MPI_THREAD_SERIALIZED 
#          (tested with Intel MPI 4.1 and higher, IBM PE 1.3 and higher)
#     -> adapt MPIFC, CC (names of the MPI compiler wrappers in Makefile) 
#
#     for free software see: 
#        http://www.open-mpi.de/
#
#  3) BLAS/LAPACK: a serial BLAS/LAPACK library
#                  (tested with Intel MKL 11.x, OpenBLAS 0.2.14)
#      -> adapt LIBLA and/or point environment variable $MKLROOT to a MKL installation 
#
#      for free software see: 
#        http://www.openblas.net/
#        http://math-atlas.sourceforge.net/
#        http://www.netlib.org/    (reference implementations)
#           
#  4) FFTW: a serial but fully thread-safe fftw3 installation or equivalent 
#           (tested with fftw 3.3 and MKL 11.1, earlier MKL versions will likely fail)
#      -> point environment variable $FFTW_HOME to a FFTW3 installation 
#      -> alternative: select  make FFTLIB=MKL (requires MKL 11.1 or later) 
#
#      for free software see: 
#        http://www.fftw.org 
#    
#  5) HDF5: an MPI-parallel HDF5 library installation
#           (tested with 1.8.9)
#      -> point environment variable $HDF5_HOME to an MPI-parallel HDF5 installation 
#      -> alternative: switch data output off with: make HDF5IO=no
#
#      for free software see: 
#        http://www.hdfgroup.org/HDF5/
#
#############################################################################

#defaults for external variables 
ARCH     ?= intel-mkl
CODE     ?= STD_CODE
HDF5IO   ?= yes
PROFLIB  ?= none
DEBUG    ?= no
#TESTS    ?= ,-DTEST1,-DTEST2

#platform specific definitions
include ARCH/make.arch.$(ARCH)

ifneq ($(DEBUG),no)
        DBGFLAGS=$(DEBUGFLAGS)
endif

ifneq ($(HDF5IO),no)
	FSRC_HDF5 = mod_hdf5io.f90
	DHDF5 = $(PPSEP)-DHDF5IO
endif

LIBPROF =
FSRC_PROF = perfdummy.f90
OBJ_PROF = perfdummy.o
ifeq ($(PROFLIB),FTIM)
	CSRC_PROF = time.c papi.c resident_set_size.c
	FSRC_PROF = ftimings_type.f90 ftimings_value.f90 ftimings.f90 \
		    mod_timings.f90 perf2ftimings.f90
#       INCLUDES += -I$(PAPI_HOME)/include
endif
ifeq ($(PROFLIB),LIKWID)
	CSRC_PROF = 
	FSRC_PROF = perf2likwid.f90
        INCLUDES += -I$(LIKWID_HOME)/include
	LIBPROF = -L$(LIKWID_HOME)/lib -llikwid $(RPATH)=$(LIKWID_HOME)/lib
endif

DOCS = doc/nsCouetteLopez2019.pdf doc/nsCouetteShi2015.pdf doc/nsCouetteUserGuide.pdf

FSOURCES = mod_getcpu.f90 wctimerclass.f90 \
           mod_fftw.f90 mod_fdInit.f90 mod_params.f90 mod_vars.f90 \
	   mod_myMpi.f90 $(FSRC_HDF5) mod_inOut.f90 mod_nonlinear.f90 \
           mod_timeStep.f90 nsCouette.f90 $(FSRC_PROF) 
CSOURCES = getcpu.c $(CSRC_PROF)


OBJECTS = $(CSOURCES:.c=.o) $(FSOURCES:.f90=.o) 
SOURCES = $(CSOURCES) $(FSOURCES)

BUILD = $(ARCH)

OBJS := $(addprefix $(BUILD)/, $(OBJECTS))
FSOURCES_PP := $(addprefix $(BUILD)/, $(FSOURCES))




INCLUDES += -I$(INCFFT)

LDFLAGS = 

GITDEV := 'GIT_REV="$(shell git rev-parse --short=8 HEAD || echo 'undef')"'
ARCH_ID := 'ARCH_ID="$(ARCH)"'
CMP_OPTS := 'CMP_OPTS="$(MPIFC) $(FFLAGS)"'

PPFLAGS = $(PPPRE)-D$(CODE)$(DHDF5)$(DEFINES)$(TESTS)
IDFLAGS = $(PPPRE)-D$(GITDEV)$(PPSEP)-D$(ARCH_ID)$(PPSEP)-D$(CMP_OPTS)

$(BUILD)/nsCouette.x: $(OBJS)
	$(MPIFC) $(FFLAGS) $(DBGFLAGS) $(LDFLAGS) -o $@ $^ $(LIBFFT) $(LIBHDF5) $(LIBLA) $(LIBPROF)

$(BUILD)/%.o: %.c
	cd $(BUILD) && $(CC) $(CFLAGS) $(INCLUDES) -c ../$< -o $*.o

$(BUILD)/%.o: %.f90
	cd $(BUILD) && $(MPIFC) $(FFLAGS) $(DBGFLAGS) $(INCLUDES) $(PPFLAGS) $(IDFLAGS) -c ../$< -o $*.o

$(BUILD)/%.f90: %.f90
	cd $(BUILD) && cpp -traditional -E -P $(PPFLAGS) ../$< > $<

$(OBJS): | $(BUILD)

$(BUILD):
	mkdir $(BUILD)

# compilation of the postprocessing utility waveSpeed
PREPROC_OBJECTS = mod_fdInit.o mod_params.o mod_vars.o mod_myMpi.o mod_preAnalysis.o waveSpeed.o
OBJS3 := $(addprefix $(BUILD)/, $(PREPROC_OBJECTS))


# compilation of the postprocessing utilities
postproc: $(BUILD)/waveSpeed.x

$(BUILD)/waveSpeed.x: $(OBJS3)
	$(MPIFC) $(FFLAGS) $(DBGFLAGS) $(LDFLAGS) -o $@ $^ $(LIBFFT) $(LIBLA)

$(BUILD)/mod_preAnalysis.o: postproc/mod_preAnalysis.f90
	@echo $@
	cd $(BUILD) && $(MPIFC) $(FFLAGS) $(DBGFLAGS) $(INCLUDES) $(PPFLAGS) -c ../$< -o ../$*.o
$(BUILD)/waveSpeed.o: postproc/waveSpeed.f90
	cd $(BUILD) && $(MPIFC) $(FFLAGS) $(DBGFLAGS) $(INCLUDES) $(PPFLAGS) -c ../$< -o ../$*.o



#dependencies created by deplist target
#include make.inc.dep
$(BUILD)/getcpu.o : getcpu.c 
$(BUILD)/perfdummy.o : perfdummy.f90 
$(BUILD)/mod_getcpu.o : mod_getcpu.f90 
$(BUILD)/mod_fftw.o : mod_fftw.f90 
$(BUILD)/mod_fdInit.o : mod_fdInit.f90 
$(BUILD)/mod_params.o : mod_params.f90 
$(BUILD)/mod_vars.o : mod_vars.f90 $(BUILD)/mod_params.o $(BUILD)/mod_fftw.o 
$(BUILD)/mod_myMpi.o : mod_myMpi.f90 $(BUILD)/mod_params.o $(BUILD)/mod_vars.o 
$(BUILD)/mod_inOut.o : mod_inOut.f90 $(BUILD)/mod_vars.o $(BUILD)/mod_fftw.o $(BUILD)/mod_myMpi.o 
$(BUILD)/mod_nonlinear.o : mod_nonlinear.f90 $(BUILD)/mod_inOut.o $(BUILD)/mod_fdInit.o $(BUILD)/mod_fftw.o $(BUILD)/mod_myMpi.o 
$(BUILD)/mod_timeStep.o : mod_timeStep.f90 $(BUILD)/mod_nonlinear.o $(BUILD)/mod_inOut.o 
$(BUILD)/nsCouette.o : nsCouette.f90 $(BUILD)/mod_timeStep.o $(BUILD)/mod_getcpu.o 
ifeq ($(HDF5IO),yes)
$(BUILD)/mod_hdf5io.o : mod_hdf5io.f90 
endif
ifeq ($(PROFLIB),FTIM)
$(BUILD)/time.o : time.c 
$(BUILD)/papi.o : papi.c 
$(BUILD)/resident_set_size.o : resident_set_size.c 
$(BUILD)/ftimings_type.o : ftimings_type.f90 
$(BUILD)/ftimings_value.o : ftimings_value.f90 $(BUILD)/ftimings_type.o 
$(BUILD)/ftimings.o : ftimings.f90 $(BUILD)/ftimings_value.o $(BUILD)/ftimings_type.o 
$(BUILD)/mod_timings.o : mod_timings.f90 $(BUILD)/ftimings.o 
$(BUILD)/perf2ftimings.o : perf2ftimings.f90 $(BUILD)/mod_timings.o 
endif


doc: doc/nsCouetteUserGuide.pdf
doc/nsCouetteUserGuide.pdf: doc/*.tex
	cd doc && pdflatex nsCouetteUserGuide;pdflatex nsCouetteUserGuide

#create list of dependencies
deplist: 
	@makedepf90 -b '$$(BUILD)' $(SOURCES) > make.inc.dep

clean:
	cd $(BUILD) && rm -f *.o *.mod *.lst *.i90 *.f90

forcheck: $(BUILD) $(FSOURCES_PP)
	cd $(BUILD) && forchk -decl -f08  -ff -ancmpl -anprg -anref -shinc -shprg -shref -l forcheck.out -moddep -refstruc -I $(INCFFT) *.f90 -include $(FCKDIR)/share/forcheck/MPI_3.flb; if [ $$? -le 2 ] ; then exit 0 ; else cat forcheck.out; exit 1;fi


ford: $(BUILD) $(FSOURCES_PP)
	cd $(BUILD) && sed "s,XXXincludeXXX,$(INCFFT)," ../ford-config.md > ford-config.md && ford ford-config.md
	mv $(BUILD)/ford-doc .

tar: Makefile ARCH $(SOURCES)
	tar cvzf nsCouette_`date "+%F"`.tgz Makefile ARCH $(SOURCES)

export: $(SOURCES) Makefile ARCH/* doc
	mkdir -p nsCouette_dist
	cp -r $(SOURCES) Makefile ARCH input_nsCouette $(DOCS) scripts nsCouette_dist/
	@./export_public.sh nsCouette_dist
	tar czf nsCouette_dist_`date "+%F"`.tgz nsCouette_dist/

check: 
	.gitlab-ci/test/validate_nscouette.sh $(ARCH) 4 4
	make ARCH=$(ARCH) clean


help:
	@echo "selected architecture: $(ARCH)"
	@echo "1) adjust platform-specific settings in ARCH/make.arch.$(ARCH)"
	@echo "2) build with: make ARCH=<architecture> HDF5IO=<yes|no> PROFLIB=<none|FTIM> DEBUG=<no|yes>"
