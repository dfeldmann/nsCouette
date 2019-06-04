# NSCOUETTE

The code NSCOUETTE is designed to simulate incompressible
Taylor-Couette flows with infinitely long cylinders (periodic in the axial 
direction) on high-performance computers. 
The code is written in modern Fortran and uses a hybrid MPI-OpenMP parallelization
strategy.

This version implements a pressure-poisson equation (PPE) formulation and a
variable timestep with automatic CFL control.

## Build status and test coverage

[![Build status](https://gitlab.mpcdf.mpg.de/mjr/nscouette/badges/master/build.svg)](https://gitlab.mpcdf.mpg.de/mjr/nscouette/commits/master)
[![Code coverage](https://gitlab.mpcdf.mpg.de/mjr/nscouette/badges/master/coverage.svg)](http://mjr.pages.mpcdf.de/nscouette/) 

## Installation

To build and install NSCOUETTE a Makefile for a standard x86_64 Linux
software stack with Intel/GNU/PGI compilers is used, for example:

```
#clone the repository
git clone ...
cd nscouette
#build the user guide
make doc
#build the code (e.g. with Intel compilers and MKL)
make ARCH=intel-mkl HDF5IO=yes
```

Dependencies and software requirements:

* Compiler: a modern Fortran compiler which is OpenMP 3 compliant
              and additionally a corresponding C-compiler
              (tested with ifort/icc 12-15 or later, GCC 4.7-4.9, PGI 14)
              
    -> select with: make COMPILER=[Intel|GNU|PGI] 

    for free software see: https://gcc.gnu.org/

* MPI: an MPI library which supports thread-level MPI_THREAD_SERIALIZED 
         (tested with Intel MPI 4.1 and higher, IBM PE 1.3 and higher)

    -> adapt MPIFC, CC (names of the MPI compiler wrappers in Makefile) 

    for free software see: http://www.open-mpi.org/

* BLAS/LAPACK: a serial BLAS/LAPACK library
                 (tested with Intel MKL 11.x, OpenBLAS 0.2.14)

     -> adapt LIBLA and/or point environment variable $MKL_ROOT to a MKL installation 

     for free software see: http://www.openblas.net/, http://math-atlas.sourceforge.net/, http://www.netlib.org/ (reference implementation)
          
* FFTW: a serial but fully thread-safe fftw3 installation or equivalent 
          (tested with fftw 3.3 and MKL 11.1, earlier MKL versions will likely fail)

     -> point environment variable $FFTW_ROOT to a FFTW3 installation
        or use make FFTLIB=MKL to select MKL (requires MKL 11.1 or later) 

     for free software see: http://www.fftw.org 
   
* HDF5 (optional): an MPI-parallel HDF5 library installation
          (tested with hdf5 1.8.9)

     -> point environment variable $HDF5_ROOT to an MPI-parallel HDF5
        installation or use make HDF5IO=no to switch off HDF5 data output

     for free software see: http://www.hdfgroup.org/HDF5/


## Contributors


Following is a running list of contributors in chronological order:

1. [Prof. Marc Avila](https://www.zarm.uni-bremen.de/en/research/fluid-dynamics/fluid-simulation-and-modeling.html), University of Bremen

2. Dr. Liang Shi, Max Planck Institute for Dynamics and Self-Organization

3. [Dr. Markus Rampp](http://home.mpcdf.mpg.de/~mjr/), Max Planck Computing and Data Facility

4. Dr. Jose Manuel Lopez, IST Austria



Specific contribution is described below:

1. Prof. Marc Avila is responsible for the numerical method/formulation and supervises the development cycle;

2. Dr. Liang Shi designed and developed the earlier MPI parallelized versions;

3. Dr. Markus Rampp developed the hybrid parallelization, parallel I/O and visualization, HPC optimization, and manages the development cycle;

4. Dr. Jose Manuel Lopez implemented a dynamic time stepper and extended the code to include temperature field.


## Documentation and References

### For the methodology, implementation and validation, please refer to 

* Shi, L.; Rampp, M.; Hof, B.; Avila, M. A hybrid MPI-OpenMP parallel implementation for pseudospectral simulations with application to Taylor–Couette flow
[Computers & Fluids, 106, 1-11 (2015)](http://www.sciencedirect.com/science/article/pii/S0045793014003582)
, preprint: [arXiv:1311.2481](http://arxiv.org/abs/1311.2481)

Please note that the version of the code which was published in Computers & Fluids has meanwhile been improved by the implementation of a pressure-poisson equation (PPE) formulation and a variable timestep with automatic CFL control.


### User's guide and documentation of the source code

* a user's guide is distributed in the doc/ subdirectory

* technical documentation of the source code, auto-generated with the
  FORD tool can be found by following [this link](http://mjr.pages.mpcdf.de/nscouette/ford-doc)


### The following works are based on simulations with NSCouette. Please drop us a note if you publish scientific results obtained with the code.

* Conductive and convective heat transfer in fluid flows between differentially heated and rotating cylinders.
Jose M. Lopez, Francisco Marques and Marc Avila.
International Journal of Heat and Mass Transfer Volume 90, November 2015, Pages 959–967

* Nonlinear waves in stratified Taylor–Couette flow. Part 1. Layer formation.
Colin Leclercq, Jamie L. Partridge, Pierre Augier, Colm-Cille P. Caulfield, Stuart B. Dalziel and Paul F. Linden
arXiv:1609.02885

* Nonlinear waves in stratified Taylor--Couette flow. Part 2. Buoyancy flux
Colin Leclercq, Jamie L. Partridge, Colm-Cille P. Caulfield, Stuart B. Dalziel, Paul F. Linden
arXiv:1609.02886

* Hydrodynamic turbulence in quasi-Keplerian rotating flows
L. Shi, B. Hof, M. Rampp and M. Avila. Physics of Fluids, 29, 044107 (2017)
arXiv:1703.01714

## Contact

If you have any questions, comments or suggestions for improvements, please
contact [Prof. Marc Avila](mailto: marc.avila@zarm.uni-bremen.de) or [Dr. Markus Rampp](mailto: markus.rampp@mpcdf.mpg.de).
General support requests and bug reports can be sent to [nscouette@zarm.uni-bremen.de](mailto: nscouette@zarm.uni-bremen.de)
