Couette

Our DNS code nsCouette is a highly scalable software tool to solve the Navier-Stokes equations for incompressible fluid flow between differentially heated and independently rotating, concentric cylinders. The governing equations for the primitive variables are discretized in a cylindrical co-ordinate system using a Fourier-Galerkin ansatz for the azimuthal and the axial direction. High-order explicit finite differences are used in the only inhomogeneous (wall-normal) direction. Periodic boundary conditions are assumed in the axial direction. nsCouette is based on a pressure-poisson equation (PPE) formulation and a dynamic time-stepping with automatic CFL control. It is implemented in modern Fortran with a hybrid MPI-OpenMP parallelization scheme and thus designed to compute turbulent flows at high Reynolds and Rayleigh numbers. An additional GPU implementation (C-CUDA) for intermediate problem sizes and a basic version for turbulent pipe flow (nsPipe) are also provided.

## Build status and test coverage

[![Build status](https://gitlab.mpcdf.mpg.de/mjr/nscouette/badges/nsCouette-1.0/build.svg)](https://gitlab.mpcdf.mpg.de/mjr/nscouette/commits/nsCouette-1.0)
[![Code coverage](https://gitlab.mpcdf.mpg.de/mjr/nscouette/badges/nsCouette-1.0/coverage.svg)](http://mjr.pages.mpcdf.de/nscouette/) 

## Installation

The source code can be easily downloaded from our git repository. Please read the included user guide for further detailed instructions.

```
# Download the repository
git clone https://github.com/dfeldmann/nsCouette
cd nsCouette
ls doc
```

To build the executable, a Makefile for a standard x86_64 Linux software stack with Intel, GNU or PGI compilers is provided. In general, the Makefile itself requires no editing by the user. Instead, all platform-specific settings can be fixed using the provided architecture files. You can find a variety of working examples for different systems in the respective directory. Have a look, copy the template which appears closest to your platform/situation and modify it accordingly.


```
ls ARCH/make.arch.*

# Build the code (e.g. with Intel compilers and MKL)
make ARCH=intel-mkl

# Build the code (more build options)
make ARCH=myPlatform CODE=<STD_CODE|TE_CODE> HDF5IO=<NO|YES> DEBUG=<NO|YES> 

# Build the user guide (requires pdflatex)
make doc
```

Dependencies and software requirements:

*  Compiler: A modern Fortran compiler which is OpenMP-3 compliant and additionally a corresponding C-compiler (for houskeeping tasks). Tested with ifort/icc 12-15 or later, PSXE2017, PSXE2018, GCC 4.7-4.9, PGI 14. For free software see: https://gcc.gnu.org/

*  MPI: An MPI library which supports thread-level MPI_THREAD_SERIALIZED. Tested with Intel MPI 4.1 and higher, IBM PE 1.3 and higher. For free software see: http://www.open-mpi.de/

*  Linear algebra: A serial BLAS/LAPACK library. Tested with Intel MKL 11.x, OpenBLAS 0.2.14. Adapt $LIBLA and/or point environment variable $MKLROOT to an MKL installation. For free software see: http://www.openblas.net/ http://math-atlas.sourceforge.net/ http://www.netlib.org/

*  Fourier transforms: A serial but fully thread-safe FFTW3 installation or equivalent. Tested with FFTW 3.3 and MKL 11.1, earlier MKL versions will likely fail. Point environment variable $FFTW_HOME to a FFTW3 installation or select MKL (11.1 or later). For free software see: http://www.fftw.org 

*  Optional data output: An MPI-parallel HDF5 library installation. Tested with 1.8.9, 1.10.0-patch1. Point environment variable $HDF5_HOME to an suitable HDF5 installation or switch data output off by setting the external variable: make HDF5IO=NO. For free software see: http://www.hdfgroup.org/HDF5/


## Contributors


Following is a running list of contributors in chronological order:

1. [Prof. Marc Avila](https://www.zarm.uni-bremen.de/en/research/fluid-dynamics/fluid-simulation-and-modeling.html), University of Bremen, ZARM

2. Dr. Liang Shi, Max Planck Institute for Dynamics and Self-Organization

3. [Dr. Markus Rampp](http://home.mpcdf.mpg.de/~mjr/), Max Planck Computing and Data Facility

4. Dr. Jose Manuel Lopez, IST Austria

5. [Dr. Daniel Feldmann](https://www.zarm.uni-bremen.de/en/research/fluid-dynamics/fluid-simulation-and-modeling.html), University of Bremen, ZARM

6. Alberto Vela-Martin, School of Aeronautics, Universidad Politecnica de Madrid



Specific contribution is described below:

1. Prof. Marc Avila is responsible for the numerical method/formulation and supervises the development cycle

2. Dr. Liang Shi designed and developed the earlier MPI parallelized versions

3. Dr. Markus Rampp developed the hybrid parallelization, parallel I/O and visualization, HPC optimization, and manages the development cycle

4. Dr. Jose Manuel Lopez implemented a dynamic time-stepper and extended the code to include the temperature field

5. Dr. Daniel Feldmann is responsible for several additional features, bug-fixes, documentation and tutorials

6. Alberto Vela-Martin designed and developed the GPU-accelerated version



## Documentation and References


### User's guide and documentation of the source code

* A user guide is distributed with the source code (nscouette/doc/nsCouetteUserGuide.pdf)

* Technical documentation of the source code, auto-generated with the FORD tool can be found [here](http://mjr.pages.mpcdf.de/nscouette/ford-doc)


### References

* If you use nsCouette, please cite

  Jose M. Lopez, Daniel Feldmann, Markus Rampp, Alberto Vela-Martin, Liang Shi, Marc Avila, nsCouette - A high-performance code for direct numerical simulations of turbulent Taylor-Couette flow, preprint: []()

* For the methodology, implementation and validation, please refer to 
  
  Liang Shi, Markus Rampp, Bjoern Hof, Marc Avila, A hybrid MPI-OpenMP parallel implementation for pseudospectral simulations with application to Taylor-Couette flow
[Computers & Fluids, 106, 1-11 (2015)](http://www.sciencedirect.com/science/article/pii/S0045793014003582), preprint: [arXiv:1311.2481](http://arxiv.org/abs/1311.2481)

 Please note that the version of the code which was published in Computers & Fluids has meanwhile been improved by the implementation of a pressure-poisson equation (PPE) formulation and a variable timestep with automatic CFL control.



### The following works are based on simulations with nsCouette. Please drop us a note if you publish scientific results obtained with our code.

* [Phase-field simulation of core-annular pipe flow](https://doi.org/10.1016/j.ijmultiphaseflow.2019.04.027), B. Song, C. Plana, J. M. Lopez, M. Avila, International Journal of Multiphase Flow 117 (2019) 14-24.

* [Dynamics of viscoelastic pipe flow at low Reynolds numbers in the maximum drag reduction limit](https://doi.org/10.1017/jfm.2019.486), J. M. Lopez, G. H. Choueiri, B. Hof, Journal of Fluid Mechanics 874 (2019) 699-719.

* [Conductive and convective heat transfer in fluid flows between differentially heated and rotating cylinders](https://doi.org/10.1016/j.ijheatmasstransfer.2015.07.026), Jose M. Lopez, Francisco Marques and Marc Avila, International Journal of Heat and Mass Transfer, Volume 90, November (2015), 959-967

* [Hydrodynamic turbulence in quasi-Keplerian rotating flows](https://doi.org/10.1063/1.4981525), L. Shi, B. Hof, M. Rampp and M. Avila. Physics of Fluids, 29, 044107 (2017)

* [Nonlinear waves in stratified Taylor–Couette flow. Part 1. Layer formation](https://arxiv.org/abs/1609.02885), Colin Leclercq, Jamie L. Partridge, Pierre Augier, Colm-Cille P. Caulfield, Stuart B. Dalziel and Paul F. Linden

* [Nonlinear waves in stratified Taylor--Couette flow. Part 2. Buoyancy flux](https://arxiv.org/abs/1609.02886v1) Colin Leclercq, Jamie L. Partridge, Colm-Cille P. Caulfield, Stuart B. Dalziel, Paul F. Linden


## Contact

If you have any questions, comments or suggestions for improvements, please contact:
[Prof. Marc Avila](mailto: marc.avila@zarm.uni-bremen.de) or
[Dr. Markus Rampp](mailto: markus.rampp@mpcdf.mpg.de).

General support requests and bug reports can be sent to:
[nsCouette@zarm.uni-bremen.de](mailto: nscouette@zarm.uni-bremen.de)

