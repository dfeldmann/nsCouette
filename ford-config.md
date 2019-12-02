project: nsCouette
src_dir: ./
include: XXXincludeXXX
exclude_dir: postproc
exclude: perf2ftimings.f90
exclude: perfdummy.f90
exclude: ftimings.f90
exclude: ftimings_type.f90
exclude: ftimings_value.f90
exclude: mod_preAnalysis.f90
exclude: waveSpeed.f90
output_dir: ./ford-doc
project_website: https://gitlab.mpcdf.mpg.de/mjr/nscouette
summary: nsCouette -- A high-performance code for direct numerical simulations of turbulent Taylor-Couette flow
author: Jose Manuel Lopez, Daniel Feldmann, Markus Rampp, Alberto Vela-Martin, Liang Shi, Marc Avila
author_description: MPG Germany, IST Austria, University Bremen, Germany
email: nsCouette@zarm.uni-bremen.de
preprocess: false
predocmark: >
media_dir: ./media
docmark_alt: #
predocmark_alt: <
display: public
         protected
         private
source: true
graph: true
search: true
license: gfdl

Our DNS code nsCouette is a highly scalable software tool to solve the
Navier-Stokes equations for incompressible fluid flow between differentially
heated and independently rotating, concentric cylinders. The governing
equations for the primitive variables are discretized in a cylindrical
co-ordinate system using a Fourier-Galerkin ansatz for the azimuthal and the
axial direction. High-order explicit finite differences are used in the only
inhomogeneous (wall-normal) direction. Periodic boundary conditions are assumed
in the axial direction. nsCouette is based on a pseudospectral spatial
discretization and dynamic time-stepping. It is implemented in modern Fortran
with a hybrid MPI-OpenMP parallelization scheme and thus designed to compute
turbulent flows at high Reynolds and Rayleigh numbers. An additional GPU
implementation (C-CUDA) for intermediate problem sizes and a basic version for
turbulent pipe flow (nsPipe) are also provided.
