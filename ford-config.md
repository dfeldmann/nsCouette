project: NSCouette
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
summary: NSCouette, a high-performance code for the direct numerical simulation of Taylor-Couette flow
author: Liang Shi, Jose-Manuel Lopez, Markus Rampp, Marc Avila
author_description: MPG Germany, IST Austria, Univ. Bremen, Germany
email: nscouette@zarm.uni-bremen.de
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

The code NSCOUETTE is designed to simulate incompressible
Taylor-Couette flows with infinitely long cylinders (periodic in the axial 
direction) on high-performance computers. 
The code is written in modern Fortran and uses a hybrid MPI-OpenMP parallelization
strategy.
