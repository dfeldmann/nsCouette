!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file is part of nsCouette -- A high-performance code for direct         !
! numerical simulations of turbulent Taylor-Couette flow                       !
!                                                                              !
! Copyright (C) 2019 Marc Avila, Bjoern Hof, Jose Manuel Lopez, Markus Rampp,  !
!                    Liang Shi, Alberto Vela-Martin, Daniel Feldmann.          !
!                                                                              !
! nsCouette is free software: you can redistribute it and/or modify it under   !
! the terms of the GNU General Public License as published by the Free         !
! Software Foundation, either version 3 of the License, or (at your option)    !
! any later version.                                                           !
!                                                                              !
! nsCouette is distributed in the hope that it will be useful, but WITHOUT ANY !
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    !
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more        !
! details.                                                                     !
!                                                                              !
! You should have received a copy of the GNU General Public License along with !
! nsCouette. If not, see <http://www.gnu.org/licenses/>.                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mod_params
!  Define all global parameters
!     th,z-directions : Fourier
!     r-direction     : Finite Difference
  
  IMPLICIT NONE
  SAVE

  !------------------------------------Mathematics constants
  REAL(KIND=8)   ,PARAMETER :: epsilon = 1D-10      ! numbers below it = 0
  REAL(KIND=8)   ,PARAMETER :: PI = ACOS(-1d0)    ! pi = 3.1415926...
  COMPLEX(KIND=8),PARAMETER :: ii = DCMPLX(0,1)    ! Complex i = sqrt(-1)

  !------------------------------------Spectral parametres
  INTEGER(KIND=4) :: m_r          
  INTEGER(KIND=4) :: m_th    
  INTEGER(KIND=4) :: m_z0    

  INTEGER(KIND=4) :: m_z 
  INTEGER(KIND=4) :: m_f  
  REAL(KIND=8)    :: k_th0
  REAL(KIND=8)    :: k_z0             

  !------------------------------------Physical parameters
  INTEGER(KIND=4) :: n_r            
  INTEGER(KIND=4) :: n_th 
  INTEGER(KIND=4) :: n_z  

  INTEGER(KIND=4) :: n_f  
  REAL(KIND=8)    :: len_r         
  REAL(KIND=8)    :: len_th 
  REAL(KIND=8)    :: len_z  
 
  REAL(KIND=8), allocatable :: r(:)!,th(:),z(:)
  
!---------------------------------------Miscellaneous

  REAL(kind=8),PARAMETER :: d_implicit= 0.5d0
  REAL(kind=8), PARAMETER   :: tolerance_dterr=5d-5
  
  !-----------------defaults for runtime parameters (can be set in input file)

  REAL(KIND=8)   :: Courant = 0.5d0
  INTEGER(KIND=4) :: print_time_screen = 250
  LOGICAL         :: variable_dt= .true.
  LOGICAL         :: const_flux= .true.
  REAL(kind=8) :: maxdt = 0.01d0
  
  !------------------------------------MPI & FFTW parameters

  INTEGER(KIND=4) :: mp_r  ! Radial points at each proc
  INTEGER(KIND=4) :: mp_f  ! Fourier points per proc
  INTEGER(KIND=4) :: mp_fmax  ! global max of Fourier points over all procs
  INTEGER(KIND=4),allocatable :: mp_f_arr(:) ! array of Fourier points over all procs    
  INTEGER(KIND=4),PARAMETER :: root = 0            ! Root processor 
  INTEGER(KIND=4),PARAMETER :: fftw_nthreads = 1
  LOGICAL        ,PARAMETER :: ifpad = .TRUE.      ! If apply '3/2' dealiasing
  
  !------------------------------------Finite Difference parameters 

  INTEGER(KIND=4),PARAMETER :: n_s = 9             ! Leading length of stencil

  !---------------defaults for runtime parameters (can be set in input file) 

  REAL(KIND=8)    :: init_dt = 1.0d-4 ! Time step
  INTEGER(KIND=4) :: dn_coeff = 100000      ! coeff per dn_coeff steps
  INTEGER(KIND=4) :: dn_ke = 10          ! energy per dn_ke steps
  INTEGER(KIND=4) :: dn_friction = 100          ! friction per dn_friction steps
  INTEGER(KIND=4) :: dn_hdf5 = 500000000       ! HDF5 output per dn_hdf5 steps
  INTEGER(KIND=4) :: numsteps = 10000    ! Number of steps 
  
  !----------------------------------- code revision identifier. automatically set by Makefile (PPFLAGS)
#ifndef GIT_REV
#define GIT_REV 'undef'
#endif
  CHARACTER(*), PARAMETER :: git_id = GIT_REV

  !----------------------------------- architecture used. automatically set by Makefile (PPFLAGS)
#ifndef ARCH_ID
#define ARCH_ID 'undef'
#endif
  CHARACTER(*), PARAMETER :: arch_id = ARCH_ID

  !----------------------------------- compiler flags used. automatically set by Makefile (PPFLAGS)
#ifndef CMP_OPTS
#define CMP_OPTS 'undef'
#endif
  CHARACTER(*), PARAMETER :: cmp_flgs = CMP_OPTS

contains
  
  subroutine init_grid
    
    implicit none

    integer :: ir

    m_z = 2*m_z0
    m_f  = (m_th+1)*m_z ! Number of fouriers modes

    n_r  = m_r          ! Number of grid points
    n_th = 2*m_th
    n_z  = m_z

    n_f  = n_th*n_z     ! points in fourier dirs
    len_r  = 1d0        ! Physical domain size
    len_th = 2*PI/k_th0
    len_z  = 2*PI/k_z0

    
    allocate(r(n_r))

  end subroutine init_grid
  


END MODULE mod_params
