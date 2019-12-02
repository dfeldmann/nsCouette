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

module mod_vars
! define global variables & types

use mod_fftw
use mod_params
use mpi
implicit none
save

private :: alloc_vec_mpi,dealloc_vec_mpi,alloc_vec_r2d,dealloc_vec_r2d

  !------------------------------------Derived Types

!!$  TYPE vspec_ptr
!!$     ! can be used to finally replace all AoS, e.g. f_hat_mp(spec)
!!$     SEQUENCE
!!$     COMPLEX(KIND=8), pointer :: ur(:,:)
!!$     COMPLEX(KIND=8), pointer :: uth(:,:)
!!$     COMPLEX(KIND=8), pointer :: uz(:,:)
!!$#ifdef TE_CODE
!!$     COMPLEX(KIND=8), pointer :: T(:,:)
!!$#endif /* TE_CODE */
!!$     COMPLEX(KIND=8), pointer :: p
!!$     REAL(KIND=8), pointer :: k_th(:,:)
!!$     REAL(KIND=8), pointer :: k_z(:,:)
!!$  END TYPE vspec_ptr


  TYPE vec_mpi
!     SEQUENCE
     COMPLEX(KIND=8), allocatable :: r(:,:)
     COMPLEX(KIND=8), allocatable :: th(:,:)
     COMPLEX(KIND=8), allocatable :: z(:,:)
#ifdef TE_CODE
     COMPLEX(KIND=8), allocatable  :: T(:,:)
#endif /* TE_CODE */
   contains
     procedure :: alloc   => alloc_vec_mpi
     procedure :: dealloc => dealloc_vec_mpi
  END TYPE vec_mpi

  TYPE vec_r2d
!     SEQUENCE
     REAL(KIND=8), allocatable :: th(:,:)
     REAL(KIND=8), allocatable :: z(:,:)
   contains
     procedure :: alloc   => alloc_vec_r2d
     procedure :: dealloc => dealloc_vec_r2d
  END TYPE vec_r2d

  !-------------------------------------Physical & spectral variables

  TYPE (vec_mpi), target  :: u_hat_mp
  COMPLEX(KIND=8),allocatable, target :: p_hat_mp(:,:)
  TYPE (vec_r2d) :: fk_mp
  TYPE (vec_mpi)  :: udu_hat_mp

  !------------------------------------MPI variables
  INTEGER(KIND=4) :: comm, numprocs, myid,counter
  INTEGER(KIND=4) :: double_complex,mpi_spec,filetype,filetype2
  INTEGER(KIND=4) :: tl_mpi           !available level of threaded MPI
  INTEGER(KIND=4) :: nomp = 1         !number of OpenMP threads

  !------------------------------------Miscellanea
  REAL(KIND=8), allocatable :: bMw_dr(:,:),bMw_drr(:,:)
  REAL(KIND=8), allocatable :: dr1(:,:),dr2(:,:) !To extrapolate to ri and ro                                                
  REAL(KIND=8), allocatable :: intrdr(:)
  REAL(KIND=8), allocatable :: row_inv_Mw_dr(:)

  INTEGER(KIND=4) :: i_time=1     ! Counter of the time steps
  REAL(KIND=8)    :: time=0.d0    ! time in viscous unit
  INTEGER(kind=4) :: iter         ! Crank-Nicholson iterations counter
  REAL(Kind=8)    :: dt           ! timestep
  LOGICAL :: new_dt
  REAL(kind=8)   ::  cfl !Courant-Friedich-Lewis condition
  INTEGER(kind=4) :: cfl_dir!critical direction for CFL (1 radial 2 azimuthal 3 axial)

  ! Non-dimensional control parameters
  real(kind=8) :: re_i, re_o    ! inner/outer cylinder Reynolds number
  real(kind=8) :: gr, pr, eps   ! Grashof and Prandtle numbers, TE_CODE only

  integer(kind=4) :: restart=0    ! initialization mode
  integer(kind=4) :: runtime=86400! maximum runtime [s] for the job
  character(40)   :: fname_ic     ! coeff file name
  character(40)   :: fbase_ic     ! base prefix for coeff file name
  integer(kind=4) :: i_start=1    ! Starting output index of coeff

  ! Set initial conditions (restart = 0)
  logical :: ic_tcbf = .true.                   ! Set Taylor-Couette base flow (T) or resting fluid (F), only when restart = 0
  logical :: ic_temp = .false.                  ! Set temperature profile (T) or zero (F), only when restart = 0, only TE_CODE
  logical :: ic_pert = .true.                   ! Add perturbation on top of base flow (T) or not (F), only when restart = 0
  integer(kind=4), parameter :: ic_np = 6       ! Number of initial perturbations
  real(kind=8), dimension(ic_np, 3) :: ic_p = 0 ! lth perturbation: amplitude and wavevector (l, (al, k_thl, k_zl))

  ! Variables of initial condition read from file
  real(kind=8)    :: re_i_ifile, re_o_ifile           ! Reynolds numbers
  real(kind=8)    :: alpha_ifile                      ! radial grid point distribution
  integer(kind=4) :: m_r_ifile, m_th_ifile, m_z_ifile ! number of points and modes
#ifdef TE_CODE
  REAL(KIND=8)    :: Gr_ifile,Pr_ifile
#endif /* TE_CODE */
  INTEGER(KIND=4) :: dn_coeff_ifile
  INTEGER(KIND=4) :: numprocs_ifile
#ifdef TEST2
  REAL(KIND=8)    :: temp_prof(m_r)
#endif /* TEST2 */

  INTEGER(kind=4) ::  m_f_ifile, mp_f_ifile, mp_fmax_ifile
  
  logical :: extended_mf_grid=.false.
 
contains 
  
  function alloc_vec_mpi(this,m,n) result(ierr)
    
    class(vec_mpi) :: this 
    integer, intent(in) :: m,n
    integer ierr

    allocate (this%r(m,n),  STAT=ierr)
    allocate (this%th(m,n), STAT=ierr)
    allocate (this%z(m,n),  STAT=ierr)
#ifdef TE_CODE
    allocate (this%T(m,n),  STAT=ierr)
#endif /* TE_CODE */
    
  end function alloc_vec_mpi


  function dealloc_vec_mpi(this) result(ierr)
    class(vec_mpi) :: this 
    integer ierr

    deallocate (this%r, STAT=ierr)
    deallocate (this%th,STAT=ierr)
    deallocate (this%z, STAT=ierr)
#ifdef TE_CODE
    deallocate (this%T, STAT=ierr)
#endif /* TE_CODE */
    
  end function dealloc_vec_mpi

  function alloc_vec_r2d(this,m,n) result(ierr)

    class(vec_r2d) :: this 
    integer, intent(in) :: m,n
    integer ierr

    allocate (this%th(m,n),STAT=ierr)
    allocate (this%z(m,n), STAT=ierr)

  end function alloc_vec_r2d


  function dealloc_vec_r2d(this) result(ierr)
    class(vec_r2d) :: this 
    integer ierr

    deallocate (this%th,STAT=ierr)
    deallocate (this%z, STAT=ierr)

  end function dealloc_vec_r2d

END MODULE mod_vars
