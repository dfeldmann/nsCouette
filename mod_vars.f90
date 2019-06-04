!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file is part of NSPipeflow, a HPC code for DNS of pipe flow          !
!                                                                           !
! Copyright (C) 2016 Marc Avila, Bjoern Hof, Jose Manuel Lopez,             !
!                    Markus Rampp, Liang Shi                                !
!                                                                           !
! NSPipeflow is free software: you can redistribute it and/or modify         !
! it under the terms of the GNU General Public License as published by      !
! the Free Software Foundation, either version 3 of the License, or         !
! (at your option) any later version.                                       !
!                                                                           !
! NSPipeflow is distributed in the hope that it will be useful,              !
! but WITHOUT ANY WARRANTY; without even the implied warranty of            !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             !
! GNU General Public License for more details.                              !
!                                                                           !
! You should have received a copy of the GNU General Public License         !
! along with NSPipeflow.  If not, see <http://www.gnu.org/licenses/>.        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=============================================
!
!     define global variables & types
!
!=============================================                                    


MODULE mod_vars
  
  USE mod_fftw
  USE mod_params
  USE mpi
  IMPLICIT NONE
  SAVE
  private :: alloc_vec_mpi,dealloc_vec_mpi,alloc_vec_r2d,dealloc_vec_r2d
  
  !------------------------------------Derived Types

  TYPE vspec_ptr
     ! can be used to finally replace all AoS, e.g. f_hat_mp(spec)
     SEQUENCE
     COMPLEX(KIND=8), pointer :: ur(:,:)
     COMPLEX(KIND=8), pointer :: uth(:,:)
     COMPLEX(KIND=8), pointer :: uz(:,:)
     COMPLEX(KIND=8), pointer :: p
     REAL(KIND=8), pointer :: k_th(:,:)
     REAL(KIND=8), pointer :: k_z(:,:)
  END TYPE vspec_ptr

  TYPE vec_mpi
     COMPLEX(KIND=8),allocatable :: r(:,:)
     COMPLEX(KIND=8),allocatable :: th(:,:)
     COMPLEX(KIND=8),allocatable :: z(:,:)
   contains
     procedure :: alloc   => alloc_vec_mpi
     procedure :: dealloc => dealloc_vec_mpi
  END TYPE vec_mpi
  TYPE vec_r2d
     REAL(KIND=8),allocatable :: th(:,:)
     REAL(KIND=8),allocatable :: z(:,:)
   contains
     procedure :: alloc   => alloc_vec_r2d
     procedure :: dealloc => dealloc_vec_r2d
  END TYPE vec_r2d

  !-------------------------------------Physical & spectral variables

  TYPE (vec_mpi), target  :: u_hat_mp
  COMPLEX(KIND=8),allocatable,DIMENSION(:,:), target :: p_hat_mp
  TYPE (vec_r2d) :: fk_mp
  TYPE (vec_mpi)  :: udu_hat_mp
  REAL(kind=8) :: vel_Pr0

  !------------------------------------MPI variables

  INTEGER(KIND=4) :: ierr, comm, np, myid,counter,numprocs
  INTEGER(KIND=4) :: double_complex,mpi_spec,filetype,filetype2 
  INTEGER(KIND=4) :: tl_mpi           !available level of threaded MPI
  INTEGER(KIND=4) :: nomp = 1         !number of OpenMP threads
  INTEGER(KIND=4) :: k0_modes
  
  !------------------------------------Miscellanea
  
  REAL(KIND=8), allocatable :: bMw_dr(:,:),bMw_drr(:,:)
  REAL(KIND=8), allocatable :: dr1(:,:),dr0(:,:) !To extrapolate to ri and ro
  REAL(KIND=8), allocatable :: intrdr(:)
  INTEGER(KIND=4),allocatable :: Sbig_map(:,:)!Parity of nodes distributed among processor
  INTEGER(KIND=4) :: i_time=1       ! Counter of the time steps
  REAL(KIND=8)    :: time=0d0         ! time (R/uc units)
  INTEGER(KIND=4) :: runtime=86400! maximum runtime [s] for the job
  INTEGER(kind=4)    :: iter    ! number of iterations in Crank-Nicholson
  REAL(kind=8)   ::  cfl !Courant-Friedich-Lewis condition
  INTEGER(kind=4) :: cfl_dir!critical direction for CFL (1 radial 2 azimuthal 3 axial)
  REAL(kind=8) :: dt  !time_step
  LOGICAL :: new_dt
  
  !------------------------------------Variable parameters

  REAL(KIND=8)    :: Re
  INTEGER(KIND=4) :: restart=0    ! initialization mode
  CHARACTER(40)   :: fName_ic     ! coeff file name
  CHARACTER(40)   :: fBase_ic     ! base prefix for coeff file name
  INTEGER(KIND=4) :: i_start=1    ! Starting output index of coeff 
  REAL(kind=8) :: dt_ifile,time_ifile
   
  !------------------------------------Variable pars of initial velocity file

  REAL(KIND=8)    :: Re_ifile
  INTEGER(KIND=4) :: m_r_ifile,m_th_ifile,m_z_ifile
  INTEGER(KIND=4) :: numprocs_ifile
  
  INTEGER(kind=4) ::  m_f_ifile,mp_f_ifile,mp_fmax_ifile
  logical :: extended_mf_grid=.false.
  
contains
  
  function alloc_vec_mpi(this,m,n) result(ierr)

    class(vec_mpi) :: this
    integer, intent(in) :: m,n
    integer ierr

    allocate (this%r(m,n),  STAT=ierr)
    allocate (this%th(m,n), STAT=ierr)
    allocate (this%z(m,n),  STAT=ierr)

  end function alloc_vec_mpi


  function dealloc_vec_mpi(this) result(ierr)
    class(vec_mpi) :: this
    integer ierr

    deallocate (this%r, STAT=ierr)
    deallocate (this%th,STAT=ierr)
    deallocate (this%z, STAT=ierr)

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
