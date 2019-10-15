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

module mod_inout
! Input & Output
!   Input: c.f. input file 'input_nsCouette'
!   Output:
!     1) spectral coefficients for velocity and pressure
!     2) modal energy and total energy
!     3) temporal evolution of velocity at middle point in r-dir
!     4) temporal evolution of torque
!     5) global parameters to file 'parameter.dat'
!     6) temporal evolution of velocity and temperature at six user
!        defined radial probe locations

use mod_mympi
use mod_fftw
use mod_vars

implicit none

  TYPE spec
     SEQUENCE
     COMPLEX(KIND=8) :: ur
     COMPLEX(KIND=8) :: uth
     COMPLEX(KIND=8) :: uz
     COMPLEX(KIND=8) :: p
#ifdef TE_CODE
     COMPLEX(KIND=8) :: T
#endif /* TE_CODE */
     REAL(KIND=8)    :: k_th
     REAL(KIND=8)    :: k_z
  END TYPE spec
  
  TYPE(spec),ALLOCATABLE,DIMENSION(:,:) :: aux_hat_mp
  logical, private :: scale = .false. ! remap, scale, interpolate inital condition to specified grid if true

  namelist /parameters_grid/ m_r,m_th,m_z0,k_th0,k_z0,eta,alpha
  namelist /parameters_physics/ re_i,re_o,gr,pr,gap,gra,nu
  namelist /parameters_timestep/ numsteps,init_dt,variable_dt,maxdt,courant
  namelist /parameters_output/ dn_coeff,dn_ke,dn_vel,dn_nu,dn_hdf5,print_time_screen,fbase_ic,dn_prbs,prl_r,prl_th,prl_z
  namelist /parameters_control/ restart,runtime
  namelist /parameters_restart/ i_time, time,dt, m_r_ifile, m_th_ifile, m_z_ifile, fbase_ic
  namelist /parameters_initialcondition/ ic_tcbf, ic_temp, ic_pert, ic_p

  private :: spec,aux_hat_mp
  private :: parameters_grid,parameters_physics,parameters_timestep,parameters_output,&
       &     parameters_control,parameters_restart

  integer(kind=4), dimension(prl_n) :: prl_rnk  ! rank of each probe location
  integer(kind=4), dimension(prl_n) :: prl_io   ! i/o unit for each probe location
  integer(kind=8), dimension(prl_n) :: prl_nr   ! radial index of each probe, global
  integer(kind=4), dimension(prl_n) :: prl_nrp  ! radial index, local/process
  integer(kind=8), dimension(prl_n) :: prl_nth  ! azimtuhal index of each probe location
  integer(kind=8), dimension(prl_n) :: prl_nz   ! axial index of each probe location

contains



subroutine read_pars()
! Purpose: Read set of parameters using predefined namelists
! Called from: nsCouette.f90 
! Content:
! + Read predefined namelists from std in, only on root
! + Compute some extra parameters, only on root
! + Broadcast final set of parameters

implicit none

if (myid .eq. root) then

 ! Read namelists from std in < nsCouette.in
 read(*, nml=parameters_grid)
 read(*, nml=parameters_physics)
 read(*, nml=parameters_timestep)
 read(*, nml=parameters_output)
 read(*, nml=parameters_control)
 read(*, nml=parameters_initialcondition)

 ! Set some extra parameters
 dt = init_dt
#ifdef TE_CODE
 eps = nu**2.0d0*gr/(gra*gap**3.0d0) ! Relative density variation
#else
 gr  = 0.0d0 ! Not used in the non-TE_CODE
 pr  = 0.0d0
 eps = 0.0d0
#endif /* TE_CODE */
#ifdef TEST2
 ! Taylor number as input parameter (turn into inner cylinder Reynolds number)
 re_i = sqrt((1.0d0-eta**2.0d0)*re_i/(2.0d0*(1.0d0-eta)**2.0d0))
#endif /* TEST2 */

end if

! Broadcast final set of parameters
call mpi_bcast(m_r, 1,mpi_integer,root,comm,ierr)
call mpi_bcast(m_th,1,mpi_integer,root,comm,ierr)
call mpi_bcast(m_z0,1,mpi_integer,root,comm,ierr)
call mpi_bcast(k_th0,1,mpi_real8,root,comm,ierr)
call mpi_bcast(k_z0,1,mpi_real8,root,comm,ierr)
call mpi_bcast(eta,1,mpi_real8,root,comm,ierr)
call mpi_bcast(alpha,1,mpi_real8,root,comm,ierr)
call mpi_bcast(re_i,1,mpi_real8,root,comm,ierr)
call mpi_bcast(re_o,1,mpi_real8,root,comm,ierr)
call mpi_bcast(restart,1,mpi_integer,root,comm,ierr)
call mpi_bcast(runtime,1,mpi_integer,root,comm,ierr)
call mpi_bcast(fbase_ic,len(fbase_ic),mpi_character,root,comm,ierr)
call mpi_bcast(init_dt,1,mpi_real8,root,comm,ierr)
call mpi_bcast(dt,1,mpi_real8,root,comm,ierr)
call mpi_bcast(numsteps,1,mpi_integer,root,comm,ierr)
call mpi_bcast(dn_coeff,1,mpi_integer,root,comm,ierr)
call mpi_bcast(dn_ke,1,mpi_integer,root,comm,ierr)
call mpi_bcast(dn_vel,1,mpi_integer,root,comm,ierr)
call mpi_bcast(dn_nu,1,mpi_integer,root,comm,ierr)
call mpi_bcast(dn_hdf5,1,mpi_integer,root,comm,ierr)
call mpi_bcast(dn_prbs,1,mpi_integer,root,comm,ierr)
call mpi_bcast(print_time_screen,1,mpi_integer,root,comm,ierr)
call mpi_bcast(variable_dt,1,mpi_logical,root,comm,ierr)
call mpi_bcast(courant,1,mpi_real8,root,comm,ierr)
call mpi_bcast(maxdt,1,mpi_real8,root,comm,ierr)
call mpi_bcast(ic_tcbf, 1, mpi_logical, root, comm, ierr)
call mpi_bcast(ic_temp, 1, mpi_logical, root, comm, ierr)
call mpi_bcast(ic_pert, 1, mpi_logical, root, comm, ierr)
call mpi_bcast(ic_p, ic_np*3, mpi_real8, root, comm, ierr)
#ifdef TE_CODE
call mpi_bcast(gr,1,mpi_real8,root,comm,ierr)
call mpi_bcast(pr,1,mpi_real8,root,comm,ierr)
call mpi_bcast(eps,1,mpi_real8,root,comm,ierr)
#endif /* TE_CODE */

end subroutine read_pars



subroutine output_pars()
! Purpose: Write final set of parameters to file
! Called from: nsCouette.f90
! + Open file, write namelists, close, only on root
 
implicit none

if (myid .eq. root) then
 open(unit=106, file='parameters.out')
 write(106, nml=parameters_grid)
 write(106, nml=parameters_physics)
 write(106, nml=parameters_timestep)
 write(106, nml=parameters_output)
 write(106, nml=parameters_control)
 write(106, nml=parameters_initialcondition)
 close(106)
end if

end subroutine output_pars



subroutine read_restart_pars()
! Purpose: Read the set of parameters for the initial velocity file when
! restarting from a checkpoint coeff file.
! Called from: nsCouette.f90
! Content:
! + Read one namelist from file restart.in
! + Overwrite restart.in to prevent loops of restart attempts in case of errors
! + Reset dt to user specified value, only in case of constant timestep size
! + Broadcast the parameters of this namelist
! + Construct file name to read initial conditions from
! + Update time and timestep according to restart mode
! + Check consistency of specified and read grid resolution -> remapping
! + Allocate new arrays in case of different grids

implicit none
integer(kind=4) :: ndims=2
integer(kind=4) :: sizes(2), subsizes(2), starts(2)
character(8) ::  suffix

if (myid .eq. root) then

 ! Read one namelist from file and erase afterwards
 open(unit=107, file='restart')
 read(107, nml=parameters_restart)
 write(107, *) ' ' ! write empty file
 close(107)

 ! Reset to specified value in case of non-variable timestep size
 if (variable_dt .eqv. .false.) then 
  dt = init_dt
  write(*,'(a, es23.15e3)') 'Constant user-specified timestep size set to dt = init_dt =', dt
 else
  write(*,'(a, es23.15e3, a)') 'Variable timestep size with initial dt =', dt, ' read from restart file'
 end if

end if

! Broadcast the parameters of this namelist
call mpi_bcast(i_time,1,mpi_integer,root,comm,ierr)
call mpi_bcast(time,1,mpi_real8,root,comm,ierr)
call mpi_bcast(dt,1,mpi_real8,root,comm,ierr)
call mpi_bcast(m_r_ifile,1,mpi_integer,root,comm,ierr)
call mpi_bcast(m_th_ifile,1,mpi_integer,root,comm,ierr)
call mpi_bcast(m_z_ifile,1,mpi_integer,root,comm,ierr)
call mpi_bcast(fbase_ic,len(fbase_ic),mpi_character,root,comm,ierr)

! Construct file name to read initial conditions from
write(suffix, '(i0.8)') i_time
fname_ic = 'coeff_'//trim(fbase_ic)//'.'//suffix

! Update time and time step according to restart mode
i_start = i_time+1
if (restart .eq. 2) then
 time = 0.0d0
 i_start = 1
end if

! Check consistency of specified and read grid resolution
if ((m_r .ne. m_r_ifile) .or. (m_z .ne. m_z_ifile) .or. (m_th .ne. m_th_ifile)) then

 scale = .true. ! Detected inconsistent resolution, report here, remap later
 if (myid .eq. root) then
  write(*, '(a)') 'Resolution of the restart file differs from the specified resolution!'
  write(*, '(a)') 'Remap the initial condition (parameters_restart -> parameters_grid):'
  write(*, '(a5, i5.5, a4, i5.5)') 'm_r  ', m_r_ifile,   ' -> ', m_r
  write(*, '(a5, i5.5, a4, i5.5)') 'm_th ', m_th_ifile,  ' -> ', m_th
  write(*, '(a5, i5.5, a4, i5.5)') 'm_z0 ', m_z_ifile/2, ' -> ', m_z0
 end if

 m_f_ifile = m_z_ifile*(m_th_ifile+1)
 mp_fmax_ifile = ceiling(real(m_f_ifile,kind=8)/real(numprocs,kind=8)) ! Fourier points per proc 

 if ((MOD(m_f_ifile,numprocs) /= 0)) then
  mp_f_ifile=(numprocs+mp_fmax_ifile*numprocs-1)/numprocs
  mp_f_ifile=max(0,min(mp_f_ifile,m_f_ifile-myid*mp_f_ifile))
  m_f_ifile=mp_fmax_ifile*numprocs
 else
  mp_f_ifile=mp_fmax_ifile
 end if
       
 sizes = (/m_r_ifile,m_f_ifile/)
 subsizes = (/m_r_ifile,mp_fmax_ifile/)
 starts = (/0,myid*mp_fmax_ifile/)
 CALL MPI_Type_create_subarray(ndims,sizes,subsizes,starts,&
    MPI_ORDER_FORTRAN,mpi_spec,filetype2,ierr)
 CALL MPI_Type_commit(filetype2,ierr)
 
 ALLOCATE(aux_hat_mp(m_r_ifile,mp_f_ifile))

end if

end subroutine read_restart_pars



subroutine init_probes()
! Purpose:
! Initialise time series output at several probe locations (in general, each
! one associated with its individual rank/task). For now, six different probe
! locations are supported. Six is hard coded, but individual locations are
! specified by user during runtime (via input_nsCouette).
! Called from:
! + nsCouette.f90 
! Contet:
! + Fix i/o units for all probe locations
! + Correct user defined relative probe locations (0<r,th,z<1) with resepct to
!   inner radius, and fundamental wave length in azimuthal and axial direction
! + Find global indices best matching user defined locations, root only
! + Find associated ranks, root only
! + Convert to task/process-local radial indices
! + Broadcast information
! + Open one file per probe location, so that each rank hosting one or more
!   probes outputs to its own file(s). I guess this is better than collecting
!   all probe data from different ranks via mpi_sent/receive/bcast and than
!   do all the output from one rank to one file, right?
! + Write header information to file in case of REWIND

implicit none
integer(kind=4) :: i ! running index over all probes
integer(kind=8) :: n ! running index over all radial locations or tasks
real(kind=8) :: d, dd ! dummy distances 
character(len=11) :: prl_fnam ! probe file name
character(len=6) :: char_position ! pointer position at i/o
integer(kind=8) :: nthfine ! number of azimtuhal points on fine grid
integer(kind=8) :: nzfine ! number of axial points on fine grid

nthfine = 3*n_th/2 ! better put all fine-grid variables to mod_vars.f90
nzfine = 3*n_z/2

if (myid .eq. root) then

 ! Fix i/o units (grep -i "unit" *.f90)
 prl_io(1) = 91
 prl_io(2) = 92
 prl_io(3) = 93
 prl_io(4) = 94
 prl_io(5) = 95
 prl_io(6) = 96

 ! Correct user defined probe locations
 prl_r = prl_r + r_i
 prl_th = prl_th * len_th
 prl_z = prl_z * len_z

 do i = 1, prl_n ! over all probes

  ! Find global indices
  d = 1.0d99 ! huge inital distance
  do n = 1, n_r ! over all r from inner to outer cylinder, n_r = m_r
   dd = abs(r(n) - prl_r(i))
   if (dd .lt. d) then
    d = dd ! update distance
    prl_nr(i) = n ! update radial index
   end if
  end do
  prl_nth(i) = int(prl_th(i)/len_th*(nthfine-1))+1 ! azimuthal index
  prl_nz(i) = int(prl_z(i)/len_z*(nzfine-1))+1 ! axial index

  ! Find associated ranks
  do n = 0, numprocs-1
   if ((r(prl_nr(i)) .ge. r(n*mp_r+1)) .and. &
    &  (r(prl_nr(i)) .le. r((n+1)*mp_r))) then
    prl_rnk(i) = n ! update rank hosting this probe location
    prl_nrp(i) = prl_nr(i) - n*mp_r ! convert to local index
   end if
  end do

 end do
end if

! Broadcast information
call mpi_bcast(prl_io, prl_n, mpi_integer4, root, comm, ierr)
call mpi_bcast(prl_r, prl_n, mpi_real8, root, comm, ierr)
call mpi_bcast(prl_th, prl_n, mpi_real8, root, comm, ierr)
call mpi_bcast(prl_z, prl_n, mpi_real8, root, comm, ierr)
call mpi_bcast(prl_nr, prl_n, mpi_integer8, root, comm, ierr)
call mpi_bcast(prl_nrp, prl_n, mpi_integer4, root, comm, ierr)
call mpi_bcast(prl_nth, prl_n, mpi_integer8, root, comm, ierr)
call mpi_bcast(prl_nz, prl_n, mpi_integer8, root, comm, ierr)
call mpi_bcast(prl_rnk, prl_n, mpi_integer4, root, comm, ierr)

do i = 1, prl_n ! over all probes
 if (myid .eq. prl_rnk(i)) then

  ! Open files
  write(prl_fnam, 1001) i ! Set probe file name
  if ((restart .eq. 0) .or. (restart .eq. 2)) then ! Access mode
   char_position = 'REWIND'
   open(unit=prl_io(i), file=prl_fnam, position=char_position)

   ! Write header
#ifndef TE_CODE
   write(prl_io(i), 1002) i, prl_rnk(i), prl_nr(i), prl_nrp(i), r(prl_nr(i)), &
   & prl_nth(i), prl_th(i), prl_nz(i), prl_z(i), ' '
#else
   write(prl_io(i), 1002) i, prl_rnk(i), prl_nr(i), prl_nrp(i), r(prl_nr(i)), &
   & prl_nth(i), prl_th(i), prl_nz(i), prl_z(i), ', t'
#endif

  else
   char_position = 'APPEND'
   open(unit=prl_io(i), file=prl_fnam, position=char_position)
  end if
 end if
end do

! Formats
1001 format('probe', i2.2, '.dat')
1002 format('# Time series data from probe ', i2.2, ' on rank ', i5.5, /, &
          & '# Radial location n_r = ', i5.5, ' n_rp = ' , i4.4, ' r =', es23.15, /, &
          & '# Azimuthal location n_th = ', i5.5, ' th =', es23.15, /, &
          & '# Axial location n_z = ', i5.5, ' z =', es23.15, /, &
          & '# Time, u_r, u_th, u_z', a)

end subroutine init_probes



subroutine final_probes()
! Purpose:
! Finalise time series output at several probes, see also init_probes()
! Called from:
! + mod_inOut.f90
! Content:
! + Close all files which were opened in init_probes()

implicit none
integer(kind=4) :: i ! running index over all probes

! Close all files
do i = 1, prl_n
 if (myid .eq. prl_rnk(i)) close(unit = prl_io(i))
end do

end subroutine final_probes



subroutine write_probes(u_r, u_th, u_z, t)
! Purpose:
! Write time series data from several probes to individual files. One file per
! probe, each one usually located on individual rank, see also init_probes().
! Flow field data in physical space is readly available from mod_nonlinear.f90
! Called from:
! + mod_nonlinear.f90
! Content:
! + Convert theta and z to single running index
! + Write to file

implicit none
integer(kind=4) :: i ! running index over all probes
integer(kind=8) :: j ! combined theta-z running index
real(kind=8), dimension(9*n_f/4, mp_r), intent(in) :: u_r ! radial velocity
real(kind=8), dimension(9*n_f/4, mp_r), intent(in) :: u_th ! azimuthal velocity
real(kind=8), dimension(9*n_f/4, mp_r), intent(in) :: u_z ! axial velocity
real(kind=8), dimension(9*n_f/4, mp_r), intent(in), optional :: t ! temperature
integer(kind=8) :: nthfine ! number of azimtuhal points on fine grid
integer(kind=8) :: nzfine ! number of axial points on fine grid

nthfine = int(3*n_th/2) ! better put all fine-grid variables to mod_vars.f90?
nzfine = int(3*n_z/2)

do i = 1, prl_n
 if (myid .eq. prl_rnk(i)) then

  ! Convert index
  j = prl_nth(i) + nthfine * (prl_nz(i) - 1)

  ! Write to file
#ifndef TE_CODE
  write(prl_io(i), 1001) time, u_r(j, prl_nrp(i)), u_th(j, prl_nrp(i)), &
  & u_z(j, prl_nrp(i))
#else
  write(prl_io(i), 1002) time, u_r(j, prl_nrp(i)), u_th(j, prl_nrp(i)), &
  & u_z(j, prl_nrp(i)), t(j, prl_nrp(i))
#endif

 end if
end do

! Formats
1001 format(4es25.13e3)
1002 format(5es25.13e3)

end subroutine write_probes



  !--------------------------------------Open files to write into
  SUBROUTINE open_files()
    IMPLICIT NONE
    
    CHARACTER(6) :: char_position   ! pointer position at I/O

    !----------------------Access mode


    IF ((restart == 0) .or. (restart == 2)) THEN
       char_position = 'REWIND'
    ELSE 
       char_position = 'APPEND'
    END IF
    
    IF (myid == root) THEN
       OPEN(unit=108,file='ke_mode',position=char_position)
       OPEN(unit=109,file='ke_th',position=char_position)
       OPEN(unit=103,file='ke_z',position=char_position)
       OPEN(unit=104,file='ke_total',position=char_position)
       OPEN(unit=105,file='torque',position=char_position)
#ifdef TE_CODE
       OPEN(unit=123,file='Temp_energy',position=char_position)
       OPEN(unit=124,file='Nusselt',position=char_position)
#endif /* TE_CODE */
    END IF

    call hdf5init

  END SUBROUTINE open_files

  !--------------------------------------Close files opened above
  SUBROUTINE close_files()
    IMPLICIT NONE
    
    IF (myid == root) THEN
       CLOSE(108)
       CLOSE(109)
       CLOSE(103)
       CLOSE(104)
       CLOSE(105)
#ifdef TE_CODE
       CLOSE(123)
       CLOSE(124)
#endif /* TE_CODE */
    end if

    call hdf5finalize

  END SUBROUTINE close_files

  subroutine write_restart_file
    
    if (myid == root) then
       !(over)write minimal info required for restart from latest coeff checkpoint
       open(unit=117,file='restart',delim='QUOTE')
       write(117,nml=parameters_restart)
       close(117)
    endif

  end subroutine write_restart_file

  !--------------------------------------Output spectral coefficients
  SUBROUTINE output_coeff()
    IMPLICIT NONE

    INTEGER(KIND=4) :: fh,status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_OFFSET_KIND) :: disp=0
    CHARACTER(8) ::  suffix
    INTEGER :: i
    TYPE(spec), DIMENSION(m_r,mp_f) :: f_hat_mp

    REAL(KIND=8) k_th0_ifile,k_z0_ifile,eta_ifile
    CHARACTER(40) git_ifile
    CHARACTER(16) arch_ifile
    CHARACTER(128) cmp_ifile

    NAMELIST /parameters_info/ k_th0_ifile,k_z0_ifile,eta_ifile,git_ifile,arch_ifile,cmp_ifile

    Call perfon(' io coeff')

    WRITE(suffix,'(I0.8)') i_time
    fName_ic='coeff_'//trim(fBase_ic)//'.'//suffix

    
    i = INDEX(fName_ic//' ',' ') - 1
    CALL Mpi_File_open(comm,fName_ic(1:i),MPI_MODE_CREATE+MPI_MODE_RDWR,MPI_INFO_NULL,fh,ierr)
    CALL MPI_File_set_view(fh,disp,mpi_spec,filetype,"native",MPI_INFO_NULL,ierr)

    f_hat_mp%ur   = u_hat_mp%r
    f_hat_mp%uth  = u_hat_mp%th
    f_hat_mp%uz   = u_hat_mp%z
    f_hat_mp%p    = p_hat_mp
#ifdef TE_CODE
    f_hat_mp%T    = u_hat_mp%T
#endif
    f_hat_mp%k_th = fk_mp%th
    f_hat_mp%k_z  = fk_mp%z
    
    counter = m_r*mp_f
    CALL MPI_File_write(fh,f_hat_mp(1,1),counter,mpi_spec,status,ierr)
    
    CALL Mpi_File_close(fh,ierr)

    call perfoff

    
!write metadata
    if (myid == root) then

       print*,'written coeff file to disk: '//fName_ic

       eta_ifile=eta
       k_th0_ifile=k_th0
       k_z0_ifile=k_z0
       m_r_ifile=m_r
       m_th_ifile=m_th
       m_z_ifile=m_z
       git_ifile=git_id
       arch_ifile=arch_id
       cmp_ifile=cmp_flgs

       ! Write metadata for each coeff file for archival purposes
       open(unit=107, file=trim(fName_ic)//'.info')
       write(107, nml=parameters_restart)
       write(107, nml=parameters_info)
       close(107)


    endif

  END SUBROUTINE output_coeff

  !------------------------------------------------------------------
  SUBROUTINE output_energy()
    !--------------------------------
    ! Output the modal kinetic enegy 
    !--------------------------------
    IMPLICIT NONE
    TYPE (vec_mpi)  :: pert_u_hat_mp
    real(kind=8)    :: eneg_pp(m_r,mp_f)
    real(kind=8)    :: eneg_mode_mp(mp_f),eneg_mode(m_th+1,m_z)
    real(kind=8)    :: eneg_th(m_th+1),eneg_z(m_z)
    real(kind=8)    :: eneg
    real(kind=8)    :: ur_base(n_r), uth_base(n_r), uz_base(n_r)
#ifdef TE_CODE
    real(kind=8)    :: enegTemp_pp(m_r,mp_f),enegTemp_mode_mp(mp_f)
    real(kind=8)    :: enegTemp_th(m_th+1),enegTemp_mode(m_th+1,m_z)
    real(kind=8)    :: t_base(n_r)
#endif
    integer(kind=4) :: i,j
    character(30)   :: fmt_th, fmt_z ! format string
    integer(kind=4),allocatable :: displ(:)

    Call perfon(' io energy')

   ! compute Taylor-Couette analytical (laminar) base flow profiles
#ifndef TE_CODE
   !CALL base_flow(uth_base)
    call base_flow(.true., .false., ur_base, uth_base, uz_base)
#else
   !CALL base_flow(uth_base,uz_base,T_base)
    call base_flow(.true., .true., ur_base, uth_base, uz_base, t_base)
#endif
    
    pert_u_hat_mp = u_hat_mp
    IF(myid == root) THEN
       pert_u_hat_mp%r(:,1)  = u_hat_mp%r(:,1)  - dcmplx(0,0)
       pert_u_hat_mp%th(:,1) = u_hat_mp%th(:,1) - dcmplx(uth_base(:),0)
#ifndef TE_CODE
       pert_u_hat_mp%z(:,1)  = u_hat_mp%z(:,1)  - dcmplx(0,0)
#else
       pert_u_hat_mp%z(:,1)  = u_hat_mp%z(:,1)  - dcmplx(uz_base(:),0)
       pert_u_hat_mp%T(:,1)  = u_hat_mp%T(:,1)  - dcmplx(t_base(:),0)
#endif
    END IF
       
    eneg_pp = pert_u_hat_mp%r*CONJG(pert_u_hat_mp%r) &
         + pert_u_hat_mp%th*CONJG(pert_u_hat_mp%th) &
         + pert_u_hat_mp%z*CONJG(pert_u_hat_mp%z)

#ifdef TE_CODE
    enegTemp_pp = pert_u_hat_mp%T*CONJG(pert_u_hat_mp%T)
#endif
    
    ! multiply the factor rdr
    DO i = 1,m_r
       eneg_pp(i,:) = eneg_pp(i,:)*r(i)
#ifdef TE_CODE
       enegTemp_pp(i,:) = enegTemp_pp(i,:)*r(i)
#endif
    END DO

    
    ! multiply the factor 1/(2*Vol) in front of integral
    eneg_pp = eneg_pp*(1.D0-eta)/(1.D0+eta)
    eneg_pp(1,:) = 0.d0
#ifdef TE_CODE
    enegTemp_pp = enegTemp_pp*(1.D0-eta)/(1.D0+eta)
#endif /* TE_CODE */

    !compute displacements for the subsequent mpi_gatherv
    allocate(displ(numprocs))
    displ(1)=0
    do i=2,numprocs
       displ(i)=displ(i-1)+mp_f_arr(i-1)
    enddo

    call DGEMV('T',m_r,mp_f,1.d0,eneg_pp,m_r,row_inv_Mw_dr,1,0.d0,eneg_mode_mp,1)
    CALL MPI_Gatherv(eneg_mode_mp,mp_f,MPI_REAL8,eneg_mode,mp_f_arr,displ, &
         MPI_REAL8,root,comm,ierr)


#ifdef TE_CODE
    enegTemp_pp(1,:) = 0.d0
    call DGEMV('T',m_r,mp_f,1.d0,enegTemp_pp,m_r,row_inv_Mw_dr,1,0.d0,enegTemp_mode_mp,1)
    CALL MPI_Gatherv(enegTemp_mode_mp,mp_f,MPI_REAL8,enegTemp_mode,mp_f_arr,displ, &
         MPI_REAL8,root,comm,ierr)
#endif
    deallocate(displ)


    
    IF (myid == root) THEN
       WRITE(fmt_th,*) m_th+1
       WRITE(fmt_z,*) m_z
       fmt_th = ADJUSTL(fmt_th)
       fmt_z  = ADJUSTL(fmt_z)
       fmt_th = "(es25.13e3,"//TRIM(fmt_th)//"es25.13e3)"
       fmt_z  = "(es25.13e3,"//TRIM(fmt_z)//"es25.13e3)"

       eneg_th = SUM(eneg_mode,2)
#ifdef TE_CODE
       enegTemp_th = SUM(enegTemp_mode,2)
#endif /* TE_CODE */

       eneg_z(:) = 2*SUM(eneg_mode,1) - eneg_mode(1,:)
       eneg = SUM(eneg_z(:))
       write(108, '(3es25.13e3)') time, eneg_mode(1,2), eneg_mode(2,2)
       write(109, TRIM(fmt_th)) time, eneg_th(:)
       write(103, TRIM(fmt_z)) time, eneg_z(1:(m_z/2+1))
       write(104, '(2es25.13e3)') time, eneg
#ifdef TE_CODE
       write(123, TRIM(fmt_th)) time, enegTemp_th(:)
#endif /* TE_CODE */
       
    END IF

    call perfoff

end subroutine output_energy



subroutine output_torque()
! Time series output of the Nusselt and friction Reynolds numbers at the inner
! and outer cylinder walls. In case of the hydrodynamic code version, this
! includes only the torque Nusselt number; i.e. the actual torque normalised by
! the torque of the analytical/laminar Taylor-Couette profile, see eq. (2) to
! (4) in Lopez et al. IJHMT 2015. In case of the code version with temperature,
! this also includes the heat transfer at both cylinder walls normalised by the
! conductive heat transfer, i.e. the regular Nusselt number.
! TODO: output to only one file e.g. nusselt.dat or friction.dat, add proper
! header information, put this to its own module e.g. mod_nusselt.f90
! add u_tau or Re_tau output so that the viscous length scale is readily
! available...
! Called from:
! nsCouette.f90
!
! t = 2*pi*mu*r_i^2*l_z*(|du_th/dr|_r = r_i - u_th/r_i)

implicit none
integer(kind=4) :: i, h ! radial running indices
integer(kind=4), parameter :: iwidth=(n_s-1)/2, idiag=iwidth+1 ! finite difference stencil
real(kind=8), dimension(n_r) :: u_th     ! azimuthal (streamwise) velocity component
real(kind=8) :: duthdr_i, duthdr_o       ! wall-normal velocity gradient at inner/outer cylinder wall
real(kind=8) :: tau_i, tau_o             ! shear stress at inner/outer cylinder wall
real(kind=8) :: nu_omega_i, nu_omega_o   ! torque Nusselt number at inner/outer wall
real(kind=8) :: u_tau_i, u_tau_o         ! friction velocity at the inner/outer wall
real(kind=8) :: re_tau_i, re_tau_o       ! friction Reynolds number at the inner/outer wall
real(kind=8) :: b                        ! for analytical/laminar Taylor-Couette profile
real(kind=8), dimension(n_r) :: tau_lam  ! analitcial/laminar shear stress profile
#ifdef TE_CODE
real(kind=8), dimension(n_r) :: tem      ! radial temperature profile
real(kind=8) :: dtdr_i, dtdr_o           ! wall-normal temperature gradient at inner/outer cylinder wall
real(kind=8), dimension(n_r) :: dtdr_lam ! analitcial/laminar temperature gradient
real(kind=8) :: nu_i, nu_o               ! heat transfer Nusselt number at inner/outer cylinder wall
#endif /* TE_CODE */

if (myid .eq. root) then

 ! analytical/circular Taylor-Couette flow
 b = eta*(re_i-eta*re_o) / ((1.0d0-eta)*(1.0d0-eta**2.0d0)) ! analytical velocity profile u_th = a*r + b/r
 tau_lam(:) = -2.0d0*b/r(:)**2                              ! shear-stress tau = du_th/dr - u_th/r
#ifdef TE_CODE
 dtdr_lam(:) = 1.0d0/(r(:)*dlog(eta))                       ! temperature gradient dT/dr = d/dr ( ln(r/r_o)/ln(eta) - 1/2 )
#endif /* TE_CODE */

 ! instantaneous mean profiles
 do i = 1, n_r                       ! from inner to outer wall
  u_th(i) = dreal(u_hat_mp%th(i, 1)) ! zeroth mode, mean velocity <u_th>_{th, z}
#ifdef TE_CODE
  tem(i) = dreal(u_hat_mp%t(i, 1))   ! zeroth mode, mean temperature <T>_{th, z}
#endif /* TE_CODE */
 end do

 ! reset gradients
 duthdr_i = 0.0d0  ! velocity
 duthdr_o = 0.0d0
#ifdef TE_CODE
 dtdr_i   = 0.0d0 ! temperature
 dtdr_o   = 0.0d0
#endif /* TE_CODE */

 ! sum up wall-normal derivatives at the inner/outer cylinder walls
 do h = 1, idiag ! loop over finite difference stencil
  duthdr_i = duthdr_i + bmw_dr(idiag-h+1, h)       * u_th(h)
  duthdr_o = duthdr_o + bmw_dr(idiag+h-1, n_r-h+1) * u_th(n_r-h+1)
#ifdef TE_CODE
  dtdr_i   = dtdr_i   + bmw_dr(idiag-h+1, h)       * tem(h)
  dtdr_o   = dtdr_o   + bmw_dr(idiag+h-1, n_r-h+1) * tem(n_r-h+1)
#endif /* TE_CODE */
 end do

 ! compute actual shear-stress
 tau_i = duthdr_i - u_th(1)/r(1)     ! at the inner wall
 tau_o = duthdr_o - u_th(n_r)/r(n_r) ! at the outer wall

 ! compute torque Nusselt numbers
 nu_omega_i = tau_i / tau_lam(1)
 nu_omega_o = tau_o / tau_lam(n_r)

 ! compute heat transfer Nusselt numbers
#ifdef TE_CODE
 nu_i = dtdr_i / dtdr_lam(1)
 nu_o = dtdr_o / dtdr_lam(n_r)
#endif /* TE_CODE */

 ! compute friction Reynolds number Re_tau = u_tau R / nu
 u_tau_i  = dsqrt(dabs(tau_i))
 u_tau_o  = dsqrt(dabs(tau_i))
 re_tau_i = u_tau_i / 2.0d0
 re_tau_o = u_tau_o / 2.0d0

 ! write time series data to file
 write(105, '(5es25.13e3)') time, nu_omega_i, nu_omega_o, re_tau_i, re_tau_o
#ifdef TE_CODE
 write(124, '(3es25.13e3)') time, nu_i, nu_o
#endif /* TE_CODE */

end if
end subroutine output_torque



subroutine setup_wavenumber(x0_k)
use mod_params, only : m_th, m_z, k_th0, k_z0
implicit none
integer(kind=4) :: j ! azimuthal running indix
integer(kind=4) :: k ! axial running indix
real(kind=8), dimension(m_f, mp_r) :: k_th ! azimuthal wave number vector
real(kind=8), dimension(m_f, mp_r) :: k_z  ! axial wave number vector
complex(kind=8), dimension(m_f, mp_r), intent(out) :: x0_k

do j = 1, m_th+1
 do k = 1, m_z
  k_th(jk(j, k), :) = (j-1) * k_th0
 end do
 do k = 1, m_z/2+1
  k_z(jk(j, k), :) = (k-1) * k_z0
 end do
 do k = 1, m_z/2-1
  k_z(jk(j, m_z/2+1+k), :) = (-m_z/2+k) * k_z0
 end do
end do
x0_k = dcmplx(k_th, k_z)

end subroutine setup_wavenumber



  !-----------------------------------------
  SUBROUTINE read_coeff()
    IMPLICIT NONE

    INTEGER(KIND=4) :: fh,status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_OFFSET_KIND) :: disp=0
    INTEGER(KIND=4) :: i
    TYPE(spec), DIMENSION(m_r,mp_f) :: f_hat_mp
    COMPLEX(KIND=8), DIMENSION(m_f,mp_r) :: x0_ur,x0_uth,x0_uz,x0_k
#ifdef TE_CODE
    COMPLEX(KIND=8), DIMENSION(m_f,mp_r) :: x0_T
#endif /* TE_CODE */
    COMPLEX(KIND=8), DIMENSION(m_r,mp_f) :: var_k

    i = INDEX(fName_ic//' ',' ') - 1
    CALL MPI_File_open(comm,fName_ic(1:i),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    if (ierr .ne. 0) then
       if (myid == root) print *,'FATAL: cannot open restart file: '//fName_ic
       stop
    endif
    if (myid == root) print*,'restarting from '//fName_ic

    if (scale .eqv. .true.) then
       
       CALL MPI_File_set_view(fh,disp,mpi_spec,filetype2,"native",&
            MPI_INFO_NULL,ierr)
       counter = m_r_ifile*mp_f_ifile
       CALL MPI_File_read(fh,aux_hat_mp(1,1),counter,mpi_spec,status,ierr)

       CALL MPI_File_close(fh,ierr)

       !extrapolate or truncate coefficients of input file to the new resolution                                                            
       !and create fk_mp%th and z (modes distributed among processors)                                                                      
       CALL readscal(aux_hat_mp%ur,u_hat_mp%r,aux_hat_mp%k_th)
       CALL readscal(aux_hat_mp%uth,u_hat_mp%th,aux_hat_mp%k_th)
       CALL readscal(aux_hat_mp%uz,u_hat_mp%z,aux_hat_mp%k_th)
#ifdef TE_CODE
       CALL readscal(aux_hat_mp%T,u_hat_mp%T,aux_hat_mp%k_th)
#endif /* TE_CODE */
 
       call setup_wavenumber(x0_k)
       call xtranspose_mpi(m_f, m_r, x0_k, var_k)
       
       fk_mp%th=DREAL(var_k)
       fk_mp%z= DIMAG(var_k)
       
       DEALLOCATE(aux_hat_mp)
       
    else
       
       CALL MPI_File_set_view(fh,disp,mpi_spec,filetype,"native",&
            MPI_INFO_NULL,ierr)
       counter = m_r*mp_f
       CALL MPI_File_read(fh,f_hat_mp(1,1),counter,mpi_spec,status,ierr)
       
       CALL MPI_File_close(fh,ierr)
       
       u_hat_mp%r  = f_hat_mp%ur
       u_hat_mp%th = f_hat_mp%uth
       u_hat_mp%z  = f_hat_mp%uz
       p_hat_mp    = f_hat_mp%p
#ifdef TE_CODE
       u_hat_mp%T  = f_hat_mp%T
#endif /* TE_CODE */
       fk_mp%th    = f_hat_mp%k_th
       fk_mp%z     = f_hat_mp%k_z

    endif

  END SUBROUTINE read_coeff
 


subroutine perturb_init()
! This subroutine name is missleading and not descriptive. Should be changed
! soon. And this subroutine can be combined with pulse_init(), whichs name
! is also not well chosen
implicit none                             
complex(kind=8), dimension(m_f, mp_r) :: x0_ur, x0_uth, x0_uz, x0_k
complex(kind=8), dimension(m_r, mp_f) :: var_k
#ifdef TE_CODE
complex(kind=8), dimension(m_f, mp_r) :: x0_t
#endif /* TE_CODE */

! should go here or even further up...   call setup_wavenumber(x0_k)

! replace the call of pulse initi with its content here. Splitting this
! up into two subroutines is not necessary, and the naming of the routines
! is especially missleading
#ifndef TE_CODE
call pulse_init(x0_ur, x0_uth, x0_uz, x0_k)
#else
call pulse_init(x0_ur, x0_uth, x0_uz, x0_k, x0_t)
#endif /* TE_CODE */



call xtranspose_mpi(m_f, m_r, x0_ur,  u_hat_mp%r)
call xtranspose_mpi(m_f, m_r, x0_uth, u_hat_mp%th)
call xtranspose_mpi(m_f, m_r, x0_uz,  u_hat_mp%z)
call xtranspose_mpi(m_f, m_r, x0_k,   var_k)
#ifdef TE_CODE
call xtranspose_mpi(m_f, m_r, x0_t,   u_hat_mp%t)
#endif /* TE_CODE */

fk_mp%th = dreal(var_k)
fk_mp%z  = dimag(var_k)
counter  = m_r*mp_f

end subroutine perturb_init


  
#ifndef TE_CODE
subroutine pulse_init(x0_ur, x0_uth, x0_uz, x0_k)
#else
subroutine pulse_init(x0_ur, x0_uth, x0_uz, x0_k, x0_t)
#endif
! Compute wave number vector, Construct inital conditions, Set base flow (resting fluid or analytical
! solution for velocty and/or temperatre), add 3d perturbation on top in
! physical space, divergence free and consistent with no-slip wall boundary
! conditions, transform initial conditions (physical space) to Fourier space
implicit none
integer(kind=4) :: i, ii, j, k ! Radial, azimuthal and axial running indices
integer(kind=4) :: l           ! Perturbation index 
integer(kind=4) :: p_th, p_z   ! lth perturbation wavevector (k_th(l), k_z(l))
real(kind=8) :: p_a, reamp     ! lth perturbation amplitude a(l)
real(kind=8) :: rr, tt, zz     ! Radial, azimuthal, axial location
real(kind=8) :: pr, prr, prth, prz ! Radius dependent perturbation terms
real(kind=8) :: psin, pcos         ! Harmonic perturbation terms
real(kind=8), dimension(n_r) :: ur_base, uth_base, uz_base             ! 1d base flow velocity profiles
real(kind=8), dimension(n_th, n_z, mp_r) :: ur_pert, uth_pert, uz_pert ! 3d velocity perturbations
real(kind=8), dimension(m_f, mp_r) :: k_th, k_z
real(kind=8), dimension(n_th, n_z) :: r_var ! Temporary variable at one radial (r) location
complex(kind=8), dimension(m_f, mp_r), intent(out) :: x0_ur, x0_uth, x0_uz, x0_k
#ifdef TE_CODE
real(kind=8), dimension(n_r):: t_base              ! 1d Base flow temperature profile
real(kind=8), dimension(n_th, n_z, mp_r) :: t_pert ! 3d temperature perturbation
complex(kind=8), dimension(m_f, mp_r), intent(out) :: x0_t
#endif

! Compute wavenumbers
call setup_wavenumber(x0_k)

! Set (analytical) base flow profiles
#ifndef TE_CODE
call base_flow(ic_tcbf, ic_temp, ur_base, uth_base, uz_base)
#else /* TE_CODE */
call base_flow(ic_tcbf, ic_temp, ur_base, uth_base, uz_base, t_base)
#endif /* TE_CODE */

! No perturbations on top of base flow by default
ur_pert(:, :, :)  = 0.0d0
uth_pert(:, :, :) = 0.0d0
uz_pert(:, :, :)  = 0.0d0
#ifdef TE_CODE
t_pert(:, :, :)   = 0.0d0
#endif

! Set finite amplitude perturbation on top of the base flow
if (ic_pert .eqv. .true.) then
 do l = 1, ic_np               ! Over all perturbations
  p_a = dabs(ic_p(l, 1))       ! amplitude
  p_th = int(dabs(ic_p(l, 2))) ! azimuthal mode
  p_z  = int(dabs(ic_p(l, 3))) ! axial mode
  if (p_a .le. 0.0d0) cycle    ! switch off, next perturbation
  reamp = dabs(re_o-re_i)      ! scale amplitude with Reynolds
  if (dabs(re_o-re_i) .eq. 0.0d0) reamp = re_i
  if (p_z .eq. 0) then
   p_a = p_a*reamp*k_th0/eta
  else
   p_a = p_a*reamp*k_z0/eta
  end if
  do i = 1, mp_r            ! Over all radial (r) points on this process
   ii = myid*mp_r+i         ! Global radial index
   rr = r(ii)               ! Radial location (r)
   pr = (rr-r_i) * (rr-r_o)
   prr = pr**2.0d0          ! Radial ansatz to fulfil wall boundary conditions
   prth = -pr*(5.0d0*rr**2.0d0 - 3.0d0*(r_i+r_o)*rr + r_i*r_o) ! -d(r*prr)/dr
   prz = prth/rr                                               ! -1/r*d(r*prr)/dr
   do k = 1, n_z                ! Over all axial (z) points
    zz = len_z/(n_z-1)*(k-1)    ! Axial location (z)
    do j = 1, n_th              ! Over all azimuthal (theta) points
     tt = len_th/(n_th-1)*(j-1) ! Azimuthal location (theta)
     psin = dsin(k_th0*p_th*tt + k_z0*p_z*zz)
     pcos = dcos(k_th0*p_th*tt + k_z0*p_z*zz)
     if (.not. ((p_th .eq. 0) .and. (p_z .eq. 0))) then
      ur_pert(j, k, i) = ur_pert(j, k, i) + p_a*prr*pcos
      if (p_z .eq. 0) then
       uth_pert(j, k, i) = uth_pert(j, k, i) + p_a*prth*psin/(k_th0*p_th)
      else
       uz_pert(j, k, i) = uz_pert(j, k, i) + p_a*prz*psin/(k_z0*p_z)
      end if
     end if
    end do
   end do
  end do
 end do
end if

! Compute spectral coefficients for initial condition
call planner_fftw(n_th, n_z, .false., fftw_nthreads)
do i = 1, mp_r ! over all radial (r) points on this process

  ! Radial component coefficients
  r_var(:, :) = ur_base(myid*mp_r+i) + ur_pert(:, :, i)
  call fwd_fft(r_var, x0_ur(:, i))

  ! Azimuthal component coefficients
  r_var(:, :) = uth_base(myid*mp_r+i) + uth_pert(:, :, i)
  call fwd_fft(r_var, x0_uth(:, i))

  ! Axial component coefficients
  r_var(:, :) = uz_base(myid*mp_r+i) + uz_pert(:, :, i)
  call fwd_fft(r_var, x0_uz(:, i))

  ! Temperature coefficients
#ifdef TE_CODE
  r_var(:, :) = t_base(myid*mp_r+i) + t_pert(:, :, i)
  call fwd_fft(r_var, x0_t(:, i))
#endif

end do
call destroyer_fftw()

end subroutine pulse_init



subroutine base_flow(tcbf, temp, ur_base, uth_base, uz_base, t_base)
! Set a base flow profile for all three velocity components and the temperature.
! Each profile is only a function of the radial location (r), and fully
! available on every process. Resting fluid or circular Taylor-Couette base
! flow is possible. Analytical solution to Navier-Stokes
! U = [u_r, u_th, u_z] = (0, c_1*r + c_2/r, 0)
! Additionally a temperature profile can be specified in combination with an
! axial velocity component.
implicit none
integer(kind=4) :: i
real(kind=8), intent(out) :: ur_base(n_r)
real(kind=8), intent(out) :: uth_base(n_r)
real(kind=8), intent(out) :: uz_base(n_r)
real(kind=8), intent(out), optional :: t_base(n_r)
real(kind=8) :: a, b, c     ! Coefficients for analytical profiles
logical, intent(in) :: tcbf ! Set Taylor-Couette profile as base flow
logical, intent(in) :: temp ! Set temperature profile and non-zero axial component

! Set resting fluid
ur_base(:)  = 0.0d0
uth_base(:) = 0.0d0
uz_base(:)  = 0.0d0

! Set Taylor-Couette analytical profile
if (tcbf .eqv. .true.) then
 a = (re_o-eta*re_i) / (1.0d0+eta)
 b = eta*(re_i-eta*re_o) / ((1.0d0-eta)*(1.0d0-eta**2.0d0))
 uth_base(:) = a*r(:) + b/r(:)
end if

! Set temperature profile and non-zero axial velocity component
#ifdef TE_CODE
if (present(t_base)) then
 if (temp .eqv. .true.) then
  t_base(:) = (dlog(r(:)/r_o)/dlog(eta))-0.5d0
  c = -(4.0d0*dlog(eta)+(1.0d0-eta**2.0d0)*(3.0d0-eta**2.0d0)) / &
  & (16.0d0*(1.0d0-eta**2.0d0)*((1.0d0+eta**2.0d0)*log(eta)+1.0d0-eta**2.0d0))
  uz_base(:) = gr*(c*(r(:)**2.0d0-r_i**2.0d0) + (c*(r_o**2.0d0-r_i**2.0d0) + &
  & 0.25d0*(r_o**2.0d0-r(:)**2.0d0))*(dlog(r(:)/r_i)/dlog(eta)))
 else
  t_base(:) = 0.0d0
 end if
end if
#endif

end subroutine base_flow

 

#ifndef HDF5IO    
! just a stub
  SUBROUTINE hdf5output(p,u_r,u_th,u_z,T)
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(9*n_f/4,mp_r),INTENT(IN) :: u_r,u_th,u_z,p
    REAL(KIND=8),DIMENSION(9*n_f/4,mp_r),INTENT(IN),OPTIONAL :: T
  end SUBROUTINE hdf5output
#else

  SUBROUTINE hdf5output(p,u_r,u_th,u_z,T)

    USE HDF5
    USE mod_hdf5io
    IMPLICIT NONE

    REAL(KIND=8),DIMENSION(9*n_f/4,mp_r),INTENT(IN) :: u_r,u_th,u_z,p
    REAL(KIND=8),DIMENSION(9*n_f/4,mp_r),INTENT(IN),OPTIONAL :: T

    type(datafile) :: file
    character*8  suffix
    character*128 filename
    integer error,i,j,k,kj


    INTEGER :: nthfine,nzfine,nfine
    REAL(KIND=8) :: dthfine,dzfine,thfine(3*n_th/2+1),zfine(3*n_z/2+1)
    INTEGER(HSIZE_T) :: dimsf(3)
    REAL(KIND=8), dimension(:,:),allocatable :: ur_pbc, uth_pbc, uz_pbc, p_pbc
#ifdef TE_CODE
    REAL(KIND=8), dimension(:,:),allocatable :: T_pbc
#endif
    
    !setup output grid
    nthfine = size(thfine)-1
    dthfine = len_th/nthfine
    thfine(1:nthfine+1) = (/(th(1)+i*dthfine, i=0,nthfine)/)

    nzfine = size(zfine)-1
    dzfine = len_z/nzfine
    zfine(1:nzfine+1) = (/(z(1)+i*dzfine, i=0,nzfine)/)

    dimsf(1:3)=(/nthfine+1,nzfine+1,m_r/)
    nfine = int(dimsf(1)*dimsf(2))

    allocate( ur_pbc(nfine,mp_r), uth_pbc(nfine,mp_r), uz_pbc(nfine,mp_r), p_pbc(nfine,mp_r))
#ifdef TE_CODE
    allocate(T_pbc(nfine,mp_r))
#endif


    !construct filename
    write(suffix, '(I0.8)') i_time
    filename='fields_'//trim(fBase_ic)//'_'//suffix
    call create_file(trim(filename)//'.h5',file)

    call perfon('    h5_1')
    !write basic setup parameters (not registered with xmdf) and grid
    call h5gcreate_f(file%id,'setup',file%current_group, error)
    call write_hdf(file,git_id,'code:git-id')
    call write_hdf(file,Re_i,'Re_i')
    call write_hdf(file,Re_o,'Re_o')
#ifdef TE_CODE
    call write_hdf(file,Gr,'Gr')
    call write_hdf(file,Pr,'Pr')
#endif /* TE_CODE */
    
    !add others here, if required 
    call h5gclose_f(file%current_group, error)

    call h5gcreate_f(file%id,'grid',file%current_group, error)
    call write_hdf(file,time,'time')
    call write_hdf(file,i_time,'step')

    call write_hdf(file,r,'r')
    call write_hdf(file,thfine,'th')
    call write_hdf(file,zfine,'z')
    call h5gclose_f(file%current_group, error)
    call perfoff

    ! Add the periodic boundary in th- & z- direction
    DO i = 1,mp_r
       p_pbc(:,i) = add_periodBC(nthfine,nzfine,p(:,i))
       ur_pbc(:,i) = add_periodBC(nthfine,nzfine,u_r(:,i))
       uth_pbc(:,i) = add_periodBC(nthfine,nzfine,u_th(:,i))
       uz_pbc(:,i) = add_periodBC(nthfine,nzfine,u_z(:,i))
#ifdef TE_CODE
       T_pbc(:,i) = add_periodBC(nthfine,nzfine,T(:,i))
#endif /* TE_CODE */
    END DO

    call perfon('    h5_2')
    !write pressure field (collective I/O)
    CALL h5gcreate_f(file%id,'fields',file%current_group, error)
    call write_hdf(file,p_pbc,dimsf,'pressure')
#ifdef TE_CODE
    call write_hdf(file,T_pbc,dimsf,'temperature')
#endif /* TE_CODE */
    CALL h5gclose_f(file%current_group, error)
    call perfoff
   
    call perfon('    h5_3')
    !write velocity vector field (collective I/O)
    CALL h5gcreate_f(file%id, "fields/velocity",file%current_group, error)
    call write_hdf(file,ur_pbc, dimsf,'u_r')
    call write_hdf(file,uth_pbc,dimsf,'u_th')
    call write_hdf(file,uz_pbc, dimsf,'u_z')
    CALL h5gclose_f(file%current_group, error)
    call perfoff

    call close_file(file)

    call perfon('    h5_4')
    if (myid==0) then
       call write_xdmf(trim(filename),time,dimsf)
       print *,'written hdf5/xdmf files to disk: ',trim(filename),'.{h5,xmf}'
    endif
    call perfoff

    deallocate(ur_pbc, uth_pbc, uz_pbc, p_pbc)
#ifdef TE_CODE
    deallocate(T_pbc)
#endif

  END SUBROUTINE hdf5output
#endif

  SUBROUTINE hdf5init
 
#ifdef HDF5IO    
    USE mod_hdf5io

    call init_io(comm,MPI_INFO_NULL,numprocs,myid)
#endif

  END SUBROUTINE hdf5init
  
  SUBROUTINE hdf5finalize

#ifdef HDF5IO    
    USE mod_hdf5io

    call finalize_io
#endif

  END SUBROUTINE hdf5finalize
  
  !---------------------------------------------
  FUNCTION add_periodBC(m,n,A) RESULT(Anew)
    !---------------------------------
    ! Add the periodic boundary
    ! M(m,n) : Double precision
    ! Mnew(m+1,n+1) : with periodic B.C.
    !---------------------------------
    IMPLICIT NONE
    INTEGER(KIND=4),INTENT(IN) :: m,n 
    REAL(KIND=8), INTENT(IN)   :: A(m*n)
    REAL(KIND=8)               :: Anew((m+1)*(n+1))
    REAL(KIND=8)               :: temp(m+1,n+1)
    INTEGER(KIND=4)            :: i,j,k

    DO j=1,n
       DO i=1,m
          temp(i,j) = A((j-1)*m+i)
       END DO
    END DO

    ! Periodic B.C.
    temp(1:m,n+1) = temp(1:m,1)
    temp(m+1,:) = temp(1,:)

    DO k=1,(m+1)*(n+1)
       i = mod(k,m+1)
       IF (i==0) i=m+1
       j = (k-i)/(m+1)+1
       Anew(k)=temp(i,j)
    END DO

  END FUNCTION add_periodBC

  !------------------------------------------------------------------------    
  !Replace a  nxn matrix A by its inverse
  !------------------------------------------------------------------------
  
  subroutine mat_inv(n,A,lda)
      integer(kind=4), intent(in) :: n, lda
      real(kind=8), intent(inout) :: A(lda,n)
      real(kind=8) :: work(n)
      integer(kind=4) :: info, ipiv(n)

      call dgetrf(n, n, A, lda, ipiv, info)
      if(info /= 0) stop 'matrix inversion error 1'

      call dgetri(n, A, lda, ipiv, work, n, info)
      if(info /= 0) stop 'matrix inversion error 2'

    end subroutine mat_inv

!----------------------------------------------------------------------
    SUBROUTINE readscal(v_input, v_output, mp_kth)                                                                                       
      
      IMPLICIT NONE                                                                                                                            
      
      COMPLEX(kind=8), intent(in), dimension(m_r_ifile,mp_f_ifile) :: v_input
      COMPLEX(kind=8), intent(out), dimension(m_r,mp_f) :: v_output 
      real(kind=8),intent(in),dimension(m_r_ifile,mp_f_ifile)  :: mp_kth                                                                       
      COMPLEX(kind=8),dimension(m_r,mp_f_ifile) :: v_int                                                                                       
      COMPLEX(kind=8),dimension(m_f_ifile,mp_r) :: v_aux                                                                                        
      COMPLEX(kind=8),dimension(m_f,mp_r) :: v_out_aux                                                                                         
      REAL(kind=8), allocatable, dimension(:) :: r_input                                                                                       
      REAL(kind=8) :: A(n_s,n_s,m_r)                                                                                                           
      complex(kind=8) :: fn(m_r_ifile)                                                                                             
      real(kind=8) ::  fo(m_r,2)                                                                                                               
      REAL(kind=8) :: par_real(1:m_r_ifile), par_imag(1:m_r_ifile)                                                         
      INTEGEr(kind=4) :: i                                                                                                                     
      
      if (m_r .ne. m_r_ifile) then                                                                                                             
         
         allocate(r_input(m_r_ifile))                                                                                              
         
         !compute radial points  (Chebyshev distributed nodes)
         r_input(1:m_r_ifile) = (/(((r_i+r_o)/2 - COS(PI*i/(m_r_ifile-1))/2), i=0,(m_r_ifile-1))/)
         
         !Get interpolation weights                                                                                                            
         
         call interp_wts(m_r_ifile, r_input, m_r ,r, A)                                                                              
         
         !Interpolate to new grid                                                                                                              
         
         do i = 1,mp_f_ifile                                                                                                                   
            
            fn(1:m_r_ifile)=v_input(:,i)                                                                                                       
            par_real = dreal(fn)                                                                                                               
            par_imag = dimag(fn) 
            
            call interp(m_r_ifile, r_input, par_real, A, m_r, r, fo(:,1))                                                            
            call interp(m_r_ifile, r_input, par_imag, A, m_r, r, fo(:,2))                                                            
            
            v_int(:,i) = dcmplx(fo(:,1),fo(:,2))
            
         end do
         deallocate(r_input)                                                                                                                   
      end if
      
      
      if ((m_th .eq. m_th_ifile) .and. (m_z .eq. m_z_ifile))  then
         
         v_output = v_int
         
      else
         
         !truncate or expand matrices
         
         !Initialise matrizx with zeros
         
         v_out_aux(:,:) = dcmplx(0d0,0d0)
         
         
         if (m_r_ifile .eq. m_r) then
            
            CALL xTranspose_mpi(m_r, m_f_ifile, v_input , v_aux,ifile = .true.)
            
         else
            
            CALL xTranspose_mpi(m_r, m_f_ifile, v_int , v_aux,ifile = .true.)
            
         end if
         
         if ((m_th .lt. m_th_ifile) .and. (m_z .eq. m_z_ifile)) then
            
            
            do i = 1,m_z
               
               v_out_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th+1,:) = &
                    v_aux((i-1)*(m_th_ifile+1)+1:(i-1)*(m_th_ifile+1)+m_th+1,:)
               
            end do
            
            v_out_aux(m_th+1:m_f:m_th+1,:)=dcmplx(0d0,0d0)
            CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output)
            
         elseif((m_th .gt. m_th_ifile) .and. (m_z .eq. m_z_ifile)) then
            
            do i = 1,m_z
               
               v_out_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th_ifile+1,:) = &
                    v_aux((i-1)*(m_th_ifile+1)+1:(i-1)*(m_th_ifile+1)+m_th_ifile+1,:)
               
            end do

            v_out_aux(m_th+1:m_f:m_th+1,:)=dcmplx(0d0,0d0)
            
            CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output)
            
         elseif((m_th .eq. m_th_ifile) .and. (m_z .lt. m_z_ifile)) then
                          
            
            do i = 2,m_z0
               v_out_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th+1,:) = &
                    v_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th+1,:)
               v_out_aux((m_z+1-i)*(m_th+1)+1:(m_z+1-i)*(m_th+1)+m_th+1,:) = &
               v_aux((m_z_ifile+1-i)*(m_th+1)+1:(m_z_ifile+1-i)*(m_th+1)+m_th+1,:)
               
            end do
            
            !mode zero and mode m_z0
            
            v_out_aux(1:m_th+1,:) = v_aux(1:m_th+1,:)
            v_out_aux(m_z0*(m_th+1)+1:m_z0*(m_th+1)+m_th+1,:) = dcmplx(0d0,0d0)
                
            
            CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output)
            
        
         elseif((m_th .eq. m_th_ifile) .and. (m_z .gt. m_z_ifile)) then
            
            
            do i = 2,m_z_ifile/2
               
               v_out_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th+1,:) = &
                    v_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th+1,:)
               v_out_aux((m_z+1-i)*(m_th+1)+1:(m_z+1-i)*(m_th+1)+m_th+1,:) = &
                    v_aux((m_z_ifile+1-i)*(m_th+1)+1:(m_z_ifile+1-i)*(m_th+1)+m_th+1,:)

            end do

            
            !mode zero and mode m_z0
            
            
            v_out_aux(1:m_th+1,:) = v_aux(1:m_th+1,:)
            v_out_aux(m_z_ifile/2*(m_th+1)+1:m_z_ifile/2*(m_th+1)+m_th+1,:) = dcmplx(0d0,0d0)
            
            
            CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output)

            
         elseif((m_th .gt. m_th_ifile) .and. (m_z .gt. m_z_ifile)) then
            
            
            do i = 2,m_z_ifile/2
               
               v_out_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th_ifile+1,:) = &
                    v_aux((i-1)*(m_th_ifile+1)+1:(i-1)*(m_th_ifile+1)+m_th_ifile+1,:)
               v_out_aux((m_z+1-i)*(m_th+1)+1:(m_z+1-i)*(m_th+1)+m_th_ifile+1,:) = &
                    v_aux((m_z_ifile+1-i)*(m_th_ifile+1)+1:(m_z_ifile+1-i)*(m_th_ifile+1)+m_th_ifile+1,:)
               
               
            end do
            
            
            !mode zero and mode m_z0
            
            
            v_out_aux(1:m_th_ifile+1,:) = v_aux(1:m_th_ifile+1,:)
            v_out_aux(m_z_ifile/2*(m_th+1)+1:m_z_ifile/2*(m_th+1)+m_th_ifile+1,:) = &
                 v_aux(m_z_ifile/2*(m_th_ifile+1)+1:m_z_ifile/2*(m_th_ifile+1)+m_th_ifile+1,:)
            
            CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output)
            
            
         elseif((m_th .lt. m_th_ifile) .and. (m_z .lt. m_z_ifile)) then
            
            
            do i = 2,m_z0
               
               v_out_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th+1,:) = &
                    v_aux((i-1)*(m_th_ifile+1)+1:(i-1)*(m_th_ifile+1)+m_th+1,:)
               v_out_aux((m_z+1-i)*(m_th+1)+1:(m_z+1-i)*(m_th+1)+m_th+1,:) = &
                    v_aux((m_z_ifile+1-i)*(m_th_ifile+1)+1:(m_z_ifile+1-i)*(m_th_ifile+1)+m_th+1,:)
               
            end do
            
            !mode zero and mode m_z0

            v_out_aux(1:m_th+1,:) = v_aux(1:m_th+1,:)
            v_out_aux(m_z0*(m_th+1)+1:m_z0*(m_th+1)+m_th+1,:) = &
                 v_aux(m_z0*(m_th_ifile+1)+1:m_z0*(m_th_ifile+1)+m_th+1,:)
            
            CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output)
            
         elseif((m_th .gt. m_th_ifile) .and. (m_z .lt. m_z_ifile)) then
            
            do i = 2,m_z0

               v_out_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th_ifile+1,:) = &
                    v_aux((i-1)*(m_th_ifile+1)+1:(i-1)*(m_th_ifile+1)+m_th_ifile+1,:)
               v_out_aux((m_z+1-i)*(m_th+1)+1:(m_z+1-i)*(m_th+1)+m_th_ifile+1,:) = &
                    v_aux((m_z_ifile+1-i)*(m_th_ifile+1)+1:(m_z_ifile+1-i)*(m_th_ifile+1)+m_th_ifile+1,:)
               
            end do
            
            !mode zero and mode m_z0
            
            v_out_aux(1:m_th_ifile+1,:) = v_aux(1:m_th_ifile+1,:)
            v_out_aux(m_z0*(m_th+1)+1:m_z0*(m_th+1)+m_th_ifile+1,:) = &
                 v_aux(m_z0*(m_th_ifile+1)+1:m_z0*(m_th_ifile+1)+m_th_ifile+1,:)
        
            CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output)
            


         elseif((m_th .lt. m_th_ifile) .and. (m_z .gt. m_z_ifile)) then
            
            do i = 2,m_z_ifile/2
               
               v_out_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th+1,:) = &
                    v_aux((i-1)*(m_th_ifile+1)+1:(i-1)*(m_th_ifile+1)+m_th+1,:)
               v_out_aux((m_z+1-i)*(m_th+1)+1:(m_z+1-i)*(m_th+1)+m_th+1,:) = &
                    v_aux((m_z_ifile+1-i)*(m_th_ifile+1)+1:(m_z_ifile+1-i)*(m_th_ifile+1)+m_th+1,:)

               
            end do
            

            !mode zero and mode m_z0
            
            
            v_out_aux(1:m_th+1,:) = v_aux(1:m_th+1,:)
            v_out_aux(m_z_ifile/2*(m_th+1)+1:m_z_ifile/2*(m_th+1)+m_th+1,:) = &
                 v_aux(m_z_ifile/2*(m_th_ifile+1)+1:m_z_ifile/2*(m_th_ifile+1)+m_th+1,:)

            
            CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output)
            
            
         end if
      end if

    END SUBROUTINE readscal



integer function jk(j, k)
! Convert 2d index tuple (j, k) to a single 1d running index (jk). Mapping from
! (1..m_th+1, 1..mz)-space to a linearised 1d coordinate in (1..m_f)-space, as
! e.g. required by subroutine xtranspose()
use mod_params, only : m_th      ! number of azimuthal Fourier modes 
implicit none
integer(kind=4), intent(in) :: j ! azimuthal running index
integer(kind=4), intent(in) :: k ! axial running index
jk = (m_th+1) * (k-1)+j
end function jk


      
!--------------------------------------------------------------------------                                                                
!  interpolate; return weights                                                                                                             
!--------------------------------------------------------------------------                                                                
                                                                                                                                           
                                                                                                                                           
Subroutine interp_wts(nr_input, r_input, nr, r, A)                                                                                         
                                                                                                                                           
IMPLICIT NONE                                                                                                                              
                                                                                                                                           
INTEGER(kind=4), intent(in) :: nr_input, nr                                                                                                
REAL(kind=8), intent(in) :: r_input(nr_input), r(nr)                                                                                       
real(kind=8) , intent(out), dimension(n_s,n_s,nr) :: A                                                                                     
integer(kind=4) :: n,nn,i,j,left,right                                                                                                     
                                                                                                                                           
do n = 1, nr                                                                                                                               
   j = 1                                                                                                                                   
   do while(r_input(j)< r(n)-1d-8 .and. j< nr_input)                                                                                       
      j = j+1                                                                                                                              
   end do                                                                                                                                  
   left = max(1,j-(n_s-1)/2)                                                                                                               
   right = min(j+(n_s-1)/2,nr_input)                                                                                                       
   nn = right-left+1                                                                                                                       
   do i = 1, nn                                                                                                                            
      A(i,1,n) = 1d0                                                                                                                       
   end do                                                                                                                                  
   do j = 2, nn                                                                                                                            
      do i = 1, nn                                                                                                                         
         A(i,j,n) = A(i,j-1,n) * (r_input(left+i-1)-r(n)) / dble(j-1)                                                                      
      end do                                                                                                                               
   end do                                                                                                                                  
   call mat_inv(nn,A(1,1,n),n_s)                                                                                                           
end do                                                                                                                                     
                                                                                                                                           
                                                                                                                                           
END Subroutine interp_wts                                                                                                                  
                                                                                                                                           
                                                                                                                                           
!--------------------------------------------------------------------------                                                                
!  interpolate, given weights from io_interp_wts()                                                                                         
!--------------------------------------------------------------------------                                                                
   subroutine interp(nr_input,r_input,fi, A, nr, r, fo)                                                                                    
      integer(kind=4), intent(in)  :: nr_input, nr                                                                                         
      REAL(kind=8), intent(in)  :: r_input(nr_input), fi(nr_input), r(nr)                                                                  
      REAL(kind=8), intent(in)  :: A(n_s,n_s,m_r)                                                                                          
      REAL(kind=8), intent(out) :: fo(nr)                                                                                                  
      INTEGER(kind=4) :: n,nn,i,j,left,right                                                                                               
                                                                                                                                           
      do n = 1, nr                                                                                                                         
         j = 1                                                                                                                             
         do while(r_input(j)< r(n)-1d-8 .and. j < nr_input)                                                                                
            j = j+1                                                                                                                         
         end do
         left = max(1,j-(n_s-1)/2)                                                                                                         
         right = min(j+(n_s-1)/2,nr_input)                                                                                                 
         nn = right-left+1                                                                                                                 
         fo(n) = dot_product(A(1,1:nn,n),fi(left:right))                                                                                   
      end do
      
end subroutine interp



#ifdef TEST1
subroutine plot_baseflow

IMPLICIT NONE
INTEGER(kind=4) :: i
REAL(KIND=8) :: A,B,C

! compute coefficients for analytical Taylor-Couette profiles
a = (re_o-eta*re_i)/(1+eta)
b = eta*(re_i-eta*re_o)/((1d0-eta)*(1d0-eta**2))
c = -(4*log(eta)+(1-eta**2)*(3-eta**2))/(16*(1-eta**2)*((1+eta**2)*log(eta)+1-eta**2))

if ((myid .eq. 0) .and. (MOD(i_time,dn_vel)== 0)) then

 open(UNIT=4,STATUS='UNKNOWN',FILE='axial_prof.out')
 do i=1,m_r
    write(4,*) r(i),dreal(u_hat_mp%z(i,1)),Gr*(C*(r(i)**2-r_i**2)+(C*(r_o**2-r_i**2)+0.25*(r_o**2-r(i)**2))*(log(r(i)/r_i)/log(eta)))
 end do
 close(4)

 open(UNIT=4,STATUS='UNKNOWN',FILE='azimuthal_prof.out')
 do i=1,m_r
    write(4,*) r(i),dreal(u_hat_mp%th(i,1)),A*r(i)+B/r(i)
 end do
 close(4)

 open(UNIT=4,STATUS='UNKNOWN',FILE='temp_prof.out')
 do i=1,m_r
    write(4,*) r(i),dreal(u_hat_mp%T(i,1)),0.5d0+log(r(i)/r_i)/log(eta)
 end do
 close(4)

end if

end subroutine plot_baseflow
#endif /* TEST1*/



end module mod_inout
