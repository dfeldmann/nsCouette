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

program nscouette
! nsCouette -- A high-performance for direct numercial simulations of turbulent
! Taylor-Couette flow written in modern Fortran.
! + Cylindrical coordinate system (r, theta, z)
! + Primitive flow field variables: Velocity (u_r, u_theta, u_z), pressure (p),
!   and optionally temperature (T).
! + Spatial approximation: Finite differences in r and Fourier-Galerkin ansatz
!   in theta and z direction
! + Temporal approximation: Predictor-corrector scheme with dynamic time
!   stepping as in openpipeflow.org
! + Nonlinear term: Pseudospectral technique
! + Hybrid parallelisation strategy: Standard one-dimensional slab decomposition
!   into Fourier modes (in theta and z), which can be treated independently of
!   each other in the solution of the linear terms. Multiple cores per MPI task
!   can be used to compute the linear therms with multiple OpenMP threads inside
!   each MPI task. To compute the non-linear terms, global data transpositions
!   based on MPI_Alltoall and task-local transposes are employed for gathering
!   all Fourier modes locally on each MPI task.
!                                                                                          
! Early versions by Liang Shi. Co-developed by Markus Rampp and Liang Shi (2011-
! 2014). Hybrid parallelisation and optimisation by Markus Rampp. Time-stepper
! and other improvements by Jose Manuel Lopez (2016). GPU version by Alberto
! Vela-Martin (2018). Minor improvements, Documentation, tutorials by Daniel
! Feldmann (2019).
!                                                                                          
! Please cite:
! [1] Jose Manuel Lopez, Daniel Feldmann, Markus Rampp, Alberto Vela-Martin,
!     Liang Shi & Marc Avila. nsCouette -- A high-performance code for direct
!     numerical simulations of turbulent Taylor-Couette flow. SoftwareX, xxx,
!     p. xx-xx, 2019.
! [2] Liang Shi, Markus Rampp, Bjoern Hof & Marc Avila. A hybrid MPI-OpenMP
!     parallel implementation for pseudospectral simulations with application to
!     Taylor-Couette flow. Computers & Fluids, 106, p. 1-11, 2015.

use mod_getcpu
use mod_timestep
use mod_params
use mod_mympi
use mod_vars
use mod_nonlinear
use mod_inout
use wctimerclass

implicit none

  TYPE(WCTimer) :: tstepTimer,jobTimer

  REAL(KIND=8) :: tstep_time,job_time
  REAL(KIND=8), allocatable :: tstp_all(:) ! total time of timestep per task
  INTEGER(KIND=4), PARAMETER :: imode=MPI_THREAD_SERIALIZED
  INTEGER(KIND=4) :: omp_get_max_threads
  LOGICAL :: terminate_job=.false.
  INTEGER(KIND=4) :: i_steps=0
  
  CALL mpi_init_thread(imode,tl_mpi,ierr)
  if (tl_mpi .lt. imode) then
     print '(A,I2,A,I2,A,I2)','mpi_init_thread failed, requested:',&
             &imode,', provided:',tl_mpi
     print '(A,I2,A)','working with provided thread level',tl_mpi,' as fallback'
  endif

!$ nomp = omp_get_max_threads()


  !----------------------------------------Initialization phase
  tstep_time=0.d0
  call jobTimer%start()

  Call perfinit
  Call perfon ('main')


  CALL pre_mpi()

  if (myid == root) then
     print '(A)','!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     print '(A)','! NSCouette, a HPC code for DNS of Taylor-Couette flow                 !'
     print '(A)','!                                                                      !'
     print '(A)','! Copyright (C) 2016 Marc Avila, Bjoern Hof, Jose Manuel Lopez,        !'
     print '(A)','!                    Markus Rampp, Liang Shi                           !'
     print '(A)','!                                                                      !'
     print '(A)','! NSCouette is free software: you can redistribute it and/or modify    !'
     print '(A)','! it under the terms of the GNU General Public License as published by !'
     print '(A)','! the Free Software Foundation, either version 3 of the License, or    !'
     print '(A)','! (at your option) any later version.                                  !'
     print '(A)','!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     print '(A)'
     print '(A,A,A,I5,A,I4,A)','starting NSCouette (',git_id, ') using',numprocs,&
          &' MPI tasks, ',nomp,' OpenMP threads/task'  
     print '(A)'
  endif


  allocate(tstp_all(numprocs))
  call mapping_info(numprocs,myid)   !requires getcpu.c mod_getcpu.f90
  
  call read_pars()    ! read set of parameters from std in
  call init_grid()    ! initialise numerical grid
  call init_mpidist() ! intitalise mpi distribution of data 

  ! set initial conditions
  if (restart .eq. 0) then   ! generate flow field
   call perturb_init() 
  else                       ! read flow field from file
   call read_restart_pars()
   call read_coeff()
  end if

 !if (dn_enfsym .gt. 0) call enfsym() ! enforce zero mode symmetry
  
  call planner_fftw(n_th, n_z, .true., fftw_nthreads)
  call open_files()
  call init_probes() ! initialise time series output at several probe locations

  call output_pars() ! write final set of parameters to file
  
  call perf_context_start('pre')
  call pre_timestep(init = .true.)
  call perf_context_end()
  
  ! main time stepping loop
  call perf_context_start('tstep')
  call perfon('tstep')
  do i_time = i_start, i_start+numsteps-1
   call tstepTimer%start()

   ! explicit predictor step 
   if (variable_dt) then           ! automatically adapt time step size dt

    call nonlinear(output=.true., cfl_eval=.true., vel_output=.true.)    ! compute non-linear terms
    call new_tstep()                                                     ! cumpute new dt
    if (new_dt) call pre_timeStep(init = .false.)                        ! recompute time stepping matrices for new dt
    call predictor()                                                     ! compute prediction for the velocity field

   else                            ! constant time step size

    call nonlinear(output=.true., cfl_eval=.false., vel_output=.true.)   ! compute non-linear terms
    call predictor()                                                     ! compute prediction for the velocity field

   end if
    
   ! iterative corrector step
   iter = 1
   do while (iter .ne. 0)                                                 ! iterate to find solution at n+1
    call nonlinear(output=.false., cfl_eval=.false., vel_output=.false.)  ! compute non-linear terms
    call corrector()                                                      ! correct velocity field at n+1
    call check_convergence()                                              ! check time step size
   end do
 
   ! enforce symmetry in k0 modes
!   if ((dn_enfsym .gt. 0) .and. (mod(i_time, dn_enfsym) .eq. 0)) call enfsym()
    
   ! update physical time
   time = time + dt
     
   call tstepTimer%stop()
   tstep_time = tstep_time + tstepTimer%elapsedTime()

#if TEST1
   call plot_baseflow()
#endif /* TEST1*/

     i_steps=i_steps+1
     
     ! Output coeff & modal kinetic energy
     IF (dn_coeff .gt. 0 .and. MOD(i_time,dn_coeff)== 0) CALL output_coeff()
     IF (MOD(i_time,dn_ke) == 0) CALL output_energy()
     IF (MOD(i_time,dn_Nu) == 0) CALL output_torque()
     
     !IF (myid == root) print '(A,2I7,1pe11.3,e11.3)','done step ',i_time,i_steps,time,jobTimer%elapsedTime()

     ! check if there is still time left, make sure that all tasks are on the same time 
     call MPI_BARRIER(comm, ierr)
     if (myid == root) then
        job_time=jobTimer%elapsedTime()
        terminate_job=( job_time .gt. runtime-2*(job_time/i_steps) )
        if (terminate_job) print '(A,I8,A)','*** wall clock time limit reached. Terminating job after ',int(job_time),' s ***'
     endif
     CALL MPI_Bcast(terminate_job,1,MPI_LOGICAL,root,comm,ierr)
     if (terminate_job) exit

  end do

  if (.not.terminate_job) i_time=i_time-1


  Call perfoff ()
  
  Call perf_context_end()
  
  if (dn_coeff .gt. 0) then
     call output_coeff()
     call write_restart_file()
  endif

  !----------------------------------------Closing phase
  call final_probes() ! finalise time series output at several probe locations
  call close_files()
  call destroyer_fftw()

  call jobTimer%stop()

  CALL MPI_Gather(tstep_time,1,MPI_DOUBLE_PRECISION,tstp_all,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ierr) 
  if (myid .eq. root) then
   print '(a)', '------------------------------------------------------------------------------------'
   print '(a, i9.8)', 'Total number of computed timesteps:', i_steps
   print '(a, f9.2, a)', 'Total elapsed WCT since start up:', jobtimer%elapsedtime(), 's'
   print '(a, 3(f9.4, a))', 'Average elapsed WCT per time step w/o coeff io (min, mean, max task):', &
         & minval(tstp_all)/real(i_steps), 's, ', tstep_time/real(i_steps), 's, ', & 
         & maxval(tstp_all)/real(i_steps), 's'
  end if

  CALL post_timeStep()
  CALL post_mpi()


  Call perfoff ()

  IF (myid == root) CALL perfout('main')


  CALL MPI_Finalize(ierr)
  !-----------------------End-------------------------------
END PROGRAM nsCouette
