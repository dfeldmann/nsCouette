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

PROGRAM nsPipeFlow
! A hybrid MPI/OpenMP code for pipe flow written in Fortran90
!
!   Coordinates system: (r,theta,z)
!   Primitive variables: (ur,uth,uz,p)
!   Spatial approximation: 
!   finite difference (r-dir) + fourier (th,z-dir)
!   Temporal approximation:
!     (see www.openpipeflow.org)
!   Nonlinear term: pseudospectral technique
!   MPI parallelisation:
!     broadcast the fourier modes into n processors
!   OpenMP threads inside each MPI task
!
!                                      Developed by Jose Manuel Lopez
!                                      built on a previous code for Taylor Couette flow
!                                      co-developed by Liang Shi and Markus Rampp. 
!                                      Some subroutines as well as the idea have been borrowed from 
!                                      openpipeflow.org. 
!                                           
!                                                               April 2015, Erlangen
!
! Please cite:  
!   A hybrid MPI-OpenMP parallel implementation for pseudospectral 
!   simulations with application to Taylor-Couette flow
!   Shi, L.; Rampp, M.; Hof, B.; Avila, M. 
!   Computers & Fluids, 106, 1-11 (2015)   
!   (arXiv:1311.2481)
!
!   The Openpipeflow Navier–Stokes solver
!   Willis, A. SoftwareX Volume 6, 2017, Pages 124-127

  USE mod_getcpu
  USE mod_timeStep
  USE mod_params
  USE WCTimerClass

  
  IMPLICIT NONE

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

  
  !MPI rank and size 
  CALL pre_mpi()


  if (myid == root) then
     print '(A)','!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     print '(A)','! NSPipeflow, a HPC code for DNS of pipe  flow                         !'
     print '(A)','!                                                                      !'
     print '(A)','! Copyright (C) 2016 Marc Avila, Bjoern Hof, Jose Manuel Lopez,        !'
     print '(A)','!                    Markus Rampp, Liang Shi                           !'
     print '(A)','!                                                                      !'
     print '(A)','! NSPipeflow is free software: you can redistribute it and/or modify    !'
     print '(A)','! it under the terms of the GNU General Public License as published by !'
     print '(A)','! the Free Software Foundation, either version 3 of the License, or    !'
     print '(A)','! (at your option) any later version.                                  !'
     print '(A)','!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     print '(A)'
     print '(A,A,A,I5,A,I4,A)','starting NSPipeflow (',git_id, ') using',numprocs,&
          &' MPI tasks, ',nomp,' OpenMP threads/task'
     print '(A)'
  endif
  


  allocate(tstp_all(numprocs))
  call mapping_info(numprocs,myid)   !requires getcpu.c mod_getcpu.f90    


  !read parameters from input file and distribute to all ranks
  CALL read_pars()

  !initialize the numerical grid
  call init_grid()
  CALL radial_grid(m_r,r)
  
  !intitalize the MPI distribution of data
  CALL init_mpidist()
  
  !Plan for FFTW
  CALL planner_fftw(n_th,n_z,.TRUE.,fftw_nthreads)
  

  ! Initialize variables (u,p)
  IF (restart ==  0) THEN
     !Perturb base flow
     CALL perturb_init()
  ELSE
     !Starting from previous file
     CALL read_restart_pars()
     CALL read_coeff()
  END IF
  
  !Open files to write out time series

  CALL open_files()
  IF (myid == root) CALL output_pars()
  
  !----------------------------------------Iteration phase
  Call perf_context_start('pre')
  CALL pre_timeStep(init = .true.)
  Call perf_context_end()

  Call perf_context_start('tstep')
  Call perfon ('tstep')

  
  DO i_time = i_start,i_start+numsteps-1

     call tstepTimer%start()
    
     if(variable_dt) then
        call nonlinear(output = .true.,cfl_eval =.true.)
        call new_tstep()!Set new dt
        if(new_dt)  call pre_timeStep(init = .false.)!If dt changes it is necessary to recompute the timestepping matrices
        !get predcition for the velocity field
        call predictor()
     else
        call nonlinear(output = .true., cfl_eval = .false.)
        CALL predictor()
     end if

     iter = 1

     do while(iter .ne. 0)!iterate to find solution at n+1
        call nonlinear(output=.false.,cfl_eval = .false.)
        call corrector()!Compute velocities at n+1
        call check_convergence()!check time step size
     end do
     
     !update time
     
     time = time + dt
    
     call tstepTimer%stop()
     tstep_time=tstep_time+tstepTimer%elapsedTime()
     
    
     i_steps=i_steps+1
     
     
     ! Output coeff and  modal kinetic energy
     
     IF ((MOD(i_time,dn_coeff)== 0) .and. (dn_coeff .gt. 0) ) then
        CALL output_coeff()
        call spectrum()
        call write_restart_file()
        call meanprof()
     end IF
     
     IF ((MOD(i_time,dn_friction) == 0) .and. (myid .eq. 0)) CALL friction()
     
     IF (MOD(i_time,dn_ke) == 0) CALL output_energy()

     ! check if there is still time left, make sure that all tasks are on the same time
     call MPI_BARRIER(comm, ierr)
     if (myid == root) then
        job_time=jobTimer%elapsedTime()
        terminate_job=( job_time .gt. runtime-2*(job_time/i_steps) )
        if (terminate_job) print '(A,I8,A)','*** wall clock time limit reached. Terminating job after ',int(job_time),' s ***'
     endif
     CALL MPI_Bcast(terminate_job,1,MPI_LOGICAL,root,comm,ierr)
     if (terminate_job) exit
   
  END DO

  Call perfoff ('tstep')
  Call perf_context_end()

  if (.not.terminate_job) i_time=i_time-1
  
  if (dn_coeff .gt. 0) then
     call output_coeff()
     call write_restart_file()
  endif
  

  !----------------------------------------Closing phase
  CALL close_files()
  CALL destroyer_fftw()

  call jobTimer%stop()

  CALL MPI_Gather(tstep_time,1,MPI_DOUBLE_PRECISION,tstp_all,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ierr)
  IF (myid == root) THEN
     PRINT '(A,I6)',    'T0: number of timesteps              :',i_steps
     PRINT '(A,F10.2)', 'T0: total time incl I/O [s]:',jobTimer%elapsedTime()
     PRINT '(A,F10.2)', 'T0:  av time per timestep excl I/O [s]:',tstep_time/real(i_steps)

     PRINT '(A,2F10.4)','ALL: min/max av time per timestep [s]:',minval(tstp_all)/real(i_steps),maxval(tstp_all)/real(i_steps)

  END IF

  CALL post_timeStep()
  CALL post_mpi()


  Call perfoff ()

  IF (myid == root) CALL perfout('main')


  CALL MPI_Finalize(ierr)
  
  !-----------------------End-------------------------------
END PROGRAM nsPipeFlow
