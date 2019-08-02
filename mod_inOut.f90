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

MODULE mod_inOut
! Input & Output
!   Input: c.f. input file 'input_nsPipeFlow'
!   Output:
!     1) spectral coefficients for velocity and pressure
!     2) modal energy and total energy
!     3) temporal evolution of wall shear, centreline velocity, friction velocity and friction coefficient
!     4) temporal evolution of torque
!     5) global parameters to file 'parameter.dat'
!     6) Mean velocity profile

  
  USE mod_myMpi
  USE mod_fftw
  USE mod_vars
  USE mod_fdInit
  
  IMPLICIT NONE

  TYPE spec
     SEQUENCE
     COMPLEX(KIND=8) :: ur
     COMPLEX(KIND=8) :: uth
     COMPLEX(KIND=8) :: uz
     COMPLEX(KIND=8) :: p
     REAL(KIND=8)    :: k_th
     REAL(KIND=8)    :: k_z
  END TYPE spec

   
  TYPE(spec), ALLOCATABLE,DIMENSION(:,:) :: aux_hat_mp
  LOGICAL, PRIVATE ::  scale = .false.
  
  
  NAMELIST /parameters_grid/ m_r,m_th,m_z0,k_th0,k_z0
  NAMELIST /parameters_physics/ Re,const_flux
  NAMELIST /parameters_timestep/ numsteps,init_dt,variable_dt,maxdt,Courant
  NAMELIST /parameters_output/ dn_coeff,dn_ke,dn_friction,dn_hdf5,print_time_screen,fBase_ic
  NAMELIST /parameters_control/ restart,runtime
  NAMELIST /parameters_restart/ i_time,time,dt,m_r_ifile,m_th_ifile,m_z_ifile,fBase_ic
  
  private :: spec,aux_hat_mp
  private :: parameters_grid,parameters_physics,parameters_timestep,parameters_output,&
       parameters_control,parameters_restart
 CONTAINS


   !--------------------------------------------------------------
   ! Read the variable parameters from the input file/Standard IO
   SUBROUTINE read_pars()
     IMPLICIT NONE
     CHARACTER(50) :: basenm,suffix


     IF (myid == root) THEN
        READ(*,nml=parameters_grid)
        READ(*,nml=parameters_physics)
        READ(*,nml=parameters_timestep)
        READ(*,nml=parameters_output)
        READ(*,nml=parameters_control)
     END IF



     CALL MPI_Bcast(m_r, 1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(m_th,1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(m_z0,1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(k_th0,1,MPI_REAL8,root,comm,ierr)
     CALL MPI_Bcast(k_z0,1,MPI_REAL8,root,comm,ierr) 
     CALL MPI_Bcast(Re,1,MPI_REAL8,root,comm,ierr)
     CALL MPI_Bcast(restart,1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(fBase_ic,len(fBase_ic),MPI_CHARACTER,root,comm,ierr)
     CALL MPI_Bcast(init_dt,1,MPI_REAL8,root,comm,ierr)
     CALL MPI_Bcast(numsteps,1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(dn_coeff,1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(dn_ke,1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(dn_friction,1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(dn_hdf5,1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(print_time_screen,1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(variable_dt,1,MPI_LOGICAL,root,comm,ierr)
     CALL MPI_Bcast(Courant,1,MPI_REAL8,root,comm,ierr)
     CALL MPI_Bcast(maxdt,1,MPI_REAL8,root,comm,ierr)
     CALL MPI_Bcast(const_flux,1,MPI_LOGICAL,root,comm,ierr)
     CALL MPI_Bcast(runtime,1,MPI_INTEGER,root,comm,ierr)

     !Initialize dt
     dt=init_dt
     
   END SUBROUTINE read_pars
   
   !--------------------------------------------------------------
   ! Read the parameters of the initial velocity file
   SUBROUTINE read_restart_pars()
     IMPLICIT NONE
     INTEGER(KIND=4) :: ndims=2
     INTEGER(KIND=4) :: sizes(2),subsizes(2),starts(2)
     CHARACTER(8) ::  suffix

     if (myid == root) then
        open(unit=107,file='restart')
        read(107,nml=parameters_restart)
        close(107)

        !erase restart inof in order to avoid loops of restart attempts
        open(unit=107,file='restart')
        write(107,*)
        close(107)
     END IF

     CALL MPI_Bcast(i_time,1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(time,1,MPI_REAL8,root,comm,ierr)
     CALL MPI_Bcast(dt,1,MPI_REAL8,root,comm,ierr)
     CALL MPI_Bcast(m_r_ifile,1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(m_th_ifile,1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(m_z_ifile,1,MPI_INTEGER,root,comm,ierr)
     CALL MPI_Bcast(fBase_ic,len(fBase_ic),MPI_CHARACTER,root,comm,ierr)
     
     
     WRITE(suffix,'(I0.8)') i_time
     fName_ic='coeff_'//trim(fBase_ic)//'.'//suffix

     i_start = i_time+1

     if (restart == 2) then
        time = 0d0 
        i_start = 1
     end if
     
     !Check if number of modes in parameters module matches the number of modes in input file

     if((m_r .ne. m_r_ifile) .or. (m_z .ne. m_z_ifile) .or.(m_th .ne. m_th_ifile))then

        if (myid == root) then
           print*,'m_r:  ',m_r,'->',m_r_ifile
           print*,'m_th: ',m_th,'->',m_th_ifile
           print*,'m_z:  ',m_z,'->',m_z_ifile
           print *,'WARNING: grid resolution changed in restart file, now remapping ...'
        endif

        m_f_ifile = m_z_ifile*(m_th_ifile+1)
        mp_fmax_ifile = ceiling(real(m_f_ifile,kind=8)/real(numprocs,kind=8)) ! Fourier points per proc


        if ((MOD(m_f_ifile,numprocs) /= 0)) then
           if (myid==numprocs-root-1) then
              mp_f_ifile=m_f_ifile-(numprocs-1)*mp_fmax_ifile
           else
              mp_f_ifile=mp_fmax_ifile
           endif
           m_f_ifile=mp_fmax_ifile*numprocs

           if (myid == root) then
              print '(A)',                     '------------------------------------------------------------------------'
              print '(A,I9,A,I5)',             'INFO: number of Fourier modes in the initial condition (m_th_ifile+1)*m_z_ifile=',&
                   &(m_th_ifile+1)*m_z_ifile,' is not divisible by number of MPI procs',numprocs
              print '(A,I9,A,I5,A,I8,A,I8,A)', '   => m_f_ifile grid enlarged to',m_f_ifile,', populated with',&
                   &(numprocs-1),'*',mp_fmax_ifile,'+',(m_th_ifile+1)*m_z_ifile-(numprocs-1)*mp_fmax_ifile,' modes'
              print '(A)',                     '------------------------------------------------------------------------'
           end if
        else
           mp_f_ifile=mp_fmax_ifile
        endif

        sizes = (/m_r_ifile,m_f_ifile/)
        subsizes = (/m_r_ifile,mp_fmax_ifile/)
        starts = (/0,myid*mp_fmax_ifile/)
        CALL MPI_Type_create_subarray(ndims,sizes,subsizes,starts,&
             MPI_ORDER_FORTRAN,mpi_spec,filetype2,ierr)
        CALL MPI_Type_commit(filetype2,ierr)
        scale=.true.
        ALLOCATE(aux_hat_mp(m_r_ifile,mp_f_ifile))        
    end if


   END SUBROUTINE read_restart_pars
  

  !--------------------------------------Open files to write into
  SUBROUTINE open_files()
    IMPLICIT NONE
    
    CHARACTER*10 :: char_position   ! pointer position at I/O
    INTEGER :: i
        
    IF ((restart == 0) .or. (restart == 2)) THEN
       char_position = 'REWIND'
    ELSE
       char_position = 'APPEND'
    END IF
    

    i = INDEX(char_position//' ',' ') - 1  ! Length of char_position
    
    IF (myid == root) THEN
       OPEN(unit=101,file='ke_mode',position=char_position(1:i))
       OPEN(unit=102,file='ke_th',position=char_position(1:i))
       OPEN(unit=103,file='ke_z',position=char_position(1:i))
       OPEN(unit=104,file='ke_total',position=char_position(1:i))
       OPEN(unit=105,file='friction',position=char_position(1:i))
    END IF
   
    call hdf5init
    
  END SUBROUTINE open_files
  
  !--------------------------------------Close files opened above
  SUBROUTINE close_files()
    IMPLICIT NONE
    
    IF (myid == root) THEN
       CLOSE(101)
       CLOSE(102)
       CLOSE(103)
       CLOSE(104)
       CLOSE(105)
    END IF

    call hdf5finalize

  END SUBROUTINE close_files


  subroutine write_restart_file
    
    if (myid == root) then
       !(over)write minimal info required for restart from latest coeff checkpoint
       open(unit=117,file='restart')
       write(117,nml=parameters_restart)
       close(117)
    endif

  end subroutine write_restart_file
 
 
  !--------------------------------------Output spectral coefficients
  SUBROUTINE output_coeff()
    IMPLICIT NONE

    INTEGER(KIND=4) :: i
    INTEGER(KIND=4) :: fh,status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_OFFSET_KIND) :: disp=0
    CHARACTER(8) ::  suffix
    TYPE(spec), DIMENSION(m_r,mp_f) :: f_hat_mp

    REAL(KIND=8) k_th0_ifile,k_z0_ifile
    CHARACTER(40) git_ifile
    CHARACTER(16) arch_ifile
    CHARACTER(128) cmp_ifile
    LOGICAL :: const_flux_ifile
    
    NAMELIST /parameters_info/ k_th0_ifile,k_z0_ifile,const_flux_ifile,git_ifile,arch_ifile,cmp_ifile

    
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
    f_hat_mp%k_th = fk_mp%th
    f_hat_mp%k_z  = fk_mp%z
    
    counter = m_r*mp_f

    CALL MPI_File_write(fh,f_hat_mp(1,1),counter,mpi_spec,status,ierr)
    
    CALL Mpi_File_close(fh,ierr)


    !write metadata
    if (myid == root) then

       print*,'written coeff file to disk: '//fName_ic

       k_th0_ifile=k_th0
       k_z0_ifile=k_z0
       const_flux_ifile=const_flux
       m_r_ifile=m_r
       m_th_ifile=m_th
       m_z_ifile=m_z
       git_ifile=git_id
       arch_ifile=arch_id
       cmp_ifile=cmp_flgs

       !write metadata for each coeff file for archival purposes
       open(unit=107,file=trim(fName_ic)//'.info')
       write(107,nml=parameters_restart)
       write(107,nml=parameters_info)
       close(107)

    endif
       
    call perfoff

  END SUBROUTINE output_coeff

  !------------------------------------------------------------------
  SUBROUTINE output_energy()
    !--------------------------------
    ! Output the modal kinetic enegy 
    !--------------------------------
    IMPLICIT NONE
    TYPE (vec_mpi)  :: pert_u_hat_mp

    REAL(KIND=8)    :: eneg_pp(m_r,mp_f)
    REAL(KIND=8)    :: eneg_mode_mp(mp_f),eneg_mode(m_th+1,m_z)
    REAL(KIND=8)    :: eneg_th(m_th+1),eneg_z(m_z),eneg_th_(m_th+1),eneg_z_(m_z)
    REAL(KIND=8)    :: eneg,uz_base(n_r),u_z_first_deriv(n_r)
    INTEGER(KIND=4) :: i,j
    CHARACTER*20    :: fmt_th,fmt_z
    integer(kind=4),allocatable :: displ(:)
    
    Call perfon(' io energy')

    pert_u_hat_mp = u_hat_mp
    IF(myid == root) THEN
       pert_u_hat_mp%r(:,1)  = u_hat_mp%r(:,1)  
       pert_u_hat_mp%z(:,1)  = u_hat_mp%z(:,1)  - DCMPLX(1d0-r*r,0)
       pert_u_hat_mp%th(:,1)  = u_hat_mp%th(:,1)
    END IF
    

    
    eneg_pp = pert_u_hat_mp%r*CONJG(pert_u_hat_mp%r) &
         + pert_u_hat_mp%th*CONJG(pert_u_hat_mp%th) &
         + pert_u_hat_mp%z*CONJG(pert_u_hat_mp%z)

    

    ! multiply the factor 4*Pi/k_z0 in front of integral (in mode 0 it is 2Pi/k_z0)
    eneg_pp = eneg_pp*(4d0*PI*PI/k_z0)
    if(k_z0 == 0d0) stop 'kinetic energy calculation: k_z=0'
    if (myid .eq. 0)  eneg_pp(:,1)=0.5d0*eneg_pp(:,1)

    
    ! Integral over r
    DO i = 1,mp_f
       eneg_mode_mp(i) = DOT_PRODUCT(intrdr(:),eneg_pp(:,i))
    END DO


    !compute displacements for the subsequent mpi_gatherv
    allocate(displ(numprocs))
    displ(1)=0
    do i=2,numprocs
       displ(i)=displ(i-1)+mp_f_arr(i-1)
    enddo

    CALL MPI_Gatherv(eneg_mode_mp,mp_f,MPI_REAL8,eneg_mode,mp_f_arr,displ, &
         MPI_REAL8,root,comm,ierr)
   

    IF (myid == root) THEN
       WRITE(fmt_th,*) m_th
       WRITE(fmt_z,*) m_z
       fmt_th = ADJUSTL(fmt_th)
       fmt_z  = ADJUSTL(fmt_z)
       fmt_th = "(f16.10,"//TRIM(fmt_th)//"e18.10)"
       fmt_z = "(f16.10,"//TRIM(fmt_z)//"e18.10)"
       eneg_th = SUM(eneg_mode,2)
       eneg_z(:) = 2*SUM(eneg_mode,1) - eneg_mode(1,:)
       eneg = SUM(eneg_z(:))
       WRITE(101,*) time, eneg_mode(2,1),eneg_mode(2,2)
       WRITE(102,TRIM(fmt_th)) time, eneg_th(1:m_th)
       WRITE(103,TRIM(fmt_z)) time, eneg_z(1:(m_z/2))
       WRITE(104,*) time, eneg
    END IF

    call perfoff

  END SUBROUTINE output_energy

  !--------------------------------------------------
  SUBROUTINE output_pars()
    IMPLICIT NONE

    open(unit=106,file='parameters')
    write(106,nml=parameters_grid)
    write(106,nml=parameters_physics)
    write(106,nml=parameters_timestep)
    write(106,nml=parameters_output)
    close(106)
    
  END SUBROUTINE output_pars
  

  !-----------------------------------------
  SUBROUTINE read_coeff()
    IMPLICIT NONE

    INTEGER(KIND=4) :: fh,status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_OFFSET_KIND) :: disp=0
    INTEGER(KIND=4) :: i
    TYPE(spec), DIMENSION(m_r,mp_f) :: f_hat_mp
    COMPLEX(KIND=8), DIMENSION(m_f,mp_r) :: x0_ur,x0_uth,x0_uz,x0_k
    COMPLEX(KIND=8), DIMENSION(m_r,mp_f) :: var_k
    
    i = INDEX(fName_ic//' ',' ') - 1
    CALL MPI_File_open(comm,fName_ic(1:i),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    if (ierr .ne. 0) then
       if (myid == root) print *,'FATAL: cannot open restart file: '//fName_ic
       stop
    endif
   
    if(scale) then
       
       
       CALL MPI_File_set_view(fh,disp,mpi_spec,filetype2,"native",&
            MPI_INFO_NULL,ierr)
       
       counter = m_r_ifile*mp_f_ifile
       
       CALL MPI_File_read(fh,aux_hat_mp(1,1),counter,mpi_spec,status,ierr)

       CALL MPI_File_close(fh,ierr)
       
       !extrapolate or truncate coefficients of input file to the new resolution  
       !and create fk_mp%th and z (modes distributed among processors)

       CALL readscal(1,aux_hat_mp%ur,u_hat_mp%r,aux_hat_mp%k_th)
       CALL readscal(1,aux_hat_mp%uth,u_hat_mp%th,aux_hat_mp%k_th)
       CALL readscal(0,aux_hat_mp%uz,u_hat_mp%z,aux_hat_mp%k_th)
       
       CALL pulse_init(.false.,x0_ur,x0_uth,x0_uz,x0_k)
       
       CALL xTranspose_mpi(m_f,m_r,x0_k,var_k)
       
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
       fk_mp%th    = f_hat_mp%k_th
       fk_mp%z     = f_hat_mp%k_z
       
    end if
    
   
  END SUBROUTINE read_coeff

  !--------------------------------------------------------------------
  SUBROUTINE perturb_init()

    use mod_params
    
    IMPLICIT NONE
    INTEGER(kind=4) :: i
    COMPLEX(KIND=8), DIMENSION(m_f,mp_r) :: x0_ur,x0_uth,x0_uz,x0_k
    COMPLEX(KIND=8), DIMENSION(m_r,mp_f) :: var_k

    !Define base flow
    if (myid .eq. 0) u_hat_mp%z(:,1)=dcmplx(1d0-r*r,0d0)

    
    CALL pulse_init(.true.,x0_ur,x0_uth,x0_uz,x0_k)

    CALL xTranspose_mpi(m_f,m_r,x0_ur,u_hat_mp%r )
    CALL xTranspose_mpi(m_f,m_r,x0_uth,u_hat_mp%th )
    CALL xTranspose_mpi(m_f,m_r,x0_k,var_k )

    fk_mp%th=DREAL(var_k)
    fk_mp%z= DIMAG(var_k)

    !Initialise some variables

    counter = m_r*mp_f
    
    
  END SUBROUTINE perturb_init

  !--------------------------------------------------------

  SUBROUTINE pulse_init(fromrest,x0_ur,x0_uth,x0_uz,x0_k)
    
    IMPLICIT NONE
    
    LOGICAL, intent(in) :: fromrest
    COMPLEX(KIND=8), DIMENSION(m_f,mp_r), INTENT(OUT) :: x0_ur,x0_uth,x0_uz,x0_k
    REAL(KIND=8), DIMENSION(m_f,mp_r) :: k_th,k_z
    INTEGER(KIND=4) :: i,j,k
    REAL(KIND=8),    DIMENSION(9*n_f/4,mp_r) :: myu_r,myu_th        
    COMPLEX(KIND=8), DIMENSION(m_f,mp_r)     :: myvar_Cr,myvar_Cth  
    REAL(KIND=8) :: myz, myth, rr                              
    integer(kind=4) :: i_Z, i_Th                                     
    
         
    if (fromrest .eqv. .true.) then    
       
      
       i_Th=3*n_th/2
       i_Z=3*n_z/2
       
       DO i = 1,mp_r                                                                        
          rr=r(i+myid*mp_r)                                                                 
          do j=1, i_Z                                                                       
             myz=len_z/i_Z*(j-1)                                                            
             do k=1, i_Th                                                                   
                
                myth=len_th/i_Th*(k-1)                                                      
             
                myu_r((j-1)*i_Th+k,i)=1d-1*(1d0-rr**2)**2*&                                   
                     exp(-10d0*(sin(PI*myz/len_z))**2)*sin(myth)                             
               
                myu_th((j-1)*i_Th+k,i) = 1d-1*((1d0-rr**2)**2 - 4d0*rr**2*(1d0-rr**2))*&          
                     exp(-10d0*(sin(PI*myz/len_z))**2) * cos(myth)                           
             
             end do
             
          end do
          
          call fwd_fft(myu_r(1,i),myvar_Cr(1,i))                                            
          
          call fwd_fft(myu_th(1,i),myvar_Cth(1,i))                                          
          
       end do
       
       x0_ur=myvar_Cr                                                                       
       
       x0_uth=myvar_Cth                                                                     
    
    end if
    
    ! Wavenumber
    DO j = 1,m_th+1
       
       DO k = 1,m_z
          k_th(jk(j,k),:) = (j-1)*k_th0
       END DO
       
       DO k = 1,m_z/2+1
          k_z(jk(j,k),:) = (k-1)*k_z0
       END DO
       DO k = 1,m_z/2-1
          k_z(jk(j,m_z/2+1+k),:) = (-m_z/2+k)*k_z0
       END DO
       
    END DO
    
    x0_k=dcmplx(k_th,k_z)
    
    return
    
     CONTAINS
       INTEGER FUNCTION jk(j,k)
         
         ! provides a mapping from 2-d (j,k) coordinates in (1..m_th+1,1..mz)-space
         ! to a linearized 1-d coordinate in (1..m_f)-space, as required, e.g. by 
         ! subroutine xTranspose()
         
         INTEGER(KIND=4), INTENT(IN) :: j,k
         jk=(m_th+1)*(k-1)+j
         
       END FUNCTION jk
       
     END SUBROUTINE pulse_init
  
#ifndef HDF5IO    
! just a stub

  SUBROUTINE hdf5output(p,u_r,u_th,u_z)
    IMPLICIT NONE
     REAL(KIND=8),DIMENSION(9*n_f/4,mp_r),INTENT(IN) :: u_r,u_th,u_z,p
   
  end SUBROUTINE hdf5output
#else

  SUBROUTINE hdf5output(p,u_r,u_th,u_z)

    USE HDF5
    USE mod_hdf5io
    IMPLICIT NONE

    REAL(KIND=8),DIMENSION(9*n_f/4,mp_r),INTENT(IN) :: u_r,u_th,u_z,p
    type(datafile) :: file
    character*8  suffix
    integer(kind=4) ::  error,i
    character*128 filename
    
    INTEGER :: nthfine,nzfine,nfine
    REAL(KIND=8) :: dthfine,dzfine,thfine(3*n_th/2+1),zfine(3*n_z/2+1)
    INTEGER(HSIZE_T) :: dimsf(3)
    REAL(KIND=8), dimension(:,:),allocatable :: ur_pbc, uth_pbc, uz_pbc, p_pbc
    
    !setup output grid
    nthfine = size(thfine)-1
    dthfine = len_th/nthfine
    thfine(1:nthfine+1) = (/(0d0+i*dthfine, i=0,nthfine)/)

    nzfine = size(zfine)-1
    dzfine = len_z/nzfine
    zfine(1:nzfine+1) = (/(0d0+i*dzfine, i=0,nzfine)/)

    dimsf(1:3)=(/nthfine+1,nzfine+1,m_r/)
    nfine = int(dimsf(1)*dimsf(2))

    allocate( ur_pbc(nfine,mp_r), uth_pbc(nfine,mp_r), uz_pbc(nfine,mp_r), p_pbc(nfine,mp_r))
    
    

    !construct filename

    write(suffix, '(I0.8)') i_time
    filename='fields_'//trim(fBase_ic)//'_'//suffix
    call create_file(trim(filename)//'.h5',file)
    
    
    call perfon('    h5_1')
    !write basic setup parameters (not registered with xmdf) and grid
    call h5gcreate_f(file%id,'setup',file%current_group, error)
    call write_hdf(file,Re,'Re')

    !add others here, if required 
    call h5gclose_f(file%current_group, error)
    call h5gcreate_f(file%id,'grid',file%current_group, error)
    call write_hdf(file,time,'time')
    call write_hdf(file,i_start+i_time-1,'step')

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
    END DO

    call perfon('    h5_2')
    !write pressure field (collective I/O)
    CALL h5gcreate_f(file%id,'fields',file%current_group, error)
    call write_hdf(file,p_pbc,dimsf,'pressure')
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
    
  END SUBROUTINE hdf5output
#endif

  SUBROUTINE hdf5init
 
#ifdef HDF5IO    
    USE mod_hdf5io

    call init_io(comm,MPI_INFO_NULL,np,myid)
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
!  norm   (1/2) \int a.a dV  for each l,m                                                                                                                                                                                                                                 
!------------------------------------------------------------------------                                                                                                                                                                                                 
   subroutine norm2(vel,Ek,Em)
      complex(kind=8),dimension(m_r,mp_f),intent(in)  :: vel
      real(kind=8), intent(out) :: Ek(m_z), Em(m_th+1)
      real(kind=8)  :: Ek_(m_z), Em_(m_th+1)
      real(kind=8) :: w, b, f(m_r)
      integer(kind=4) :: i

      Ek = 0d0
      Em = 0d0
      w  = 4d0 * PI* PI / k_z0
      if(k_z0 == 0d0) stop 'k_z0 = 0 :cannot compute kinetic energy'
      do i=1,mp_f
         f =  w * vel(:,i)*conjg(vel(:,i))
         b = dot_product(f,intrdr)
         if (abs(fk_mp%th(1,i)) <= epsilon .and. abs(fk_mp%z(1,i)) <= epsilon) b = 0.5d0*b
         Em(int(fk_mp%th(1,i))+1)      = Em(int(fk_mp%th(1,i))+1)      + b
         Ek(abs(int(fk_mp%z(1,i)))+1) = Ek(abs(int(fk_mp%z(1,i)))+1) + b
      end do

      call mpi_allreduce( Ek, Ek_, m_z, mpi_real8, mpi_sum, comm, ierr)
      Ek = Ek_
      call mpi_allreduce( Em, Em_,m_th+1, mpi_real8, mpi_sum, comm, ierr)
      Em = Em_


   end subroutine norm2

!*************************************************************************        


!--------------------------------------------------------------------------                                                                                                                                                                                               
!  write:  1, time;  2,  bulk vel / excess pressure fraction if fixed flux;                                                                                                                                                                                               
!          3, mean uz at r=0;  4, friction velocity, 5 friction coefficient                                                                                                                                                                                                                      
!--------------------------------------------------------------------------                                                                                                                                                                                               
   subroutine friction()

     real(kind=8) :: Ub, Uc, Ufr, velz_real(m_r), cf

     
     velz_real= dreal(u_hat_mp%z(:,1))
     
     if(const_flux)  then
        Ub = vel_Pr0  
     else
        Ub = 2d0*dot_product(velz_real,intrdr)
     end if
     
     Uc = dot_product(velz_real(1:1+(n_s-1)/2),dr0(:,0))
     
     Ufr = dot_product(velz_real(m_r-(n_s-1)/2:m_r),dr1(:,1))
     Ufr = dsqrt(dabs(Ufr)/Re)
     
     
     if (const_flux)  then

        cf = 2d0*(Ufr*Ufr)/(0.25)
     else
        cf = 2d0*(Ufr*Ufr)/(Ub*Ub)
     end if
     if(myid .eq. 0) write(105,'(5e20.12)') time, Ub, Uc, Ufr, cf
   end subroutine friction


!--------------------------------------------------------------------------                                                                                                                                                                                               
!  save axial mean flow profile                                                                                                                                                                                                                                           
!--------------------------------------------------------------------------                                                                                                                                                                                               
   subroutine meanprof()

     character(8) :: cnum
     integer(kind=4) :: n

     if (myid .eq. 0) then
        write(cnum,'(I0.8)') i_time
        open(11, status='unknown', file='meanprof'//cnum//'.dat')
        write(11,*) '# t = ', time
        write(11,*) '# r  uz(r)'
        do n = 1, m_r
           write(11,'(2e20.12)')  r(n), dreal(u_hat_mp%z(n,1)) 
        end do
        close(11)
     end if
     
    

   end subroutine meanprof



   SUBROUTINE readscal(S, v_input, v_output, mp_kth)

  IMPLICIT NONE
  INTEGER(kind=4),intent(in) :: S
  COMPLEX(kind=8), intent(in), dimension(m_r_ifile,mp_f_ifile) :: v_input
  COMPLEX(kind=8), intent(out), dimension(m_r,mp_f) :: v_output
  real(kind=8),intent(in),dimension(m_r_ifile,mp_f_ifile)  :: mp_kth
  COMPLEX(kind=8),dimension(m_r,mp_f_ifile) :: v_int
  COMPLEX(kind=8),dimension(m_f_ifile,mp_r) :: v_aux
  COMPLEX(kind=8),dimension(m_f,mp_r) :: v_out_aux
  REAL(kind=8), allocatable, dimension(:) :: r_input
  REAL(kind=8) :: A(n_s,n_s,m_r)
  complex(kind=8) :: fn(1-(n_s-1)/2:m_r_ifile)
  real(kind=8) ::  fo(m_r,2)
  REAL(kind=8) :: par_real(1-(n_s-1)/2:m_r_ifile), par_imag(1-(n_s-1)/2:m_r_ifile)
  INTEGEr(kind=4) :: i
  
  if (m_r .ne. m_r_ifile) then

     allocate(r_input(1-(n_s-1)/2:m_r_ifile))
     call radial_grid(m_r_ifile,r_input(1:m_r_ifile))
     
     !extrapolate to r=0                                                                                                                                                                                                                                                  
     r_input(1-(n_s-1)/2:0) = -r_input((n_s-1)/2:1:-1)
     

     !Get interpolation weights                                                                                                                                                                                                                                           
     call interp_wts((n_s-1)/2+m_r_ifile, r_input, m_r ,r, A)

     !Interpolate to new grid                                                                                                                                                                                                                                             

     do i = 1,mp_f_ifile
        
        fn(1:m_r_ifile)=v_input(:,i)
        fn(:0) = fn((n_s-1)/2:1:-1)
        if(modulo(int(mp_kth(1,i))+S,2) .eq. 1)  fn(:0) = -fn(:0)
        par_real = dreal(fn)
        par_imag = dimag(fn)
        call interp((n_s-1)/2+m_r_ifile, r_input, par_real, A, m_r, r, fo(:,1))
        call interp((n_s-1)/2+m_r_ifile, r_input, par_imag, A, m_r, r, fo(:,2))
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
        
        CALL xTranspose_mpi(m_r, m_f_ifile, v_input , v_aux,ifile=.true.)
        
     else
        
        CALL xTranspose_mpi(m_r, m_f_ifile, v_int , v_aux,ifile=.true.)
        
     end if
     
     if ((m_th .lt. m_th_ifile) .and. (m_z .eq. m_z_ifile)) then
        

        do i = 1,m_z

           v_out_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th+1,:) = &
                v_aux((i-1)*(m_th_ifile+1)+1:(i-1)*(m_th_ifile+1)+m_th+1,:)

        end do

        v_out_aux(m_th+1:m_f:m_th+1,:)=dcmplx(0d0,0d0)
        CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output )

     elseif((m_th .gt. m_th_ifile) .and. (m_z .eq. m_z_ifile)) then


        do i = 1,m_z

           v_out_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th_ifile+1,:) = &
                v_aux((i-1)*(m_th_ifile+1)+1:(i-1)*(m_th_ifile+1)+m_th_ifile+1,:)
           
        end do
        
        v_out_aux(m_th+1:m_f:m_th+1,:)=dcmplx(0d0,0d0)
        CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output )

     elseif((m_th .eq. m_th_ifile) .and. (m_z .lt. m_z_ifile)) then

        
        do i = 2,m_z0
           
           v_out_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th+1,:) = &
                v_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th+1,:)
           v_out_aux((m_z+1-i)*(m_th+1)+1:(m_z+1-i)*(m_th+1)+m_th+1,:) = &
                v_aux((m_z_ifile+1-i)*(m_th+1)+1:(m_z_ifile+1-i)*(m_th+1)+m_th+1,:)
           
        end do
        
        !mode zero and mode m_z0

        v_out_aux(1:m_th+1,:) = v_aux(1:m_th+1,:)
        v_out_aux(m_z0*(m_th+1)+1:m_z0*(m_th+1)+m_th+1,:) = &
             v_aux(m_z0*(m_th+1)+1:m_z0*(m_th+1)+m_th+1,:)
        
       
        CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output )


     elseif((m_th .eq. m_th_ifile) .and. (m_z .gt. m_z_ifile)) then
        
        
        
        do i = 2,m_z_ifile/2

           v_out_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th+1,:) = &
                v_aux((i-1)*(m_th+1)+1:(i-1)*(m_th+1)+m_th+1,:)
           v_out_aux((m_z+1-i)*(m_th+1)+1:(m_z+1-i)*(m_th+1)+m_th+1,:) = &
                v_aux((m_z_ifile+1-i)*(m_th+1)+1:(m_z_ifile+1-i)*(m_th+1)+m_th+1,:)
           

        end do

        
        !mode zero and mode m_z0

        
        v_out_aux(1:m_th+1,:) = v_aux(1:m_th+1,:)
        v_out_aux(m_z_ifile/2*(m_th+1)+1:m_z_ifile/2*(m_th+1)+m_th+1,:) = &
             v_aux(m_z_ifile/2*(m_th+1)+1:m_z_ifile/2*(m_th+1)+m_th+1,:)
        
      
        CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output )

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
        
        CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output )
        
        
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
        
        CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output )

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
                
        CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output )
         


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
              
        
        CALL xTranspose_mpi(m_f,m_r, v_out_aux , v_output )
       

     end if
  end if

  

END SUBROUTINE readscal
                                                                   

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
         A(i,j,n) = A(i,j-1,n) * (r_input(left+i-1)-r(n)) / real(j-1)
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

!------------------------------------------------------------------------                                                                                                                                                                                                 
!  Replace nxn matrix A by its inverse                                                                                                                                                                                                                                    
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
    


      SUBROUTINE CHEBY (M, ND, X, NP, Y)                                                  
!     Last change: 20 Jan 1990                                                            
!     Purpose    : to find the values Y(1),...Y(NP)                                       
!                  of the NDthe derivative of the Cheby poly                              
!                  of degree M at the points                                              
!                  X(1),...X(NP)                                                          
      INTEGER ONE,M,NP,IA,MM,J,ND,NJ                                                      
      DOUBLE PRECISION A(500),AX(500),X(NP),Y(NP)                                         
      IA = 500                                                                            
      MM = M + 1                                                                          
      ONE = 1                                                                             
      IF (MM .GT. IA) GO TO 70                                                            
      IF (M .EQ. 0) GO TO 50                                                              
      DO 10 J = 1, M                                                                      
        A(J) = 0D0                                                                        
   10 CONTINUE                                                                            
      A(MM) = 1D0                                                                         
      IF (ND .EQ. 0) GO TO 40                                                             
      DO 30 NJ = 1, ND                                                                    
        CALL DERIV2(A, AX, MM)                                                            
        DO 20 J = 1, MM                                                                   
          A(J) = AX(J)                                                                    
   20   CONTINUE                                                                          
   30 CONTINUE                                                                            
   40 CALL EVAL2(A, MM, X, NP, Y)                                                         
      RETURN                                                                              
   50 DO 60 J = 1, NP                                                                     
        IF (ND .EQ. 0) Y(J) = 1D0                                                         
        IF (ND .GT. 0) Y(J) = 0D0                                                         
   60 CONTINUE                                                                            
      RETURN                                                                              
   70 WRITE (6,80)                                                                        
   80 FORMAT (' Sub CHEBY, M too big')                                                    
      RETURN                                                                              
      END                                                                                 
!                                                                                         
      SUBROUTINE EVAL2 (A, M, X, NX, Y)                                                   
!     Last change: 20 Jan 1990                                                            
!     Purpose    : to compute the value of the Cheby series                               
!                  A(1),...A(M) of length M, where A(1)....A(M)                           
!                  are respectively the coeffs of T(0)...T(M-1),                          
!                  at the points X(j) j=1,...NX.                                          
!                  The output is Y(j) j=1,...NX, which is the                             
!                  value of the series at these NX points.                                
      DOUBLE PRECISION A(M),X(NX),Y(NX),B(0:502),XA,XB,XX                                 
      INTEGER M,NX,IA,M1,M2,J,I                                                           
      IA = 502                                                                            
      M1 = M + 1                                                                          
      M2 = M + 2                                                                          
      IF (M2 .GT. IA) GO TO 30                                                            
      DO 20 I = 1, NX                                                                     
        XX = X(I)                                                                         
        XA = 2D0 * XX - 1D0                                                               
        XB = 2D0 * XA                                                                     
        B(M2) = 0D0                                                                       
        B(M1) = 0D0                                                                       
        DO 10 J = M - 1, 0, -1                                                            
          B(J + 1) = XB * B(J + 2) - B(J + 3) + A(J + 1)                                  
   10   CONTINUE                                                                          
        Y(I) = B(1) - XA * B(2)                                                           
   20 CONTINUE                                                                            
      RETURN                                                                              
   30 WRITE (6,40)                                                                        
   40 FORMAT (' M too big in Sub EVAL2')                                                  
      RETURN                                                                              
      END                                                                                 
!                                                                                         
      SUBROUTINE DERIV2 (A, AX, M)                                                        
!     Last change: 20 Jan 1990                                                            
!     Purpose    :   given the Cheby series A(1),...A(M)                                  
!                    which represents the function                                        
!                    f=Sum(j=1..M) A(j)*T(j-1)                                            
!                    the subroutines computes the coeffs                                  
!                    AX(1).......AX(M) of the function df/dx                              
!                    Note that AX(M)=0                                                    
!                                                                                         
      DOUBLE PRECISION A(M),AX(M)                                                         
      INTEGER          M,J                                                                
      AX(M) = 0D0                                                                         
      IF (M .LE. 0) GO TO 70                                                              
      IF (M .GT. 4) GO TO 50                                                              
      IF (M .EQ. 1) GO TO 10                                                              
      IF (M .EQ. 2) GO TO 20                                                              
      IF (M .EQ. 3) GO TO 30                                                              
      IF (M .EQ. 4) GO TO 40                                                              
   10 AX(1) = 0D0                                                                         
      RETURN                                                                              
   20 AX(1) = 2D0 * A(2)                                                                  
      RETURN                                                                              
   30 AX(1) = 2D0 * A(2)                                                                  
      AX(2) = 8D0 * A(3)                                                                  
      RETURN                                                                              
   40 AX(1) = 2D0 * A(2) + 6D0 * A(4)                                                     
      AX(2) = 8D0 * A(3)                                                                  
      AX(3) = 12D0 * A(4)                                                                 
      RETURN                                                                              
   50 AX(M - 1) = - 4D0 * (1D0 - M) * A(M)                                                
      AX(M - 2) = - 4D0 * (2D0 - M) * A(M - 1)                                            
      DO 60 J = M - 3, 2, -1                                                              
        AX(J) = 4D0 * J * A(J + 1) + AX(J + 2)                                            
   60 CONTINUE                                                                            
      AX(1) = 2D0 * A(2) + AX(3) / 2D0                                                    
      RETURN                                                                              
   70 WRITE (6,80)                                                                        
   80 FORMAT (' M .LE. 0 in Sub DERIV2')                                                  
      RETURN                                                                              
      END                                                                                 
!                                                                                         
      SUBROUTINE CMESH (N, X, NKIND)                                                      
!     Last changed: 20 Jan 1990                                                           
!     Purpose :     Chebyshev mesh:                                                       
!                   it finds N points X(j) j=1,...N which are the                         
!                   zeroes of the Cheby polynomial of first kind Tn                       
!                   or of the Cheby polynomial of second kind Un,                         
!                   over the interval 0 < x < 1                                           
      DOUBLE PRECISION X(N),PI,y                                                          
      INTEGER N,IV,I,NKIND                                                                
      PI = 4D0 * DATAN (1D0)                                                              
      IF (NKIND .EQ. 1) GO TO 10                                                          
      IF (NKIND .EQ. 2) GO TO 30                                                          
!     Zeroes of the Chebyshev polynomials of the first kind                               
   10 DO 20 I = 1, N                                                                      
        IV = N + 1 - I                                                                    
        X(I) = 0.5D0 * (1D0+DCOS((2*IV - 1)*PI/(2D0*N)))                                  
   20 CONTINUE                                                                            
      RETURN                                                                              
!     Zeroes of the Chebyshev polynomials of the second kind                              
   30 DO 40 I = 1, N                                                                      
        Y = DCOS(PI*I/(N + 1))                                                            
        X(I) = 0.5D0 * (Y + 1D0)                                                          
   40 CONTINUE                                                                            
      RETURN                                                                              
   END                                                                                    
                                                                                          
                                                                                          
  subroutine spectrum                                                                     
                                                                                          
    IMPLICIT NONE                                                                         
                                                                                          
                                                                                          
    REAL(kind=8),allocatable,SAVE :: TM(:,:)                                              
    REAL(kind=8) :: r_(m_r), z_(0:m_z0-1), th_(0:m_th),x(m_r)                             
    REAL(kind=8) :: r__(m_r), z__(0:m_z0-1), th__(0:m_th)                                 
    integer(kind=4) :: i, n,kp,m                                                          
    REAL(kind=8) :: d(m_r), dRe(m_r), dIm(m_r)                                            
    character(8) :: cnum                                                                  
                                                                                          
10  format(i4,1e20.12)                                                                    
    !Compute chebyshev tranformation matrix (only the first time the subroutine is called)
                                                                                          
    if(.not.allocated(TM)) then                                                           
       ALLOCATE(TM(m_r,m_r))                                                              
       do n = 0, m_r-1                                                                    
          x(n+1) = 0.5d0 * ( 1d0 + dcos(PI*(m_r-n)/dble(m_r)) )                           
       end do                                                                             
       do n = 1, m_r                                                                      
          call cheby(n-1, 0, x, m_r, TM(1,n))                                             
       end do                                                                             
       call mat_inv(m_r,TM,m_r)                                                           
       TM = transpose(TM)                                                                 
    end if                                                                                
                                                                                          
    r_ = 0d0                                                                              
    z_ = 0d0                                                                              
    th_ = 0d0                                                                             
                                                                                          
    do i=1,mp_f                                                                           
       dRe = matmul(dreal(u_hat_mp%r(:,i)), TM)                                           
       dIm = matmul(dimag(u_hat_mp%r(:,i)), TM)                                           
       d = dRe*dRe+dIm*dIm                                                                
       dRe = matmul(dreal(u_hat_mp%th(:,i)), TM)                                          
       dIm = matmul(dimag(u_hat_mp%th(:,i)), TM)                                          
       d = max(d, dRe*dRe+dIm*dIm)                                                        
       dRe = matmul(dreal(u_hat_mp%z(:,i)), TM)                                           
       dIm = matmul(dimag(u_hat_mp%z(:,i)), TM)                                           
       d = max(d, dRe*dRe+dIm*dIm)                                                        
       d = dsqrt(d)                                                                       
       kp = int(abs(fk_mp%z(1,i))/k_z0)                                                   
       m = int(fk_mp%th(1,i)/k_th0)                                                       
       do n = 1, m_r                                                                      
          r_(n)  = max(d(n), r_(n))                                                       
          z_(kp) = max(d(n), z_(kp))                                                      
          th_(m)  = max(d(n), th_(m))                                                     
       end do                                                                             
    end do                                                                                
                                                                                          
                                                                                          
    call mpi_allreduce(r_, r__, m_r, MPI_REAL8,  &                                        
         MPI_MAX, COMM, ierr)                                                             
    r_ = r__                                                                              
    call mpi_allreduce(z_, z__, m_z0, MPI_REAL8,  &                                       
         MPI_MAX, COMM, ierr)                                                             
    z_ = z__                                                                              
    call mpi_allreduce(th_, th__, m_th+1, MPI_REAL8,  &                                   
         MPI_MAX, COMM, ierr)                                                             
    th_ = th__                                                                            
                                                                                          
    if (myid .eq. 0) then                                                                 
       write(cnum,'(I0.8)') i_time                                                        
       open(11, status='unknown', file='velocity_spectrum'//cnum//'.dat')                 
       write(11,*) '# t = ', time                                                         
       write(11,*) '# r'                                                                  
       do i = 1, m_r                                                                      
          write(11,10) i, r_(i)/maxval(r_)                                                
       end do                                                                             
       write(11,*) '&'                                                                    
       write(11,*) '# z'                                                                  
       do i = 0, m_z0-1                                                                   
          write(11,10) i, z_(i)/maxval(z_)                                                
       end do                                                                             
       write(11,*) '&'                                                                    
       write(11,*) '# th, k_th0=',k_th0                                                   
       do i = 0, m_th-1                                                                   
          write(11,10) i, th_(i)/maxval(th_)                                              
       end do                                                                             
       write(11,*) '&'                                                                    
       close(11)                                                                          
                                                                                          
    end if                                                                                
                                                                                          
                                                                                          
  end subroutine spectrum

  !-------------end module------------------

END MODULE mod_inOut
