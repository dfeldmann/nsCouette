MODULE mod_preAnalysis

  USE mod_myMpi
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
  TYPE(spec), DIMENSION(:,:), allocatable :: f_hat_mp,f_hat_mp_dt

  CHARACTER*40    :: fName_dt     ! coefficient file name at t+dt
  REAL(KIND=8)    :: delta_t      ! time shift of two coeff files
  
  private :: spec,f_hat_mp,f_hat_mp_dt

CONTAINS
  !--------------------------------------------------------------
  ! Read the variable parameters from the input file/Standard IO
  SUBROUTINE read_pars()
    IMPLICIT NONE
    INTEGER(KIND=4) :: i

    IF (myid == root) THEN
       READ*,m_r
       READ*,m_th
       READ*,m_z0
       READ*,k_th0
       READ*,k_z0
       READ*,eta

       READ*, Re_i
       READ*, Re_o
      !mu = Re_o/Re_i
       
       READ*, fName_ic
       READ*, fName_dt
       i=INDEX(fName_ic//' ',' ')-1
    END IF

    CALL MPI_Bcast(m_r, 1,MPI_INTEGER,root,comm,ierr)
    CALL MPI_Bcast(m_th,1,MPI_INTEGER,root,comm,ierr)
    CALL MPI_Bcast(m_z0,1,MPI_INTEGER,root,comm,ierr)
    CALL MPI_Bcast(k_th0,1,MPI_REAL8,root,comm,ierr)
    CALL MPI_Bcast(k_z0,1,MPI_REAL8,root,comm,ierr)
    CALL MPI_Bcast(eta,1,MPI_REAL8,root,comm,ierr)

    CALL MPI_Bcast(Re_i,1,MPI_REAL8,root,comm,ierr)
    CALL MPI_Bcast(Re_o,1,MPI_REAL8,root,comm,ierr)
   !CALL MPI_Bcast(mu,1,MPI_REAL8,root,comm,ierr)
    CALL MPI_Bcast(delta_t,1,MPI_REAL8,root,comm,ierr)
    CALL MPI_Bcast(i,1,MPI_INTEGER,root,comm,ierr)
    CALL MPI_Bcast(fName_ic,i,MPI_CHARACTER,root,comm,ierr)
    CALL MPI_Bcast(fName_dt,i,MPI_CHARACTER,root,comm,ierr)

  END SUBROUTINE read_pars

  !-----------------------------------------
  SUBROUTINE read_coeff()
    IMPLICIT NONE

    INTEGER(KIND=4) :: fh,fh2,status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_OFFSET_KIND) :: disp=0
    INTEGER(KIND=4) :: i
    INTEGER(KIND=4) :: i_time,m_r_ifile,m_th_ifile,m_z_ifile
    REAL(KIND=8)    :: time,dt
    CHARACTER(40)   :: fBase_ic 
    NAMELIST /parameters_restart/ i_time,time,dt,m_r_ifile,m_th_ifile,m_z_ifile,fBase_ic

    allocate(f_hat_mp(m_r,mp_f))
    allocate(f_hat_mp_dt(m_r,mp_f))


    i = INDEX(fName_ic//' ',' ') - 1
    counter = m_r*mp_f
    ! Read the first coeff file
    CALL MPI_File_open(comm,fName_ic(1:i),MPI_MODE_RDONLY,&
         MPI_INFO_NULL,fh,ierr)
    if (ierr .ne. 0) then
       PRINT*, fName_ic(1:i)//': error opening file'
       CALL MPI_Finalize(ierr)
       stop
    endif

    CALL MPI_File_set_view(fh,disp,mpi_spec,filetype,"native",&
         MPI_INFO_NULL,ierr)
    CALL MPI_File_read(fh,f_hat_mp(1,1),counter,mpi_spec,status,ierr)

    ! Read the second coeff file
    CALL MPI_File_open(comm,fName_dt(1:i),MPI_MODE_RDONLY,&
         MPI_INFO_NULL,fh2,ierr)
    if (ierr .ne. 0) then
       PRINT*, fName_dt(1:i)//': error opening file'
       CALL MPI_Finalize(ierr)
       stop
    endif
    CALL MPI_File_set_view(fh2,disp,mpi_spec,filetype,"native",&
         MPI_INFO_NULL,ierr)
    CALL MPI_File_read(fh2,f_hat_mp_dt(1,1),counter,mpi_spec,status,ierr)

    ! Close both coeff files
    CALL MPI_File_close(fh,ierr)
    CALL MPI_File_close(fh2,ierr)


    !read metadata
    if (myid == root) then
       open(unit=107,file=trim(fName_ic)//'.info')
       read(107,nml=parameters_restart)
       close(107)
       delta_t=time

       open(unit=107,file=trim(fName_dt)//'.info')
       read(107,nml=parameters_restart)
       close(107)
       delta_t=time-delta_t
    endif
    CALL MPI_Bcast(delta_t,1,MPI_REAL8,root,comm,ierr)

    print *,'# ',trim(fName_ic),' -- ',trim(fName_dt)
    print *,'# delta_t:',delta_t

  END SUBROUTINE read_coeff

  !----------------------------------------
  SUBROUTINE norm_du(phi,fPhi)
    !--------------------------------
    ! evaluate the f(phi)=sqrt{|int_{volume} {u(t)-exp(i*phi)*u(t+dt)}|}
    ! but in spectral space
    !--------------------------------
    USE mod_fdInit
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: phi
    REAL(KIND=8), INTENT(OUT) :: fPhi

    COMPLEX(KIND=8) :: eiPhi ! exp(i*phi)
    TYPE(spec), DIMENSION(m_r,mp_f)  :: fPhi_hat_mp
    REAL(KIND=8)    :: eneg_pp(m_r,mp_f)
    REAL(KIND=8)    :: eneg_mode_mp(mp_f),eneg_mode(m_th+1,m_z)
    REAL(KIND=8)    :: eneg_th(m_th+1),eneg_z(m_z), eneg
    INTEGER(KIND=4) :: h,i,j
    REAL(KIND=8)    :: row_inv_Mw_dr(m_r),inv_Mw_dr(m_r,m_r)
    REAL(KIND=8)    :: bA(n_s,m_r),bMw_dr(n_s,m_r),bMw_drr(n_s,m_r)
    REAL(KIND=8)    :: dr1((n_s-1)/2+1,0:(n_s-1)/2),dr2((n_s-1)/2+1,0:(n_s-1)/2),intrdr(m_r)

    INTEGER(KIND=4), PARAMETER :: iwidth=(n_s-1)/2, idiag=iwidth+1


    ! Calculate \hat{f(phi)} in nodes
    DO j = 1,mp_f
       eiPhi = DCMPLX(DCOS(f_hat_mp(1,j)%k_th*phi),DSIN(f_hat_mp(1,j)%k_th*phi))
       DO i = 1,m_r
          fPhi_hat_mp(i,j)%ur = f_hat_mp(i,j)%ur - eiPhi*f_hat_mp_dt(i,j)%ur
          fPhi_hat_mp(i,j)%uth = f_hat_mp(i,j)%uth - eiPhi*f_hat_mp_dt(i,j)%uth
          fPhi_hat_mp(i,j)%uz = f_hat_mp(i,j)%uz - eiPhi*f_hat_mp_dt(i,j)%uz
       END DO
    END DO

    !eneg_pp = fPhi_hat_mp*CONJG(fPhi_hat_mp)
    eneg_pp = fPhi_hat_mp%ur*CONJG(fPhi_hat_mp%ur) &
         + fPhi_hat_mp%uth*CONJG(fPhi_hat_mp%uth) &
         + fPhi_hat_mp%uz*CONJG(fPhi_hat_mp%uz)

    ! multiply the factor rdr
    DO i = 1,m_r
       eneg_pp(i,:) = eneg_pp(i,:)*r(i)
    END DO

    ! Radial derivatives 
    Call fd_Matrix(n_s-1,m_r-1,r,bMw_dr,bMw_drr,dr1,dr2,intrdr)
       

    ! Precompute the inverse of bMw_dr (-> row_inv_Mw_dr in output_energy)
    bA(:,:)=bMw_dr
    do h=1,idiag
       bA(idiag-h+1,h) = 0.d0     
    enddo
    bA(idiag,1) = 1.d0                 
    call bmatinv(m_r,(n_s-1)/2,bA,inv_Mw_dr)
    row_inv_Mw_dr(:)=inv_Mw_dr(m_r,:)


    ! multiply the factor 1/(2*Vol) in front of integral
    eneg_pp = eneg_pp*(1.D0-eta)/(1.D0+eta)
    eneg_pp(1,:) = 0.d0
    call DGEMV('T',m_r,mp_f,1.d0,eneg_pp,m_r,row_inv_Mw_dr,1,0.d0,eneg_mode_mp,1)
    CALL MPI_Gather(eneg_mode_mp,mp_f,MPI_REAL8,eneg_mode,mp_f,&
         MPI_REAL8,root,comm,ierr)


    IF (myid == root) THEN
       eneg_th = SUM(eneg_mode,2)
       eneg_z(:) = 2*SUM(eneg_mode,1) - eneg_mode(1,:)
       eneg = SUM(eneg_z(:))
       fPhi = SQRT(eneg)
    END IF

  END SUBROUTINE norm_du

  !---------------------------------------------------------------------
  subroutine bmatinv(n,k,A,B)
  !---------------------------------------------------------------------
  ! Compute inverse B of a real nxn band matrix A 
  !
  ! input:  A   banded matrix of dimension n (LAPACK storage)
  ! input:  n   dimension of A
  ! input:  k   number of sub/super-diagonals within the band of A
  ! output: B   inverse of B in regular storage
  !---------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN)  :: n,k
    REAL(KIND=8), INTENT(IN)     :: A(2*k+1,n)
    REAL(KIND=8), INTENT(OUT)    :: B(n,n)
    INTEGER(KIND=4)              :: i,j

    type LUdcmp
       sequence
       real(kind=8), allocatable :: LU(:,:)
       integer(kind=4), allocatable :: ipiv(:)
       integer(kind=4) :: n,k,lda,ldb,info
    end type LUdcmp
    TYPE(LUdcmp) :: LU

    ! identity matrix
    B(:,:) = 0.d0
    do i=1,n
       B(i,i)=1.d0
    enddo

    allocate(LU%LU(3*k+1,n),LU%ipiv(n))

    LU%n=n 
    LU%k=k 
    LU%lda=3*k +1
    LU%ldb=n

    DO j = 1,n
       DO i = 1,2*k+1
          LU%LU(k+i,j) = A(i,j)
       END DO
    END DO

    CALL dgbtrf(LU%n,LU%n,LU%k,LU%k,LU%LU,LU%lda,LU%ipiv,LU%info)
    if (LU%info /= 0) STOP 'bmatinv: dgbtrf'

    CALL dgbtrs('N',LU%n,LU%k,LU%k,n,LU%LU,LU%lda,LU%ipiv,B,LU%ldb,LU%info)
    if (LU%info /= 0) STOP 'bmatinv: dgbtrs'

    deallocate(LU%LU,LU%ipiv)

    RETURN
  end subroutine bmatinv

END MODULE mod_preAnalysis
