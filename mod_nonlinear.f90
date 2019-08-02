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

MODULE mod_nonlinear

!   nonlinear term (udu) in NS equation
!                by
!   pseudospectral method with 3/2 rule

  USE mod_myMpi
  USE mod_fftw
  USE mod_fdInit
  USE mod_inOut
  IMPLICIT NONE

  
  
CONTAINS
!----------------------------------------------------------------------
  SUBROUTINE nonlinear(output, cfl_eval)
    
    IMPLICIT NONE

    logical, intent(in) :: output,cfl_eval
    TYPE (vec_mpi)                           :: drdu_hat_mp
    COMPLEX(KIND=8), DIMENSION(m_f,mp_r)     :: var_Cr,var_Cth,var_Cz,var_Cp
    COMPLEX(KIND=8), DIMENSION(m_f,mp_r)     :: dvar_Cr,dvar_Cth,dvar_Cz
    REAL(KIND=8),    DIMENSION(9*n_f/4,mp_r) :: u_r,u_th,u_z,u_thOr,p
    REAL(KIND=8),    DIMENSION(9*n_f/4,mp_r) :: dthdu_r,dthdu_th,dthdu_z
    REAL(KIND=8),    DIMENSION(9*n_f/4,mp_r) :: dzdu_r,dzdu_th,dzdu_z
    REAL(KIND=8),    DIMENSION(9*n_f/4,mp_r) :: drdu_r,drdu_th,drdu_z
    REAL(KIND=8),    DIMENSION(9*n_f/4,mp_r) :: udu_th,udu_z,udu_r
    INTEGER(KIND=4)                          :: i,j,iloop
    real(kind=8)                             :: vel_z_pert(m_r)

    call perfon('  nonlin')

    ierr=drdu_hat_mp%alloc(m_r,mp_f)
    
    call perfon('   matvec')
    
!$OMP PARALLEL DO &      
!$OMP  SCHEDULE (RUNTIME) &
!$OMP  DEFAULT(NONE) &
!$OMP  PRIVATE(i) &
!$OMP  SHARED(mp_f,drdu_hat_mp,bMw_dr,u_hat_mp)
    DO i = 1,mp_f
       call bmatvec(bMw_dr,u_hat_mp%r(:,i), drdu_hat_mp%r(:,i))
       call bmatvec(bMw_dr,u_hat_mp%th(:,i),drdu_hat_mp%th(:,i))
       call bmatvec(bMw_dr,u_hat_mp%z(:,i), drdu_hat_mp%z(:,i))
    END DO

    call perfoff('   matvec')
    
    call perfon('tranf')                                                


    
    ! Transpose velocities to perform Fourier transformation                          
    CALL xTranspose_mpi(m_r,m_f,u_hat_mp%r,var_Cr )                         
    CALL xTranspose_mpi(m_r,m_f,u_hat_mp%th,var_Cth )
    CALL xTranspose_mpi(m_r,m_f,u_hat_mp%z,var_Cz )                         
                      
    if ((MOD(i_time,dn_hdf5) == 0) .and. (output)) then     
       CALL xTranspose_mpi(m_r,m_f,p_hat_mp,var_Cp )                        
    endif                                                                  
    call perfoff()                                                         


  
    if (tl_mpi .lt. MPI_THREAD_SERIALIZED) then
       !fallback if threaded MPI is not available
       call perfon('   tranf')
       CALL xTranspose_mpi(m_r,m_f,drdu_hat_mp%r,dvar_Cr )
       CALL xTranspose_mpi(m_r,m_f,drdu_hat_mp%th,dvar_Cth )
       CALL xTranspose_mpi(m_r,m_f,drdu_hat_mp%z,dvar_Cz )
       call perfoff()
    endif

    call perfon('   deriv')

!$OMP PARALLEL &      
!$OMP  DEFAULT(SHARED) 
!$OMP MASTER
    if (tl_mpi .ge. MPI_THREAD_SERIALIZED) then
       CALL xTranspose_mpi(m_r,m_f,drdu_hat_mp%r,dvar_Cr )
       CALL xTranspose_mpi(m_r,m_f,drdu_hat_mp%th,dvar_Cth )
       CALL xTranspose_mpi(m_r,m_f,drdu_hat_mp%z,dvar_Cz )
    endif
!$OMP END MASTER

!$OMP SINGLE
    DO i = 1,mp_r

!$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(i)
       CALL bwd_fft(var_Cr(1,i),u_r(1,i))
!$OMP END TASK


!$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(i)
       CALL bwd_fft(var_Cth(1,i),u_th(1,i))
       u_thOr(:,i) = u_th(:,i)/r(i+myid*mp_r)
!$OMP END TASK

!$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(i)
       CALL bwd_fft(var_Cz(1,i),u_z(1,i))
!$OMP END TASK

!$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(i)
       CALL deriv_x(k_th0,var_Cr(1,i),dthdu_r(1,i))
!$OMP END TASK

!$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(i)
       CALL deriv_y(k_z0,var_Cr(1,i),dzdu_r(1,i))
!$OMP END TASK

!$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(i)
       CALL deriv_x(k_th0,var_Cth(1,i),dthdu_th(1,i))
!$OMP END TASK

!$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(i)
       CALL deriv_y(k_z0,var_Cth(1,i),dzdu_th(1,i))
!$OMP END TASK

!$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(i)
       CALL deriv_x(k_th0,var_Cz(1,i),dthdu_z(1,i))
!$OMP END TASK

!$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(i)
       CALL deriv_y(k_z0,var_Cz(1,i),dzdu_z(1,i))
!$OMP END TASK

          END DO
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL

          
                    
    call perfoff('   deriv')
          

    call perfon('   fft')



    !$OMP PARALLEL DO &
    !$OMP  DEFAULT(SHARED) &
    !$OMP  PRIVATE(i,iloop) &
    !$OMP  SCHEDULE (STATIC)
    do iloop = 0,3*mp_r-1
       i=mod(iloop,mp_r)+1
       
       select case (iloop/mp_r)
          
       case(0)
        
          CALL bwd_fft(dvar_Cr(1,i),drdu_r(1,i))
          udu_r(:,i)  = u_r(:,i)*drdu_r(:,i)  + u_thOr(:,i)*dthdu_r(:,i)&
               + u_z(:,i)*dzdu_r(:,i)  - u_th(:,i)*u_thOr(:,i)
          CALL fwd_fft(udu_r(1,i),var_Cr(1,i))
          
       case(1)
          
          CALL bwd_fft(dvar_Cth(1,i),drdu_th(1,i))
          
          udu_th(:,i) = u_r(:,i)*drdu_th(:,i) + u_thOr(:,i)*dthdu_th(:,i)&
               + u_z(:,i)*dzdu_th(:,i) + u_r(:,i)*u_thOr(:,i)
          
          CALL fwd_fft(udu_th(1,i),var_Cth(1,i))
          
       case(2)
          
          CALL bwd_fft(dvar_Cz(1,i),drdu_z(1,i))
          udu_z(:,i)  = u_r(:,i)*drdu_z(:,i)  + u_thOr(:,i)*dthdu_z(:,i)&
               + u_z(:,i)*dzdu_z(:,i)
          CALL fwd_fft(udu_z(1,i),var_Cz(1,i))
       case default
          print *,i,iloop/mp_r
          stop 'this should never happen'
       end select

    END DO
    !$OMP END PARALLEL DO  
         
    
    call perfoff('   fft')
    
!HDF5 output
#ifdef HDF5IO    
    IF((MOD(i_time,dn_hdf5) == 0) .and. (output)) THEN
       call perfon('   hdf5')

!$OMP PARALLEL DO &
!$OMP  DEFAULT(NONE) &
!$OMP  SHARED(var_Cp,p,mp_r) &
!$OMP  PRIVATE(i) &
!$OMP  SCHEDULE (STATIC)
       do i = 1,mp_r
          call bwd_fft(var_Cp(1,i),p(1,i))
       enddo
       
       call hdf5output(p,u_r,u_th,u_z)

       call perfoff()
    ENDIF
#endif /* HDF5IO */    

    call perfon('   tranb')

   
    ! Tranpose back to solve the linear algebraic equations
    CALL xTranspose_mpi(m_f,m_r,var_Cr,udu_hat_mp%r )
    CALL xTranspose_mpi(m_f,m_r,var_Cth,udu_hat_mp%th )
    CALL xTranspose_mpi(m_f,m_r,var_Cz,udu_hat_mp%z )

    call perfoff('   tranb')
    
    !change sign (rhs of equation)
    udu_hat_mp%r=-1d0*udu_hat_mp%r
    udu_hat_mp%th=-1d0*udu_hat_mp%th
    udu_hat_mp%z=-1d0*udu_hat_mp%z

    ierr=drdu_hat_mp%dealloc()
    call perfoff('  nonlin')

    !Introducing forcing (either constant mass flux or constant pressure) 
    
    
    if (const_flux) then 

       if (myid .eq. 0) then
           
          
          !Compute friction and adjust flux to constant mass flow
          vel_z_pert = dreal(u_hat_mp%z(:,1)) - (1d0-r*r)
          vel_Pr0 = -0.5d0*dot_product(vel_z_pert(m_r-(n_s-1)/2:m_r),dr1(:,1))
          udu_hat_mp%z(:,1)=udu_hat_mp%z(:,1)+ dcmplx(4d0*(1d0+vel_Pr0)/Re,0d0)
          
       end if
       
    else
       !Add extra term (fixed Pressure) in z equation 
       if(myid .eq. root) udu_hat_mp%z(:,1)=udu_hat_mp%z(:,1)+dcmplx((4d0/Re),0d0)
    end if

    
 
    !Evaluate CFL if needed and change timestep size                                                                                                                                                                                                                 

    IF ( cfl_eval ) call maxstep(u_r,u_th,u_z)


  END SUBROUTINE nonlinear


!----------------------------------------------------------------------
  SUBROUTINE bmatvec(A,x,y)
    !--------------------------------------------
    ! Matrix-Vector Multiplication: y = Ax
    ! x,y : double complex vectors
    ! A:    real matrix in BLAS-like band storage 
    !--------------------------------------------
    IMPLICIT NONE
    REAL(KIND=8),    INTENT(IN)        :: A(:,:)
    COMPLEX(KIND=8), INTENT(IN)        :: x(:)
    COMPLEX(KIND=8), INTENT(OUT)       :: y(:)
    REAL(KIND=8), DIMENSION(SIZE(x),2) :: par_o,par_i
    INTEGER(KIND=4)  :: n,iw,lda
    
    n=SIZE(A,dim=2)
    lda=SIZE(A,dim=1)
    iw=(lda-1)/2
    par_i(:,1) = DREAL(x)
    par_i(:,2) = DIMAG(x)

    call dgbmv('N',n,n,iw,iw,1.d0,A,lda,par_i(:,1),1,0.d0,par_o(:,1),1)
    call dgbmv('N',n,n,iw,iw,1.d0,A,lda,par_i(:,2),1,0.d0,par_o(:,2),1)

    y(:) = DCMPLX(par_o(:,1),par_o(:,2))

  END SUBROUTINE bmatvec

!------------------------------------------------------------------------                                              
!  Obtain Courant-Friedrich-Lewis condition                                                                            
!------------------------------------------------------------------------                                              
  subroutine maxstep(u_r,u_th,u_z)                                                                                     
                                                                                                                       
    REAL(kind=8),dimension(9*n_f/4,mp_r),intent(in)  :: u_r,u_th,u_z                                                   
    INTEGER :: nthfine, nzfine                                                           
    REAL(KIND=8) :: dthfine, dzfine                                           
    real(kind=8)  :: drfine, mx, deltat(3),deltat_(3)
    real(kind=8), allocatable :: r_(:)
    integer(kind=4) :: N_,n, n__                                                                                       

    nthfine = 3*n_th/2
    nzfine = 3*n_z/2
    dthfine = len_th/nthfine
    dzfine = len_z/nzfine
    allocate(r_(0:mp_r+1))
    
    
    N_ = mp_r                                                                                                          
    do n = 0, N_+1                                                                                                     
       n__ = n+myid*mp_r                                                                                               
       if(n__>= 1 .and. n__<= m_r)  r_(n) = r(n__)                                                                     
    end do                                                                                                             
                                                                                                                       
    deltat = 1d99                                                                                                      
                                                                                                                       
    do n = 1, N_                                                                                                       
       
       
       !radial direction                                                                                               
       if(n .eq. 1 .and. myid .eq. 0) then                                                                             
          drfine = r_(2) - r_(1)                                                                                       
       else if(n .eq. N_ .and. myid .eq. numprocs-1) then                                                              
          drfine = r_(N_) - r_(N_-1)                                                                                   
       else                                                                                                            
          drfine = min( r_(n)-r_(n-1), r_(n+1)-r_(n) )                                                                 
       end if
       mx = maxval( dabs(u_r(:,n)) )                                                                                   
       if (mx .ne. 0d0) deltat(1) = min( deltat(1), drfine/mx )                                                        
       
       !azimuthal direction                                                                                            
       mx = maxval( dabs(u_th(:,n)) )                                                                                  
       if(mx .ne. 0d0) deltat(2) = min( deltat(2), dthfine/mx )                                                        
       
       !Axial direction                                                                                                
       
       mx = maxval( dabs(u_z(:,n))) !+ (1d0- r_(n)*r_(n))) ! (This is only if the parabolic profile is imposed)           
       if(mx .ne. 0d0) deltat(3) = min( deltat(3), dzfine/mx )                                                         
    end do
                                                                                                                       
    
    call mpi_allreduce(deltat, deltat_, 3, MPI_REAL8,  &                                                               
         MPI_MIN, COMM, ierr)                                                                                          
                                                                                                                       
    deltat = deltat_                                                                                                   
                                                                                                                       
    cfl = minval(deltat)                                                                                               
    if(cfl .eq. deltat(1))  cfl_dir=1                                                                                  
    if(cfl .eq. deltat(2))  cfl_dir=2                                                                                  
    if(cfl .eq. deltat(3))  cfl_dir=3                                                                                  

    deallocate(r_)
    
  end subroutine maxstep
  

END MODULE mod_nonlinear

