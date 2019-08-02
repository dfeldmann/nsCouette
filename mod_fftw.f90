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

MODULE mod_fftw
! Do 2d-DFT by using fftw including dealiased '3/2' rule:
!       REAL var_R(m,n) <-> COMPLEX var_C(m/2+1,n)
!
! SVN $HeadURL$
! SVN $Id$
! SVN $LastChangedDate$

  IMPLICIT NONE
  PRIVATE
  
  !-------------------------------------------------------
  ! Work arrays
!DEC$ ATTRIBUTES ALIGN:16 :: var_C
  COMPLEX(KIND=8), SAVE, ALLOCATABLE :: var_C(:,:)
!DEC$ ATTRIBUTES ALIGN:16 :: var_R
  REAL(KIND=8), SAVE,    ALLOCATABLE :: var_R(:,:)
!$OMP threadprivate(var_C,var_R)

  ! Physical and logical sizes for 2d-DFT
  INTEGER(KIND=4), SAVE        :: n(2), n_log(2)
  ! FFTW plans
  INTEGER(KIND=8), SAVE        :: plan_r2c, plan_c2r
  ! Normalization scale for fftw
  REAL(KIND=8), SAVE           :: scale,scale_log
  ! Indicator of dealiasing
  LOGICAL , SAVE               :: ifpad

  PUBLIC planner_fftw, fwd_fft, bwd_fft, destroyer_fftw, deriv_x, deriv_y

CONTAINS

  !-------------------------------------------------------
  SUBROUTINE planner_fftw(n1,n2,ifpad_in,fftw_nthreads)

    IMPLICIT NONE
    include 'fftw3.f'

    ! Flag for fftw
    !INTEGER(KIND=4), PARAMETER  :: flag_fftw = FFTW_ESTIMATE
    INTEGER(KIND=4), PARAMETER  :: flag_fftw = FFTW_PATIENT ! EXHAUSTIVE
    
    INTEGER(KIND=4), INTENT(IN) :: n1,n2
    LOGICAL,         INTENT(IN) :: ifpad_in
    INTEGER(KIND=4), INTENT(IN) :: fftw_nthreads
    INTEGER(KIND=4) ierr,omp_get_max_threads
    LOGICAL omp_get_nested

#ifndef BLUEGENE
    if (fftw_nthreads .gt. 1) then
! Set up multi-threaded fftw only if compiled with -openmp
!$       if (omp_get_nested()) then
!$         call dfftw_init_threads(ierr)
!$         if (ierr .eq. 0) then
!$            print *, 'dfftw_init_threads failed: ',ierr
!$            stop 'planner_fftw'
!$         endif
!$         call dfftw_plan_with_nthreads(fftw_nthreads)
!$         print *, 'dfftw_plan_with_nthreads: ',fftw_nthreads
!$       endif
    endif
#endif

    n(1) = n1
    n(2) = n2
    ifpad = ifpad_in
      
    ! Set logical sizes for 2d FFTW
    IF (ifpad) THEN
       ! Applying '3/2' rule for dealiasing
       n_log = 3*n/2  
    ELSE
       n_log = n
    END IF

    ! Set scale parameters for 2d FFTW
    scale = 1.0D0/DBLE(n(1)*n(2))
    scale_log = 1.0D0/DBLE(n_log(1)*n_log(2))

    ! Allocate work arrays

!$OMP PARALLEL
    ALLOCATE(var_R(n_log(1)    , n_log(2)))
    ALLOCATE(var_C(n_log(1)/2+1, n_log(2)))
!$OMP END PARALLEL


    ! Set up fftw plans
    CALL dfftw_plan_dft_r2c_2d(plan_r2c, n_log(1), n_log(2), &
         var_R, var_C, flag_fftw)
    CALL dfftw_plan_dft_c2r_2d(plan_c2r, n_log(1), n_log(2), &
         var_C, var_R, flag_fftw)
    
! check output whether SSE2/AVX kernels are used
    !call dfftw_print_plan(plan_r2c)
    !call dfftw_print_plan(plan_c2r)

  END SUBROUTINE planner_fftw

  !-------------------------------------------------------
  ! Forward fourier transformation
  SUBROUTINE fwd_fft(x,x_hat)
    
    IMPLICIT NONE

    REAL(KIND=8),    INTENT(IN)  :: x(n_log(1),n_log(2))
    COMPLEX(KIND=8), INTENT(OUT) :: x_hat(n(1)/2+1,n(2))
    INTEGER(KIND=4) :: i

 
    var_R = x
    CALL dfftw_execute_dft_r2c(plan_r2c,var_R,var_C)
    
    ! Chop high wavenumbers' coefficients if padded
    IF (ifpad) THEN
       x_hat = chopper(var_C)
       x_hat = scale_log*x_hat
    ELSE
       x_hat = var_C
       x_hat = scale*x_hat
       ! Zeroing Nyquist wavenumbers
       x_hat(n_log(1)/2+1,:) = 0.0D0
       x_hat(:,n_log(2)/2+1) = 0.0D0
    END IF

    !Enforce symmetry (0,k_z) = (0,-k_z)

    do i=2,n(2)/2
       x_hat(1,n(2)+2-i) = conjg(x_hat(1,i))
    end do
    
  END SUBROUTINE fwd_fft

  !-------------------------------------------------------
  ! Backward fourier transformation
  SUBROUTINE bwd_fft(x_hat,x)
    
    IMPLICIT NONE

    REAL(KIND=8), INTENT(OUT)   :: x(n_log(1),n_log(2))
    COMPLEX(KIND=8), INTENT(IN) :: x_hat(n(1)/2+1,n(2))

    if (ifpad) THEN
       var_C = padding(x_hat)
    ELSE
       var_C = x_hat
       ! Zeroing Nyquist wavenumbers
       var_C(n_log(1)/2+1,:) = 0.0D0
       var_C(:,n_log(2)/2+1) = 0.0D0
    END if
    CALL dfftw_execute_dft_c2r(plan_c2r,var_C,var_R)
    x = var_R

  END SUBROUTINE bwd_fft

  !-------------------------------------------------------
  SUBROUTINE destroyer_fftw()

    ! Destroy FFTW plans
    CALL dfftw_destroy_plan(plan_r2c)
    CALL dfftw_destroy_plan(plan_c2r)

    ! Deallocate work arrays
!$OMP PARALLEL
    DEALLOCATE(var_R, var_C)
!$OMP END PARALLEL

  END SUBROUTINE destroyer_fftw

  !------------------------------------------------------
  ! Dealiasing, in spectral space
  ! Padding zeros to the high wavenumber coefficients
  FUNCTION padding(x) RESULT(x_pad)
      
    IMPLICIT NONE

    COMPLEX(KIND=8) :: x(n(1)/2+1,n(2))
    COMPLEX(KIND=8) :: x_pad(n_log(1)/2+1,n_log(2))
    INTEGER(KIND=4) :: i,j

    x_pad = (0.0D0,0.0D0)
    DO j=1,n(2)/2+1
       DO i=1,n(1)/2+1
          x_pad(i,j) = x(i,j)
       END DO
    END DO
    DO j=n_log(2)-n(2)/2+2,n_log(2)
       DO i=1,n(1)/2+1
          x_pad(i,j) = x(i,j-n_log(2)+n(2))
       END DO
    END DO

    RETURN
  END FUNCTION padding
  
  !------------------------------------------------------
  ! Do after the pseudo-technique for the nonlinearity
  ! Chop the high wavenumber coefficients
  FUNCTION chopper(x) RESULT(x_chop)
    
    IMPLICIT NONE

    COMPLEX(KIND=8) :: x(n_log(1)/2+1,n_log(2))
    COMPLEX(KIND=8) :: x_chop(n(1)/2+1,n(2))
    INTEGER(KIND=4) :: i,j

    x_chop = (0.0D0,0.0D0)
    DO j=1,n(2)/2+1
       DO i=1,n(1)/2+1
          x_chop(i,j) = x(i,j)
       END DO
    END DO
    DO j=n(2)/2+2,n(2)
       DO i=1,n(1)/2+1
          x_chop(i,j) = x(i,j+n_log(2)-n(2))
       END DO
    END DO

  END FUNCTION chopper
  
!!$  !---------------------------------------------------------
!!$  ! Calculate the derivatives of a scalar function f(x,y)
!!$  ! Giving f_hat(:,:) & the wavenumbers k_{x,y}
!!$  SUBROUTINE derivation(k_x,k_y,f_hat,f,dfdx,dfdy)
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    REAL(KIND=8),    INTENT(IN)  :: k_x,k_y
!!$    COMPLEX(KIND=8), INTENT(IN)  :: f_hat(n(1)/2+1,n(2))
!!$    COMPLEX(KIND=8)              :: dfdx_hat(n(1)/2+1,n(2))
!!$    COMPLEX(KIND=8)              :: dfdy_hat(n(1)/2+1,n(2))
!!$    REAL(KIND=8),    INTENT(OUT) :: f(n_log(1),n_log(2))
!!$    REAL(KIND=8),    INTENT(OUT) :: dfdx(n_log(1),n_log(2))
!!$    REAL(KIND=8),    INTENT(OUT) :: dfdy(n_log(1),n_log(2))
!!$    COMPLEX(KIND=8)              :: c_x,c_y
!!$    INTEGER(KIND=4)              :: i,j
!!$
!!$    CALL bwd_fft(f_hat,f)
!!$    
!!$    DO j=1,n(2)/2+1
!!$       DO i=1,n(1)/2
!!$          c_x = DCMPLX(0.D0,(i-1)*k_x)
!!$          c_y = DCMPLX(0.D0,(j-1)*k_y)
!!$          dfdx_hat(i,j) = c_x*f_hat(i,j)
!!$          dfdy_hat(i,j) = c_y*f_hat(i,j)
!!$       END DO
!!$    END DO
!!$    dfdx_hat(n(1)/2+1,:) = 0.D0
!!$    dfdy_hat(n(1)/2+1,:) = 0.D0
!!$    dfdx_hat(:,n(2)/2+1) = 0.D0
!!$    dfdy_hat(:,n(2)/2+1) = 0.D0
!!$    DO j=n(2)/2+2,n(2)
!!$       DO i=1,n(1)/2
!!$          c_x = DCMPLX(0.D0,(i-1)*k_x)
!!$          c_y = DCMPLX(0.D0,(j-n(2)-1)*k_y)
!!$          dfdx_hat(i,j) = c_x*f_hat(i,j)
!!$          dfdy_hat(i,j) = c_y*f_hat(i,j)
!!$       END DO
!!$    END DO
!!$
!!$    CALL bwd_fft(dfdx_hat,dfdx)
!!$    CALL bwd_fft(dfdy_hat,dfdy)
!!$
!!$    RETURN
!!$  END SUBROUTINE derivation

  !----------------------------------------------------------

  SUBROUTINE deriv_x(k_x,f_hat,df)
    
    IMPLICIT NONE
    
    REAL(KIND=8),    INTENT(IN)  :: k_x
    COMPLEX(KIND=8), INTENT(IN)  :: f_hat(n(1)/2+1,n(2))
    COMPLEX(KIND=8)              :: df_hat(n(1)/2+1,n(2))
    REAL(KIND=8),    INTENT(OUT) :: df(n_log(1),n_log(2))
    COMPLEX(KIND=8)              :: c_x
    INTEGER(KIND=4)              :: i,j

    df_hat(n(1)/2+1,:) = (0.0D0,0.0D0)
    df_hat(:,n(2)/2+1) = (0.0D0,0.0D0)

    DO j=1,n(2)/2
       DO i=1,n(1)/2
          c_x = DCMPLX(0.D0,(i-1)*k_x)
          df_hat(i,j) = c_x*f_hat(i,j)
       END DO
    END DO

    DO j=n(2)/2+2,n(2)
       DO i=1,n(1)/2
          c_x = DCMPLX(0.D0,(i-1)*k_x)
          df_hat(i,j) = c_x*f_hat(i,j)
       END DO
    END DO

    CALL bwd_fft(df_hat,df)

  END SUBROUTINE deriv_x
  !----------------------------------------------------------

  SUBROUTINE deriv_y(k_y,f_hat,df)
    
    IMPLICIT NONE
    
    REAL(KIND=8),    INTENT(IN)  :: k_y
    COMPLEX(KIND=8), INTENT(IN)  :: f_hat(n(1)/2+1,n(2))
    COMPLEX(KIND=8)              :: df_hat(n(1)/2+1,n(2))
    REAL(KIND=8),    INTENT(OUT) :: df(n_log(1),n_log(2))
    COMPLEX(KIND=8)              :: c_y
    INTEGER(KIND=4)              :: i,j

    df_hat(n(1)/2+1,:) = (0.0D0,0.0D0)
    df_hat(:,n(2)/2+1) = (0.0D0,0.0D0)

    DO j=1,n(2)/2
       c_y = DCMPLX(0.D0,(j-1)*k_y)
       DO i=1,n(1)/2
          df_hat(i,j) = c_y*f_hat(i,j)
       END DO
    END DO

    DO j=n(2)/2+2,n(2)
       c_y = DCMPLX(0.D0,(j-n(2)-1)*k_y)
       DO i=1,n(1)/2
          df_hat(i,j) = c_y*f_hat(i,j)
       END DO
    END DO

    CALL bwd_fft(df_hat,df)

  END SUBROUTINE deriv_y

  !----------------------------------------------------------
END MODULE mod_fftw
