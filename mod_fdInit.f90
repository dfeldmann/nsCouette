!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file is part of NSCouette, a HPC code for DNS of Taylor-Couette flow !
!                                                                           !
! Copyright (C) 2016 Marc Avila, Bjoern Hof, Jose Manuel Lopez,             !
!                    Markus Rampp, Liang Shi                                !
!                                                                           !
! NSCouette is free software: you can redistribute it and/or modify         !
! it under the terms of the GNU General Public License as published by      !
! the Free Software Foundation, either version 3 of the License, or         !
! (at your option) any later version.                                       !
!                                                                           !
! NSCouette is distributed in the hope that it will be useful,              !
! but WITHOUT ANY WARRANTY; without even the implied warranty of            !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             !
! GNU General Public License for more details.                              !
!                                                                           !
! You should have received a copy of the GNU General Public License         !
! along with NSCouette.  If not, see <http://www.gnu.org/licenses/>.        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=============================================
!
!      A module to compute 
!      the finite difference weights
!      for the derivatives
!
!=============================================

MODULE mod_fdInit
  
  IMPLICIT NONE
  private

  public :: fd_Matrix

CONTAINS

  !------------------------------------------------------------------
  SUBROUTINE fd_Matrix(m,n,x,Mw_dx,Mw_dxx,dr1,dr2,intrdr)
  !--------------------------------------------------------------
  ! FD weights matrix M such that df = M*f, for arbitrary f(x)
  ! Input:
  !      x       Descretized points in x-direction
  !      n       Length(x)-1
  !      m       Length(stencil)-1
  ! Output:
  !      Mw_dx   Weights matrix for the 1st derivative
  !      Mw_dxx  Weights matrix for the 2nd derivative
  !      dr1     Weights to extrapolate to r=r_i
  !      dr2     Weights to extrapolate to r=r_o
  !      intrdr  Weights for integration
  !      Matrices in BLAS-like band storage
  !--------------------------------------------------------------

    IMPLICIT NONE
    
    INTEGER(KIND=4), INTENT(IN)  :: m,n
    REAL(KIND=8), INTENT(IN)  :: x(0:n)
    REAL(KIND=8), INTENT(OUT) :: Mw_dx(0:m,0:n)
    REAL(KIND=8), INTENT(OUT) :: Mw_dxx(0:m,0:n)
    REAL(KIND=8), INTENT(OUT) :: dr1(0:m/2,0:m/2),dr2(0:m/2,0:m/2),intrdr(0:n)
    REAL(KIND=8), ALLOCATABLE :: s(:),c(:,:)
    INTEGER(KIND=4)           :: i,j,mt,k,iw,left,right
    REAL(kind=8)              :: e1,e2

    Mw_dx  = 0.0D0
    Mw_dxx = 0.0D0


    iw=m/2


    ALLOCATE(s(0:m),c(0:m,0:m))
    DO i = m/2,n-m/2
       s = x(i-m/2:i+m/2)
       CALL fd_Weigths(x(i),s,m,2,c)


       do j = i-m/2,i+m/2

          k = iw - j
          Mw_dx(k+i,j)  = c(j-(i-m/2),1)
          Mw_dxx(k+i,j) = c(j-(i-m/2),2) 
           

       end do
          
    END DO
    DEALLOCATE(s,c)

    DO i = 0,m/2-1
       mt = m/2+i
       ALLOCATE(s(0:mt),c(0:mt,0:mt))

       s = x(0:mt)
       CALL fd_Weigths(x(i),s,mt,2,c)

       
       do j = 0,mt

          k = iw - j
          Mw_dx(k+i,j)  = c(j,1)
          Mw_dxx(k+i,j) = c(j,2) 
           

       end do
       s = x(n-mt:n)
       CALL fd_Weigths(x(n-i),s,mt,2,c)

       do j = n-mt,n

          k = iw - j
          Mw_dx(k+n-i,j)  = c(j-(n-mt),1)
          Mw_dxx(k+n-i,j) = c(j-(n-mt),2) 
           

       end do
       DEALLOCATE(s,c)
    END DO
        
    !Compute weights to extrapolate to r_i and r_o                                                                                          

    CALL fd_Weigths(x(0),x(0:m/2),m/2,m/2,dr1)
    mt = m/2
    ALLOCATE(s(0:mt),c(0:mt,0:mt))
    s = x(n-m/2:n)
    CALL fd_Weigths(x(n),s,mt,mt,c)
    dr2=c(:,:)
    DEALLOCATE(s,c)

       ! weights for integration  r dr                                                                                                          

   intrdr = 0d0
   do i = 0, n
      left = max(0,i-iw)
      right = min(i+iw,n)
      mt = right-left
      ALLOCATE(s(0:mt),c(0:mt,0:mt))
      s(0:mt)=x(left:right)
      CALL fd_Weigths(x(i),s,mt,mt,c)
      e1 = 1d0
      e2 = 1d0
      do j = 0, mt
         e1 = e1 * (x(min(i+1,n))-x(i)) / dble(j+1)
         e2 = e2 * (x(max(0,i-1))-x(i)) / dble(j+1)
         intrdr(left:left+mt) = intrdr(left:left+mt) + 0.5d0*(e1-e2)*c(0:mt,j)*x(i)
      end do
      DEALLOCATE(s,c)
   end do
    

    RETURN
  END SUBROUTINE fd_Matrix

  !-------------------------------------------------------------------
  SUBROUTINE fd_Weigths(x0,x,n,m,c)
  !------------------------------------------------------------------
  ! A general algorithm of FD weights based on Lagrange interpolation
  !
  ! Input:
  !     x0  Evaluation point
  !     x   Interpolation points or stencil
  !     n   Length(x)-1
  !     m   Highest order of derivatives to be approximated
  !
  ! Output
  !     c   Weights of derivative approximation
  !------------------------------------------------------------------

    IMPLICIT NONE
    
    INTEGER(KIND=4), INTENT(IN) :: n,m
    REAL(KIND=8), INTENT(IN)    :: x0,x(0:n)
    REAL(KIND=8), INTENT(OUT)   :: c(0:n,0:m)
    INTEGER(KIND=4)             :: i,j,k,mn
    REAL(KIND=8)                :: c1,c2,c3,c4,c5

    c1 = 1.0D0
    c4 = x(0)-x0
    c  = 0.0D0
    c(0,0) = 1.0D0
    DO i = 1,n
       mn = MIN(i,m)
       c2 = 1.0D0
       c5 = c4
       c4 = x(i) - x0
       DO j = 0,i-1
          c3 = x(i) - x(j)
          c2 = c2*c3
          if (j .EQ. i-1) THEN
             DO k = mn,1,-1
                c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
             END DO
             c(i,0) = -c1*c5*c(i-1,0)/c2
          END if
          DO k = mn,1,-1
             c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
          END DO
          c(j,0) = c4*c(j,0)/c3
       END DO
       c1 = c2
    END DO

    RETURN
  END SUBROUTINE fd_Weigths

!------------------------------------------------------------------
END MODULE mod_fdInit
