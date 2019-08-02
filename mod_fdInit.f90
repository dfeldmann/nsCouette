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

MODULE mod_fdInit
!      A module to compute 
!      the finite difference weights
!      for the derivatives

  IMPLICIT NONE
  
CONTAINS




  SUBROUTINE radial_grid(m,rad)

    !Compute radial points as in Ashleys code.
    
    IMPLICIT NONE 
    
    INTEGER(kind=4),intent(in) :: m
    REAL(kind=8),intent(out) :: rad(m)
    integer(kind=4) :: n, N_
    real(kind=8) :: dr
    REAL(KIND=8)   ,PARAMETER :: PI = ACOS(-1d0) 
    ! T_{N-1}(x) has N extrema
    ! Numerical Recipes 5.8.5
      
      
    !Chebyshev collocation points (!zero cannot be a point of the domain)
    N_ = m+int(sqrt(real(m)))
       
    do n = N_-m+1, N_
       rad(n-N_+m) = 0.5d0*( 1d0+dcos(PI*(N_-n)/real(N_)) )
          
    end do
       
    do n = 1, 10
       dr = 1.5d0*rad(1)-0.5d0*rad(2)
       rad(:) = rad(:)*(1d0+dr) - dr
    end do
        
    
  end subroutine radial_grid
                                                                                                                         


  !------------------------------------------------------------------                         
  SUBROUTINE fd_Matrix(m,n,x,Mw_dx,Mw_dxx,dr0,dr1,intrdr)                                     
  !--------------------------------------------------------------                             
  ! FD weights matrix M such that df = M*f, for arbitrary f(x)                                
  ! Input:                                                                                    
  !      x       Descretized points in x-direction                                            
  !      n       Length(x)-1                                                                  
  !      m       Length(stencil)-1                                                            
  ! Output:                                                                                   
  !      Mw_dx   Weights matrix for the 1st derivative                                        
  !      Mw_dxx  Weights matrix for the 2nd derivative                                        
  !      dr0     Weights matrix to extrapolate to r=0                                         
  !      Matrices in BLAS-like band storage                                                   
  !--------------------------------------------------------------                             
                                                                                              
    IMPLICIT NONE                                                                             
                                                                                              
    INTEGER(KIND=4), INTENT(IN)  :: m,n                                                       
    REAL(KIND=8), INTENT(IN)  :: x(0:n)                                                       
    REAL(KIND=8), INTENT(OUT) :: Mw_dx(0:m,0:n)                                               
    REAL(KIND=8), INTENT(OUT) :: Mw_dxx(0:m,0:n)                                              
    REAL(KIND=8), INTENT(OUT) :: dr0(0:m/2,0:1),dr1(0:m/2,0:1),intrdr(0:n)                    
    REAL(KIND=8), ALLOCATABLE :: s(:),c(:,:)                                                  
    INTEGER(KIND=4)           :: i,j,mt,k,iw,left,right                                       
    REAL(kind=8)              :: e1,e2,r_(-1:n)                                               
                                                                                              
    Mw_dx  = 0.0D0                                                                            
    Mw_dxx = 0.0D0                                                                            
                                                                                              
    dr0=0d0                                                                                   
                                                                                              
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
                                                                                              
                                                                                              
    !Compute weights to extrapolate to r=0 and r=1                                            
                                                                                              
   CALL fd_Weigths(0d0,x(0:m/2),m/2,1,dr0)                                                    
   mt = m/2                                                                                   
   ALLOCATE(s(0:mt),c(0:mt,0:mt))                                                             
   s = x(n-m/2:n)                                                                             
   CALL fd_Weigths(x(n),s,mt,1,c)                                                             
   dr1=c(:,0:1)                                                                               
   DEALLOCATE(s,c)                                                                            
                                                                                              
                                                                                              
   !radial points                                                                             
                                                                                              
   r_(0:n)=x(0:n)                                                                             
   r_(-1)= 0d0                                                                                
                                                                                              
                                                                                              
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
         e1 = e1 * (r_(min(i+1,n))-r_(i)) / real(j+1)                                         
         e2 = e2 * (r_(max(-1,i-1))-r_(i)) / real(j+1)                                        
         intrdr(left:left+mt) = intrdr(left:left+mt) + 0.5d0*(e1-e2)*c(0:mt,j)*r_(i)          
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

END MODULE mod_fdInit

