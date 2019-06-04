!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file is part of NSPipeflow, a HPC code for DNS of pipe flow !
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

!===========================================================================
! TIME STEP
! Time splitting scheme based on a predictor-corrector method
! Influence matrices are used to impose the boundary conditions
! For further details see documentation of openpipeflow.org
! Main idea:
!  1) Compute predictor of velocity field
!  2) Corrector step: recompute non linear terms and iterate up to velocity field converges
!  whithin the tolerance requested in mod_params.f90
! Note that all boundary conditions are imposed on the velocity field and thus
! no artificial boundary conditions for the pressure are needed.
! Linear equations
!      (A*x=b)
!      are solved by LU decomposition method and in spectral space
!
!===========================================================================           

MODULE mod_timeStep

  USE mod_inOut
  USE mod_nonlinear


  IMPLICIT NONE
  
  ! Variables for linear equations A*x = b, L is laplacian operator

  COMPLEX(KIND=8), PRIVATE, allocatable :: b(:)

  ! Velocity u_old,u,u _new and nonlinear term udu
  TYPE (vec_mpi) , PRIVATE :: u_hat_mp_old, u_hat_mp_int,udu_hat_mp_old
  COMPLEX(KIND=8), PRIVATE, allocatable, dimension(:,:) :: u_plus_old,u_minus_old
  COMPLEX(KIND=8), PRIVATE, allocatable, dimension(:,:) :: udu_plus_old,udu_minus_old
 
  !variables to check convergence 
  REAL(kind=8),PRIVATE :: corr_dt,dterr
  !Other parameters
  INTEGER(KIND=4), private, PARAMETER :: iwidth=(n_s-1)/2, idiag=iwidth+1
  
  !Influence matrices
  REAL(KIND=8),allocatable, dimension(:,:,:),private :: u_functions
  REAL(kind=8), dimension(:,:,:), allocatable, private  :: inf_matrix
  REAL(kind=8),private    :: BRe(4),BIm(4), aR(4),aI(4)
    
  ! see http://software.intel.com/en-us/forums/topic/293641
  !DEC$ OPTIONS /NOWARN
  type LUdcmp
     sequence
     real(kind=8),allocatable  :: LU(:,:)
     integer(kind=4),allocatable :: ipiv(:)
     integer(kind=4) :: n,k,lda,ldb,info
  end type LUdcmp
  !DEC$ END OPTIONS
  type(LUdcmp), allocatable, dimension(:), private :: LU_p,LU_uz,LU_up,LU_um
  !Laplacian matrices for the RHS of the NS equations
  REAL(KIND=8), dimension(:,:,:), allocatable, private ::  Lap_uz, Lap_up, Lap_um
  !-----------variables to adjust axial so that net mass flux is zero

  REAL(kind=8), dimension(:), allocatable, private :: adj_flux
  REAL(kind=8), private :: c_adj
  
CONTAINS
  !--------------------------------------------------------------------
  SUBROUTINE pre_timeStep(init)
    
    IMPLICIT NONE
    LOGICAL, intent(in),OPTIONAL :: init
    REAL(KIND=8)    :: bA(n_s,m_r),bL(n_s,m_r)
    INTEGER(KIND=4) :: i,j,k,h
    INTEGER(KIND=4) :: S_big(m_f,2)
    Complex(kind=8) :: aux1(m_r),aux2(m_r),aux3(m_r),auxp(m_r),auxm(m_r),x(m_r)
    
    call perfon('pretstp')

    if (init) then

       !allocate variables

       ierr=u_hat_mp_old%alloc(m_r,mp_f)
       ierr=u_hat_mp_int%alloc(m_r,mp_f)
       ierr=udu_hat_mp_old%alloc(m_r,mp_f)

       allocate(u_plus_old(m_r,mp_f))
       allocate(u_minus_old(m_r,mp_f))
       allocate(udu_plus_old(m_r,mp_f))
       allocate(udu_minus_old(m_r,mp_f))

       allocate(Sbig_map(mp_f,2))
       allocate(LU_p(mp_f),LU_uz(mp_f),LU_up(mp_f),LU_um(mp_f))
       allocate(Lap_uz(n_s,m_r,mp_f), Lap_up(n_s,m_r,mp_f), Lap_um(n_s,m_r,mp_f))
       
       allocate(inf_matrix(4,4,mp_f))
       allocate(u_functions(m_r,mp_f,6))

       allocate(adj_flux(m_r))
       allocate(b(m_r))
       allocate(bMw_dr(n_s,m_r),bMw_drr(n_s,m_r))
       allocate(dr1((n_s-1)/2+1,0:1),dr0((n_s-1)/2+1,0:1))
       allocate(intrdr(m_r))
       

       
       !Initialize variables
       
       corr_dt = 0d0
       dterr = 0d0
       new_dt = .false.

       ! Radial derivatives 
       Call fd_Matrix(n_s-1,m_r-1,r,bMw_dr,bMw_drr,dr0,dr1,intrdr)
       

       ! Precompute parity
       
       !S index, 1(Pressure and axial velocity )                                                                                                                                                                                                 
       !2(u_pos and u_neg).                                                                                                                                                                                                                                                 
       !symmetry on axis            ! S==1 fn odd => f=0                                                                                                                                                                                                          
       ! S==0 fn even => dr(f)==0                                                                                                                                                                                                                                             
       
       do j=1,m_th+1
          do k=1,m_z
             S_big((m_th+1)*(k-1)+j,1)=modulo((j-1)*int(k_th0),2)
             S_big((m_th+1)*(k-1)+j,2)=modulo((j-1)*int(k_th0)+1,2)
             
          end do
       end do
       
       !Distribute in each processor                                                                                                                                                                                                                                           
       
       Sbig_map(:,1)=S_big(1+(myid*mp_f):mp_f+(myid*mp_f),1)
       Sbig_map(:,2)=S_big(1+(myid*mp_f):mp_f+(myid*mp_f),2)
       
    end if
       
   

    ! calculate Laplacian differentiation Matrix
    ! and precompute the LU decompositions in both
    ! left and right hand side of N-S 
  
!$OMP PARALLEL DO &                                                               
!$OMP  SCHEDULE (RUNTIME) &                                                                                
!$OMP  DEFAULT(NONE) &                                                                                                                                       
!$OMP  PRIVATE(h,i,j,k,bA,bL,BRe,BIm,aux1,aux2,aux3,auxp,auxm,x) &                
!$OMP  SHARED(init,m_r,const_flux,Re,mp_f,dt,bMw_dr,bMw_drr,fk_mp,LU_p,LU_uz,LU_up,LU_um,dr0,Sbig_map,r& 
!$OMP         ,u_functions,inf_matrix,Lap_uz,Lap_um,Lap_up,adj_flux,c_adj,intrdr)                       

    DO k=1,mp_f
       DO j=1,m_r
          h = iwidth + 1 - j
          do i = max(1,j-iwidth), min(m_r,j+iwidth)
             bL(h+i,j) = bMw_dr(h+i,j)/r(i) + bMw_drr(h+i,j)
          enddo
          bL(idiag,j) = bL(idiag,j) - (fk_mp%th(j,k)/r(j))**2 - (fk_mp%z(j,k))**2
       END DO

       if (init) then 
          !For pressure p                                                                                                                                                                                                                                                         
          
          bA(:,:) = bL(:,:)
          
          !Introduce proper parity at r=0                                                                                                                                                                                                                                      
          
          do h=1,idiag
             bA(idiag-h+1,h)=dr0(h,1-Sbig_map(k,1))
          end do
          
          
          ! Neumann B.C. at r=1                                                                                                                                                                                                                                                
          do h=1,idiag
             bA(idiag+h-1,m_r-h+1) = bMw_dr(idiag+h-1,m_r-h+1)!  A(m_r,:) = Mw_dr(m_r,:)                                                                                                                                                                                       
          enddo
          
          ! Replace the Neumann BC by Dirichlet at r = r_o                                                                                                                                                                                                                        
          if (abs(fk_mp%th(1,k)) <= epsilon .and. abs(fk_mp%z(1,k)) <= epsilon) then
             do h=1,idiag
                bA(idiag+h-1,m_r-h+1) = 0.d0                  !  A(m_r,:) = 0.d0                                                                                                                                                                                    
             enddo
             bA(idiag,m_r) = 1.d0                             !  A(m_r,m_r) = 1.d0                                                                                                                                                                                                
          end if
          
          call luDecmp(m_r,(n_s-1)/2,bA,LU_p(k))

       end if
            
       !predictor_v                                                                                                                                                                                                                                                            
       
       !uz                                                                                                                                                                                                                                                                     
       
       !RHS equation                                                                                                                                                                                                                                                           
       
       Lap_uz(:,:,k) = (1d0-d_implicit)/Re*bL(:,:)
       Lap_uz(idiag,:,k) =  1.0D0/dt + Lap_uz(idiag,:,k)
       
       !LHS equation                                                                                                                                                                                                                                                           
       bA(:,:) = -(d_implicit/Re)*bL(:,:)
       bA(idiag,:) = 1.D0/dt + bA(idiag,:)
       
       !Introduce proper parity at r=0                                                                                                                                                                                                                                         
       
       
       do h=1,idiag
          bA(idiag-h+1,h)=dr0(h,1-Sbig_map(k,1))
       end do
       
       
       
       !Dirichlet boundary condition at r=R                                                                                                                                                                                                                                    
       do h=1,idiag
          
          bA(idiag+h-1,m_r-h+1) = 0.d0      !  A(m_r,:) = 0.d0                                                                                                                                                                                                                 
       enddo
       
       bA(idiag,m_r) = 1.d0                 !  A(m_r,m_r) = 1.d0                                                                                                                                                                                                               
       
       call luDecmp(m_r,(n_s-1)/2,bA,LU_uz(k))
     

       if ((const_flux) .and. (abs(fk_mp%th(1,k)) <= epsilon .and. abs(fk_mp%z(1,k)) <= epsilon)) then
          
          aux1(:) = dcmplx(1d0,0d0)
          aux1(1) = dcmplx(0d0,0d0)
          aux1(m_r) = dcmplx(0d0,0d0)


          call luSolve(LU_uz(k),aux1(:), x)

          adj_flux = dreal(x)                                                                                                    
          c_adj = 2d0*dot_product(adj_flux, intrdr)                                                                                        

       end if   

       
       !u_+                                                                                                                                                                                                                                                               
       
       bA(:,:) = bL(:,:)
       
       do j=1,m_r
          bA(idiag,j) =  bA(idiag,j) - (1d0 + 2d0*fk_mp%th(j,k))/r(j)**2
       end do
       
       !RHS equation                                                                                                                                                                                                                                                      
       Lap_up(:,:,k) = (1d0-d_implicit)/Re*bA(:,:)
       Lap_up(idiag,:,k) =  1.0D0/dt + Lap_up(idiag,:,k)
       
       !LHS equation                                                                                                                                                                                                                                                      
       bA(:,:) = -(d_implicit/Re)*bA(:,:)
       bA(idiag,:) = 1.0D0/dt + bA(idiag,:)
       
       
       !Introduce proper parity at r=0                                                                                                                                                                                                                                    
          
       
       do h=1,idiag
          bA(idiag-h+1,h)=dr0(h,1-Sbig_map(k,2))
       end do
       
       
       !Dirichlet boundary condition at r=R                                                                                                                                                                                                                               
       
       do h=1,idiag
          bA(idiag+h-1,m_r-h+1) = 0.d0      !  A(m_r,:) = 0.d0                                                                                                                                                                                                            
       enddo
       
       bA(idiag,m_r) = 1.d0                 !  A(m_r,m_r) = 1.d0                                                                                                                                                                                                          
       
       
       call luDecmp(m_r,(n_s-1)/2,bA,LU_up(k))
       
       
       !u_-                                                                                                                                                                                                                                                               
       
       
       bA(:,:) = bL(:,:)

       do j=1,m_r
          bA(idiag,j) = bA(idiag,j) + (-1d0 + 2d0*fk_mp%th(j,k))/r(j)**2
       end do
       
       !RHS equation                                                                                                                                                                                                                                                      
       Lap_um(:,:,k) = (1d0-d_implicit)/Re*bA(:,:)
       Lap_um(idiag,:,k) =  1.0D0/dt + Lap_um(idiag,:,k)
       
       
       bA(:,:) = -(d_implicit/Re)*bA(:,:)
       bA(idiag,:) = 1.0D0/dt + bA(idiag,:)
       
       
       !LHS equation                                                                                                                                                                                                                                                      
       
       !Introduce proper parity at r=0                                                                                                                                                                                                                                    
       
       
       do h=1,idiag
          bA(idiag-h+1,h)=dr0(h,1-Sbig_map(k,2))
       end do
       
       
       
       !Dirichlet boundary condition at r=R                                                                                                                                                                                                                               
       
       do h=1,idiag
          bA(idiag+h-1,m_r-h+1) = 0.d0      !  A(m_r,:) = 0.d0                                                                                                                                                                                                            
       enddo
       
       bA(idiag,m_r) = 1.d0                 !  A(m_r,m_r) = 1.d0                                                                                                                                                                                                          
       
          
       call luDecmp(m_r,(n_s-1)/2,bA,LU_um(k))
       
       
       !Get influence matrices                                                                                                                                                                                                                                            
       if (abs(fk_mp%th(1,k)) <= epsilon .and. abs(fk_mp%z(1,k)) <= epsilon) cycle

       aux1(:) = dcmplx(0d0,0d0)
       aux2(:) = dcmplx(0d0,0d0)
       aux3(:) = dcmplx(0d0,0d0)
       
       
       !Get u function when u+ = 1 is the boundary condition                                                                                                                                                                                                                
       aux1(m_r) = dcmplx(1d0,0d0)
       call luSolve(LU_up(k),aux1(:), x)
       u_functions(:,k,1) = DREAL(x)
       call evalBC(k, x, aux2, aux3, BRe, BIm)
       inf_matrix(:,1,k) = BRe(:)
       aux1(m_r) = dcmplx(0d0,0d0)
       
       !Get u function when u-=1 is the boundary condition                                                                                                                                                                                                                
       
       aux2(m_r) = dcmplx(1d0,0d0)
       call luSolve(LU_um(k),aux2(:), x)
       u_functions(:,k,2) = DREAL(x)
       call evalBC(k, aux1, x ,aux3 , BRe ,BIm )
       inf_matrix(:,2,k) = BRe(:)
       aux2(m_r) = dcmplx(0d0,0d0)
       
       
       !Get u function when uz=i is the boundary condition                                                                                                                                                                                                                
       
       aux3(m_r) = dcmplx(0d0,1d0)
       call luSolve(LU_uz(k),aux3(:), x)
       u_functions(:,k,3) = DIMAG(x)
       call evalBC(k, aux1, aux2, x, BRe, BIm)
       inf_matrix(:,3,k) = BRe(:)
       aux3(m_r) = dcmplx(0d0,0d0)
       
       !Get functions corresponding to the pressure                                                                                                                                                                                                                       
       
       aux1(m_r) = dcmplx(-1d0,0d0)
       call luSolve(LU_p(k),aux1(:), x)
       
       !Obtain gradient of pressure                                                                                                                                                                                                                                       
       
       
       aux1(:) = DZGBMV(bMw_dr,x(:),iwidth,iwidth) !radial derivative                                                                                                                                                                                                                  
       aux2(:)= ii*fk_mp%th(:,k)/r(:)*x(:)!azimuthal derivative                                                                                                                                                                                                           
       aux3(:) = ii*fk_mp%z(:,k)*x(:) !axial derivative                                                                                                                                                                                                                   
       
       !change of variables to decouple radial and azimuthal equations                                                                                                                                                                                                    
       auxp = aux1 + ii*aux2
       auxm = aux1 - ii*aux2
       
       
       u_functions(:,k,4) = DREAL(auxp)
       u_functions(:,k,5) = DREAL(auxm)
       u_functions(:,k,6) = DIMAG(aux3)
       
       
       
       call evalBC(k, auxp, auxm, aux3, BRe,BIm)
       
       inf_matrix(:,4,k) =  BRe(:)

       
       call mat_inv(4,inf_matrix(:,:,k),4)
       
       
         
    ENDDO

    !$OMP END PARALLEL DO
    call perfoff()
  
        
  END SUBROUTINE pre_timeStep



  !--------------------------------------------------------------------
  SUBROUTINE post_timeStep()

    IMPLICIT NONE

    !deallocate variables
    ierr=u_hat_mp_old%dealloc()
    ierr=u_hat_mp_int%dealloc()
    ierr=udu_hat_mp_old%dealloc()
    deallocate(u_plus_old,u_minus_old)
    deallocate(udu_plus_old,udu_minus_old)


    deallocate(LU_p,LU_uz,LU_up,LU_um)
    deallocate(Lap_uz,Lap_up,Lap_um)
    deallocate(inf_matrix,u_functions,p_hat_mp)

    deallocate(adj_flux)
    deallocate(b,r)
    deallocate(bMw_dr,bMw_drr)
    deallocate(dr1,dr0)
    deallocate(intrdr)
    
    ierr=u_hat_mp%dealloc()
    ierr=udu_hat_mp%dealloc()
    ierr=fk_mp%dealloc()
    
    

  END SUBROUTINE post_timeStep
  
  
  Subroutine predictor
    
    IMPLICIT NONE
    
    Integer(kind=4) :: k,j
    COMPLEX(KIND=8) :: rhs_p(m_r),rhs_m(m_r),rhs_z(m_r)
    COMPLEX(kind=8) :: rhs_r(m_r),rhs_th(m_r),dotp(m_r)
    COMPLEX(KIND=8) :: u_plus(m_r),u_minus(m_r)
    COMPLEX(KIND=8) :: udu_plus(m_r),udu_minus(m_r)
    
    call perfon('predic')
    !$OMP PARALLEL DO &
    !$OMP  SCHEDULE (RUNTIME) &
    !$OMP  DEFAULT(NONE) &
    !$OMP  PRIVATE(k,j,b,dotp,aR,aI,BRe,BIm,rhs_p,rhs_m,rhs_z,rhs_r,rhs_th,u_plus,u_minus,udu_plus,udu_minus) &
    !$OMP  SHARED(mp_f,m_r,Re,LU_p,Lu_up,Lu_um,Lu_uz,u_hat_mp,udu_hat_mp,fk_mp,p_hat_mp,bMw_dr,&
    !$OMP  inf_matrix,u_functions,Lap_uz,Lap_up,Lap_um,r,dr1,udu_hat_mp_old,u_plus_old,u_minus_old,u_hat_mp_old)   

    do k=1,mp_f
    
       ! Do change of variables to decouple equations for u_r & u_th                                                                                                                                                                                                             
       u_plus(:)    = u_hat_mp%r(:,k) + ii*u_hat_mp%th(:,k)
       u_minus(:)   = u_hat_mp%r(:,k) - ii*u_hat_mp%th(:,k)
       
       udu_plus(:)  = udu_hat_mp%r(:,k) + ii*udu_hat_mp%th(:,k)
       udu_minus(:) = udu_hat_mp%r(:,k) - ii*udu_hat_mp%th(:,k)
       
       
       
       ! Save velocities and nonlinear terms at step n                                                                                                                                                                                                                           
       
       udu_hat_mp_old%r(:,k) = udu_hat_mp%r(:,k)
       udu_hat_mp_old%th(:,k) = udu_hat_mp%th(:,k)
       udu_hat_mp_old%z(:,k) = udu_hat_mp%z(:,k)
       
       u_hat_mp_old%r(:,k) = u_hat_mp%r(:,k)
       u_hat_mp_old%th(:,k) = u_hat_mp%th(:,k)
       u_hat_mp_old%z(:,k) = u_hat_mp%z(:,k)
       
       
       u_plus_old(:,k) = u_plus(:)
       u_minus_old(:,k) = u_minus(:)
    
   
       
    !check if DZGBMV is faster or not
       
       rhs_p(:)=DZGBMV(Lap_up(:,:,k),u_plus(:),iwidth,iwidth)
       rhs_m(:)=DZGBMV(Lap_um(:,:,k),u_minus(:),iwidth,iwidth)
       rhs_z(:)=DZGBMV(Lap_uz(:,:,k),u_hat_mp%z(:,k),iwidth,iwidth)
       
       
       !Adding non-linear terms                                                                                                                                                                                                                                             
       
       rhs_p(:)=rhs_p(:)+udu_plus(:)
       rhs_m(:)=rhs_m(:)+udu_minus(:)
       rhs_z(:)=rhs_z(:)+udu_hat_mp%z(:,k)
       
       
       !Get pressure (project RHS)                                                                                                                                                                                                                                          
       
       rhs_r(:)=  0.5D0 * (rhs_p(:) + rhs_m(:))
       rhs_th(:)=  -0.5D0*ii*(rhs_p(:) - rhs_m(:))
       
       !Compute the divergence of RHS                                                                                                                                                                                                                                       
       dotp(:) = DZGBMV(bMw_dr,rhs_r(:),iwidth,iwidth)
       
       DO j=2,m_r-1
          
          b(j) = dotp(j) &
               + rhs_r(j)/r(j) &
               +ii*fk_mp%th(j,k)*rhs_th(j)/r(j)&
               +ii*fk_mp%z(j,k)*rhs_z(j)
          
       END DO
       
       
       !set boundary conditions to zero                                                                                                                                                                            
       
       b(1)=dcmplx(0d0,0d0)
       b(m_r)=dcmplx(0d0,0d0)
       
       if (abs(fk_mp%th(1,k)) <= epsilon .and. abs(fk_mp%z(1,k)) <= epsilon) then
          p_hat_mp(:,k) = dcmplx(0d0,0d0)
       else
          ! solve system using LU decomposition precomputed in pre_timeStep()
          call luSolve(LU_p(k),b,p_hat_mp(:,k))
       end if
       
       !Compute pressure gradient and substract from rhs (written directly as rhs_p and rhs_m)
       
       
       dotp(:)=DZGBMV(bMw_dr,p_hat_mp(:,k),iwidth,iwidth)
       
       
     rhs_p(:)=rhs_p(:)-dotp(:)+fk_mp%th(:,k)*p_hat_mp(:,k)/r(:)
     rhs_m(:)=rhs_m(:)-dotp(:)-fk_mp%th(:,k)*p_hat_mp(:,k)/r(:)
     rhs_z(:)=rhs_z(:)-ii*fk_mp%z(:,k)*p_hat_mp(:,k)
     
     
     
     !Set boundary conditions to zero
     
     
     rhs_p(1)=dcmplx(0d0,0d0)
     rhs_p(m_r)=dcmplx(0d0,0d0)
     rhs_m(1)=dcmplx(0d0,0d0)
     rhs_m(m_r)=dcmplx(0d0,0d0)
     rhs_z(1)=dcmplx(0d0,0d0)
     rhs_z(m_r)=dcmplx(0d0,0d0)
     
     
     !Solve equations
     
     
     call luSolve(LU_up(k),rhs_p(:),u_plus(:))
     call luSolve(LU_um(k),rhs_m(:),u_minus(:))
     call luSolve(LU_uz(k),rhs_z(:),u_hat_mp%z(:,k))
             
    
     !correct boundary conditions with influence matrix

     if (abs(fk_mp%th(1,k)) <= epsilon .and. abs(fk_mp%z(1,k)) <= epsilon) cycle
     
     call evalBC(k, u_plus(:),u_minus(:),u_hat_mp%z(:,k),BRe,BIm)
     
   
     aR(:) = -matmul(inf_matrix(:,:,k),BRe(:))
     aI(:) = -matmul(inf_matrix(:,:,k),BIm(:))

     u_plus(:) = dcmplx(dreal(u_plus(:)) + aR(1)*u_functions(:,k,1)+ aR(4)*u_functions(:,k,4),&
          dimag(u_plus(:)) + aI(1)*u_functions(:,k,1) + aI(4)*u_functions(:,k,4))
     u_minus(:) = dcmplx(dreal(u_minus(:)) + aR(2)*u_functions(:,k,2)+ aR(4)*u_functions(:,k,5),&
          dimag(u_minus(:)) + aI(2)*u_functions(:,k,2) + aI(4)*u_functions(:,k,5))                    
     u_hat_mp%z(:,k) = dcmplx(dreal(u_hat_mp%z(:,k)) - aI(3)*u_functions(:,k,3) - aI(4)*u_functions(:,k,6),&
          dimag(u_hat_mp%z(:,k)) + aR(3)*u_functions(:,k,3) + aR(4)*u_functions(:,k,6))                  


     !undo change of variables
     u_hat_mp%r(:,k)  =  0.5D0   *(u_plus(:) + u_minus(:))
     u_hat_mp%th(:,k) = -0.5D0*ii*(u_plus(:) - u_minus(:))
     
  end do
       
  !$OMP END PARALLEL DO                                                                                  
       
  
      
  !free divergence conditions for the mode 0                                                             
  
       
  if(myid .eq. root) then                                                                                
     u_hat_mp%r(:,1)=dcmplx(0d0,0d0)                                                                     
     u_hat_mp%th(:,1)=dcmplx(dreal(u_hat_mp%th(:,1)),0d0)                                                 
     u_hat_mp%z(:,1)=dcmplx(dreal(u_hat_mp%z(:,1)),0d0)
     if (const_flux) call adjust_flux
  end if
    
  call perfoff()
       
end Subroutine predictor


Subroutine corrector
  
  IMPLICIT NONE
  INTEGER(kind=4) :: j,k
  COMPLEX(KIND=8) :: rhs_p(m_r),rhs_m(m_r),rhs_z(m_r)
  COMPLEX(kind=8) :: rhs_r(m_r),rhs_th(m_r),dotp(m_r)
  COMPLEX(KIND=8) :: u_plus(m_r),u_minus(m_r)
  COMPLEX(KIND=8) :: udu_plus(m_r),udu_minus(m_r)
  
  
  call perfon (' correc')
  
  !$OMP PARALLEL DO &
  !$OMP  SCHEDULE (RUNTIME) &
  !$OMP  DEFAULT(NONE) &
  !$OMP  PRIVATE(k,j,b,dotp,aR,aI,BRe,BIm,rhs_p,rhs_m,rhs_z,rhs_r,rhs_th, u_plus,u_minus,udu_plus,udu_minus) &
  !$OMP  SHARED(mp_f,m_r,Re,LU_p,Lu_up,Lu_um,Lu_uz,u_hat_mp,udu_hat_mp,fk_mp,p_hat_mp,bMw_dr,u_minus_old,&
  !$OMP  inf_matrix,u_functions,Lap_uz,Lap_up,Lap_um,r,dr1,udu_hat_mp_old,u_hat_mp_int,u_hat_mp_old,u_plus_old)   
  


  do k=1,mp_f
     
     !Save velocity at intermediate step
     
     u_hat_mp_int%r(:,k)   = u_hat_mp%r(:,k)
     u_hat_mp_int%th(:,k) = u_hat_mp%th(:,k)
     u_hat_mp_int%z(:,k)  = u_hat_mp%z(:,k)
     
     !multiply b by time-stepping coefficient (c NLq+1 + (1-c) NLq)
     
     udu_hat_mp%r(:,k)  = d_implicit*udu_hat_mp%r(:,k) + (1d0-d_implicit)* udu_hat_mp_old%r(:,k)
     udu_hat_mp%th(:,k) = d_implicit*udu_hat_mp%th(:,k) + (1d0-d_implicit)* udu_hat_mp_old%th(:,k)
     udu_hat_mp%z(:,k)  = d_implicit*udu_hat_mp%z(:,k) + (1d0-d_implicit)* udu_hat_mp_old%z(:,k)
       
       
     !Change of variables
       
     udu_plus(:)   = udu_hat_mp%r(:,k) + ii*udu_hat_mp%th(:,k)
     udu_minus(:)  = udu_hat_mp%r(:,k) - ii*udu_hat_mp%th(:,k)
       
     
       
     rhs_p(:)=DZGBMV(Lap_up(:,:,k),u_plus_old(:,k),iwidth,iwidth)
     rhs_m(:)=DZGBMV(Lap_um(:,:,k),u_minus_old(:,k),iwidth,iwidth)
     rhs_z(:)=DZGBMV(Lap_uz(:,:,k),u_hat_mp_old%z(:,k),iwidth,iwidth)
     
     
     !Adding non-linear terms
       
     
     rhs_p(:)=rhs_p(:)+udu_plus(:)
     rhs_m(:)=rhs_m(:)+udu_minus(:)
     rhs_z(:)=rhs_z(:)+udu_hat_mp%z(:,k)
       
     
     !Get pressure (project RHS)
     
          
     
     rhs_r(:)=  0.5D0 * (rhs_p(:) + rhs_m(:))
     rhs_th(:)=  -0.5D0*ii*(rhs_p(:) - rhs_m(:))
     
     !Compute the divergence of RHS
     
     
     
     dotp(:) = DZGBMV(bMw_dr,rhs_r(:),iwidth,iwidth)
     
     DO j=2,m_r-1
        
        b(j) = dotp(j) &
             + rhs_r(j)/r(j) &
             +ii*fk_mp%th(j,k)*rhs_th(j)/r(j)&
             +ii*fk_mp%z(j,k)*rhs_z(j)
        
     END DO
     
     !set boundary conditions to zero
     
     
     b(1)=dcmplx(0d0,0d0)
     b(m_r)=dcmplx(0d0,0d0)
     
     if (abs(fk_mp%th(1,k)) <= epsilon .and. abs(fk_mp%z(1,k)) <= epsilon) then
        
        p_hat_mp(:,k) = dcmplx(0d0,0d0)
     else
        
        ! solve system using LU decomposition precomputed in pre_timeStep()
        
        call luSolve(LU_p(k),b,p_hat_mp(:,k))
        
     end if
     
     !Compute pressure gradient and substract from rhs (written directly as rhs_p and rhs_m)
     
     
     
     
     dotp(:)=DZGBMV(bMw_dr,p_hat_mp(:,k),iwidth,iwidth)
     rhs_p(:)=rhs_p(:)-dotp(:)+fk_mp%th(:,k)*p_hat_mp(:,k)/r(:)
     rhs_m(:)=rhs_m(:)-dotp(:)-fk_mp%th(:,k)*p_hat_mp(:,k)/r(:)
     rhs_z(:)=rhs_z(:)-ii*fk_mp%z(:,k)*p_hat_mp(:,k)
     
     
     !Set boundary conditions to zero
     
     
     rhs_p(1)=dcmplx(0d0,0d0)
     rhs_p(m_r)=dcmplx(0d0,0d0)
     rhs_m(1)=dcmplx(0d0,0d0)
     rhs_m(m_r)=dcmplx(0d0,0d0)
     rhs_z(1)=dcmplx(0d0,0d0)
     rhs_z(m_r)=dcmplx(0d0,0d0)
     
     
     !Solve equations                                                                                                                                                                                                                                                     
     
     call luSolve(LU_up(k),rhs_p(:),u_plus(:))
     call luSolve(LU_um(k),rhs_m(:),u_minus(:))
     call luSolve(LU_uz(k),rhs_z(:),u_hat_mp%z(:,k))
     
     
     
     !correct boundary conditions with influence matrix                                                                                                                                                                                                                  
     if (abs(fk_mp%th(1,k)) <= epsilon .and. abs(fk_mp%z(1,k)) <= epsilon) cycle
     
     call evalBC(k, u_plus(:),u_minus(:),u_hat_mp%z(:,k),BRe,BIm)
     
     aR(:) = -matmul(inf_matrix(:,:,k),BRe(:))
     aI(:) = -matmul(inf_matrix(:,:,k),BIm(:))
     
     u_plus(:) = dcmplx(dreal(u_plus(:)) + aR(1)*u_functions(:,k,1)+ aR(4)*u_functions(:,k,4),&                            
          dimag(u_plus(:)) + aI(1)*u_functions(:,k,1) + aI(4)*u_functions(:,k,4))                                           
     u_minus(:) = dcmplx(dreal(u_minus(:)) + aR(2)*u_functions(:,k,2)+ aR(4)*u_functions(:,k,5),&                          
          dimag(u_minus(:)) + aI(2)*u_functions(:,k,2) + aI(4)*u_functions(:,k,5))                                          
     u_hat_mp%z(:,k) = dcmplx(dreal(u_hat_mp%z(:,k)) - aI(3)*u_functions(:,k,3) - aI(4)*u_functions(:,k,6),&                   
          dimag(u_hat_mp%z(:,k)) + aR(3)*u_functions(:,k,3) + aR(4)*u_functions(:,k,6))                                       


     !undo change of variables
     u_hat_mp%r(:,k)  =  0.5D0   *(u_plus(:) + u_minus(:))
     u_hat_mp%th(:,k) = -0.5D0*ii*(u_plus(:) - u_minus(:))
     
     
  end do
  
    !$OMP END PARALLEL DO                                                                                                  
  
                                                                                
                                                                   
    !free divergence conditions for the mode 0                                                                             
    
    
    if(myid .eq. root) then                                                                                                
       u_hat_mp%r(:,1)=dcmplx(0d0,0d0)                                                                                     
       u_hat_mp%th(:,1)=dcmplx(dreal(u_hat_mp%th(:,1)),0d0)                                                                 
       u_hat_mp%z(:,1)=dcmplx(dreal(u_hat_mp%z(:,1)),0d0)
       if (const_flux) call adjust_flux
    end if
    
    
    !Check whether the difference is within the requested tolerance                                                                                                                                                                                                         
    call measurecorr(u_hat_mp%r, u_hat_mp%th, u_hat_mp%z, u_hat_mp_int%r, u_hat_mp_int%th, u_hat_mp_int%z)
    call perfoff()
    
end Subroutine corrector



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
                                                                                                            

  !---------------------------------------------------------------------                                                                                                                                                                                                       
  subroutine luDecmp(n,k,A,LU)
  !---------------------------------------------------------------------                                                                                                                                                                                                       
  ! Compute LU decomposition for nxn band matrix A                                                                                                                                                                                                                             
  ! variant for calling LA routines DGBTRF, DGBTRS of the LAPACK library                                                                                                                                                                                                       
  !                                                                                                                                                                                                                                                                            
  ! input:  A   banded matrix of dimension n (LAPACK storage)                                                                                                                                                                                                                  
  ! input:  n   dimension of A                                                                                                                                                                                                                                                 
  ! input:  k   number of sub/super-diagonals within the band of A                                                                                                                                                                                                             
  ! output: LU  structure encapsulating LU decomposition to be passed to                                                                                                                                                                                                       
  !             solver routine luSolve                                                                                                                                                                                                                                         
  !---------------------------------------------------------------------                                                                                                                                                                                                       
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN)  :: n,k
    REAL(KIND=8), INTENT(IN)     :: A(2*k+1,n)
    INTEGER(KIND=4)              :: i,j
    TYPE(LUdcmp),  INTENT(OUT)   :: LU

    if (.not. allocated(LU%LU)) allocate(LU%LU(3*k+1,n))                                                                                                                                                                                                                                
    if (.not. allocated(LU%ipiv)) allocate(LU%ipiv(n))                                                                                                                                                                                                                        

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
    
   
    if (LU%info /= 0) STOP 'luDecmp: dgbtrf'
    

    RETURN
  end subroutine luDecmp

 !---------------------------------------------------------------------                                                                                                                                                                                                  
  subroutine luSolve(LU,b,x)
  !---------------------------------------------------------------------                                                                                                                                                                                                  
  ! solve the linear system A*x = b using a precomputed LU decomposition                                                                                                                                                                                                  
  ! from routine luDecmp                                                                                                                                                                                                                                                  
  ! variant for calling LA routines DGBTRF, DGBTRS of the LAPACK library                                                                                                                                                                                                  
  !                                                                                                                                                                                                                                                                       
  ! input: LU  structure encapsulating LU decomposition computed by                                                                                                                                                                                                       
  !            routine luDecmp                                                                                                                                                                                                                                            
  ! input: b   right-hand side vector (complex)                                                                                                                                                                                                                           
  ! input: x   solution vector (complex)                                                                                                                                                                                                                                  
  !                                                                                                                                                                                                                                                                       
  !---------------------------------------------------------------------                                                                                                                                                                                                  
    IMPLICIT NONE
    COMPLEX(KIND=8), INTENT(IN)  :: b(:)
    COMPLEX(KIND=8), INTENT(OUT) :: x(:)
    TYPE(LUdcmp),  INTENT(IN)    :: LU
    REAL(KIND=8)                 :: par(LU%n,2)
    INTEGER(KIND=4)              :: i,info


    par(:,1) = DREAL(b)
    par(:,2) = DIMAG(b)
    CALL dgbtrs('N',LU%n,LU%k,LU%k,2,LU%LU,LU%lda,LU%ipiv,par,LU%ldb,info)
    if (info /= 0) STOP 'luSolve: dgbtrs at REAL or COMPLEX part'

    DO i=1,LU%n
       x(i) = DCMPLX(par(i,1),par(i,2))
    END DO

    RETURN
  end subroutine luSolve
  
 SUBROUTINE evalBC(k, vr, vth, vz, BRe, BIm)

    IMPLICIT NONE
    INTEGER(kind=4), intent(in) :: k
    complex(kind=8), dimension(m_r), intent(in) :: vr, vth, vz 
    real(kind=8), intent(out)   :: BRe(4),BIm(4)
    real(kind=8) :: drRe,drIm, urRe,urIm, utRe,utIm, uzRe,uzIm
    real(kind=8) :: d
    
    
    urRe =  0.5d0 * (dreal(vr(m_r)) + dreal(vth(m_r)))
    urIm =  0.5d0 * (dimag(vr(m_r)) + dimag(vth(m_r)))
    utRe =  0.5d0 * (dimag(vr(m_r)) - dimag(vth(m_r)))
    utIm = -0.5d0 * (dreal(vr(m_r)) - dreal(vth(m_r)))
    uzRe =  dreal(vz(m_r))
    uzIm =  dimag(vz(m_r))
    drRe =  0.5d0 * dot_product(dr1(:,1), dreal(vr(m_r-(n_s-1)/2:m_r))+dreal(vth(m_r-(n_s-1)/2:m_r)))
    drIm =  0.5d0 * dot_product(dr1(:,1), dimag(vr(m_r-(n_s-1)/2:m_r))+dimag(vth(m_r-(n_s-1)/2:m_r)))
    
    !evaluate no-slip b.c.                                                                                                                                                                                                                                              
    BRe(1) = dreal(vr(m_r))
    BIm(1) = dimag(vr(m_r))
    BRe(2) = dreal(vth(m_r))
    BIm(2) = dimag(vth(m_r))
    BRe(3) = -uzIm
    BIm(3) =  uzRe
    d = 1d0/r(m_r)
    
    !Evaluate continuity equation                                                                                                                                                                                                                                        
    BRe(4) = urRe*d + drRe - d*fk_mp%th(m_r,k)*utIm - fk_mp%z(m_r,k)*uzIm
    BIm(4) = urIm*d + drIm + d*fk_mp%th(m_r,k)*utRe + fk_mp%z(m_r,k)*uzRe
    

    
  END SUBROUTINE evalBC
  
  !DEC$ ATTRIBUTES FORCEINLINE :: DZGBMV                                                                                                                                                                                                                                    
 FUNCTION DZGBMV(A,x,kl,ku) RESULT(v)

!------------------------------------------------------------------                                                                                                                                                                                                       
! Computes a matrix-vector product of a real, banded Matrix with a                                                                                                                                                                                                        
! complex vector.                                                                                                                                                                                                                                                         
!                                                                                                                                                                                                                                                                         
! Input:  A(kl+ku+1,:) band matrix (real) in BLAS-style band storage                                                                                                                                                                                                      
!         x(:)         dense vector (complex)                                                                                                                                                                                                                             
!         kl,ku        number of lower and upper diagonals                                                                                                                                                                                                                
! Output: v(:)         result (complex) of v=A*x                                                                                                                                                                                                                          
!                                                                                                                                                                                                                                                                         
! Implementation note:                                                                                                                                                                                                                                                    
!  Although this BLAS-2 operation is naively implemented, we are                                                                                                                                                                                                          
!  usually faster (n~1000) as compared with promoting A to double                                                                                                                                                                                                         
!  complex and calling a highly optimized ZGBMV                                                                                                                                                                                                                           
!------------------------------------------------------------------                                                                                                                                                                                                       

    implicit none

    real(kind=8), intent(in) :: A(:,:)
    complex(kind=8), intent(in) :: x(:)
    integer(kind=4), intent(in) :: kl,ku 
    
    complex(kind=8) :: v(SIZE(x))
    
    integer(kind=4) :: n,i,j,k,ilo,ihi

    n=size(x)
    
    v(:)=dcmplx(0.d0,0.d0)
    do j=1,n
       k=ku+1-j
       ilo=max(j-ku,1)
       ihi=min(j+kl,n)
       
       
       
     DO i=ilo,ihi
        v(i) = v(i) + A(k+i,j)*x(j)
     END DO
  enddo
  
END FUNCTION DZGBMV

!  measure the magnitude of the correction: max diff                                                                                                                                                                                                                      
subroutine measurecorr(c1,c2,c3, c1_,c2_,c3_)
  implicit none
  complex(kind=8), dimension(m_r,mp_f), intent(in)  :: c1,c2,c3,c1_,c2_,c3_
  real(kind=8), dimension(2) :: d, d_
  integer(kind=4) :: nh

  call perfon('mcorr')

  d(:) = -1d9

  !$OMP PARALLEL &
  !$OMP  DEFAULT(NONE) &
  !$OMP  PRIVATE(nh) &
  !$OMP  SHARED(mp_f,d,c1,c2,c3,c1_,c2_,c3_)

  !$OMP DO SCHEDULE (RUNTIME) &
  !$OMP  REDUCTION(max:d)
  do nh = 1,mp_f

     d(1) = max(d(1),  &
          maxval( dreal(c1(:,nh))*dreal(c1(:,nh))+dimag(c1(:,nh))*dimag(c1(:,nh))),  &
          maxval( dreal(c2(:,nh))*dreal(c2(:,nh))+dimag(c2(:,nh))*dimag(c2(:,nh))),  &
          maxval( dreal(c3(:,nh))*dreal(c3(:,nh))+dimag(c3(:,nh))*dimag(c3(:,nh))))

     d(2) = max(d(2),  &
          maxval( (dreal(c1(:,nh))-dreal(c1_(:,nh)))*(dreal(c1(:,nh))-dreal(c1_(:,nh)))&
          +(dimag(c1(:,nh))-dimag(c1_(:,nh)))*(dimag(c1(:,nh))-dimag(c1_(:,nh)))),  &
          maxval( (dreal(c2(:,nh))-dreal(c2_(:,nh)))*(dreal(c2(:,nh))-dreal(c2_(:,nh)))&
          +(dimag(c2(:,nh))-dimag(c2_(:,nh)))*(dimag(c2(:,nh))-dimag(c2_(:,nh)))),  &
          maxval( (dreal(c3(:,nh))-dreal(c3_(:,nh)))*(dreal(c3(:,nh))-dreal(c3_(:,nh)))&
          +(dimag(c3(:,nh))-dimag(c3_(:,nh)))*(dimag(c3(:,nh))-dimag(c3_(:,nh)))))

  end do
  !$OMP END DO
  !$OMP END PARALLEL


  call mpi_allreduce(d, d_, 2, MPI_REAL8, MPI_MAX, COMM, ierr)
  d = d_

  if (d(1) .ne. 0) then
     dterr = max(dterr,dsqrt(d(2)/d(1)))
  else
     dterr=0d0
  endif

  call perfoff()

  return
end subroutine measurecorr
                                                                                      

!-------------------------------------------------------------------------                                                                                                                                                                                                
!  check for convergence via the 2-norm of the correction                                                                                                                                                                                                                 
!  -- see tim_addcorr()                                                                                                                                                                                                                                                   
!-------------------------------------------------------------------------                                                                                                                                                                                                
  subroutine check_convergence()
    implicit NONE
    real(kind=8), save :: lasterr
    
    !Measure correction of dt                                                                                                                                                                                                                                           
    
    if(iter .eq. 1) then
       if (dterr .eq. 0.0d0) then 
          corr_dt = 0d0
       else
          corr_dt =  dt * dsqrt( tolerance_dterr/dterr )
       end if
       lasterr = 1d99
    end if
      
    !Abort simulation if dt is very small                                                                                                                                                                                                                               
    
    if(dt .lt. 1d-9 .and. i_time .gt. 30) then
       if(myid .eq. root) print*, 'check_convergence: dt --> 0 !!!?'
       CALL MPI_ABORT(comm,102,ierr)
       
       !Abort simulation if number of iterations is large                                                                                                                                                                                                               
    else if(iter .gt. 10) then
       if(myid .eq. root) print*, 'check_convergence: too many its!!!'
       CALL MPI_ABORT(comm,103,ierr)
       !Abort simulation if error is increasing  or change correction time?                                                                                                                                                                                             
    else if(dterr .gt.  lasterr) then
       if(myid .eq. root)  print*, 'check_convergence: increasing error!!!'
       if(dterr .gt. 2d0*tolerance_dterr) CALL MPI_ABORT(comm,104,ierr)
       if(dterr .lt. 2d0*tolerance_dterr) corr_dt = dt/(1d1*Courant)
       iter = 0
    else if(dterr .gt. tolerance_dterr) then
       lasterr = dterr
       iter = iter + 1
    else
       if(myid .eq. 0 .and. modulo(i_time,print_time_screen) .eq. 0) then

          if(variable_dt) then
             print*,' step=',i_time,' dt=', dt
          else
             print*,' step=',i_time,' its=',iter
          endif
       end if
       iter = 0
    end if
    dterr = 0d0
    
  end subroutine check_convergence

!-------------------------------------------------------------------------                                                                                                                                                                                                
!  set a new timestep size, if necessary.                                                                                                                                                                                                                                 
!-------------------------------------------------------------------------                                                                                                                                                                                                
    subroutine new_tstep()
      implicit none
      real(kind=8)  :: deltat
      integer(kind=4), save :: ind = 0

      deltat = min(dt*1.11d0, maxdt)
      if(i_time .eq. 1 .and. corr_dt .eq. 0d0) deltat = min(deltat, cfl*0.1d0)
      if(cfl .gt. 0d0)  deltat = min(deltat, cfl*Courant)
      if(corr_dt .gt. 0d0)  deltat = min(deltat, corr_dt*0.95d0)

      ind = ind - 1
      new_dt = (deltat .lt. dt*0.95d0 .or. (ind .lt. 0 .and. deltat .gt.  dt*1.10d0))
      if(new_dt) then
         if(deltat .gt. dt*1.10d0) ind = 100
         dt = deltat
      end if
      
    end subroutine new_tstep

    Subroutine adjust_flux()

      IMPLICIT NONE


      REAL(kind=8) :: d,d_,par_r(m_r)

      !Mean flow for perturbation must be zero for fixed mass flux

      ! compute flow rate
      d_ = 2d0*dot_product(dreal(u_hat_mp%z(:,1)), intrdr)
      d = (0.5d0-d_)/c_adj

      
      par_r = dreal(u_hat_mp%z(:,1)) + d*adj_flux

      !Add laminar flow


      u_hat_mp%z(:,1)= dcmplx(par_r,0d0)


    end Subroutine adjust_flux
    

                                                                                                                                                                                        
END MODULE mod_timeStep

