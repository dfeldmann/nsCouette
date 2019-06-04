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

!=========================================
!         Tranpose matrix by MPI:
!          in(m,n) -> out(n,m)
!
! SVN $HeadURL$
! SVN $Id$
! SVN $LastChangedDate$
!=========================================


MODULE mod_myMpi

  USE mpi
  USE mod_vars
  USE mod_params
  IMPLICIT NONE

  INTEGER(KIND=4), PRIVATE :: block_lengths(2),types(2),extent,np
  INTEGER(KIND=MPI_ADDRESS_KIND), PRIVATE :: displacements(2)

CONTAINS
  !------------------------------------------------------------
  SUBROUTINE pre_mpi()
    IMPLICIT NONE

    comm = MPI_COMM_WORLD
    CALL MPI_COMM_SIZE (comm, numprocs, ierr)
    CALL MPI_COMM_RANK (comm, myid, ierr)
    
   
  end SUBROUTINE pre_mpi

  subroutine init_mpidist
    IMPLICIT NONE
    INTEGER(KIND=4) :: ndims=2
    INTEGER(KIND=4) :: sizes(2),subsizes(2),starts(2),i
#ifdef TE_CODE
    INTEGER(KIND=4), PARAMETER :: blocks=5
#else
    INTEGER(KIND=4), PARAMETER :: blocks=4
#endif

    

    ! Check the divisibility of dimensions
    IF ((MOD(m_r,numprocs) /= 0)) THEN
       if (myid==root) PRINT*, 'number of radial modes',m_r,' must be divisible by number of MPI procs',numprocs
       CALL MPI_ABORT(comm,1,ierr)
    END IF

    mp_r=m_r/numprocs ! Radial points at each proc
    mp_fmax=ceiling(real(m_f,kind=8)/real(numprocs,kind=8)) ! Fourier points per proc
    np=numprocs

    if ((MOD(m_f,numprocs) /= 0)) then

       extended_mf_grid=.true.
       mp_f=(numprocs+mp_fmax*numprocs-1)/numprocs
       !this is the OpenMP SCHEDULE(STATIC)-like way without zero load for any task     (e.g. 11=2+3+3+3)
       !we cannot use it since we have to rely on a 'compact' population of the extended grid here
       !if (myid .lt. mp_f*numprocs-m_f) mp_f=mp_f-1
       !this is the canonical MPI-like way that can imply zero load for the 'last' task (e.g. 11=4+4+3+0)
       mp_f=max(0,min(mp_f,m_f-myid*mp_f))
       m_f=mp_fmax*numprocs

       if (myid==root) then
          print '(A)',                     '------------------------------------------------------------------------'
          print '(A,I9,A,I5)',             'INFO: chosen number of Fourier modes (m_th+1)*m_z=',&
               &(m_th+1)*m_z,' is not divisible by number of MPI procs',numprocs
          print '(A,I9)',                  '   => m_f grid enlarged to',m_f
       endif
    else
       extended_mf_grid=.false.
       mp_f=mp_fmax
    endif

!a global view on the distribution of fourier modes
    allocate(mp_f_arr(numprocs)) 
    call mpi_allgather(mp_f,1,MPI_INTEGER,mp_f_arr,1,MPI_INTEGER,comm,ierr)


    if (extended_mf_grid) then
       if (myid==root) then
          do i=1,numprocs
             print '(A,I5,A,I5,A,I5,A,I6)', 'MPI proc=',i,', mp_f=',mp_f_arr(i),', mp_fmax=',mp_fmax,', m_f=',m_f
          enddo
          print '(A,f7.2,A)','waste due to imbalance: ',100.d0*(1.d0-SUM(mp_f_arr)/real(m_f,kind=8)),' %'
          print '(A)',                     '------------------------------------------------------------------------'
       endif
    endif

    ! final consistency check
    if (SUM(mp_f_arr) .ne. (m_th+1)*m_z) then
       if (myid==root) print '(A,I6,A,I6)', 'mismatch in total number of active Fourier modes, is:',&
            &SUM(mp_f_arr),'expected: ',(m_th+1)*m_z
       CALL MPI_ABORT(comm,1,ierr)
    endif


    ierr=u_hat_mp%alloc(m_r,mp_f)
    ierr=udu_hat_mp%alloc(m_r,mp_f)
    ierr=fk_mp%alloc(m_r,mp_f)
    allocate(p_hat_mp(m_r,mp_f))



    ! Define new MPI type
    CALL MPI_TYPE_CONTIGUOUS (2, MPI_DOUBLE_PRECISION, & 
         double_complex, ierr)
    CALL MPI_TYPE_EXTENT(double_complex,extent,ierr)
    
    block_lengths = (/blocks,2/)
    displacements = (/0,blocks*extent/)

    types = (/double_complex,MPI_REAL8/)
    CALL MPI_Type_create_struct(2,block_lengths,displacements,types,mpi_spec,ierr)
    CALL MPI_TYPE_COMMIT(mpi_spec, ierr)

    sizes = (/m_r,m_f/)
    subsizes = (/m_r,mp_fmax/)
    starts = (/0,myid*mp_fmax/) ! Originally from (1, myid*mp_f)
    CALL MPI_Type_create_subarray(ndims,sizes,subsizes,starts,&
         MPI_ORDER_FORTRAN,mpi_spec,filetype,ierr)
    CALL MPI_Type_commit(filetype,ierr)

  END SUBROUTINE init_mpidist


  !------------------------------------------------------------
  subroutine xTranspose_mpi(m,n,in,out,ifile)

    implicit none
  ! preliminary wrapper routine for cases where MOD(m_f,numprocs)!=0
  ! extends the local arrays and calls the original xTranspose_mpi implementation

    integer(KIND=4), intent(in)  :: m,n
    complex(KIND=8), intent(in)  :: in(:,:)
    complex(KIND=8), intent(out) :: out(:,:)
    logical,intent(in), optional :: ifile
    complex(KIND=8) :: tmp_i(m,n/numprocs),tmp_o(n,m/numprocs)
    integer(KIND=4) :: shape_i(2),shape_o(2)
    logical :: extend

    if (present(ifile)) then
       if (ifile) then
          extend=(mp_f_ifile .ne. mp_fmax_ifile)
       else
          extend=(mp_f .ne. mp_fmax)
       endif
    else
       extend=(mp_f .ne. mp_fmax)
    endif

    if (extend) then
       shape_i=shape(in)
       shape_o=shape(out)
                   
       tmp_i(1:shape_i(1),1:shape_i(2))=in(1:shape_i(1),1:shape_i(2))
       call xTranspose_mpi_impl(m,n,tmp_i,tmp_o)
       out(1:shape_o(1),1:shape_o(2))=tmp_o(1:shape_o(1),1:shape_o(2))   
    else
       call xTranspose_mpi_impl(m,n,in,out)
    endif


    
  end subroutine xTranspose_mpi

  !------------------------------------------------------------
  SUBROUTINE xTranspose_mpi_impl(M,N,in,out)

  ! written by Florian Merz (IBM,RZG)

    implicit none

    ! Declarations
    integer(KIND=4), intent(in)  :: M,N
    complex(KIND=8), intent(in)  :: in(M,N/np)
    complex(KIND=8), intent(out) :: out(N,M/np)
    complex(KIND=8) :: temp(N/np,M),temp2(N/np,M/np,np)
    integer:: blocksize, n_bl, ierr, pe, i

    blocksize=(N/np)*(M/np)
    n_bl=N/np

    if (np.eq.1) then
       !only local transpose necessary
       call local_transp(M,n_bl,in,out)
    elseif (np.eq.N) then
       !processor transpose (=alltoall) and reordering
       !local transpose is not needed because N/np=1
       call mpi_alltoall(in, blocksize, MPI_DOUBLE_COMPLEX,&
            temp2, blocksize, MPI_DOUBLE_COMPLEX,&
            MPI_COMM_WORLD, ierr)
       call local_transp(M/np,np,temp2,out)
    elseif (np.eq.M) then
       !local transpose, processor transpose (=alltoall)
       !reordering is not needed because M/np=1
       call local_transp(M,n_bl,in,temp)
       call mpi_alltoall(temp, blocksize, MPI_DOUBLE_COMPLEX,&
            out, blocksize, MPI_DOUBLE_COMPLEX,&
            MPI_COMM_WORLD, ierr)
    else
       !local transpose, processor transpose (=alltoall) and reordering
       call local_transp(M,n_bl,in,temp)
       call mpi_alltoall(temp, blocksize, MPI_DOUBLE_COMPLEX,&
            temp2, blocksize, MPI_DOUBLE_COMPLEX,&
            MPI_COMM_WORLD, ierr)
       do i=1,M/np
          do pe=1,np
             out((pe-1)*n_bl+1:pe*n_bl,i)=temp2(:,i,pe)
          enddo
       enddo
    end if

  END SUBROUTINE xTranspose_mpi_impl

  !------------------------------------------------------------
  SUBROUTINE local_transp(d1,d2,loc_in,loc_out)

    integer:: d1,d2
    complex(KIND=8):: loc_in(d1,d2), loc_out(d2,d1), scale=(1.d0,0.d0)
  
#if defined(WITHESSL)
    call zgetmo(loc_in,d1,d1,d2,loc_out,d2)
#elif defined(WITHMKL)
!mjr: do _not_ use mkl_zomatcopy with MKL versions prior to 11.1.1 
!     because of its inferiour performance wrt. F90 intrinsic transpose()
!     the issue is adressed first in MKL 11.1.1 Product Build 20131010
!     but there is still no clear performance benefit of MKL vs. transpose()
    call mkl_zomatcopy('C','T',d1,d2,scale,loc_in,d1,loc_out,d2)
#else
    loc_out=transpose(loc_in)
#endif
    
  END SUBROUTINE local_transp

  !---------------------------------------------------
  SUBROUTINE post_mpi()
    IMPLICIT NONE

    ! Free data types
    CALL MPI_Type_free(double_complex, ierr)
    CALL MPI_Type_free(mpi_spec,ierr)
    CALL MPI_Type_free(filetype,ierr)

  END SUBROUTINE post_mpi
  !----------------------------------------------------END MODULE---------

END MODULE mod_myMpi
