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

module mod_getcpu

  interface
    function threadGetProcessorId( ) bind( c, name="threadGetProcessorId" )
      use, intrinsic :: iso_c_binding
      integer( kind = c_int )        :: threadGetProcessorId
    end function threadGetProcessorId
  end interface

  interface
    function where_am_i_running( ) bind( c, name="where_am_i_running" )
      use, intrinsic :: iso_c_binding
      integer( kind = c_int )        :: where_am_i_running
    end function where_am_i_running
  end interface

contains

  subroutine mapping_info(numprocs,myrank)

    IMPLICIT NONE

    integer, intent(in) :: numprocs,myrank 

    integer, allocatable :: procids(:)
    integer numthreads,threadnum,ic
!$  integer omp_get_max_threads,omp_get_thread_num
    character(16) hostname
    character(3) fmt
    numthreads=1 !fallback
    threadnum=0 !fallback

!$ numthreads=omp_get_max_threads()
    call hostnm(hostname)
    allocate(procids(0:numthreads-1))
    procids(:)=-1
!$OMP PARALLEL DO ORDERED SCHEDULE(static,1) &
!$OMP   PRIVATE(ic,threadnum) SHARED(procids)
    do ic=1,numthreads
!$     threadnum=omp_get_thread_num()
       procids(threadnum)=threadGetProcessorId()
    enddo


    write(fmt,'(I3)') numthreads
!    write(*,*)
!    write(*,*)
!    write(*,*) '-----------------------------------------------------------------------------'
    write(*,'(A9,I4,A9,I2,A9,I4,A10,A10,A7,'//fmt//'I3)') 'TASKS:',numprocs,&
         &'THREADS:',numthreads,'THIS:',myrank,' Host:',hostname,' Cores:',procids
!    write(*,*) '-----------------------------------------------------------------------------'

    deallocate(procids)

  end subroutine mapping_info

end module mod_getcpu
