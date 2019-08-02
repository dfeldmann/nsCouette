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
! function c_to_f_string and subroutine hostname are adapted from
! https://www.rosettacode.org/wiki/Hostname#Fortran

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

  interface !to function: int gethostname(char *name, size_t namelen);
     integer(c_int) function gethostname(name, namelen) bind(c)
       use, intrinsic  :: iso_c_binding, only: c_char, c_int, c_size_t
       integer(c_size_t), value, intent(in) :: namelen
       character(len=1,kind=c_char), dimension(namelen),  intent(inout) ::  name
     end function gethostname
  end interface

contains

  subroutine mapping_info(numprocs,myrank)

    IMPLICIT NONE

    integer, intent(in) :: numprocs,myrank 

    integer, allocatable :: procids(:)
    integer numthreads,threadnum,ic
!$  integer omp_get_max_threads,omp_get_thread_num
    character(16) hostnam
    character(3) fmt
    numthreads=1 !fallback
    threadnum=0 !fallback

!$ numthreads=omp_get_max_threads()
    call hostname(hostnam)
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
    write(*,'(A9,I4,A9,I2,A9,I4,A10,A10,A7,'//fmt//'I4)') 'TASKS:',numprocs,&
         &'THREADS:',numthreads,'THIS:',myrank,' Host:',hostnam,' Cores:',procids
!    write(*,*) '-----------------------------------------------------------------------------'

    deallocate(procids)

  end subroutine mapping_info

  pure function c_to_f_string(c_string) result(f_string)
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char
    character(kind=c_char,len=1), intent(in) :: c_string(:)
    character(len=:), allocatable :: f_string
    integer i, n
    i = 1
    do
       if (c_string(i) == c_null_char) exit
       i = i + 1
    end do
    n = i - 1  ! exclude c_null_char
    allocate(character(len=n) :: f_string)
    f_string = transfer(c_string(1:n), f_string)
  end function c_to_f_string

  subroutine hostname(name)

    use, intrinsic  :: iso_c_binding, only: c_char, c_int, c_size_t
    character(*), intent(out) :: name
    integer(c_int) :: status
    integer,parameter :: HOST_NAME_MAX=255
    integer(c_size_t) :: lenstr
    character(kind=c_char,len=1),dimension(HOST_NAME_MAX) :: cstr_hostname
    lenstr=HOST_NAME_MAX
    status = gethostname(cstr_hostname,lenstr)
    name = c_to_f_string(cstr_hostname)

  end subroutine hostname

end module mod_getcpu
