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
