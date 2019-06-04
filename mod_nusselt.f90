
subroutine init_nusselt()
! Purpose:
! Initialise time series output of the Nusselt numbers at the inner and outer
! cylinder walls. In case of the hydrodynamic code version, this includes only
! the torque Nusselt number (actual torque normalised by the torque of the
! analytical Couette profile), whereas in case of the code with temperature
! this also includes the heat transfer at both cylinder walls normalised by the
! conductive heat transfer. 
! Called from:
! + nsCouette.f90 
! Contet:
! + Fix i/o units for all probe locations
! + Broadcast information
! + Open one file per probe location, so that each rank hosting one or more
!   probes outputs to its own file(s). I guess this is better than collecting
!   all probe data from different ranks via mpi_sent/receive/bcast and than
!   do all the output from one rank to one file, right?
! + Write header information to file in case of REWIND

implicit none
character(len=11) :: nus_fnam ! nusselt file name
character(len=6) :: char_position ! pointer position at i/o

if (myid .eq. root) then

 ! Fix i/o units (grep -i "unit" *.f90)
 nus_io = 105

end if

! Broadcast information
call mpi_bcast(nus_io, prl_n, mpi_integer4, root, comm, ierr)

if (myid .eq. root) then

  ! Open files
  write(nus_fnam, 1001) ! Set probe file name
  if ((restart .eq. 0) .or. (restart .eq. 2)) then ! Access mode
   char_position = 'REWIND'
   open(unit=nus_io(i), file=nus_fnam, position=char_position)

   ! Write header
#ifndef TE_CODE
   write(nus_io, 1002) nus_rnk
#else
#endif

  else
   char_position = 'APPEND'
   open(unit=nus_io(i), file=nus_fnam, position=char_position)
