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

! perflib wrapper for likwid marker-API calls
! see 
! https://github.com/RRZE-HPC/likwid
! https://github.com/RRZE-HPC/likwid/wiki/TutorialMarkerF90

subroutine perfinit
  use likwid
  call likwid_markerInit()

!$OMP PARALLEL
!$   call likwid_markerthreadInit()
!$OMP END PARALLEL

end subroutine perfinit

subroutine perfon(region)
  use likwid
  character(*), intent(in) :: region
  call likwid_markerStartRegion(region)
end subroutine perfon

subroutine perfoff(region)
  use likwid
  character(*), intent(in) :: region
  call likwid_markerStopRegion(region)
end subroutine perfoff

subroutine perfout()
  use likwid
  call likwid_markerClose()
end subroutine perfout


subroutine perf_context_start()

end subroutine perf_context_start

subroutine perf_context_end()

end subroutine perf_context_end
