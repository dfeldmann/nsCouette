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

! simple wrappers for mapping perflib calls to the ftimings
! library (by L. Huedepohl, RZG) 
subroutine perfinit
  use timings

  call timer%enable()
end subroutine perfinit

subroutine perfon(label)
  use timings
  character(*), intent(in) :: label

  call timer%start(trim(adjustl(label)))
end subroutine perfon

subroutine perfoff(label)
  use timings
  character(*), intent(in) :: label

  call timer%stop()
end subroutine perfoff

subroutine perfout(label)
  use timings
  character(*), intent(in) :: label

  call timer%print()
end subroutine perfout

subroutine perf_context_start(label)
  character(*), intent(in) :: label

end subroutine perf_context_start

subroutine perf_context_end()

end subroutine perf_context_end

