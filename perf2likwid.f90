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
