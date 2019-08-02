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

MODULE WCTimerClass 
!      Defines a class for timing Fortran program  execution
!
!      Based on FTTimerClass.f90 by NocturnalAviationSoftware 
!
!      Original code downloaded by mjr 2017-10-19 from 
!      http://www.nocturnalaviationsoftware.com/blog/using-type-bound-procedures.html
!      the webpages says: "... You can download the source code for the Class here. You are free to use it as you like."
!
!      exchanged (largely useless) cpu_time measurement by a wall_clock timer + further adaptations
!      
!      Usage
!
!         * Starting the timer: CALL timer%start()
!         * Stopping the timer: CALL timer%stop()
!         * Reading the time:   time = timer%elapsedTime([units])
!             units (optional) = TC_SECONDS or TC_MINUTES or TC_HOURS
!             timer%elapsedTime() does not stop the timer
      IMPLICIT NONE
      PRIVATE

      INTEGER, PARAMETER, PRIVATE :: d = 15

      INTEGER, PARAMETER, PUBLIC  :: TP = SELECTED_REAL_KIND(d)
      INTEGER, PARAMETER, PUBLIC  :: TC_SECONDS = 0, TC_MINUTES = 1, TC_HOURS = 2

      TYPE, PUBLIC :: WCTimer
         LOGICAL      , PRIVATE :: started = .FALSE.,stopped = .FALSE.,initialized= .FALSE.
         REAL(KIND=TP), PRIVATE :: startTime = 0.0_TP
         REAL(KIND=TP), PRIVATE :: finishTime = 0.0_TP
         INTEGER, PRIVATE :: clock_rate,clock_max


         CONTAINS

         PROCEDURE, PASS :: start => startTimer
         PROCEDURE, PASS :: stop  => stopTimer
         PROCEDURE, PASS :: elapsedTime
         
      END TYPE WCTimer

      CONTAINS

!////////////////////////////////////////////////////////////////////////  
!  
      SUBROUTINE initTimer(self)  
         IMPLICIT NONE

         CLASS(WCTimer) :: self

         call system_clock (count_rate=self%clock_rate, count_max=self%clock_max)
         self%initialized = .TRUE.
       END SUBROUTINE initTimer


!////////////////////////////////////////////////////////////////////////  
!  
      FUNCTION queryTimer() result(clock)
         IMPLICIT NONE
         integer :: clock

         call system_clock (count=clock)

      END FUNCTION queryTimer

!////////////////////////////////////////////////////////////////////////  
!  
      SUBROUTINE startTimer(self)  
         IMPLICIT NONE

         CLASS(WCTimer) :: self

         if (.not.self%initialized) call initTimer(self)
         self%started = .TRUE.
         self%startTime=queryTimer()
      END SUBROUTINE startTimer
!
!////////////////////////////////////////////////////////////////////////  
!  
      SUBROUTINE stopTimer(self)  
         IMPLICIT NONE
         CLASS(WCTimer) :: self

         self%finishTime=queryTimer()
         self%stopped = .TRUE.
      END SUBROUTINE stopTimer
!
!//////////////////////////////////////////////////////////////////////// 
! 
       FUNCTION elapsedTime(self,units)  
         IMPLICIT NONE

         CLASS(WCTimer)    :: self
         INTEGER, OPTIONAL :: units
         REAL(KIND=TP)  :: elapsedTime
         INTEGER :: ticks
!        ------------------------------------------
!        Return zero if the timer was never started
!        ------------------------------------------
         IF ( .NOT.self%started )     THEN
            elapsedTime = 0.0_TP
            RETURN
         END IF 

! query time if timer is not yet stopped
         IF ( .NOT.self%stopped )     THEN
            ticks =  queryTimer() - self%startTime
         else
            ticks =  self%finishTime - self%startTime
         END IF
         if (ticks<0) ticks = ticks+self%clock_max
         elapsedTime = ticks/real(self%clock_rate,kind=TP)

!        -------------------------------------
!        Convert to requested units if present
!        -------------------------------------
         IF ( PRESENT(units) )     THEN
         
            SELECT CASE ( units )
               CASE( TC_MINUTES ) 
                  elapsedTime = elapsedTime/60.0_TP
               CASE( TC_HOURS )
                  elapsedTime = elapsedTime/3600.0_TP
               CASE DEFAULT 
               
            END SELECT 
         END IF 
      
      END FUNCTION elapsedTime
      
      END MODULE WCTimerClass
