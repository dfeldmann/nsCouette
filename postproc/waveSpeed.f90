!=======================================================
!
!  A program to calculate the wave speed of wavy vortices
!by using the Brent's method (minimal of a scalar function)
!
!                                        Liang Shi
!                                        28 Sep 2012
!=======================================================

PROGRAM waveSpeed
  
  USE mod_preAnalysis
  IMPLICIT NONE
  REAL(KIND=8)    :: t_start, t_end  ! Start & end time of MPI
  REAL(KIND=8)    :: f_phi,dPhi, x_min, x_max
  INTEGER (KIND=4) :: status
  INTEGER :: i
  integer, parameter :: maxiter = 500

  CALL Mpi_init(ierr)
  !-------------------------------Initialization phase
  t_start = Mpi_wtime()
  CALL pre_mpi()
  CALL read_pars()
  call init_grid()                   
  CALL init_mpidist()                
  CALL read_coeff()

!  CALL MPI_Gather(f_hat_mp(1,1),counter,mpi_spec,f_hat(1,1,1),counter,&
!       mpi_spec,root,comm,ierr)

  IF (myid == root) THEN
     status = 0
     iter=0
     READ*, x_min     !=0D0
     READ*, x_max     !=5D0
     CALL local_min_rc(x_min,x_max,dPhi,status,f_phi)
     DO WHILE (iter< maxiter .and. status .NE. 0)
        CALL norm_du(dPhi, f_phi) ! this function call constrains the number of MPI task to be 1
        CALL local_min_rc(x_min,x_max,dPhi,status,f_phi)
        iter=iter+1
     END DO
     PRINT*, '# number of iterations: ', iter,' max: ',maxiter
     PRINT*, '# min(f(d_phi)) = ', f_phi,' at',dPhi
     PRINT*, '# non-dimensional wavespeed d_phi*r_i/dt/Re_i = ', dPhi*0.868/0.132/delta_t/Re_i
  END IF

  DO i = 1,1001
     dPhi = (i-1)*0.01D0
     CALL norm_du(dPhi, f_phi)
     IF (myid == root) PRINT*, dPhi,f_phi
  END DO
  
  t_end = Mpi_wtime()
!  IF (myid == root) PRINT*,t_end-t_start

  CALL post_mpi()
  CALL MPI_Finalize(ierr)

END PROGRAM waveSpeed

!----------------------------------------------------------------------------------
subroutine local_min_rc ( a, b, arg, status, value )

!*****************************************************************************80
!
!! LOCAL_MIN_RC seeks a minimizer of a scalar function of a scalar variable.
!
!  Discussion:
!
!    This routine seeks an approximation to the point where a function
!    F attains a minimum on the interval (A,B).
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much
!    slower than that for a Fibonacci search.  If F has a continuous
!    second derivative which is positive at the minimum (which is not
!    at A or B), then convergence is superlinear, and usually of the
!    order of about 1.324...
!
!    The routine is a revised version of the Brent local minimization 
!    algorithm, using reverse communication.
!
!    It is worth stating explicitly that this routine will NOT be
!    able to detect a minimizer that occurs at either initial endpoint
!    A or B.  If this is a concern to the user, then the user must
!    either ensure that the initial interval is larger, or to check
!    the function value at the returned minimizer against the values
!    at either endpoint.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters
!
!    Input/output, real ( kind = 8 ) A, B.  On input, the left and right
!    endpoints of the initial interval.  On output, the lower and upper
!    bounds for an interval containing the minimizer.  It is required
!    that A < B.
!
!    Output, real ( kind = 8 ) ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS zero, ARG is the routine's
!    estimate for the function minimizer.
!
!    Input/output, integer ( kind = 4 ) STATUS, used to communicate between 
!    the user and the routine.  The user only sets STATUS to zero on the first 
!    call, to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated minimizer.
!
!    Input, real ( kind = 8 ) VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
!  Local parameters:
!
!    C is the squared inverse of the golden ratio.
!
!    EPS is the square root of the relative machine precision.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ), save :: c
  real ( kind = 8 ), save :: d
  real ( kind = 8 ), save :: e
  real ( kind = 8 ), save :: eps
  real ( kind = 8 ), save :: fu
  real ( kind = 8 ), save :: fv
  real ( kind = 8 ), save :: fw
  real ( kind = 8 ), save :: fx
  real ( kind = 8 ), save :: midpoint
  real ( kind = 8 ), save :: p
  real ( kind = 8 ), save :: q
  real ( kind = 8 ), save :: r
  integer ( kind = 4 ) status
  real ( kind = 8 ), save :: tol
  real ( kind = 8 ), save :: tol1
  real ( kind = 8 ), save :: tol2
  real ( kind = 8 ), save :: u
  real ( kind = 8 ), save :: v
  real ( kind = 8 ) value
  real ( kind = 8 ), save :: w
  real ( kind = 8 ), save :: x
!
!  STATUS (INPUT) = 0, startup.
!
  if ( status == 0 ) then

    if ( b <= a ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LOCAL_MIN_RC - Fatal error!'
      write ( *, '(a)' ) '  A < B is required, but'
      write ( *, '(a,g14.6)' ) '  A = ', a
      write ( *, '(a,g14.6)' ) '  B = ', b
      status = -1
      stop
    end if

    c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

    eps = sqrt ( epsilon ( eps ) )
    tol = epsilon ( tol )

    v = a + c * ( b - a )
    w = v
    x = v
    e = 0.0D+00

    status = 1
    arg = x

    return
!
!  STATUS (INPUT) = 1, return with initial function value of FX.
!
  else if ( status == 1 ) then

    fx = value
    fv = fx
    fw = fx
!
!  STATUS (INPUT) = 2 or more, update the data.
!
  else if ( 2 <= status ) then

    fu = value

    if ( fu <= fx ) then

      if ( x <= u ) then
        a = x
      else
        b = x
      end if

      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu

    else

      if ( u < x ) then
        a = u
      else
        b = u
      end if

      if ( fu <= fw .or. w == x ) then
        v = w
        fv = fw
        w = u
        fw = fu
      else if ( fu <= fv .or. v == x .or. v == w ) then
        v = u
        fv = fu
      end if

    end if

  end if
!
!  Take the next step.
!
  midpoint = 0.5D+00 * ( a + b )
  tol1 = eps * abs ( x ) + tol / 3.0D+00
  tol2 = 2.0D+00 * tol1
!
!  If the stopping criterion is satisfied, we can exit.
!
  if ( abs ( x - midpoint ) <= ( tol2 - 0.5D+00 * ( b - a ) ) ) then
    status = 0
    return
  end if
!
!  Is golden-section necessary?
!
  if ( abs ( e ) <= tol1 ) then
    if ( midpoint <= x ) then
      e = a - x
    else
      e = b - x
    end if

    d = c * e
!
!  Consider fitting a parabola.
!
  else

    r = ( x - w ) * ( fx - fv )
    q = ( x - v ) * ( fx - fw )
    p = ( x - v ) * q - ( x - w ) * r
    q = 2.0D+00 * ( q - r )
    if ( 0.0D+00 < q ) then
      p = - p
    end if
    q = abs ( q )
    r = e
    e = d
!
!  Choose a golden-section step if the parabola is not advised.
!
    if ( &
      ( abs ( 0.5D+00 * q * r ) <= abs ( p ) ) .or. &
      ( p <= q * ( a - x ) ) .or. &
      ( q * ( b - x ) <= p ) ) then

      if ( midpoint <= x ) then
        e = a - x
      else
        e = b - x
      end if

      d = c * e
!
!  Choose a parabolic interpolation step.
!
    else

      d = p / q
      u = x + d

      if ( ( u - a ) < tol2 ) then
        d = sign ( tol1, midpoint - x )
      end if

      if ( ( b - u ) < tol2 ) then
        d = sign ( tol1, midpoint - x )
      end if

    end if

  end if
!
!  F must not be evaluated too close to X.
!
  if ( tol1 <= abs ( d ) ) then
    u = x + d
  end if

  if ( abs ( d ) < tol1 ) then
    u = x + sign ( tol1, d )
  end if
!
!  Request value of F(U).
!
  arg = u
  status = status + 1

  return
end subroutine local_min_rc


