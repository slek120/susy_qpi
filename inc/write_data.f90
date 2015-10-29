external sG0
external dcuhre

! Double precision
! integer, parameter :: dp=kind(0.d0)

! sqlite interface
type(SQLITE_DATABASE)                      :: db
type(SQLITE_COLUMN), dimension(:), pointer :: column

! Parameters
!    omega = om + i * del
real(dp), intent(in) :: om, del

! Global variables
!     omega  complex.
!            Frequency (energy) used for G
common /gvars/ Pi, t, mu, x0, epsf, V, Uc, Uf, Gkinv, Gkqinv
common /freq/  omega, qx, qy
common /date/  date, time
real(dp)    :: Pi, t, mu, x0, epsf, V, Uc, Uf, qx, qy
complex(dp) :: omega, Gkinv(16), Gkqinv(16)
character(len=8) :: date
character(len=4) :: time
character(len=12):: filename
character(len=64):: title

! Local variables
!     steps  integer.
!            Number of iterations to do
!            steps*steps/2 total iterations
real(dp)         :: qstep
integer          :: steps, iqx, iqy, log=10
real             :: start, end

! Variables for dcuhre
!     ndim   integer.
!            number of variables. 1 < ndim <=  15.
!     nfun   integer.
!            number of components of the integral.
!     a      real array of dimension ndim.
!            lower limits of integration.
!     b      real array of dimension ndim.
!            upper limits of integration.
!     minpts integer.
!            minimum number of function evaluations.
!     maxpts integer.
!            maximum number of function evaluations.
!     abserr real.
!            requested absolute error.
!     relerr real.
!            requested relative error.
!     nwork  integer.
!            defines the length of the working array work.
!     result real array of dimension nfun.
!            approximations to all components of the integral.
!     absest real array of dimension nfun.
!            estimates of absolute errors.
!     neval  integer.
!            number of function evaluations used by dcuhre.
!     ifail  integer.
!            ifail = 0 for normal exit, when abserr(k) <=  epsabs or
!              abserr(k) <=  abs(result(k))*epsrel with maxpts or less
!              function evaluations for all values of k,
!              1 <= k <= nfun .
!            ifail = 1 if maxpts was too small for dcuhre
!              to obtain the required accuracy. in this case dcuhre
!              returns values of result with estimated absolute
!              errors abserr.
!            ifail = 2 if key is less than 0 or key greater than 4.
!            ifail = 3 if ndim is less than 2 or ndim greater than 15.
!            ifail = 4 if key = 1 and ndim not equal to 2.
!            ifail = 5 if key = 2 and ndim not equal to 3.
!            ifail = 6 if nfun is less than 1.
!            ifail = 7 if volume of region of integration is zero.
!            ifail = 8 if maxpts is less than 3*num.
!            ifail = 9 if maxpts is less than minpts.
!            ifail = 10 if epsabs < 0 and epsrel < 0.
!            ifail = 11 if nw is too small.
!            ifail = 12 if unlegal restar.
!     work   real array of dimension nwork.
!            used as working storage.
integer            :: minpts, maxpts, neval, ifail
integer, parameter :: ndim = 2, nfun = 1, nwork = 6000000
real(dp)           :: a(ndim), b(ndim), abserr, relerr, &
                      result(nfun), absest(nfun), work(nwork)

! Interface for called procedures
! TODO: fix segmentation fault caused by non-declared array sizes
!   interface
!     subroutine dcuhre(ndim, nfun, a, b, minpts, maxpts, sG0, &
!                     abserr, relerr, key, nwork, restar, result,&
!                     absest, neval, ifail, work)
!       implicit none
!       integer, parameter      :: dp=kind(0.d0)
!       integer, intent(in)     :: ndim, nfun, minpts, maxpts, key, nwork, restar
!       integer, intent(out)    :: neval, ifail
!       real (dp), intent(in)   :: a(:), b(:), abserr, relerr
!       real (dp), intent(out)  :: result(:), absest(:)
!       real (dp), intent(in out)  :: work(:)

!       interface
!         subroutine sG0(ndim, z, nfun, f)
!           implicit none
!           integer, parameter     :: dp=kind(0.d0)
!           integer, intent(in)    :: ndim, nfun
!           real (dp), intent(in)  :: z(:)
!           real (dp), intent(out) :: f(:)
!         end subroutine sG0
!       end interface
!     end subroutine dcuhre
!   end interface