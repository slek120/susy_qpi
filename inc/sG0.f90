! Double precision
integer, parameter:: dp=kind(0.d0)

! Parameters
integer, intent(in)   :: ndim, nfun
real(dp), intent(in)  :: z(ndim)
real(dp), intent(out) :: f(nfun)

! Globar variables
common /freq/  omega, qx, qy
real(dp)    :: qx, qy, t1,t2,t3,t4,mu
complex(dp) :: omega

! Local variables
integer     :: status, log=10
real(dp)    :: kx, ky, kqx, kqy, &
               epsk, epskq
complex(dp) :: Gk, Gkq

! Interface for called procedures
interface
  subroutine invert4x4(m, invOut, fail)
    implicit none
    complex(kind=8), dimension(16), intent(in) :: m
    complex(kind=8), dimension(16), intent(out) :: invOut
    integer, intent(out) :: fail
  end subroutine invert4x4
end interface