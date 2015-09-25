! Double precision
integer, parameter:: dp=kind(0.d0)

! Parameters
integer, intent(in)   :: ndim, nfun
real(dp), intent(in)  :: z(ndim)
real(dp), intent(out) :: f(nfun)

! Globar variables
common /gvars/ Pi, t, mu, x0, epsf, V, Uc, Uf, Gkinv, Gkqinv
common /freq/  omega, qx, qy
real(dp)    :: Pi, t, mu, x0, epsf, V, Uc, Uf, qx, qy
complex(dp) :: omega, Gkinv(16), Gkqinv(16)

! Local variables
integer     :: status, log=10
real(dp)    :: kx, ky, kqx, kqy, &
               epsk, epskpQ, epskq, epskqpQ, &
               chik, chikq, chikpQ, chikqpQ
complex(dp) :: Gk(16), Gkq(16)

! Interface for called procedures
interface
  subroutine invert4x4(m, invOut, fail)
    implicit none
    complex(kind=8), dimension(16), intent(in) :: m
    complex(kind=8), dimension(16), intent(out) :: invOut
    integer, intent(out) :: fail
  end subroutine invert4x4
end interface