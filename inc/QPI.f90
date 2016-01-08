! Double precision
integer, parameter:: dp=kind(0.d0)

! Parameters
integer, intent(in)   :: ndim, nfun
real(dp), intent(in)  :: z(ndim)
real(dp), intent(out) :: f(nfun)

! Globar variables
common /freq/  omega, qx, qy
real(dp)    :: qx, qy
complex(dp) :: omega

! Local variables
real(dp)    :: kx, ky, kqx, kqy
complex(dp) :: Gk(16), Gkq(16)

interface
  subroutine Gmatrix(kx, ky, Gk)
    implicit none
    real(kind=8), intent(in) :: kx, ky
    complex(kind=8), intent(out) :: Gk(16)
  end subroutine Gmatrix
end interface