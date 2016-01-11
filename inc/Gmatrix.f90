! Double precision
integer, parameter :: dp=kind(0.d0)

! Parameters
real(dp), intent(in)     :: kx, ky
complex(dp), intent(out) :: Gk(16)

! Globar variables
common /freq/ omega, qx, qy
real(dp)    :: qx, qy
complex(dp) :: omega

! Local variables
integer     :: status, log=10
complex(dp) :: a, b, c, d, Gkinv(16)
real(dp)    :: Pi, kQx, kQy
real(dp)    :: epsk, chik, epskQ, chikQ
real(dp), parameter :: t10 = -185._dp
real(dp), parameter :: mu  = 3.8_dp
real(dp), parameter :: x10 = -5._dp
real(dp), parameter :: x11 =  0.35_dp
real(dp), parameter :: x30 =  0.035_dp
real(dp), parameter :: x31 = -0.3_dp
real(dp), parameter :: epsf= 0.8_dp
real(dp), parameter :: s   = 9._dp
real(dp), parameter :: Uc  = 0._dp
real(dp), parameter :: Uf  = 0._dp

! Interface for called procedures
interface
  subroutine invert4x4(m, invOut, fail)
    implicit none
    complex(kind=8), intent(in)  :: m(16)
    complex(kind=8), intent(out) :: invOut(16)
    integer, intent(out) :: fail
  end subroutine invert4x4
end interface