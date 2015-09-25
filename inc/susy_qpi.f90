! Double precision
integer, parameter :: dp=kind(0.d0)

! Global variables
!   epsq = -2 * t ( cos(kx) + cos(ky) ) - mu
!     t      real.
!     mu     real.
!   chiq = -2 * x0 ( cos(kx) + cos(ky) ) - epsf
!     x0     real.
!     epsf   real.
!   Inv[G(k,omega)] = (omega)*I - H_k
!     Gkinv  real array of dimension 16 (4x4 matrix)
!            Inverse of Green function at k
!     Gkqinv real array of dimension 16 (4x4 matrix)
!            Inverse of Green function at k+q
common /gvars/ Pi, t, mu, x0, epsf, V, Uc, Uf, Gkinv, Gkqinv
real(dp)    :: Pi, t, mu, x0, epsf, V, Uc, Uf
complex(dp) :: Gkinv(16), Gkqinv(16)

! Local variables
integer  :: i, j

! Interface for called procedures
interface
  subroutine write_data(om, del)
    implicit none
    integer, parameter   :: dp=kind(0.d0)
    real(dp), intent(in) :: om, del
  end subroutine write_data
end interface