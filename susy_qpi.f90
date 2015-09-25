program susy_qpi
  implicit none
  integer, parameter :: dp=kind(0.d0)
  integer :: i,j

  do i=1,100
    do j=1,10
      call write_data(i*0.5_dp,j*0.2_dp)
    end do
  end do

end program susy_qpi

subroutine write_data(om, del)
  implicit none
  external sG0

! Double precision
  integer, parameter :: dp=kind(0.d0)

! Parameters
!    omega = om + i * del
  real(dp), intent(in) :: om, del

! Global variables
!   epsq = -2 * t ( cos(kx) + cos(ky) ) - mu
!     t      real.
!     mu     real.
!   chiq = -2 * x0 ( cos(kx) + cos(ky) ) - epsf
!     x0     real.
!     epsf   real.
!   Inv[G(k,omega)] = (omega)*I - H_k
!     Gkinv  real array of dimension 16 (4x4 matrix)
!            Inverse of greens function
!     Gkqinv real array of dimension 16 (4x4 matrix)
!            Inverse of greens function at k+q
!     omega  complex.
!            Frequency (energy) used for G

  common /gvars/      t, mu, x0, epsf, qx, qy, Pi, Gkinv, Gkqinv
  real(dp)         :: t, mu, x0, epsf, qx, qy, Pi, Uc, Uf, V
  complex(dp)      :: Gkinv(16), Gkqinv(16)
  common /freq/       omega
  complex(dp)      :: omega

! Local variables
!     steps  integer.
!            Number of iterations to do
!            steps*steps/2 total iterations
  real(dp)         :: qstep
  integer          :: steps, iqx, iqy, log=10, dat=20
  real             :: start, end
  character(len=8) :: date
  character(len=4) :: time
  character(len=16):: somega
  character(len=64):: filename
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

! Start timer
  call cpu_time(start)
  call date_and_time(date,time)

! Save results to file
  open(log, file="susy_qpi.log", position="append", status="unknown")
  write(somega, '(f6.2,"+",f6.2,"i")') om, del
  filename = date//"_"//time//"_w="//trim(somega)//"_susy_qpi.dat"
  write(log,*)
  write(log,*) "======================================================="
  write(log,*)
  write(log,*) "START: "//filename
  close(log)
  open(dat, file=filename, status="unknown")

! Set constants
  Pi=4.0_dp*datan(1.0_dp)

! Settings for dcuhre
  a(1) = -1.0_dp*Pi
  b(1) = Pi
  a(2) = -1.0_dp*Pi
  b(2) = Pi
  minpts = 2000
  maxpts = 60000000
  abserr = 1e-6_dp
  relerr = 1e-3_dp

! Set experimental data
  t   = 1_dp
  mu  =-0.9_dp*t
  x0  = 1e-2_dp
  epsf=-0.9_dp*x0
  V   = 1e-2_dp
  Uc  = 1_dp
  Uf  = 1e-3_dp

! frequency
  omega = dcmplx(om,del)

! Set up matrix as an array of size 16
  Gkinv = (/ complex(dp) :: &
       0,  V, -Uc,  0, &
       V,  0,   0, Uf, &
     -Uc,  0,   0,  V, &
       0, Uf,   V,  0 /)
  Gkqinv = (/ complex(dp) :: &
       0,  V, -Uc,  0, &
       V,  0,   0, Uf, &
     -Uc,  0,   0,  V, &
       0, Uf,   V,  0 /)

  open(log, file="susy_qpi.log", position="append", status="old")
  write(log, *) "t:     ", t
  write(log, *) "mu:    ", mu
  write(log, *) "x0:    ", x0
  write(log, *) "epsf:  ", epsf
  write(log, *) "V:     ", V
  write(log, *) "Uc:    ", Uc
  write(log, *) "Uf:    ", Uf
  write(log, *) "omega: ", omega
  close(log)

! Step through 1/8 triangle of BZ
  steps=100
  qstep=Pi/(steps*1.0_dp)

  do iqx=0,steps
    qx=iqx*qstep
    do iqy = 0,iqx
      qy=iqy*qstep
!     Integrate sG0 subroutine and write results to file
      call dcuhre(ndim, nfun, a, b, minpts, maxpts, sG0, &
                  abserr, relerr, 0, nwork, 0, result,&
                  absest, neval, ifail, work)
!     Log results
      if (ifail>0) then
        open(log, file="susy_qpi.log", position="append", status="old")
        write (log,*) "ERROR: DCUHRE exit code ", ifail
        close(log)
      end if
      write(dat,*)  qx,  qy, result(1)
      write(dat,*)  qy,  qx, result(1)
      write(dat,*) -qy,  qx, result(1)
      write(dat,*) -qx,  qy, result(1)
      write(dat,*) -qx, -qy, result(1)
      write(dat,*) -qy, -qx, result(1)      
      write(dat,*)  qy, -qx, result(1)
      write(dat,*)  qx, -qy, result(1)
      call flush(dat)
    end do
  end do

  close(dat)
  open(log, file="susy_qpi.log", position="append", status="old")
  write(log,*) "END: "//date//"_"//time//"_susy_qpi.dat"
  call cpu_time(end)
  write(log,*) "TOTAL TIME: ", end-start
  close(log)

end subroutine write_data

!==============================================================

subroutine sG0(ndim, z, nfun, f)
  implicit none

! Double precision
  integer, parameter:: dp=kind(0.d0)

! Parameters
  integer, intent(in)   :: ndim, nfun
  real(dp), intent(in)  :: z(ndim)
  real(dp), intent(out) :: f(nfun)

! Globar variables
  common /gvars/    t, mu, x0, epsf, qx, qy, Pi, Gkinv, Gkqinv
  real(dp)       :: t, mu, x0, epsf, qx, qy, Pi
  complex(dp)    :: Gk(16), Gkq(16), Gkinv(16), Gkqinv(16)
  common /freq/     omega
  complex(dp)    :: omega

! Local variables
  integer        :: status, log=10
  real(dp)       :: kx, ky, kqx, kqy
  real(dp)       :: epsk, epskpQ, epskq, epskqpQ, chik, chikq, chikpQ, chikqpQ

  kx = z(1)
  ky = z(2)
  kqx=kx+qx
  kqy=ky+qy

  epsk    = -2.0_dp * t * (dcos(kx)     + dcos(ky))     - mu
  epskpQ  = -2.0_dp * t * (dcos(kx+Pi)  + dcos(ky+Pi))  - mu
  epskq   = -2.0_dp * t * (dcos(kqx)    + dcos(kqy))    - mu
  epskqpQ = -2.0_dp * t * (dcos(kqx+Pi) + dcos(kqy+Pi)) - mu

  chik    = -2.0_dp * x0 * (dcos(kx)     + dcos(ky))     - epsf
  chikpQ  = -2.0_dp * x0 * (dcos(kx+Pi)  + dcos(ky+Pi))  - epsf
  chikq   = -2.0_dp * x0 * (dcos(kqx)    + dcos(kqy))    - epsf
  chikqpQ = -2.0_dp * x0 * (dcos(kqx+Pi) + dcos(kqy+Pi)) - epsf

  Gkinv(1)   = omega - epsk
  Gkinv(6)   = omega - chik
  Gkinv(11)  = omega - epskpQ
  Gkinv(16)  = omega - chikpQ
  Gkqinv(1)  = omega - epskq
  Gkqinv(6)  = omega - chikq
  Gkqinv(11) = omega - epskqpQ
  Gkqinv(16) = omega - chikqpQ

  call invert4x4(Gkinv, Gk, status)
  if ( status>0 ) then
    open(log, file="susy_qpi.log", position="append", status="old")
    write (log,*) "ERROR: Singular matrix "
    write (log,*) Gkinv
    close(log)
  end if
  call invert4x4(Gkqinv, Gkq, status)
  if ( status>0 ) then
    open(log, file="susy_qpi.log", position="append", status="old")
    write (log,*) "ERROR: Singular matrix "
    write (log,*) Gkqinv
    close(log)
  end if

  f(1)=dimag(Gk(1)*Gkq(1))

end subroutine sG0
