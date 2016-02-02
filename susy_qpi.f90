! Susy QPI
! This program outputs data used to plot
! the QPI spectrum of the Kondo Lattice Model
! using a supersymmetric representation
! 
! by Eric Mascot
! 
!==============================================================
! susy_qpi    Set experimental data and frequency
!==============================================================

program susy_qpi
  implicit none
  include "inc/susy_qpi.f90"

  call date_and_time(date,time)
  open(10, file="susy_qpi.log", position="append", status="unknown")
    write(10,*)
    write(10,*) "******************START: susy_qpi********************"
    write(10,*) "DATE: "//date//" TIME: "//time
  close(10)

!   do i=1,4
!     do j=1,4
!       do k=0,4
!         call write_data( -15._dp, 0.1_dp)
!         call write_data( -10._dp, 0.1_dp)
!         call write_data( -5._dp, 0.1_dp)
        call write_data(  0._dp, 0.1_dp)
!       end do
!     end do
!   end do

  open(10, file="susy_qpi.log", position="append", status="old")
    write(10,*) "********************END: susy_qpi********************"
    write(10,*)
  close(10)
end program susy_qpi


!==============================================================
! write_data  Writes data to OMEGA_susy_qpi.dat
!   om        real part of frequency
!   del       imaginary part of frequency
!==============================================================

subroutine write_data(om, del)
  implicit none
  include "inc/write_data.f90"

! Start timer
  call system_clock(start,rate)
  call date_and_time(date,time)
  
! frequency
  omega = dcmplx(om,del)
! Set constants
  Pi=4.0_dp*datan(1.0_dp)

! Start log
  filename = "data/"//date//time//".dat"
  open(dat, file=filename, position="append", status="unknown")
    write(dat, '("w=", f0.2, "+", f0.2, "i")') om, del
  close(dat)

  open(log, file="susy_qpi.log", position="append", status="old")
    write(log,*)
    write(log,*) "================="
    write(log,*) "omega: ", omega
    write(log,*) "date:  ", date
    write(log,*) "time:  ", time
  close(log)

! Settings for dcuhre
  a(1) = -Pi
  b(1) = Pi
  a(2) = -Pi
  b(2) = Pi
  minpts = 2000
  maxpts = 60000000
  abserr = 1e-6_dp
  relerr = 1e-3_dp

! Step through 1/8 triangle of BZ
  steps=100
  qstep=Pi/(steps*1.0_dp)
  i=0

  do iqx=0,steps
    do iqy = 0,iqx
      qx=iqx*qstep
      qy=iqy*qstep

!     Integrate QPI subroutine
      call dcuhre(ndim, nfun, a, b, minpts, maxpts, QPI, &
                  abserr, relerr, 0, nwork, 0, result,&
                  absest, neval, ifail, work)
!     Log results
      if (ifail>0) then
        open(log, file="susy_qpi.log", position="append", status="old")
        write(log,*) "ERROR: DCUHRE exit code ", ifail, qx, qy
        close(log)
      end if

      open(dat, file=filename, position="append", status="old")
        write(dat,*) qx, qy, result(1), result(2), result(3), result(4)
      close(dat)

      i=i+1
      call progress(i/5151.0, start)
    end do
  end do

  open(log, file="susy_qpi.log", position="append", status="old")
    call system_clock(end)
    write(log,*) "CALCULATION TIME: ", real(end-start)/real(rate)
  close(log)

end subroutine write_data


!==============================================================
! QPI         Integrand for dcuhre
!   ndim      number of variables
!   z         variables for integration
!   nfun      number of components of integral
!   f         value of integrand
!==============================================================

subroutine QPI(ndim, z, nfun, f)
  implicit none
  include "inc/QPI.f90"

  kx = z(1)
  ky = z(2)
  kqx=kx+qx
  kqy=ky+qy

  call Gmatrix(kx, ky, Gk)
  call Gmatrix(kqx,kqy,Gkq)

! c-band
! G(k,k)G(k+q,k+q)
  f(1)=dimag(Gk(1)*Gkq(1))
! G(k,k+Q)G(k+Q+q,k+q)
  f(2)=dimag(Gk(3)*Gkq(9))

! f-band
! G(k,k)G(k+q,k+q)
  f(3)=dimag(Gk(6)*Gkq(6))
! G(k,k+Q)G(k+Q+q,k+q)
  f(4)=dimag(Gk(4)*Gkq(13))

end subroutine QPI


!==============================================================
! Gmatrix     Returns Greens function matrix
!   kx        x component of k vector
!   ky        x component of k vector
!   Gk        Resultant matrix
!==============================================================

subroutine Gmatrix(kx,ky,Gk)
  implicit none
  include "inc/Gmatrix.f90"

! Set constants
  Pi=4.0_dp*datan(1.0_dp)
  kQx = kx+Pi
  kQy = ky+Pi

! Energy dispersion of c band
  epsk  = t10*( -2._dp * (dcos(kx) + dcos(ky)) - mu)
  epskQ = t10*( -2._dp * (dcos(kQx) + dcos(kQy)) - mu)

! Energy dispersion of f band
  chik  = x10*( -2._dp * (dcos(kx) + dcos(ky))  &
          -4._dp * x11 * dcos(kx) * dcos(ky)  &
          -2._dp * x30 * (dcos(3*kx) + dcos(3*ky))  &
          -2._dp * x31 * (dcos(3*kx) * dcos(ky) + dcos(kx) * dcos(3*ky))  &
          -epsf)
  chikQ = x10*( -2._dp * (dcos(kQx) + dcos(kQy))  &
          -4._dp * x11 * dcos(kQx) * dcos(kQy)  &
          -2._dp * x30 * (dcos(3*kQx) + dcos(3*kQy))  &
          -2._dp * x31 * (dcos(3*kQx) * dcos(kQy) + dcos(kQx) * dcos(3*kQy))  &
          -epsf)

  a = omega - epsk
  b = omega - chik
  c = omega - epskQ
  d = omega - chikQ

  Gkinv = (/ complex(dp) :: &
             a,  s, -Uc,  0, &
             s,  b,   0, Uf, &
           -Uc,  0,   c,  s, &
             0, Uf,   s,  d /)

  call invert4x4(Gkinv, Gk, status)
  if ( status>0 ) then
    open(log, file="susy_qpi.log", position="append", status="old")
    write(log,*) "ERROR: Singular matrix "
    write(log,*) Gkinv
    close(log)
  end if
end subroutine Gmatrix


!==============================================================
! progress    Show progress bar
!   percent   Percent completed
!   start     Time of start from system_clock
!==============================================================

subroutine progress(percent, start)
  implicit none
  integer, intent(in) :: start
  real, intent(in)    :: percent
  integer             :: ticks, end, rate, elapsed, remaining
  character(len=52)   :: bar

  ticks = int(percent*50)
  if (ticks>50) ticks=50
  if (ticks<0) ticks=0

  call system_clock(end, rate)
  elapsed   = int(real(end-start)/real(rate))
  remaining = int(elapsed*(1.0/percent-1.0))
  bar  = "["//repeat("=",ticks)//repeat(" ",50-ticks)//"]"
  
  write(*,"(A,I3,'% ',I3,':',I2.2,' elapsed',I3,':',I2.2,' remaining')") &
    bar, int(percent*100), elapsed/3600, mod(elapsed/60,60), elapsed/3600, mod(remaining/60,60)
  call flush()
end subroutine progress