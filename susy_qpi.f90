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
        call write_data( -8._dp, 0.1_dp)
        call write_data( -6._dp, 0.1_dp)
        call write_data(  0._dp, 0.1_dp)
        call write_data(  6._dp, 0.1_dp)
        call write_data( 14._dp, 0.1_dp)
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
        write(dat,*) qx, qy, result(1), dabs(result(1))
      close(dat)

!     Calculate percentage complete and estimated time remaining
      i=i+1
      call system_clock(end)
      est = real(end-start)/real(rate)*(5150.0/i-1.0)
      write(*,"(I3,'% ',I4,':',I2,' remaining')") int(i/51.5), est/60, mod(est,60)
!     Flush the stdout (for nohup)
      call flush()
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

  f(1)=dimag(Gk(1)*Gkq(1))

end subroutine QPI

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