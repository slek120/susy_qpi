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


! Set constants
  Pi=4.0_dp*datan(1.0_dp)

!   do i=1,4
!     do j=1,4
!       do k=0,4
! Set experimental data
        t1 = -15._dp
        t2 =  0.35_dp *t1
        t3 =  0.035_dp*t1
        t4 = -0.15_dp *t1
        mu =  0.8_dp  *t1

        call write_data(-20._dp, 0.01_dp)
        call write_data(-15._dp, 0.01_dp)
        call write_data(-10._dp, 0.01_dp)
        call write_data( -8._dp, 0.01_dp)
        call write_data( -6._dp, 0.01_dp)
        call write_data( -4._dp, 0.01_dp)
        call write_data( -2._dp, 0.01_dp)
        call write_data(  0._dp, 0.01_dp)
        call write_data(  2._dp, 0.01_dp)
        call write_data(  4._dp, 0.01_dp)
        call write_data(  6._dp, 0.01_dp)
        call write_data(  8._dp, 0.01_dp)
        call write_data( 10._dp, 0.01_dp)
        call write_data( 12._dp, 0.01_dp)
        call write_data( 14._dp, 0.01_dp)
        call write_data( 16._dp, 0.01_dp)
        call write_data( 19._dp, 0.01_dp)
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
  use sqlite
  implicit none
  include "inc/write_data.f90"

! Start timer
  call system_clock(start,rate)
  call date_and_time(date,time)
  
! frequency
  omega = dcmplx(om,del)

! Start log
  filename = "data/"//date//time//".dat"
  open(dat, file=filename, position="append", status="unknown")
    write(dat, '("w=", f0.2, "+", f0.2, "i")') om, del
  close(dat)

  open(log, file="susy_qpi.log", position="append", status="old")
    write(log,*)
    write(log,*) "======================================"
    write(log,*) "t1:    ", t1
    write(log,*) "t2:    ", t2
    write(log,*) "t3:    ", t3
    write(log,*) "t4:    ", t4
    write(log,*) "mu:    ", mu
    write(log,*) "omega: ", omega
    write(log,*) "date:  ", date
    write(log,*) "time:  ", time
  close(log)

! Settings for dcuhre
  a(1) = -1.0_dp*Pi
  b(1) = Pi
  a(2) = -1.0_dp*Pi
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
    qx=iqx*qstep
    do iqy = 0,iqx
      i=i+1
      qy=iqy*qstep
!     Integrate sG0 subroutine
      call dcuhre(ndim, nfun, a, b, minpts, maxpts, sG0, &
                  abserr, relerr, 0, nwork, 0, result,&
                  absest, neval, ifail, work)
!     Log results
      if (ifail>0) then
        open(log, file="susy_qpi.log", position="append", status="old")
        write(log,*) "ERROR: DCUHRE exit code ", ifail
        close(log)
      end if

      open(dat, file=filename, position="append", status="old")
        write(dat,*) qx, qy, result(1), dabs(result(1))
      close(dat)

<<<<<<< HEAD
      call cpu_time(end)
      print *, 0.000194175*(0.5*(iqx*(iqx+1))+iqy), end-start
=======
!     Calculate percentage complete and estimated time remaining
      call system_clock(end)
      est = real(end-start)/real(rate)*(5150.0/i-1.0)
      write(*,"(I3,'% ',I4,':',I2,' remaining')") int(i/51.5), est/60, mod(est,60)
!     Flush the stdout (for nohup)
>>>>>>> single_band
      call flush()
    end do
  end do

  open(log, file="susy_qpi.log", position="append", status="old")
    call system_clock(end)
    write(log,*) "CALCULATION TIME: ", real(end-start)/real(rate)
  close(log)

end subroutine write_data

!==============================================================
! sG0         Integrand for dcuhre
!   ndim      number of variables
!   z         variables for integration
!   nfun      number of components of integral
!   f         value of integrand
!==============================================================

subroutine sG0(ndim, z, nfun, f)
  implicit none
  include "inc/sG0.f90"

  kx = z(1)
  ky = z(2)
  kqx=kx+qx
  kqy=ky+qy

  epsk = -2_dp * t1 * (dcos(kx) + dcos(ky))  &
         -4._dp * t2 * dcos(kx) * dcos(ky)  &
         -2_dp * t3 * (dcos(3*kx) + dcos(3*ky))  &
         -4._dp * t4 * (dcos(3*kx) * dcos(ky) + dcos(kx) * dcos(3*ky))  &
         -mu
  epskq= -2_dp * t1 * (dcos(kqx) + dcos(kqy))  &
         -4._dp * t2 * dcos(kqx) * dcos(kqy)  &
         -2_dp * t3 * (dcos(3*kqx) + dcos(3*kqy))  &
         -4._dp * t4 * (dcos(3*kqx) * dcos(kqy) + dcos(kqx) * dcos(3*kqy))  &
         -mu

  Gk  = 1_dp/(omega - epsk)
  Gkq = 1_dp/(omega - epskq)

  f(1)=dimag(Gk*Gkq)

end subroutine sG0
