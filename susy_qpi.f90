! Susy QPI
! This program outputs data used to chart
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

  open(10, file="susy_qpi.log", position="append", status="unknown")
  write(10,*) "******************START: susy_qpi********************"
  close(10)


! Set constants
  Pi=4.0_dp*datan(1.0_dp)

! Set experimental data
  t   = 1_dp
  mu  =-0.9_dp*t
  x0  = 1e-2_dp
  epsf=-0.9_dp*x0
  V   = 1e-2_dp
  Uc  = 1_dp
  Uf  = 1e-3_dp

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

  do i=1,100
    do j=1,10
      call write_data(i*0.5_dp,j*0.2_dp)
    end do
  end do

  open(10, file="susy_qpi.log", position="append", status="old")
  write(10,*) "********************END: susy_qpi********************"
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
  call cpu_time(start)
  call date_and_time(date,time)

! frequency
  omega = dcmplx(om,del)

! Save results to file
  open(log, file="susy_qpi.log", position="append", status="old")
  write(somega, '(f0.2,"+",f0.2,"i")') om, del
!   filename = date//"_"//time//"_w="//trim(somega)//"_susy_qpi.dat"
  filename = trim(somega)//"_susy_qpi.dat"
  write(log,*)
  write(log,*) "======================================================="
  write(log,*)
  write(log,*) "START: "//filename
  write(log, *) "t:     ", t
  write(log, *) "mu:    ", mu
  write(log, *) "x0:    ", x0
  write(log, *) "epsf:  ", epsf
  write(log, *) "V:     ", V
  write(log, *) "Uc:    ", Uc
  write(log, *) "Uf:    ", Uf
  write(log, *) "omega: ", omega
  close(log)
  open(dat, file=filename, status="unknown")

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
        write(log,*) "ERROR: DCUHRE exit code ", ifail
        close(log)
      end if
      write(dat,*)  qx,  qy, result(1)
      write(dat,*)  qx, -qy, result(1)
      write(dat,*) -qx,  qy, result(1)
      write(dat,*) -qx, -qy, result(1)
      write(dat,*)  qy,  qx, result(1)
      write(dat,*)  qy, -qx, result(1)
      write(dat,*) -qy,  qx, result(1)
      write(dat,*) -qy, -qx, result(1)      
      call flush(dat)
    end do
  end do

  close(dat)
  open(log, file="susy_qpi.log", position="append", status="old")
  write(log,*) "END: "//filename
  call cpu_time(end)
  write(log,*) "TOTAL TIME: ", end-start
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
    write(log,*) "ERROR: Singular matrix "
    write(log,*) Gkinv
    close(log)
  end if
  call invert4x4(Gkqinv, Gkq, status)
  if ( status>0 ) then
    open(log, file="susy_qpi.log", position="append", status="old")
    write(log,*) "ERROR: Singular matrix "
    write(log,*) Gkqinv
    close(log)
  end if

  f(1)=dimag(Gk(1)*Gkq(1))

end subroutine sG0
