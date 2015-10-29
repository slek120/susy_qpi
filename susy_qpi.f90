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

  do i=0,4
    do j=0,4
      do k=0,4

! Set experimental data
        t   = 1_dp
        mu  = 0.1_dp*t
        x0  = 0.02_dp
        epsf= 0.1_dp*x0
        Uc  = 0.25_dp*i
        Uf  = 0.25_dp*j
        V   = 0.25_dp*k

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

        call write_data(0.0_dp,0.3_dp)
      end do
    end do
  end do

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
  call cpu_time(start)
  call date_and_time(date,time)
  filename = date//time
  write(title, '("w=", f0.2, "+", f0.2, "i", &
    &" V=", f0.2, " U_c=", f0.2, " U_f=", f0.2)') om, del, V, Uc, Uf

! frequency
  omega = dcmplx(om,del)

! Set column properties for sqlite
  call sqlite3_open( 'data.db', db )
  allocate( column(15) )
  call sqlite3_column_props( column(1), "t", SQLITE_REAL )
  call sqlite3_column_props( column(2), "mu", SQLITE_REAL )
  call sqlite3_column_props( column(3), "x0", SQLITE_REAL )
  call sqlite3_column_props( column(4), "epsf", SQLITE_REAL )
  call sqlite3_column_props( column(5), "V", SQLITE_REAL )
  call sqlite3_column_props( column(6), "Uc", SQLITE_REAL )
  call sqlite3_column_props( column(7), "Uf", SQLITE_REAL )
  call sqlite3_column_props( column(8), "omega", SQLITE_REAL )
  call sqlite3_column_props( column(9), "delta", SQLITE_REAL )
  call sqlite3_column_props( column(10), "qx", SQLITE_REAL )
  call sqlite3_column_props( column(11), "qy", SQLITE_REAL )
  call sqlite3_column_props( column(12), "result", SQLITE_REAL )
  call sqlite3_column_props( column(13), "absresult", SQLITE_REAL )
  call sqlite3_column_props( column(14), "date", SQLITE_CHAR, 8 )
  call sqlite3_column_props( column(15), "time", SQLITE_CHAR, 4 )
  call sqlite3_create_table( db, "susy_qpi", column)
  call sqlite3_create_table( db, "runs", column)

  call sqlite3_set_column( column(1), t )
  call sqlite3_set_column( column(2), mu )
  call sqlite3_set_column( column(3), x0 )
  call sqlite3_set_column( column(4), epsf )
  call sqlite3_set_column( column(5), V )
  call sqlite3_set_column( column(6), Uc )
  call sqlite3_set_column( column(7), Uf )
  call sqlite3_set_column( column(8), om )
  call sqlite3_set_column( column(9), del )
  call sqlite3_set_column( column(14), date )
  call sqlite3_set_column( column(15), time )

! Start log
  open(log, file="susy_qpi.log", position="append", status="old")
    write(log,*)
    write(log,*) "======================================"
    write(log,*) "t:     ", t
    write(log,*) "mu:    ", mu
    write(log,*) "x0:    ", x0
    write(log,*) "epsf:  ", epsf
    write(log,*) "V:     ", V
    write(log,*) "Uc:    ", Uc
    write(log,*) "Uf:    ", Uf
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

  do iqx=0,steps
    qx=iqx*qstep
    do iqy = 0,iqx
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

      call sqlite3_begin( db )

      call sqlite3_set_column( column(10), qx )
      call sqlite3_set_column( column(11), qy )
      call sqlite3_set_column( column(12), result(1) )
      call sqlite3_set_column( column(13), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi', column )

      call sqlite3_set_column( column(10), -qx )
      call sqlite3_set_column( column(11), qy )
      call sqlite3_set_column( column(12), result(1) )
      call sqlite3_set_column( column(13), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi', column )

      call sqlite3_set_column( column(10), qx )
      call sqlite3_set_column( column(11), -qy )
      call sqlite3_set_column( column(12), result(1) )
      call sqlite3_set_column( column(13), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi', column )

      call sqlite3_set_column( column(10), -qx )
      call sqlite3_set_column( column(11), -qy )
      call sqlite3_set_column( column(12), result(1) )
      call sqlite3_set_column( column(13), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi', column )

      call sqlite3_set_column( column(10), qy )
      call sqlite3_set_column( column(11), qx )
      call sqlite3_set_column( column(12), result(1) )
      call sqlite3_set_column( column(13), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi', column )

      call sqlite3_set_column( column(10), qy )
      call sqlite3_set_column( column(11), -qx )
      call sqlite3_set_column( column(12), result(1) )
      call sqlite3_set_column( column(13), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi', column )

      call sqlite3_set_column( column(10), -qy )
      call sqlite3_set_column( column(11), qx )
      call sqlite3_set_column( column(12), result(1) )
      call sqlite3_set_column( column(13), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi', column )

      call sqlite3_set_column( column(10), -qy )
      call sqlite3_set_column( column(11), -qx )
      call sqlite3_set_column( column(12), result(1) )
      call sqlite3_set_column( column(13), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi', column )

      call sqlite3_commit( db )
    end do
  end do

  call sqlite3_begin( db )
  call sqlite3_insert( db, 'runs', column)
  call sqlite3_commit( db )
  call sqlite3_close( db )

  open(log, file="susy_qpi.log", position="append", status="old")
    call cpu_time(end)
    write(log,*) "CALCULATION TIME: ", end-start
  close(log)

! Export plot to png
  call system('sqlite3 -column data.db "select distinct qx, qy, absresult from susy_qpi where date='&
    //date//' and time='//time//' order by qx, qy;" | awk -f add_blanks.awk > data/'//filename//".dat")
  call system('gnuplot -e ''filename="'//filename//'"; name="'//trim(title)//'"'' plot.gp')
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
