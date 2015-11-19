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

        call write_data(-25._dp, 0.01_dp)
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
  call cpu_time(start)
  call date_and_time(date,time)
  filename = date//time//".dat"
  write(title, '("w=", f0.2, "+", f0.2, "i")') om, del

! frequency
  omega = dcmplx(om,del)

! Set column properties for sqlite
  call sqlite3_open( 'data.db', db )
  allocate( column(13) )
  call sqlite3_column_props( column(1), "t1", SQLITE_REAL )
  call sqlite3_column_props( column(2), "t2", SQLITE_REAL )
  call sqlite3_column_props( column(3), "t3", SQLITE_REAL )
  call sqlite3_column_props( column(4), "t4", SQLITE_REAL )
  call sqlite3_column_props( column(5), "mu", SQLITE_REAL )
  call sqlite3_column_props( column(6), "omega", SQLITE_REAL )
  call sqlite3_column_props( column(7), "delta", SQLITE_REAL )
  call sqlite3_column_props( column(8), "qx", SQLITE_REAL )
  call sqlite3_column_props( column(9), "qy", SQLITE_REAL )
  call sqlite3_column_props( column(10), "result", SQLITE_REAL )
  call sqlite3_column_props( column(11), "absresult", SQLITE_REAL )
  call sqlite3_column_props( column(12), "date", SQLITE_CHAR, 8 )
  call sqlite3_column_props( column(13), "time", SQLITE_CHAR, 4 )
  call sqlite3_create_table( db, "susy_qpi_single_band", column)
  call sqlite3_create_table( db, "runs_single_band", column)

  call sqlite3_set_column( column(1), t1 )
  call sqlite3_set_column( column(2), t2 )
  call sqlite3_set_column( column(3), t3 )
  call sqlite3_set_column( column(4), t4 )
  call sqlite3_set_column( column(5), mu )
  call sqlite3_set_column( column(6), om )
  call sqlite3_set_column( column(7), del )
  call sqlite3_set_column( column(12), date )
  call sqlite3_set_column( column(13), time )

! Start log
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

      call sqlite3_set_column( column(8), qx )
      call sqlite3_set_column( column(9), qy )
      call sqlite3_set_column( column(10), result(1) )
      call sqlite3_set_column( column(11), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi_single_band', column )

      call sqlite3_set_column( column(8), -qx )
      call sqlite3_set_column( column(9), qy )
      call sqlite3_set_column( column(10), result(1) )
      call sqlite3_set_column( column(11), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi_single_band', column )

      call sqlite3_set_column( column(8), qx )
      call sqlite3_set_column( column(9), -qy )
      call sqlite3_set_column( column(10), result(1) )
      call sqlite3_set_column( column(11), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi_single_band', column )

      call sqlite3_set_column( column(8), -qx )
      call sqlite3_set_column( column(9), -qy )
      call sqlite3_set_column( column(10), result(1) )
      call sqlite3_set_column( column(11), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi_single_band', column )

      call sqlite3_set_column( column(8), qy )
      call sqlite3_set_column( column(9), qx )
      call sqlite3_set_column( column(10), result(1) )
      call sqlite3_set_column( column(11), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi_single_band', column )

      call sqlite3_set_column( column(8), qy )
      call sqlite3_set_column( column(9), -qx )
      call sqlite3_set_column( column(10), result(1) )
      call sqlite3_set_column( column(11), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi_single_band', column )

      call sqlite3_set_column( column(8), -qy )
      call sqlite3_set_column( column(9), qx )
      call sqlite3_set_column( column(10), result(1) )
      call sqlite3_set_column( column(11), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi_single_band', column )

      call sqlite3_set_column( column(8), -qy )
      call sqlite3_set_column( column(9), -qx )
      call sqlite3_set_column( column(10), result(1) )
      call sqlite3_set_column( column(11), dabs(result(1)) )
      call sqlite3_insert( db, 'susy_qpi_single_band', column )

      call sqlite3_commit( db )

!     Calculate percentage complete and time elapsed and flush stdout
      call cpu_time(end)
      print *, 0.000194175*(0.5*(iqx*(iqx+1))+iqy), end-start
      call flush()
    end do
  end do

  call sqlite3_begin( db )
  call sqlite3_insert( db, 'runs_single_band', column)
  call sqlite3_commit( db )
  call sqlite3_close( db )

  open(log, file="susy_qpi.log", position="append", status="old")
    call cpu_time(end)
    write(log,*) "CALCULATION TIME: ", end-start
  close(log)

! Export plot to png
  call system('sqlite3 -column data.db "select distinct qx, qy, absresult from susy_qpi_single_band where date='&
    //date//' and time='//time//' order by qx, qy;" | awk -f add_blanks.awk > data/'//filename)
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
