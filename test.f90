program test
  implicit none
  integer :: start, i
  call system_clock(start)
  do i=1,97
    call sleep(1)
    call progress(i/97.0,start)
  end do
end program test

subroutine progress(percent, start)
  implicit none
  integer, intent(in) :: start
  real, intent(in)    :: percent
  integer             :: end, rate, elapsed, remaining
  character(len=102)  :: bar

  call system_clock(end, rate)
  elapsed   = real(end-start)/real(rate)
  remaining = elapsed*(1.0/percent-1.0)
  bar  = "["//repeat("=",int(percent*100))//repeat(" ",int((1.0-percent)*100))//"]"
  write(*,"(A,I3,'% ',I4,':',I2.2,' elapsed',I4,':',I2.2,' remaining')") &
    bar, int(percent*100), elapsed/60, mod(elapsed,60), remaining/60, mod(remaining,60)
end subroutine progress