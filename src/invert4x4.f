c      FORTRAN implementation of gluInvertMatrix by Eric Mascot
c      http://www.mesa3d.org/
c
c      Inverts a 4 by 4 matrix
c
c      m - matrix to invert
c        array of size 16
c
c      invOut - inverted matrix
c        array of size 16
c
c      fail - status of result
c        0 - Normal exit
c        1 - Singular matrix
c
c============================================================
c
c      program example
c        implicit none
c        complex, dimension(16) :: a,b
c        integer :: status
c
c        complex, parameter :: ci = dcmplx(0.d0,1.d0)
c
c        a = (/ complex :: &
c          ci, 0, 0, 0, &
c           0,-1, 0, 0, &
c           0, 0, 1, 0, &
c           0, 0, 0,-1 /)
c
c        call invert4x4(a,b,status)
c
c        if (status.eq.0) then
c          print *, b
c        else
c          print *, "Inversion failed"
c        end if
c
c      end program example

      subroutine invert4x4(m, invOut, fail)
      implicit none

      complex(kind=8), dimension(16), intent(in) :: m
      complex(kind=8), dimension(16), intent(out) :: invOut
      integer, intent(out) :: fail

      integer :: i
      complex(kind=8) :: det
      complex(kind=8), dimension(16) :: inv

      inv(1) = m(6)  * m(11) * m(16) -
     +         m(6)  * m(12) * m(15) -
     +         m(10) * m(7)  * m(16) +
     +         m(10) * m(8)  * m(15) +
     +         m(14) * m(7)  * m(12) -
     +         m(14) * m(8)  * m(11)

      inv(5) = -m(5)  * m(11) * m(16) +
     +          m(5)  * m(12) * m(15) +
     +          m(9)  * m(7)  * m(16) -
     +          m(9)  * m(8)  * m(15) -
     +          m(13) * m(7)  * m(12) +
     +          m(13) * m(8)  * m(11)

      inv(9) = m(5)  * m(10) * m(16) -
     +         m(5)  * m(12) * m(14) -
     +         m(9)  * m(6)  * m(16) +
     +         m(9)  * m(8)  * m(14) +
     +         m(13) * m(6)  * m(12) -
     +         m(13) * m(8)  * m(10)

      inv(13) = -m(5)  * m(10) * m(15) +
     +           m(5)  * m(11) * m(14) +
     +           m(9)  * m(6)  * m(15) -
     +           m(9)  * m(7)  * m(14) -
     +           m(13) * m(6)  * m(11) +
     +           m(13) * m(7)  * m(10)

      inv(2) = -m(2)  * m(11) * m(16) +
     +          m(2)  * m(12) * m(15) +
     +          m(10) * m(3)  * m(16) -
     +          m(10) * m(4)  * m(15) -
     +          m(14) * m(3)  * m(12) +
     +          m(14) * m(4)  * m(11)

      inv(6) = m(1)  * m(11) * m(16) -
     +         m(1)  * m(12) * m(15) -
     +         m(9)  * m(3)  * m(16) +
     +         m(9)  * m(4)  * m(15) +
     +         m(13) * m(3)  * m(12) -
     +         m(13) * m(4)  * m(11)

      inv(10) = -m(1) * m(10) * m(16) +
     +          m(1)  * m(12) * m(14) +
     +          m(9)  * m(2)  * m(16) -
     +          m(9)  * m(4)  * m(14) -
     +          m(13) * m(2)  * m(12) +
     +          m(13) * m(4)  * m(10)

      inv(14) = m(1)  * m(10) * m(15) -
     +          m(1)  * m(11) * m(14) -
     +          m(9)  * m(2)  * m(15) +
     +          m(9)  * m(3)  * m(14) +
     +          m(13) * m(2)  * m(11) -
     +          m(13) * m(3)  * m(10)

      inv(3) = m(2)  * m(7) * m(16) -
     +         m(2)  * m(8) * m(15) -
     +         m(6)  * m(3) * m(16) +
     +         m(6)  * m(4) * m(15) +
     +         m(14) * m(3) * m(8)  -
     +         m(14) * m(4) * m(7)

      inv(7) = -m(1)  * m(7) * m(16) +
     +          m(1)  * m(8) * m(15) +
     +          m(5)  * m(3) * m(16) -
     +          m(5)  * m(4) * m(15) -
     +          m(13) * m(3) * m(8)  +
     +          m(13) * m(4) * m(7)

      inv(11) = m(1)  * m(6) * m(16) -
     +          m(1)  * m(8) * m(14) -
     +          m(5)  * m(2) * m(16) +
     +          m(5)  * m(4) * m(14) +
     +          m(13) * m(2) * m(8)  -
     +          m(13) * m(4) * m(6)

      inv(15) = -m(1)  * m(6) * m(15) +
     +           m(1)  * m(7) * m(14) +
     +           m(5)  * m(2) * m(15) -
     +           m(5)  * m(3) * m(14) -
     +           m(13) * m(2) * m(7)  +
     +           m(13) * m(3) * m(6)

      inv(4) = -m(2)  * m(7) * m(12) +
     +          m(2)  * m(8) * m(11) +
     +          m(6)  * m(3) * m(12) -
     +          m(6)  * m(4) * m(11) -
     +          m(10) * m(3) * m(8)  +
     +          m(10) * m(4) * m(7)

      inv(8) = m(1) * m(7) * m(12) -
     +         m(1) * m(8) * m(11) -
     +         m(5) * m(3) * m(12) +
     +         m(5) * m(4) * m(11) +
     +         m(9) * m(3) * m(8)  -
     +         m(9) * m(4) * m(7)

      inv(12) = -m(1) * m(6) * m(12) +
     +           m(1) * m(8) * m(10) +
     +           m(5) * m(2) * m(12) -
     +           m(5) * m(4) * m(10) -
     +           m(9) * m(2) * m(8)  +
     +           m(9) * m(4) * m(6)

      inv(16) = m(1) * m(6) * m(11) -
     +          m(1) * m(7) * m(10) -
     +          m(5) * m(2) * m(11) +
     +          m(5) * m(3) * m(10) +
     +          m(9) * m(2) * m(7)  -
     +          m(9) * m(3) * m(6)

      det = m(1) * inv(1) +
     +      m(2) * inv(5) +
     +      m(3) * inv(9) +
     +      m(4) * inv(13)

      if (det.eq.0) then
          fail = 1
          return
      end if

      det = 1 / det

      do 10 i = 1,16
          invOut(i) = inv(i) * det
10    continue

      fail = 0
      end subroutine invert4x4