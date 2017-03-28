! Use a cubic spline to approximate a function.

module constants

  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: ndp = kind(1.0d0)

  integer, parameter :: npts = 20
  integer, parameter :: mpts = 200

  real(ndp), parameter :: zero = 0.0d0
  real(ndp), parameter :: one = 1.0d0
  real(ndp), parameter :: two = 2.0d0

end module constants

module common

  use constants

  implicit real(ndp) (a-h,o-z)
  implicit integer(i4b) (i-n)

end module common

module procedures

  use constants

  implicit real(ndp) (a-h,o-z)
  implicit integer(i4b) (i-n)

contains

  subroutine dodp(fpts,flevel,fdp,npts)
    intent(in) :: fpts, flevel, npts 
    intent(out) :: fdp !!  second derivative

    dimension vdp(npts),a(npts),b(npts),c(npts),r(npts),fpts(npts),flevel(npts), &
         fdp(npts)

    an = -one
    c1 = -one
    ! take care of first and last second derivatives
    a(1) = zero
    a(npts) = an
    b(1) = one
    b(npts) = one
    c(1) = c1
    c(npts) = zero
    r(1) = zero
    r(npts) = zero

    do i = 2,npts-1
       a(i) = (fpts(i)-fpts(i-1))/6.0d0
       b(i) = (fpts(i+1)-fpts(i-1))/3.0d0
       c(i) = (fpts(i+1)-fpts(i))/6.0d0
       r(i) = (flevel(i+1)-flevel(i))/ &
            (fpts(i+1)-fpts(i)) - &
            (flevel(i)-flevel(i-1))/ &
            (fpts(i)-fpts(i-1))
    enddo

    call tridag(a,b,c,r,vdp,npts)

    fdp = vdp

  end subroutine dodp

  subroutine interp(x,fpts,flevel,fdp,y,yp,ny,nyp,npts)

    intent(in) :: x, fpts, flevel, fdp, ny, nyp, npts
    intent(out) :: y, yp
    
    ! x at which i want to interpolate
    ! fdp ... second derivatives
    ! y ... interpolated level
    ! yp ... interpolated derivative
    ! ny, nyp ... flags if y or yp should be computed
    
    dimension fpts(npts),flevel(npts),fdp(npts)

    klo = 1
    khi = npts
1   if ((khi-klo) > 1) then
       k = (khi+klo)/2
       if (fpts(k) > x) then
          khi = k
       else
          klo = k
       endif
       goto 1
    endif

    h = fpts(khi) - fpts(klo)
    a = (fpts(khi) - x)/h
    b = (x - fpts(klo))/h
    asq = a*a
    bsq = b*b

    if (ny == 1)  &
         y = a*flevel(klo) + b*flevel(khi) + &
         ((asq*a-a)*fdp(klo) &
         +(bsq*b-b)*fdp(khi))*(h*h)/6.0D+00

    if (nyp == 1) &
         yp = (flevel(khi)-flevel(klo))/h - &
         (3.0d0*asq-one)/6.0d0*h*fdp(klo) + &
         (3.0d0*bsq-one)/6.0d0*h*fdp(khi)

    return

  end subroutine interp

  subroutine tridag(a,b,c,r,u,n)

    real(ndp), parameter :: toler = 1.0d-12

    dimension gam(n),a(n),b(n),c(n),r(n),u(n)

    bet = b(1)
    u(1) = r(1)/bet
    do j = 2,n
       gam(j) = c(j-1)/bet
       bet = b(j) - a(j)*gam(j)
       if (dabs(bet) <= toler) then
          write(6,"(' Failure in subroutine tridag')")
          call wait
       endif
       u(j) = (r(j)-a(j)*u(j-1))/bet
    enddo
    do j = n-1,1,-1
       u(j) = u(j) - gam(j+1)*u(j+1)
    enddo

  end subroutine tridag

  subroutine dogrid(xpts,xlow,xhigh,xinc,npts)

    dimension xpts(npts)

    xinc = (xhigh-xlow)/dble(real(npts-1))

    xpts(1) = xlow
    do i = 2,npts
       xpts(i) = xpts(i-1) + xinc
    enddo

    return

  end subroutine dogrid

  double precision function fct(x)

    fct = dlog(x)

    return

  end function fct

  subroutine wait

    write(6,"(' Waiting...')")
    read *

    return

  end subroutine wait

end module procedures

program main

  use constants
  use common
  use procedures

  implicit real(ndp) (a-h,o-z)
  implicit integer(i4b) (i-n)

  dimension x(npts),y(npts),ydp(npts),xfine(mpts)

  xlow = 0.1d0
  xhigh = 2.0d0

  call dogrid(x,xlow,xhigh,xinc,npts)

  do i = 1,npts   
     y(i) = fct(x(i))
  enddo

  call dodp(x,y,ydp,npts)

  do i = 1,npts
     call interp(x(i),x,y,ydp,yval,yvalp,1,1,npts)
     write(6,"(i6,7f15.8)") i,x(i),y(i),yval,one/x(i),yvalp,-one/(x(i)*x(i)),ydp(i)
  enddo

  call dogrid(xfine,xlow,xhigh,xinc,mpts)

  open(1,file='spline1.dat')

  do i = 1,mpts
     xval = xfine(i)
     yval = fct(xval)
     call interp(xval,x,y,ydp,yval2,yvalp,1,0,npts)
     diff = yval2 - yval
     write(1,"(i6,4f15.8)") i,xval,yval,yval2,diff
  enddo

  close(1)

end program main

