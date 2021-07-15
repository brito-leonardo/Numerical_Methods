 !code name: 'newton_FD_method.f90'
!
!description: This code aims to apply the Newton's method,
!              to find roots.
!
!             problem: we want to find x that satisfies cos(x)=x**3
!
!             we define:
!                   f(x)=cos(x)-x**3
!             suppose we don't know its first derivative,
!             we will use finite differences (FD) method
!             to calculate the first derivative
!
!reference: https://en.wikipedia.org/wiki/Newton%27s_method
!
!date: July 9th, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!_______________________________________________________________________
program newton_FD_code
implicit none
real::f,fd !fucntion and its first derivative
real::err !error
real::epsilon !convergence criterion
real::x0,xold,x !solutions

!solution guess x0
x0=0.5
!convergence criterion
epsilon=1.e-6
!set starting variables
xold=x0
x=x0
err=1.0
!Newton's scheme
do while(err>epsilon)
   x=x-f(x)/fd(x)
   err=abs(x-xold)
   xold=x
end do

write(*,*)'the solution is x: ',x


end program newton_FD_code
!*********************************
function f(x0)
implicit none
real::f,x0

f=cos(x0)-x0**3

end function f
!*********************************
! here we calculate the firts derivative by finite differences
! it is useful to see: https://en.wikipedia.org/wiki/Finite_difference_method
!
function fd(x)
implicit none
real::f,x
real::dx
real::fd,x1,x2

!set an arbitrary step
dx=0.0001
x1=x
x2=x+dx
fd=(f(x2)-f(x1))/(dx)

end function fd
