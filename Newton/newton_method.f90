!code name: 'newton_method.f90'
!
!description: This code aims to apply the Newton's method,
!              to find roots.
!
!             problem: we want to find x that satisfies cos(x)=x**3
!
!             we define:
!                   f(x)=cos(x)-x**3
!             so, we have the first derivative:
!                   f'(x)=-sin(x)-3x**2
!
!
!reference: https://en.wikipedia.org/wiki/Newton%27s_method
!
!date: July 9th, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!_______________________________________________________________________
program newton_code
implicit none
real::f,fd !fucntion and its first derivative
real::err !error
real::epsilon !convergence criterium
real::x0,xold,x !solutions

!solution guess x0
x0=0.5
!convergence criterium
epsilon=1.e-8
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


end program newton_code
!*********************************
function f(x0)
implicit none
real::f,x0

f=cos(x0)-x0**3

end function f
!*********************************
function fd(x0)
implicit none
real::fd,x0

fd=-sin(x0)-3*x0**2

end function fd
