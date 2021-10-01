!code name: 'euler_method4.f90'
!
!description: This code aims to apply the Euler's method
!             to solve a second order ordinary differential
!             equation (ODE) as a system.
!
!             problem: we have an second order ODE :
!                   y''-4y'=12x, y(0)=4, y'(0)=1
!             exact solution: y(x)=exp(-x/2)*(-cos(x)+(3/4)*sin(2*x))
!             we handle this problem rewriting it as a ODE system
!                 y1'(x)=y2(x) , y1(0)=4
!                 y2'(x)=y''(x)=4y+12x=4y1+12x, y2(0)=1
!
!
!reference: Differential Equation with Boundary-Value Problems (2018)
!                                 (Dennis G.Zill)
!           [Example 2. section 4.1.1, p.120]
!
!date: October 1st, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!_______________________________________________________________________
program euler_code4
implicit none
real::dx !x-step size
real,dimension(1:2)::y!solution
real::yex! exact solution
real,dimension(1:2)::f !function
real::x0,x !time variable
real::xf !final-time
integer::Nx !total number of steps
integer::i !auxiliary variable

dx=0.01
x0=0.0
y(1)=4.0
y(2)=1.0
xf=2.0
Nx=int((xf-x0)/dx)

open(1,file='n-exact-yex.txt')
open(2,file='n-euler-y.txt')
!write the initial values
write(1,*)x0,yex(x0)
write(2,*)x0,y(1)
!evolve the system by Euler's method
do i=1,Nx
   x=x0+i*dx
   call func(x,y,f)
   y=y+dx*f
   write(1,*)x,yex(x)
   write(2,*)x,y(1)
end do
read*


end program euler_code4
!***************************************
!function y'(t,y)=f(t,y) , f(t,y)=y(t)-t**2+1
subroutine func(x,y,f)
implicit none
real,dimension(1:2)::y,f
real::x

f(1)=y(2)
f(2)=4*y(1)+12*x

end subroutine func
!******************************************
!exact solution yex
function yex(x)
implicit none
real::x,yex

yex=3*exp(2*x)+exp(-2*x)-3*x

end function yex
