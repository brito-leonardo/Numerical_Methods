!code name: 'euler_improved.f90'
!
!description: This code aims to apply the Improved Euler's
!             method, to solve a ordinary differential
!             equation (ODE).
!
!             problem: we want to find y(1.0)
!
!             we have the ODE:
!                   y'=x-2y, y(0)=1
!
!
!reference: Differential Equation with Boundary-Value Problems (2018)
!                                 (Dennis G.Zill)
!           [Exercise 15, section 9.1, p.373]
!
!date: September 23rd, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!_______________________________________________________________________
program euleri_code
implicit none
real::dx !step size
real::yaux,y !solution
real::yex !exact solution
real::f !function
real::x0,xold,x !x-variable
real::xf !final-point
integer::Nx !total number of steps
integer::i !auxiliary variable

dx=0.01
x0=0.0
x=x0
xf=1.0
y=1.0
Nx=int(xf/dx)
!write the initial values
!evolve the system by the Improved Euler's method
do i=1,Nx
   xold=x
   yaux=y+dx*f(x,y)
   x=x0+i*dx
   y=y+dx*(f(xold,y)+f(x,yaux))/2.0
end do

write(*,*)'Exact solution yex(1.5): ',yex(xf)
write(*,*)'Improved Euler method y(1.5): ',y

read*

end program euleri_code
!***************************************
!function y'(t,y)=f(t,y) , f(t,y)=y
function f(x,y)
implicit none
real::x,y,f

f=x-2*y

end function f
!****************************************
!exact solution yex
function yex(x)
implicit none
real::x,yex

yex=(1.0/2)*x-(1.0/4)+(5.0/4)*exp(-2*x)

end function yex
