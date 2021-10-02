!code name: 'rk4_method.f90'
!
!description: This code aims to apply the 4th order Runge-Kutta (RK4)
!             method to solve a first order ordinary differential
!             equation (ODE).
!
!             problem: we have an second order ODE :
!                   y'=2xy, y(1)=1
!             exact solution: y(x)=exp(x**2)/e , where 'e' is the Euler's number, e=2,71828...
!
!
!reference: Differential Equation with Boundary-Value Problems (2018)
!                                 (Dennis G. Zill)
!           [Example 1, section 9.2, p.375-376]
!
!date: October 1st, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!_______________________________________________________________________
program rk4_code
implicit none
real::dx !step-size
real::x0,xf,x !x-variable
real::y !rk4  approximate solution
real::yex !exact solution
integer::Nx !number of steps
real::f !function
real::k1,k2,k3,k4 !auxiliary variables
integer::i


dx=0.1
x0=1.0
x=x0
xf=1.5
Nx=int((xf-x0)/dx)
!initial value
y=1.0

!evolve the system by the RK4  method
do i=1,Nx+1
   k1=f(x,y)
   k2=f(x+0.5*dx,y+0.5*dx*k1)
   k3=f(x+0.5*dx,y+0.5*dx*k2)
   k4=f(x+dx,y+dx*k3)
   y=y+dx*(k1+2*k2+2*k3+k4)/6.0
   x=x0+i*dx
end do

write(*,*)'exact solution yex(xf): ',yex(xf)
write(*,*)'RK4 aproximation y(xf): ',x,y

read*

end program rk4_code
!*********************
!function f(x,y)=2xy
function f(x,y)
implicit none
real::x,y,f
f=2*x*y
end function f
!*********************
!exact solution yex(x,y)=2xy
function yex(x)
implicit none
real,parameter::e=exp(1.0)
real::x,yex
yex=exp(x**2)/e
end function yex
