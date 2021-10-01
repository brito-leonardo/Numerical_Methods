!code name: 'euler_method3.f90'
!
!description: This code aims to apply the Euler's method
!             to solve an ordinary differential
!             equation (ODE) system.
!
!             problem: we have 1st order ODE system:
!                   x1'=-(2/25)*x1+(1/50)*x2, x1(0)=25
!                   x2'=(2/25)*x1-(2/25)*x2, x2(0)=0
!             exact solutions: x1(t)=(25/2)*exp(-t/25)+(25/2)*exp(-3*t/25)
!                              x2(t)=25*exp(-t/25)-25*exp(-3*t/25
!
!
!reference: Differential Equation with Boundary-Value Problems (2018)
!                                 (Dennis G.Zill)
!           [Example 3. section 4.9, p.186-187]
!
!date: October 1st, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!_______________________________________________________________________
program euler_code3
implicit none
real::dt !x-step size
real,dimension(1:2)::x!solution
real::x1_ex,x2_ex! exact solution
real,dimension(1:2)::f !function
real::t0,t !time variable
real::tf !final-time
integer::Nt !total number of steps
integer::i !auxiliary variable

dt=0.01
t0=0.0
x(1)=25.0
x(2)=0.0
tf=100.0
Nt=int((tf-t0)/dt)

open(1,file='x1-exact.txt')
open(2,file='x1-euler.txt')
open(3,file='x2-exact.txt')
open(4,file='x2-euler.txt')

!write the initial values
write(1,*)t0,x1_ex(t0)
write(2,*)t0,x(1)
write(3,*)t0,x2_ex(t0)
write(4,*)t0,x(2)
!evolve the system by Euler's method
do i=1,Nt
   t=t0+i*dt
   call func(t,x,f)
   x=x+dt*f
   write(1,*)t,x1_ex(t)
   write(2,*)t,x(1)
   write(3,*)t,x2_ex(t)
   write(4,*)t,x(2)
end do

read*


end program euler_code3
!***************************************
!function y'(t,y)=f(t,y) , f(t,y)=y(t)-t**2+1
subroutine func(t,x,f)
implicit none
real,dimension(1:2)::x,f
real::t

f(1)=-(2.0/25)*x(1)+(1.0/50)*x(2)
f(2)=(2.0/25)*x(1)-(2.0/25)*x(2)

end subroutine func
!******************************************
!exact solution x1_ex
function x1_ex(t)
implicit none
real::t,x1_ex

x1_ex=(25.0/2)*exp(-t/25)+(25.0/2)*exp(-3*t/25)

end function x1_ex
!******************************************
!exact solution x2_ex
function x2_ex(t)
implicit none
real::t,x2_ex

x2_ex=25.0*exp(-t/25)-25.0*exp(-3*t/25)

end function x2_ex
