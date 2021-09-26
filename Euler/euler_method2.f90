!code name: 'euler_method2.f90'
!
!description: This code aims to apply the Euler's method,
!             to solve a ordinary differential
!             equation (ODE).
!
!             problem: we have the ODE:
!                   y'=f(t,y), f(t,y)=y-t**2+1 , y(0)=0.5
!             exact solution: y(t)=1+2*t+t**2-0.5*exp(t)
!
!
!reference: Materials from 'Computational Physics Course'
!           Physics Department (DF)
!           Federal university of SÆo carlos (UFSCAR)
!
!date: September 26th, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!_______________________________________________________________________
program euler_code2
implicit none
real::dt !time-step size
real::y !solution
real::f !function
real::f0 !exact solution
real::t !time variable
real::tf !final-time
integer::Nt !total number of steps
integer::i !auxiliary variable

dt=0.025
t=0.0
y=0.5
tf=0.5
Nt=int(tf/dt)

open(1,file='exact2-fort.txt')
open(2,file='euler2-fort.txt')
!write the initial values
write(1,*)t,f0(t)
write(2,*)t,y
!evolve the system by Euler's method
do i=1,Nt+1
   t=i*dt
   y=y+dt*f(t,y)
   write(1,*)t,f0(t)
   write(2,*)t,y
end do
read*


end program euler_code2
!***************************************
!function y'(t,y)=f(t,y) , f(t,y)=y(t)-t**2+1
function f(t,y)
implicit none
real::y,t,f

f=y-t**2+1.0

end function f
!****************************************
!exact solution
function f0(t)
implicit none
real::t,f0

f0=1.0+2*t+t**2-0.5*exp(t)

end function f0
