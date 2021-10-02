!code name: 'rk4_method3.f90'
!
!description: This code aims to apply the 4th order Runge-Kutta (RK4)
!             method to solve a system of first order ordinary differential
!             equations (ODE).
!
!             problem: we have an second order ODE :
!                   y''-4y'+4y=0, y(0)=-2, y'(0)=1
!
!             exact solution: y(t)=(-2+5t)exp(2t)
!
!             we want to find the solutions y(0.2)=x1(0.6) and y(0.6)=x2(0.6)
!             we conveniently rewrite the ODE as asystem with variables y1=y and y'=y2
!                   y1'=y2, y1(0)=-2
!                   y2'=4y2-4y1, y2(0)=1
!
!
!reference: Differential Equation with Boundary-Value Problems (2018)
!                                 (Dennis G. Zill)
!           [problem 3, section 9.4, p.385]
!
!date: October 2nd, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!_______________________________________________________________________
program  rk4_code3
implicit none
real::dt !step-size
real::t0,tf,t !t-variable
real,dimension(1:2)::y0,y !rk4  approximate solution
!we store y=y1 at y(1) and y'=y2=y2 at y(2)
real::y_ex !exact solution
integer::Nt ! number of steps
integer::i


t0=0.0
tf=0.2
!initial value
y0(1)=-2.0 !y(0)
y0(2)=1.0  !y'(0)

write(*,*)'exact solution y_ex(tf): ',y_ex(tf)
write(*,*)
!aproximate solution for dt=0.2
t=t0
y=y0
dt=0.1
Nt=int((tf-t0)/dt)
write(*,*)'KK4 aproximation when dt is ',dt
write(*,*)'Columns t and y'
do i=1,Nt+1
   write(*,992)t,y(1)
   call RK4_rout(dt,t,y)
   t=t0+i*dt
end do
write(*,*)
!aproximate solution for dt=0.2
t=t0
y=y0
dt=0.05
Nt=int((tf-t0)/dt)
write(*,*)'KK4 aproximation when dt is ',dt
write(*,*)'Columns t and y'
do i=1,Nt+1
   write(*,992)t,y(1)
   call RK4_rout(dt,t,y)
   t=t0+i*dt
end do
write(*,*)

992 format(2f10.4)

read*

end program rk4_code3
!*******************************************
subroutine RK4_rout(dt,t,y)
implicit none
real,dimension(1:2)::y !RK4 aproximate solution
real,dimension(1:2)::f !function
real,dimension(1:2)::k1,k2,k3,k4 !auxiliary variables
real::dt,t

!evolve the system by the RK4  method
call func(t,y,f)
k1=f
call func(t+0.5*dt,y+0.5*dt*k1,f)
k2=f
call func(t+0.5*dt,y+0.5*dt*k2,f)
k3=f
call func(t+dt,y+dt*k3,f)
k4=f
y=y+dt*(k1+2*k2+2*k3+k4)/6.0

end subroutine RK4_rout
!*******************************************
!function f(x)=f(x)
subroutine func(t,y,f)
implicit none
real,dimension(1:2)::y,f
real::t
f(1)=y(2)
f(2)=4*y(2)-4*y(1)
end subroutine func
!*******************************
!exact solution y_ex(t)=(-2+5t)exp(2t)
function y_ex(t)
implicit none
real::t,y_ex
y_ex=(-2.0+5*t)*exp(2*t)
end function y_ex

