!code name: 'rk4_method2.f90'
!
!description: This code aims to apply the 4th order Runge-Kutta (RK4)
!             method to solve a system of first order ordinary differential
!             equations (ODE).
!
!             problem: we have an ODE system:
!                   x'=2x+4y, x(0)=-1
!                   y'=-x+6y, y(0)=6
!             exact solutions: x(t)=(26t-1)exp(4t)
!                              y(t)=(13t+6)exp(4t)
!
!             we want to find the solutions x(0.6)=x1(0.6) and y(0.6)=x2(0.6)
!             we conveniently rewrite the system with variables x1 and x2
!                   x1'=2x1+4x2, x1(0)=-1
!                   x2'=-x1+6x2, x2(0)=6
!
!
!reference: Differential Equation with Boundary-Value Problems (2018)
!                                 (Dennis G. Zill)
!           [Example 3, section 9.4, p.383-384]
!
!date: October 2nd, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!_______________________________________________________________________
program rk4_code2
implicit none
real::dt !step-size
real::t0,tf,t !t-variable
real,dimension(1:2)::x0,x !rk4  approximate solution
!we store x1 at x(1) and x2 at x(2)
real::x1_ex,x2_ex !exact solution
integer::Nt ! number of steps
integer::i


t0=0.0
tf=0.6
!initial value
x0(1)=-1.0 !x1(0)
x0(2)=6.0  !x2(0)

write(*,*)'exact solution x1_ex(tf): ',x1_ex(tf)
write(*,*)'exact solution x2_ex(tf): ',x2_ex(tf)
write(*,*)
!aproximate solution for dt=0.2
! see table 9.4.1 form reference
t=t0
x=x0
dt=0.2
Nt=int((tf-t0)/dt)
write(*,*)'KK4 aproximation when dt is ',dt
write(*,*)'Columns t, x1 and x2 '
do i=1,Nt+1
   write(*,993)t,x(1),x(2)
   call RK4_rout(dt,t,x)
   t=t0+i*dt
end do
write(*,*)
!aproximate solution for dt=0.2
! see table 9.4.2 form reference
t=t0
x=x0
dt=0.1
Nt=int((tf-t0)/dt)
write(*,*)'KK4 aproximation when dt is ',dt
write(*,*)'Columns t, x1 and x2 '
do i=1,Nt+1
   write(*,993)t,x(1),x(2)
   call RK4_rout(dt,t,x)
   t=t0+i*dt
end do
write(*,*)


993 format(3f10.4)

read*

end program rk4_code2
!*********************
subroutine RK4_rout(dt,t,x)
implicit none
real,dimension(1:2)::x !RK4 aproximate solution
real,dimension(1:2)::f !function
real,dimension(1:2)::k1,k2,k3,k4 !auxiliary variables
real::dt,t

!evolve the system by the RK4  method
call func(t,x,f)
k1=f
call func(t+0.5*dt,x+0.5*dt*k1,f)
k2=f
call func(t+0.5*dt,x+0.5*dt*k2,f)
k3=f
call func(t+dt,x+dt*k3,f)
k4=f
x=x+dt*(k1+2*k2+2*k3+k4)/6.0

end subroutine RK4_rout
!**********************
!function f(x,y)=f(x1,x2)
subroutine func(t,x,f)
implicit none
real,dimension(1:2)::x,f
real::t
f(1)=2*x(1)+4*x(2)
f(2)=-x(1)+6*x(2)
end subroutine func
!*************************************
!exact solution x1_ex(t)=(26t-1)exp(4t)
function x1_ex(t)
implicit none
real::t,x1_ex
x1_ex=(26*t-1.0)*exp(4*t)
end function x1_ex
!*************************************
!exact solution x2_ex(t)=(13t+6)exp(4t)
function x2_ex(t)
implicit none
real::t,x2_ex
x2_ex=(13*t+6.0)*exp(4*t)
end function x2_ex
