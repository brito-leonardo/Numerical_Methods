!code name: 'rkf45_method.f90'
!
!description: This code aims to apply the 4-5th order 
!             Runge-Kutta-Fehlberg (RKF45)
!             method to solve a first order ordinary differential
!             equation (ODE) with adaptative steps.
!
!             problem: we have an second order ODE :
!                   y'=2xy, y(1)=1
!             exact solution: y(x)=exp(x**2-1) 
!
!method reference RKF45 method: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method#cite_note-1
!
!
!problem reference: Differential Equation with Boundary-Value Problems (2018)
!                                 (Dennis G. Zill)
!                   [Example 1, section 9.2, p.375-376]
!
!
!date: December 17th, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of São Paulo (USP)
!_______________________________________________________________________
program rkf45_code
implicit none
real,dimension(1:6)::a,c,ch,ct
real,dimension(1:6,1:5)::b
real::h,hnew !step-size
real::x0,xf,x !x-variable
real::y0,y !rk4  approximate solution
real::yex !exact solution
integer::Nx !number of steps
real::f !function
real::k1,k2,k3,k4,k5,k6 !auxiliary variables
real::TE
real::epsilon
integer::i,step


open(1,file='data.txt')

call rkf45_coef(a,b,c,ch,ct)

epsilon=1.e-5
h=0.05
x0=1.0
x=x0
xf=1.5
Nx=int((xf-x0)/h)
!initial value
y0=1.0

!evolve the system by the RKF45  method
step=0
write(1,*)step,x,y,yex(x)
do while(x<=xf)
100   k1=h*f(x+a(1)*h,y0)
      k2=h*f(x+a(2)*h,y0+b(2,1)*k1)
      k3=h*f(x+a(3)*h,y0+b(3,1)*k1+b(3,2)*k2)
      k4=h*f(x+a(4)*h,y0+b(4,1)*k1+b(4,2)*k2+b(4,3)*k3)
      k5=h*f(x+a(5)*h,y0+b(5,1)*k1+b(5,2)*k2+b(5,3)*k3+b(5,4)*k4)
      k6=h*f(x+a(6)*h,y0+b(6,1)*k1+b(6,2)*k2+b(6,3)*k3+b(6,4)*k4+b(6,5)*k5)
      y=y0+ch(1)*k1+ch(2)*k2+ch(3)*k3+ch(4)*k4+ch(5)*k5+ch(6)*k6
      TE=ct(1)*k1+ct(2)*k2+ct(3)*k3+ct(4)*k4+ct(5)*k5+ct(6)*k6
      hnew=0.9*h*(epsilon/TE)**(1.0/5.0)
      if(TE>epsilon)then
         h=hnew
         go to 100
      end if
      x=x+h
      y0=y
      h=hnew
      step=step+1
      !the columns follow: step number, x, y(x) and exact solution yex(x).
      write(1,*)step,x,y,yex(x)
end do


end program rkf45_code
!***********************************
subroutine rkf45_coef(a,b,c,ch,ct)
implicit none
real,dimension(1:6)::a,c,ch,ct
real,dimension(1:6,1:5)::b

a(1)=0.0;a(2)=2.0/9.0;a(3)=1.0/3.0;a(4)=3.0/4.0;a(5)=1.0;a(6)=5.0/6.0
c(1)=1.0/9.0;c(2)=0.0;c(3)=9.0/20.0;c(4)=16.0/45.0;c(5)=1.0/12.0
ch(1)=47.0/450.0;ch(2)=0.0;ch(3)=12.0/25.0;ch(4)=32.0/225.0;ch(5)=1.0/30.0;ch(6)=6.0/25.0
ct(1)=-1.0/150.0;ct(2)=0.0;ct(3)=3.0/100.0;ct(4)=-16.0/75.0;ct(5)=-1.0/20.0;ct(6)=6.0/25.0

b(2,1)=2.0/9.0
b(3,1)=1.0/12.0;b(3,2)=1.0/4.0
b(4,1)=69.0/128.0;b(4,2)=-243.0/128.0;b(4,3)=135.0/64.0
b(5,1)=-17.0/12.0;b(5,2)=27.0/4.0;b(5,3)=-27.0/5.0;b(5,4)=16.0/15.0
b(6,1)=5.0/6.0;b(6,2)=-5.0/16.0;b(6,3)=13.0/16.0;b(6,4)=4.0/27.0;b(6,5)=5.0/144.0

end subroutine rkf45_coef
!*************************************
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
yex=exp(x**2-1.0d0)
end function yex
