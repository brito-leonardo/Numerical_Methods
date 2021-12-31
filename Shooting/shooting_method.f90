!code name: 'shooting_method.f90'
!
!description:   We want to solve an ODE when we only have
!               the solution at the first and the final point
!               and we guess the first derivative at the starting point
!               by the 'shooting method', i.e, we try two arbitrary guesses
!               and then we estimate a third better guess.
!               To evolve our guesses, we use the euler method.
!
!             problem: we want to solve the ODE
!
!                   x''-(1-t/5)x=t, x(1)=2, x(3)=-1
!
!             x'=v
!             v'=x''=t+(1-t/5)x
!             x(1)=2 , x'(1)=v(1)=? (here we use a guess, the better guess
!                                    will find x(3)=-1)
!
!             we map
!                    x(1)=x   and   f(1)=x(2)
!                    x(2)=v         f(2)=t+(1-t/5)x(1)
!
!reference:    Applied Bumerical Analysis, 1989
!          (Curtis F. Gerald, Patrick O. Wheatley)
!                    [Section 6.1, p.412]
!
!date: December 31th, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of São Paulo (USP)
!_______________________________________________________________________
program shooting_code
implicit none
real::dt !step size
real::t0,tf !starting and final points
real::x0,xf !solution at the first and the last points
real::D !desired solution, i.e, the final point
real::G1,G2,G !fisrt derivative guesses
real::R1,R2,R !final solution by euler method

!enter the parameters
dt=0.2
t0=1.0
tf=3.0
x0=2.0
xf=-1.0

!desired final value
D=xf

!guess 1 for the first derivative x'(t0)
write(*,*)'Guess 1'
write(*,*)
G1=-1.5
call euler(dt,t0,tf,x0,G1,R1)
!guess2
write(*,*)'Guess 2'
write(*,*)
G2=-3.0
call euler(dt,t0,tf,x0,G2,R2)
!estimate a better guess
write(*,*)'The better Guess'
write(*,*)
G=G1+(G2-G1)*(D-R1)/(R2-R1)
call euler(dt,t0,tf,x0,G,R)

write(*,*)'The better guess got the endpoint x:', R
write(*,*)
write(*,*)"The program 'shooting_method.f90' is finished!"
read*

end program shooting_code
!***************************************
subroutine euler(dt,t0,tf,x0,xd,xr)
implicit none
real::dt,t,t0,tf
real::x0,
real::xd !first derivative guess at t=t0
real::xr !final point by euler evolution
real,dimension(1:2)::x !solution
real,dimension(1:2)::f !function

t=t0
x(1)=x0
x(2)=xd

write(*,*)
write(*,*)'-----------------------------------------'
write(*,*)'first derivative guess xd:',xd
write(*,*)
write(*,993)t,x(1),x(2)
do while(t<tf)
   call func(t,x,f)
   x=x+dt*f
   t=t+dt
   write(*,993)t,x(1),x(2)
end do
write(*,*)
write(*,*)'-----------------------------------------'
write(*,*)
write(*,*)

xr=x(1)

993 format(3f10.3)

end subroutine
!*****************************************
!function x''(t,x)=f(t,x) , f(t,x)=t+(1-t/5)*x
subroutine func(t,x,f)
implicit none
real,dimension(1:2)::x,f
real::t

f(1)=x(2)
f(2)=t+(1.0-t/5.0)*x(1)

end subroutine func
