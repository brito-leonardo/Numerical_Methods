!code name: 'simpson_1d_method.f90'
!
!description: This code aims to apply the Simpson's rule to integrate
!             a one-dimensional function f(x)
!            (we advise that a choice of an odd number Nx
!            of points can be more suitable)
!
!     scheme:
!               result=dx*[f(x1)+4*f(x2)+f(x3)]/3
!                        + dx*[f(x3)+4*f(x4)+f(x5)]/3
!                        + dx*[f(x5)+4*f(x6)+f(x7)]/3
!                        + ...
!                        + dx*[f(xN-2)+4*f(xN-1)+f(xN)]/3
!
!     problem: we want to find integral of the function
!     defined by
!              f(x)=sin(x)**2
!     in the interval x=[0,pi].
!     We know the exact answer: Result=pi/2=1.570796327...
!
!reference: https://en.wikipedia.org/wiki/Simpson%27s_rule
!
!date: July 27th, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!_______________________________________________________________________
program simpson
implicit none
integer,parameter::Nx=100!arrays size
real,dimension(1:Nx)::x !spatial array
real,dimension(1:Nx)::f !function to be integrated
real,parameter::pi=dacos(-1.0d0) !usual way to define the number pi=3,1415...
real::a,b !initial and final points of the integration
real::dx !spatial-step length
real:: result !integral result

!set the integration interval x=[a,b]
a=0.0
b=pi
!define a function to be integrated
!and the spatial interval where
!the function shall be integrated
call func(Nx,a,b,dx,x,f)

!call the integration routine
call simpson_rout(Nx,dx,f,result)

write(*,*)"The Simpson's rule returns the integration"
write(*,*)'of the function f(x)=sin(x)**2'
write(*,*)'in the interval x=[a,b]=[0,pi] as'
write(*,*)'the result: ', result
write(*,*)'while the exact solution is: ', pi/2


read* !command to stop the terminal

end program simpson
!***********************************************************************
subroutine func(Nx,a,b,dx,x,f)
implicit none
integer::Nx !arrays size
real,dimension(1:Nx)::x !spatial array
real,dimension(1:Nx)::f !function to be integrated
real::dx !spatial-step length
real::a,b
integer::i

!define a interval into a spatial array
!interval is x=[0,pi], then we set Lx=pi

!set the spatial-step size
!(valid for an even or odd number Nx)
dx=(b-a)/(Nx-1)

x(1)=a
do i=2,Nx
   x(i)=(i-1)*dx
end do

!define the fucntion in the spatial grid
do i=1,Nx
   f(i)=sin(x(i))**2
end do


end subroutine func
!***********************************************************************
subroutine simpson_rout(Nx,dx,f,result)
implicit none
integer::Nx ! size of the arrays
real,dimension(1:Nx)::f !function to be integrated
real::dx !spatial-step length
real,dimension(1:Nx)::fx !array of factors
real::result !result of the integration
integer::i !auxiliary integer

!simpson factors
!valid for an even or odd number of points Nx
!odd
if(mod(Nx,2)==0.0)then !even, the error is a little bigger
   fx(1)=1
   do i=2,Nx,2
      fx(i)=4
   end do
   do i=3,Nx-1,2
      fx(i)=2
   end do
   result=0.0d0
   do i=1,Nx
      result=result+(fx(i)*dx*f(i))/3.0d0
   end do
else  !odd
   fx(1)=1
   do i=2,Nx-1,2
      fx(i)=4
   end do
   do i=3,Nx-2,2
      fx(i)=2
   end do
   fx(Nx)=1
   result=0.0d0
   do i=1,Nx
      result=result+(fx(i)*dx*f(i))/3.0d0
   end do
end if
!sum of the trapezium rule

end subroutine simpson_rout
