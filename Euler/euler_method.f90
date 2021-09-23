!code name: 'euler_method.f90'
!
!description: This code aims to apply the Euler's method,
!             to solve a ordinary differential
!             equation (ODE).
!
!             problem: we want to find y(4)
!
!             we have the ODE:
!                   y'=y, y(0)=1
!
!
!reference: Wikipedia contributors. (2021, September 19). Euler method.
!           In Wikipedia, The Free Encyclopedia. Retrieved 15:59,
!           September 22, 2021,
!           from https://en.wikipedia.org/w/index.php?title=Euler_method&oldid=1045231921
!
!date: September 23rd, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!_______________________________________________________________________
program euler_code
implicit none
real::h !step size
real::y !solution
real::f !function
real::t !time variable
real::tf !final-time
integer::Nt !total number of steps
integer::i !auxiliary variable

h=1.0
t=0.0
y=1.0
tf=4.0
Nt=int(tf/h)
!write the initial values
write(*,992)t,y
!evolve the system by Euler's method
do i=1,Nt
   t=i*h
   y=y+h*f(y)
   write(*,992)t,y
end do
read*

992 format(2f10.2)

end program euler_code
!***************************************
!function y'(t,y)=f(t,y) , f(t,y)=y
function f(y)
implicit none
real::y,f

f=y

end function f
