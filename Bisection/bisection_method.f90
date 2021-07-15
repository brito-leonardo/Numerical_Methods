!code name: 'bisection_method.f90'
!
!description: This code aims to apply the Bisection method,
!              to find the root of a fucntion.
!
!             problem: f(x)=x**3-x-2
!
!
!reference: https://en.wikipedia.org/wiki/Bisection_method
!
!date: July 10th, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!_______________________________________________________________________
program bisec_code
implicit none
real::a,b,c !auxiliary variables
real::x !roots
real::f !function
real::epsilon !convergence criterion
real::err !error



!enter an initial interval [a,b] (we should have some insight of a good interval)
a=1.0
b=2.0
!this choice can be understood if we note
!f(a=1)=-2 and f(b=2)=+4
!there is a root within the interval [1,2]
!because the function change its sign.

!set a convergence criterion
epsilon=1.e-6
!set an arbitrary error
err=1.0
!bisection algorithm (see the reference)
do while(err>epsilon)
   !calculate the middle point
   c=(a+b)/2.0
   if(f(a)*f(c)<0.0)then !if 'true' there a change of sign in the interval [a,c]
      b=c   !if 'true' we make the interval closer to 'a'
      else
      a=c !if 'false' we make the interval closer to 'b'
   end if
   err=abs(f(c))
end do
!the root is the last middle point 'c' found
x=c

write(*,*)'The root is: ', x




end program bisec_code
!***********************************************
function f(x)
implicit none
real::f !function
real::x!argument

f=x**3-x-2.0

end function f
