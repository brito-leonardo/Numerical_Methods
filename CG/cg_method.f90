!code name: 'cg_method.f90'
!
!description: This code aims to apply the Conjugate gradient (CG) method,
!              to solve a linear matrix problem.
!
!             problem: A*x=b
!
!     A=(4  1)  ; x=(x1) ; b=(1)
!       (1  3)      (x2)     (2)
!
!reference: https://en.wikipedia.org/wiki/Conjugate_gradient_method
!
!date: July 9th, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of São Paulo (USP)
!_______________________________________________________________________
program cg_code
implicit none
integer,parameter::N=2
real,dimension(1:N,1:N)::A
real,dimension(1:N)::x0,x !solutions
real,dimension(1:N)::b !rght side of the equation
real::epsilon
integer::i,j

!define A,B and the guess solution x0
call initial(N,A,x0,b)
write(*,*)'The initial conditions are:'
write(*,*)
write(*,*)'A:'
do i=1,N
   write(*,*)(A(i,j), j=1,N)
end do
write(*,*)
write(*,*)'b:'
write(*,*)b
write(*,*)
write(*,*)'x0:'
write(*,*)x0
write(*,*)
write(*,*)

!set the initial guess
x=x0
!set the convergence criterium
epsilon=1.e-6
!call the CG method
call CG(N,A,x,b,epsilon)

write(*,*)"The solution is x: "
write(*,*)x

pause

end program cg_code
!*********************************
subroutine initial(N,A,x0,b)
implicit none
integer::N
real,dimension(1:N,1:N)::A
real,dimension(1:N)::x0,b

!define matrix A
A(1,1)=4.0;A(1,2)=1.0
A(2,1)=1.0;A(2,2)=3.0
!define right side of the equation
b(1)=1.0;b(2)=2.0
!set the solution initial guess
x0=0.0
!sometimes we need a more accurated guess

end subroutine initial
!*********************************
subroutine CG(N,A,x,b,epsilon)
implicit none
integer::N
real,dimension(1:N,1:N)::A
real,dimension(1:N)::x,b
real::epsilon
real,dimension(1:N)::r0,p,r
real::alpha,beta,err

r0=b-matmul(A,x)
p=r0

!set an initial error
err=1.0
do while(err>epsilon)
   alpha=dot_product(r0,r0)/dot_product(p,matmul(A,p))
   x=x+alpha*p
   r=r0-alpha*matmul(A,p)
   !calculate the error
   err=maxval(abs(r))
   !update variables
   beta=dot_product(r,r)/dot_product(r0,r0)
   p=r+beta*p
   r0=r
end do

end subroutine CG
