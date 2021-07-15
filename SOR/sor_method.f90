!code name: 'sor_method.f90'
!
!description: This code aims to solve a linear matrix problem by
!             Sucessive-Over relaxation (SOR) method
!
!             problem:
!                      A*x=b
!
!     consider the linear system
!      4x1 - x2 - 6x3 + 0x4 = 2
!     -5x1 - 4x2 + 10x3 + 8x4 = 21
!      0x1 + 9x2 + 4x3 - 2x4 = -12
!      1x1 + 0x2 - 7x3 + 5x4 = -6
!
!     in the matrix form, it becomes
!
!       ( 4 -1  -6  0)      (x1)     (2)
!       (-5 -4  10  8)      (x2)     (21)
!    A= ( 0  9   4 -2)  ; x=(x3) ; b=(-12)
!       ( 1  0  -7  5)      (x4)     (-6)
!
!reference: https://en.wikipedia.org/wiki/Successive_over-relaxation
!
!date: July 8th, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!_______________________________________________________________________
program sor_code
implicit none
integer,parameter::N=4
real,dimension(1:N)::b,x
real,dimension(1:N,1:N)::A
real::omega !relaxation factor
real::epsilon !convergence criterion
integer::iter !number of SOR iterations
integer::i,j,k,Nk

!set the initial matrix A and vectors x and b
call input(N,A,x,b)
!print initial conditions
write(*,*)'The initial conditions are:'
write(*,*)
call print_matrix(N,A,x,b)
write(*,*)
write(*,*)
!set SOR parameters
!set omega
omega=0.5
!set the convergence criterion
epsilon=1.e-8
!call the SOR scheme
call sor_rout(N,A,x,b,omega,epsilon,iter)
!print the final matrix A and vectors x and b
write(*,*)'The final conditions are:'
write(*,*)
call print_matrix(N,A,x,b)
write(*,*)
write(*,*)
write(*,*)"The number of iterations was: ",iter



end program sor_code
!*******************************
subroutine input(N0,A0,x0,b0)
implicit none
integer::N0
real,dimension(1:N0)::x0,b0
real,dimension(1:N0,1:N0)::A0

!define matrix A
A0(1,1)=4.0;A0(1,2)=-1.0;A0(1,3)=-6.0;A0(1,4)=0.0
A0(2,1)=-5.0;A0(2,2)=-4.0;A0(2,3)=10.0;A0(2,4)=8.0
A0(3,1)=0.0;A0(3,2)=9.0;A0(3,3)=4.0;A0(3,4)=-2.0
A0(4,1)=1.0;A0(4,2)=0.0;A0(4,3)=-7.0;A0(4,4)=5.0

!define vector b
b0(1)=2.0;b0(2)=21.0;b0(3)=-12.0;b0(4)=-6.0

!initial guess for x
x0=0.0

end subroutine input
!*********************************
subroutine print_matrix(N0,A0,x0,b0)
implicit none
integer::N0
real,dimension(1:N0)::x0,b0
real,dimension(1:N0,1:N0)::A0
integer::i,j

write(*,*)"matrix A:"
write(*,*)
do i=1,N0
   write(*,*)(A0(i,j), j=1,N0)
end do
write(*,*)
write(*,*)"vector b:"
write(*,*)
do i=1,N0
   write(*,*)b0(i)
end do
write(*,*)"vector x"
write(*,*)
do i=1,N0
   write(*,*)x0(i)
end do

end subroutine print_matrix
!*****************************************************
!SOR scheme based on the algorith provided by:
!https://en.wikipedia.org/wiki/Successive_over-relaxation
!
subroutine sor_rout(N0,A0,x0,b0,omega0,epsilon0,iter0)
implicit none
integer::N0
real,dimension(1:N0)::b0,x0,err_vec
real,dimension(1:N0,1:N0)::A0
real::omega0 !relaxation factor
real::epsilon0 !convergence criterion
real::err !residual error
real::sigma !auxiliaru variable
integer::i,j,k,Nk,iter0


!set an arbitrary initial error
err=1.0

iter0=0
do while(err>epsilon0)
   do i=1,N0
      sigma=0.0
      do j=1,N0
         if(j/=i)then
            sigma=sigma+A0(i,j)*x0(j)
         end if
      end do
      x0(i)=(1.0-omega0)*x0(i)+omega0*(b0(i)-sigma)/A0(i,i)
   end do
      err_vec=abs(Matmul(A0,x0)-b0)
      err=maxval(err_vec)
   iter0=iter0+1
end do


end subroutine sor_rout
