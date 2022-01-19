!code name: 'combinatorial_function.f90'
!
!description: This code provides some routines to calculate the
!             factorial and combinatory functions.
!
!date: January 18th, 2022
!author: Leonardo Brito
!Ph.D. Physics Student
!Institute of Physics (IF)
!University of SÆo Paulo (USP)
!SÆo Paulo-SP, Brazil
!_______________________________________________________________________
program main
implicit none
integer::i,j,fact,comb

!calculate some factorials
write(*,*)'5! =',fact(5)
write(*,*)'10! =',fact(10)

!calculate some combinatorial functions
write(*,*)'C(10,2) =',comb(10,2)
write(*,*)'C(0,0) =',comb(0,0)
write(*,*)'C(1,1) =',comb(1,1)
write(*,*)'C(1,2) =',comb(1,2)
write(*,*)'C(2,1) =',comb(2,1)
write(*,*)'C(50,45) =',comb(50,45)

read*

end program main
!***********************************************************************
function comb(s,m)
!reference: https://en.wikipedia.org/wiki/Combination
!parameters:
!     input,integer: s !top number
!     input,integer: m !low number
!     output,integer: comb !combinatory function result
implicit none
integer::comb
integer::s,m!,comb2
integer::i,p
integer::N_2,fact
real::Combr,aux,aux0


!calculate the number of configurations Nc
! case s>m
if((s>=0).and.(m>=0).and.(s>m))then
   aux0=s
   aux=aux0
   if((s-m)<m)then
      N_2=s-m
      p=s-m-1
   else
      N_2=m
      p=m-1
   end if
   !combination
   do i=1,p
      aux=aux*(aux0-1.0*i)
   end do
   Combr=aux/fact(N_2)
end if

!especial cases
!case m=s
if((s>=0).and.(m>=0).and.(m==s))then
   Combr=1.0
end if

!case s<m (prohibited)
if(s<0)then
   Combr=0.0
end if

if(s<m)then
   Combr=0.0
end if

!case m=0
if(m==0)then
   Combr=1.0
end if
!case m<0
if(m<0)then
   Combr=0.0
end if
!case s=1
if((s==1).and.(m/=1))then
   Combr=0.0
end if
!case s=0
if(s==0)then
   Combr=0.0
end if
!case s=m=0
if((s==0).and.(m==0))then
   Combr=1.0
end if


 comb=int(Combr)

end function comb
!***********************************************************************
!factorial function
!            fact(k)=k!
function fact(k)
!reference: https://en.wikipedia.org/wiki/Factorial
!parameters:
!     input,integer: k !desired factorial number
!     output,integer: fact !factorial function result
implicit none
integer::i,k,fact

fact=1
do i=1,k
   fact=i*fact
end do

end function fact
!***********************************************************************
