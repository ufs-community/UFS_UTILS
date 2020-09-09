!-------------------------------------------------------------------------------
!$$$  Subprogram documentation block
!
! Subprogram:    ludcmp      lower and upper triangular decomposition
!   Prgmmr: Iredell    Org: W/NP23      Date: 2008-08-01
!
! Abstract: This subprogram decomposes a matrix into a product of
!   lower and upper triangular matrices.
!
! Program history log:
!   2008-08-01  Mark Iredell
!
! Usage:    call ludcmp(a,n,np,indx,d)
!   Input argument list:
!     a        real(np,np) matrix (will be overwritten)
!     n        integer order of matrix
!     np       integer dimension of matrix
!
!   Output argument list:
!     a        real(np,np) LU-decomposed matrix
!              (U is upper part of A, including diagonal;
!               L is lower part of A, with 1 as diagonal;
!               L*U equals original A after permuting)
!     indx     integer(n) pivot indices
!              (original A rows are permuted in order i with indx(i))
!     d        real determinant permutation (1 or -1, or 0 if singular)
!              (determinant is output diagonal product times d)
!
! Attributes:
!   Language: Fortran 90
!
!$$$
subroutine ludcmp(a,n,np,indx,d)
  implicit none
  integer,intent(in):: n,np
  real,intent(inout):: a(np,np)
  integer,intent(out):: indx(n)
  real,intent(out):: d
  integer i,j,k,imax
  real aamax,sum,dum
  real vv(n)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  d=1
  do i=1,n
    aamax=0
    do j=1,n
      if(abs(a(i,j))>aamax) aamax=abs(a(i,j))
    enddo
    if(aamax==0) then
      d=0
      return
    endif
    vv(i)=1/aamax
  enddo
  do j=1,n
    do i=1,j-1
      sum=a(i,j)
      do k=1,i-1
        sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
    enddo
    aamax=0.
    do i=j,n
      sum=a(i,j)
      do k=1,j-1
        sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
      dum=vv(i)*abs(sum)
      if(dum>=aamax) then
        imax=i
        aamax=dum
      endif
    enddo
    if (j/=imax)then
      do k=1,n
        dum=a(imax,k)
        a(imax,k)=a(j,k)
        a(j,k)=dum
      enddo
      d=-d
      vv(imax)=vv(j)
    endif
    indx(j)=imax
    if(a(j,j)==0) then
      d=0
      return
    endif
    if(j/=n)then
      dum=1/a(j,j)
      do i=j+1,n
        a(i,j)=a(i,j)*dum
      enddo
    endif
  enddo
end subroutine
!-------------------------------------------------------------------------------
!$$$  Subprogram documentation block
!
! Subprogram:    lubksb      lower and upper triangular back substitution
!   Prgmmr: Iredell    Org: W/NP23      Date: 2008-08-01
!
! Abstract: This subprogram back substitutes to solve decomposed
!   lower and upper triangular matrices as outputted by ludcmp.
!
! Program history log:
!   2008-08-01  Mark Iredell
!
! Usage:    call lubksb(a,n,np,indx,b)
!   Input argument list:
!     a        real(np,np) LU-decomposed matrix
!              (from ludcmp)
!     n        integer order of matrix
!     np       integer dimension of matrix
!     indx     integer(n) pivot indices
!              (from ludcmp)
!     b        real(n) rhs vector of linear problem (will be overwritten)
!
!   Output argument list:
!     b        real(n) solution of linear problem 
!              (original A times output B equals original B)
!
! Attributes:
!   Language: Fortran 90
!
!$$$
subroutine lubksb(a,n,np,indx,b)
  implicit none
  integer,intent(in):: n,np
  real,intent(in):: a(np,np)
  integer,intent(in):: indx(n)
  real,intent(inout):: b(n)
  integer i,j,ii,ll
  real sum
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ii=0
  do i=1,n
    ll=indx(i)
    sum=b(ll)
    b(ll)=b(i)
    if (ii/=0)then
      do j=ii,i-1
        sum=sum-a(i,j)*b(j)
      enddo
    elseif(sum/=0) then
      ii=i
    endif
    b(i)=sum
  enddo
  do i=n,1,-1
    sum=b(i)
    do j=i+1,n
      sum=sum-a(i,j)*b(j)
    enddo
    b(i)=sum/a(i,i)
  enddo
end subroutine
