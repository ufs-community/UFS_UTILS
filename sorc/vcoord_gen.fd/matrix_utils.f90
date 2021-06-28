!> @file
!! Lower and upper triangular decomposition.
!! @author Mark Iredell @date 2008-08-01

!> This subprogram decomposes a matrix into a product of
!! lower and upper triangular matrices.
!!
!! @param[inout] a - input: a real(np,np) matrix (will be overwritten) output:
!!                 - output real(np,np) LU-decomposed matrix
!!                   (U is upper part of A, including diagonal;
!!                   L is lower part of A, with 1 as diagonal;
!!                   L*U equals original A after permuting)
!! @param[in] n integer order of matrix
!! @param[in] np integer dimension of matrix
!! @param[out] indx integer(n) pivot indices
!!              (original A rows are permuted in order i with indx(i))
!! @param[out] d real determinant permutation (1 or -1, or 0 if singular)
!!              (determinant is output diagonal product times d)
!! @author Mark Iredell @date 2008-08-01
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

!> Lower and upper triangular back substitution
!! @author Iredell @date 2008-08-01
!!
!! This subprogram back substitutes to solve decomposed
!! lower and upper triangular matrices as outputted by ludcmp.
!!
!! @param[in] a real(np,np) LU-decomposed matrix (from ludcmp)
!! @param[in] n integer order of matrix
!! @param[in] np integer dimension of matrix
!! @param[in] indx integer(n) pivot indices (from ludcmp)
!! @param[inout] b  - input real(n) rhs vector of linear problem (will be overwritten)
!!                  - output real(n) solution of linear problem 
!!                    (original A times output B equals original B)
!! @author Mark Iredell @date 2008-08-01
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
