!!
!! Luke McCulloch
!! Incorporation of block tri diag solver to module form
!! Nov. 2011
!! rename to .F90 to support modules   
!! I did not code this!  goto? not me   

module blocktri

use precise, only : defaultp

integer, parameter, private :: WP=defaultp

contains

!* sub, main, super,  *!
subroutine solve (e,d,f,b,min,max)
!     parameter(ni=51,nj=11,ni1=ni+1,nj1=nj+1)
!     parameter (neq=3,no=neq*(nj-1))
      parameter (nq=3,ni=52)
      parameter (no=3)
!     parameter (neq=3,no=neq*nit)
!
!   This subroutine solves block tridiagonal systems where the
!   e's are the lower diagonal blocks, d's are the diagonal blocks,
!   f's are the upper diagonal blocks, and b is the right-hand-side
!   vector.  The solution is returned in the vector b.  The nomenclature
!   followed is that on page 111 of Golub and van Loan.  However, to
!   save storage, their L matrices are written back in e, and their U
!   matrices are written back in d.  At this point in time, the blocks
!   are square, are order no, and there are ((max+1)-min) of them.
!   Go til ya hear glass!
!
!    Swaff's addendum to the above:
!
!
!        nq is the size of each block (i.e., nq x nq)
!
!        d=diagonal matrices (NOTE: the blocks, and I mean ALL
!                            blocks, are dimensioned "i,j,k", where
!                            i=row, j=column, and k=position in
!                            overall block matrix -- see notes below)
!
!        f=upper diagonal matrices (start indexing them at i,j,1)
!          (i.e., if the number of diagonal blocks is "max", then the
!          upper diagonal blocks are indexed from 1 -> max-1)
!
!        e=lower diagonal matrices (start indexing them at i,j,2)
!          (i.e., if the number of diagonal blocks is "max", then the
!          lower diagonal blocks are indexed from 2 -> max)
!
!        Note that all blocks (and I mean ALL blocks) are dimensioned
!        as (nq x nq x max).  In this example, max=ni-2.
!
!        b=rhs vectors
!
      Real(wp), dimension(no,no,ni-2) :: e,d,f
      Real(wp), dimension(no,ni-2) :: b
      !dimension e(no,no,ni-2),d(no,no,ni-2),f(no,no,ni-2),b(no,ni-2)
      dimension dt(no,no),y(no),z(no),yy(no,0:ni-2),xx(no,ni)
      dimension mm(ni-2),nn(ni-2)
!
      do 4 i=min,max
      mm(i)=no
    4 nn(i)=no
!
      do 80 i=min+1,max
      do 75 m=1,nn(i)
      do 75 n=1,nn(i)
   75 dt(m,n)=d(n,m,i-1)
      do 10 m=1,nn(i)
      if (m.eq.1) go to 60
      do 20 n=m,nn(i)
      sum=0.
      do 30 ip=1,m-1
   30 sum=sum+dt(m,ip)*dt(ip,n)
   20 dt(m,n)=dt(m,n)-sum
   60 continue
      do 10 n=m+1,nn(i)
      sum=0.
      do 50 ip=1,m-1
   50 sum=sum+dt(n,ip)*dt(ip,m)
   10 dt(n,m)=(dt(n,m)-sum)/dt(m,m)
      do 80 irow=1,nn(i)
      do 81 icol=1,nn(i)
   81 z(icol)=e(irow,icol,i)
      y(1)=z(1)
      do 100 m=2,mm(i)
      sum=0.
      do 101 n=1,m-1
  101 sum=sum+dt(m,n)*y(n)
  100 y(m)=z(m)-sum
      z(nn(i))=y(nn(i))/dt(nn(i),nn(i))
      do 102 m=mm(i)-1,1,-1
      sum=0.
      do 103 n=nn(i),m+1,-1
  103 sum=sum+dt(m,n)*z(n)
  102 z(m)=(y(m)-sum)/dt(m,m)
      do 82 icol=1,nn(i)
   82 e(irow,icol,i)=z(icol)
      do 80 icol=1,nn(i)
      sum=0.
      do 83 j=1,nn(i)
   83 sum=sum+e(irow,j,i)*f(j,icol,i-1)
   80 d(irow,icol,i)=d(irow,icol,i)-sum
      do 301 i=min,max
      do 301 j=1,nn(i)
      sum=0.
      do 300 k=1,nn(i)
  300 sum=sum+e(j,k,i)*yy(k,i-1)
  301 yy(j,i)=b(j,i)-sum
      call aeqlu (d)
      do 210 i=max,min,-1
      do 303 j=1,nn(i)
      sum=0.
      do 302 k=1,nn(i)
  302 sum=sum+f(j,k,i)*xx(k,i+1)
  303 z(j)=yy(j,i)-sum
      y(1)=z(1)
      do 220 m=2,mm(i)
      sum=0.
      do 221 n=1,m-1
  221 sum=sum+d(m,n,i)*y(n)
  220 y(m)=z(m)-sum
      z(nn(i))=y(nn(i))/d(nn(i),nn(i),i)
      do 222 m=mm(i)-1,1,-1
      sum=0.
      do 223 n=nn(i),m+1,-1
  223 sum=sum+d(m,n,i)*z(n)
  222 z(m)=(y(m)-sum)/d(m,m,i)
      do 210 j=1,nn(i)
      b(j,i)=z(j)
  210 xx(j,i)=z(j)
      return
      end subroutine solve

!---------------------------------------------------------------------------
      subroutine aeqlu (a)
!     PARAMETER(NI=21,NJ=11,NI1=NI+1,NJ1=NJ+1)
!     parameter (neq=3,no=neq*(nj-1))
      parameter (nq=3,no=nq,ni=52)
!   This subroutine performs a Doolittle decomposition of the ni-1
!   noxno matrices brought-in in a, and then stores the results in a.
!   How do ya like all of these good comments?
      !dimension a(no,no,ni-2)
      Real(wp), dimension(no,no,ni-2) :: a
      do 10 i=1,ni-2
      do 10 m=1,no
      if (m.eq.1) go to 60
      do 20 n=m,no
      sum=0.
      do 30 ip=1,m-1
   30 sum=sum+a(m,ip,i)*a(ip,n,i)
   20 a(m,n,i)=a(m,n,i)-sum
   60 continue
      do 10 n=m+1,no
      sum=0.
      do 50 ip=1,m-1
   50 sum=sum+a(n,ip,i)*a(ip,m,i)
   10 a(n,m,i)=(a(n,m,i)-sum)/a(m,m,i)
      return
      end subroutine aeqlu

End Module blocktri
