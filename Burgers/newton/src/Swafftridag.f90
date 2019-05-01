!!
!! Luke McCulloch
!! Incorporation of the trid subroutine into module form
!! September 2011
!! Rename to F90 file to support modules!
!! 
MODULE Swafftridag

  USE precise, only : defaultp
 ! Implicit None
 ! Integer  :: n
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp

Contains

      subroutine trid (a,b,c,r,u,nn)
      integer :: nn
      !Integer, parameter :: nn=n
      Integer, parameter :: nmax=500
      
!!$c
!!$c  routine to solve scalar tridiagonal matrics using the thomas
!!$c  algorithm.  a is the lower diagonal, b is the main diagonal,
!!$c  c is the upper diagonal, and r is the rhs.  solution returned in u.
!!$c  n is the number of elements in the main diagonal.
!!$c  Reference:  
!!$c     Press, William H., Blannery Brian P., Teukolsky, Saul A.,
!!$c     and Vetterling, William T., Numerical Recipes: The Art of Scientific
!!$c     Computing, Cambridge University Press, Cambridge, MA, 1986, p. 40.
!!$
!!$c   Note that:
!!$c   indices for the lower diagonal (a) should range from 2 -> n (contains n-1 elements)
!!$c   indices for the main diagonal (b) should range from 1 -> n (contains n elements)
!!$c   indices for the upper diagonal (c) should range from 1 -> n-1 (contains n-1 elements)
!!$c   indices for the rhs should range from 1 -> n (contains n elements)
!!$c   note that a(1) and c(n) are undefined and are not referenced by the routine.

      Real(wp), dimension(nn) :: gam,a,b,c,r,u

  

!!$    Do i=j,nn
!!$       write(*,*) 'rhs', r(j)
!!$    End do
      print*
      if(b(1).eq.0.0) stop  'stopping...need to modify matrix'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,nn
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j)*gam(j)
      if(bet.eq.0.0) stop 'stopping....zero on the diagonal'
      u(j)=(r(j)-a(j)*u(j-1))/bet
   11 continue
      do 12 j=nn-1,1,-1
      u(j)=u(j)-gam(j+1)*u(j+1)
   12 continue
      return

!!$    Do i=j,N
!!$       write(*,*) 'U(j) in the swaff inv', U(j)
!!$    End do

 End Subroutine trid

End Module SwaffTridag
