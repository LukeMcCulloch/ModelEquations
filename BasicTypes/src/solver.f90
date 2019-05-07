!!
!! Luke McCulloch
!! solver.f90
!! 8-25-2011
!! Serial Solution Various equations - This will be Inviscid Burgers.
!!
!!   Variable Listing
!! Comp   - Compare Acutal to Computed Solutions
!! Compo  - Old Comparison (used to max Comp as a scalar)
!! diffo  - intial difference between T(i,j)new and T(i,j)old
!! diff   - latest difference between T(i,j)new and T(i,j)old
!! lx     - domain space in x
!! ly     - domain space in y
!! dx     - element size in x
!! dy     - element size in y
!! NI     - # points in x
!! NJ     - # points in y
!! n      - dummy vbl for points (mabye points is more machinery than needed)
!! nipt   - Total number of points on the plate = NI*NJ
!! K      - Iteration counter  ????Never used????
!! i      - ith position in x
!! j      - jth position in y
!! ipt    - dummy points index "the ipt'th position of nipt total"
!! C      - Cnst term in the discrete second differences
!! mn     - max norm of all "diffs"  //T(i,j)-To(i,j)//, and max norm of the solution variance (loosley speaking)
!! ux, uy - Unit vectors
!! X, Y   - Actual Geometry Position maps of i,j (Scaler values)
!! Points - Vector Valued (X(i,j),Y(i,j))  - Again more machinery than needed
!! T      - Temperature at T(i,j) at K iterations - i.e. the solution at the end
!! S      - Given Source F(i,j)->F(x,y)  S(x,y) = xe^y
!! To     - T(i,j) at K-1 iterations ( i.e. before update )
!! Treal  - Mathematical solution to the problem
!! w      - Relaxation factor {: (0,2)  <1 is slower but more convergent.  >1 is faster but less convergent...



PROGRAM solver

  use precise, only : defaultp   ! Module to handle precision
  use GEO                        ! Module to write the Points on the surface
  !!use Source                     ! Module to compute an arbitrary source term.
  use inv                        ! My tri dag inversion module from Spring 2011
  use Swafftridag                ! Swaffs version

  IMPLICIT NONE  

  integer, parameter :: WP=defaultp
  real, parameter :: e=2.718281828459045
  INTEGER, DIMENSION(5), PARAMETER :: it1 = (/ 1, 5, 10, 15, 40 /)   ! GNU output times

  integer::count_0, count1, count2, count_rate, count_max, walltime, tstart, tend

  Integer :: NI, i, t, NT, m, mmax
  Real(WP) :: dx, L, dt, Umx
  Real(WP), DIMENSION(2) :: ux           ! Unit Vectors
  REAL(wp), Allocatable, DIMENSION(:) :: X, ratiodtdx  ! X locations
  REAL(wp), Allocatable, DIMENSION(:) :: am, bm, cm ! Diagonal Matrices
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: RHS, Delta, U, Uexact, Aplus, Aminus, L1, L2, fplus, fminus, F ! Velocity, scheme components
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: Unm  

  ! Simple start up message
  WRITE(6,'(A)') 'CFD Burgers Schemes Code'
  Write(6,*)     ' 10-09-2011 Implementation By Luke McCulloch'
  write(6,*) ''
  Write(6,*)     ' Implicit Solver FV Inviscid Burgers Equation in 1D'
  write(6,*) '--------------------------------------------------------'


  Umx=-1.0e10
  dx=0.05
  L=3.
  NI=int(real(L/dx))
  NI=NI+1
  dt = 0.03333333  !!Swaffs Recommendation-may try others to see what it does.
  NT=40

  Allocate( X(NI) )
  Allocate( U(NI,0:NT) )
  Allocate( Uexact(NI,0:NT) )
  Allocate( am(NI-3) )
  Allocate( bm(NI-2) )
  Allocate( cm(NI-3) )

  Allocate( L1(NI,0:NT) )
  Allocate( L2(NI,0:NT) )
  Allocate( Aplus(NI,0:NT) )
  Allocate( Aminus(NI,0:NT) )

  Allocate( Fplus(NI-2,0:NT) )
  Allocate( Fminus(NI-2,0:NT) )
  Allocate( F(NI-2,0:NT) )

  Allocate( Delta(NI-2,0:NT) ) 
  Allocate( RHS(NI,0:NT) )
  Allocate( Unm(NI,0:NT) )

  write(*,*) 'memory allocation complete'

 
  write(*,*) 'compute points on a line'
  Call D1(NI, dx, X )   ! Subroutine to compute the vector of Points.

  !Check Points - Comment this in the full Version
!  Write(6,'(AA)') 'Points'

!!$  do i=1,NI
!!$     !write(6,'(2D16.8)') X(i)
!!$     write(*,*) 'I, X(i)', i, X(i)
!!$  end do
!!$  Write(6,*) ''
  !Get rid of this in the full Version, obv

!  write(6,*) NI

  write(*,*) 'establish B.C.s'
  !Establish the B.C.'s-------------------------------------------------
  Do i=1,NI
     U(i,0)=0.0  !Going to need to specifdy U(0.,t)=0.0, U(3.,t)=0.0
     !Delta(i,0)=0.0
  End Do
  !Do i=1,NI
     !Going to need to specifdy U(0.,t)=0.0, U(3.,t)=0.0
   !  Delta(i,0)=0.0
  !End Do

  write(*,*) '----------------Boundary Conditions--------------'
  write(*,*) '                  0. <= x <= 1.0'
  count1=1
  count2=int(real(1.0/dx))+1
  Do i=count1,count2
     u(i,0)=1.5
     Unm(i,0)=1.5
     write(*,*) 'i, x(0), u(x,0)', i, X(i), U(i,0)
  End do
  write(*,*) '                  1.0 < x <= 3.0'
  count1=int(real(1.0/dx))+2
  count2=int(real(3.0/dx))+1
  Do i=count1,count2
     u(i,0)=0.5
     Unm(i,0)=0.5
     write(*,*) 'i, x(0), u(x,0)', i, X(i), U(i,0)
  End do
  write(*,*) ''  
 
  !End of the B.C. portion------------------------------------------------

  !!c=1.0  !! Non constant



  !! Subroutine to establish the initial equation terms
!!$  do i=1, NI
!!$  ratiodtdx= 0.8/max(Umx, U(i,:))
!!$  end do
!!$  write(*,*) 'ratiodtdx', ratiodtdx
 !! Umx=1.5


!!--------------Simple Upwind Difference (backward space difference) - Working
!!$  Write (*,*) '------------------Begin Simple Upwind Difference------------------------'
!!$  Do t=0,NT
!!$     Write(*,*) 'Begin Time step t =', t
!!$     do i=2,NI
!!$        U(1,t+1)=U(1,t)
!!$        U(i,t+1)=U(i,t)-0.5*dt*((U(i,t)**2.)-(U(i-1,t)**2.))/dx  !! Watch your signs
!!$     write(*,*) '   u(x,t)   ',i, x(i), U(i,t+1)
!!$     end do
!!$     write(*,*) 'End timestep t'
!!$  end do

!!--------------Lax-Friedrichs (Central Space Difference)
  
!!$  Write (*,*) '------------------Begin Lax-Friedrichs Method---------------------------'
!!$  Do t=1,NT
!!$     Write(*,*) 'Begin Time step t'
!!$     do i=2,NI-1
!!$        U(1,t)=U(1,t-1)
!!$        U(NI,t)=U(NI,t-1)
!!$        U(i,t) =0.5*(U(i+1,t-1)+U(i-1,t-1))-0.25*(dt/dx)*((U(i+1,t-1))**2.-(U(i-1,t-1))**2.)
!!$        write(*,*) ' i, X(i)   u(x,t)   ', i, X(i), U(i,t)
!!$     end do
!!$     write(*,*) 'End timestep t'
!!$  end do


!!--------------Lax-Wendroff (Central Space Difference) 
  
!!$  Write (*,*) '------------------Begin Lax-Wendroff Method---------------------------'
!!$  Do t=1,NT
!!$     Write(*,*) 'Begin Time step t'
!!$     do i=2,NI-1
!!$        U(1,t)=U(1,t-1)
!!$        U(NI,t)=U(NI,t-1)
!!$
!!$        L1(i,t-1)=0.25*(dt/dx)*((U(i+1,t-1))**2.-(U(i-1,t-1))**2.)
!!$        !write(*,*) 'L1(i,t-1)', L1(i,t-1)
!!$
!!$        Aplus(i,t-1)=0.5*( U(i,t-1) + U(i+1,t-1)   )
!!$        !write(*,*) 'Aplus', Aplus(i,t-1)
!!$
!!$        Aminus(i,t-1)=0.5*( U(i,t-1) + U(i-1,t-1)   )
!!$        !write(*,*) 'Aminus', Aminus(i,t-1)
!!$
!!$        Fplus(i,t-1)=(U(i+1,t-1))**2.-(U(i,t-1))**2.  !* remember to include the 1/2 's *!
!!$        Fminus(i,t-1)=(U(i,t-1))**2.-(U(i-1,t-1))**2.
!!$
!!$        L2(i,t-1)=0.25*((dt/dx)**2.)*((Aplus(i,t-1)*Fplus(i,t-1))-(Aminus(i,t-1)*Fminus(i,t-1)))
!!$        write(*,*) 'L2', L2(i,t-1)
!!$        U(i,t) =U(i,t-1) - L1(i,t-1) + L2(i,t-1)
!!$        write(*,*) ' X(i)   u(x,t)   ', X(i), U(i,t)
!!$     end do
!!$     write(*,*) 'End timestep t'
!!$    ! pause
!!$  end do


  Write (*,*) '------------------Begin Implicit Euler Method--------------------------'
  ! NI is the number of elements in the diagonal
  ! am is the subdiagonal
  ! bm is the diagonal
  ! cm is the superdiagonal
  ! F(:,t-1) is the RHS
  ! delta(:,t) is the Solution at time t
  ! t is indexed from 0 to 30.  NT=30



  Do t=1,NT

  Do i=1,NI-2  !!! Make the Delta matrix which excludes the bc's
     Delta(i,t-1)=0.0
  End do
  



     ! Compute the Subdiagonal, a, Diagonal, b, And Superdiagonal, c, Terms
     !Delta(1,t-1)=1.5!U(1,t-2)
     !Delta(NI,t-1)=0.5!U(NI,t-2)
     Do i=2,NI
        write(*,*) 'i, t, Unm(i,t-1)', i, t, Unm(i,t-1)
        write(*,*) 'i,t,Unm(i-1,t-1)', i, t, Unm(i-1,t-1)
        write(*,*) ''  
     End do
     
     !!Do the Newton Iterations from 1-10
     mmax=10
     do m=1,mmax
        write(*,*)  'Begin m loop'
        !pause!----------------------------------------
        do i=1,NI-3
           !!am(i)=(-dt/(2.*dx))*Delta(i-1,t-1)
           am(i)=(dt/dx)*Unm(i+2,t-1)   !! Range = 3,NI-1
        End do
     
        do i=2,NI-2
           !!cm(i)=(dt/(2.*dx))*Delta(i+1,t-1)
           cm(i)=0.0    !! Range = 2,NI-2
        End do



        
        do i=1,NI-2               
           Fplus(i,t-1) =((Unm(i+1,t-1))**2.)  !! Range = 2,NI-1
           Fminus(i,t-1)=((Unm(i,t-1))**2.)    !! Range = 1,NI-2
           F(i,t-1)=Unm(i+1,t-1)-U(i+1,t-1)+0.5*dt*(Fplus(i,t-1)-Fminus(i,t-1))/dx
           
           bm(i)=1.0+(dt/dx)*Unm(i+1,t-1)   !! Range = 2,NI-1
        
        end do
        
!!$        do i=1,NI-2 
!!$
!!$        end do



      !   write(*,*) 'a', am(:)
      !   write(*,*) ''
      !   write(*,*) 'b', bm(:)
      !   write(*,*) ''
      !   write(*,*) 'c', cm(:)
      !   write(*,*) ''    
        !! Solve for the delta U(:,t)
        
        !! Call my solver
        !Call tlmtri (am(:),bm(:),cm(:),-.5*(dt/dx)*(Fplus(:,t-1)-Fminus(:,t-1)),Delta(:,t),NI-2)
        
        !! Call Swaff Solver
        Call trid (am,bm,cm,-F(:,t-1),Delta(:,t),NI-2)

        !Do i=1,NI-2
        !   write(*,*) 'am, bm, cm, Fplus, Fminus, F', am(i), bm(i), cm(i), Fplus(i,t-1), Fminus(i,t-1), F(i,t-1)
        !End Do

        
        Do i=1,NI-2
           Unm(i+1,t-1)=Delta(i,t)+Unm(i+1,t-1)  !!Caareful, Unm range = [1,NI]
           !write(*,*) 'Delta(i,t,m+1), Delta(i,t,m), Unm(i,t-1)', Delta(i,t), Delta(i,t-1), Unm(i,t-1)
        End Do

!!$        !!Convert to U(:,t-1)
!!$        Do i=1,NI-2
!!$           U(i+1,t)=Delta(i,t)+U(i+1,t-1)
!!$        End Do
        
     End do  !! End the iterative m loop for implicit NM.
     
     
     !!Convert to U(:,t)
     Do i=1,NI-2
        U(i+1,t)=Unm(i+1,t-1) !!Updating inner range = [2,NI-1] could try whole domain but this should work...
        Unm(i+1,t)=Unm(i+1,t-1)  
     End Do
     
!! Should be able to include these immediately above since they never change...
     U(1,t)=1.5
     U(NI,t)=0.5
     Unm(1,t)=1.5
     Unm(NI,t)=0.5    
     
     write(6,*) 'Timestep', T
     Do i=1,NI
        write(*,*) 'new U(i,t)', U(i,t) 
     End do
  End do    !! End Timestep

!!----------------------------End Euler Scheme----------------------------------------





!!----------------------------Exact Solution----------------------------------------

 !Do t=1,NT
 
 t=40*dt
    Do i=1,NI
       Uexact(i,40)=U(i,0)
    End Do
 !End Do
 print*
 Write(6,*) 'exact solution at t=40'
 Do i=1,NI
    write(*,*) ' Uexact ', Uexact(i,40)
 End do



 Write(6,*) 'end exact solution at t=40'
!!----------------------------End Exact Solution----------------------------------------

  Write(6,*) 'output GNUplot file'
  ! output Gnuplot file
  OPEN(UNIT=1, FILE='EulerImplicit.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  !IF (ios /= 0) THEN
  !  PRINT *, "Error opening Gnu plot file "
  !  STOP
  !END IF

  Write(6,*) 'output x,U'
  Do t=1, 5
     DO i=1, NI
        !WRITE(UNIT=1,FMT='(E19.10," ",E19.10," 0.0")') x(i,it1(t)),u(i,it1(t))
        WRITE(UNIT=1,FMT='(E19.10," ",E19.10)') x(i),U(i,it1(t))
     END DO
     WRITE(UNIT=1,FMT='(" ")')
  End Do

  Write(6,*) 'output x,Uexact'
  Do i=1,NI
     WRITE(UNIT=1,FMT='(E19.10," ",E19.10)') x(i)+dt*Nt*0.5*(1.5+.5),Uexact(i,40) 
  End do
  CLOSE(UNIT=1)




!! Plotting in gnu plot
!! once in the proper directory in terminal, open a second tab
!! type:  gnuplot   to open gnuplot
!! type:  plot "Gnu.dat" w l        to plot







  Write(6,*) 'deallocating'

!!$  IF (ALLOCATED(Points)) DEALLOCATE(Points)
!!$  IF (ALLOCATED(S)) DEALLOCATE(S)
!!$  IF (ALLOCATED(T)) DEALLOCATE(T) 
!!$  IF (ALLOCATED(To)) DEALLOCATE(To) 
!!$  IF (ALLOCATED(Treal)) DEALLOCATE(Treal) 
  
  Write(6,*) 'deallocating 1'
  IF (ALLOCATED(X)) DEALLOCATE(X)
  IF (ALLOCATED(U)) DEALLOCATE(U) 
  IF (ALLOCATED(Uexact)) DEALLOCATE(Uexact) 
  IF (ALLOCATED(am)) DEALLOCATE(am) 
  IF (ALLOCATED(bm)) DEALLOCATE(bm) 
  IF (ALLOCATED(cm)) DEALLOCATE(cm)  

  Write(6,*) 'deallocating 2'
  IF (ALLOCATED(L1)) DEALLOCATE(L1) 
  IF (ALLOCATED(L2)) DEALLOCATE(L2) 
  IF (ALLOCATED(Aplus)) DEALLOCATE(Aplus)  
  IF (ALLOCATED(Aminus)) DEALLOCATE(Aminus) 


  Write(6,*) 'deallocating 4'
  IF (ALLOCATED(Delta)) DEALLOCATE(Delta) 
  IF (ALLOCATED(RHS)) DEALLOCATE(RHS) 
  if (allocated(Unm)) deallocate(Unm)

  Write(6,*) 'deallocating Fplus'
  IF (ALLOCATED(Fplus)) DEALLOCATE(Fplus) 
  Write(6,*) 'deallocating Fminus'
  IF (ALLOCATED(Fminus)) DEALLOCATE(Fminus)  
  Write(6,*) 'deallocating F'
  IF (ALLOCATED(F)) DEALLOCATE(F) 


END PROGRAM solver


