!!
!! Luke McCulloch
!! solver.f90
!! 8-25-2011
!! Finite Volume Solver.
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
  !use inv                        ! My tri dag inversion module from Spring 2011
  use Swafftridag                ! Swaffs version
 ! use newtri

  IMPLICIT NONE  

  integer, parameter :: WP=defaultp
  real, parameter :: e=2.718281828459045
  INTEGER, DIMENSION(5), PARAMETER :: it1 = (/ 1, 5, 10, 15, 40 /)   ! GNU output times

  integer::count_0, count1, count2, count_rate, count_max, walltime, tstart, tend

  Integer :: NI, i, t, NT, NV
  Real(WP) :: dx, L, dt, Umx
  Real(WP), DIMENSION(2) :: ux           ! Unit Vectors
  REAL(wp), Allocatable, DIMENSION(:) :: X, ratiodtdx, V  ! X locations
  REAL(wp), Allocatable, DIMENSION(:) :: am, bm, cm ! Diagonal Matrices
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: RHS, Delta, vU, xU, Uexact, Aplus, Aminus, L1, L2, fplus, fminus ! Velocity, scheme components

  ! Simple start up message
  WRITE(6,'(A)') 'CFD Burgers Schemes Code'
  Write(6,*)     ' 9-27-2011 Implementation By Luke McCulloch'
  write(6,*) ''
  Write(6,*)     ' Soving the Inviscid Burgers Equation in 1D'
  write(6,*) '-----------------------------------------------------'


  Umx=-1.0e10
  dx=0.05
  L=3.
  NI=int(real(L/dx))   !! Gives 60 Spaces.  Add 1 to get Points on each side of the space
  NI=NI+1        !!Change to Finite Volue index.  Add 2 Boundary Grid points - 1 Phantom on each side
  NV=NI+1        !!Volume centroid index.  One on each side of a boundary grid point.
  dt = 0.026667  !!Swaffs Recommendation
  NT=40

  Allocate( X(NI) )
  Allocate( xU(NI,0:NT) )   !!U at each volumetric boundary
  Allocate( vU(NV,0:NT) )   !!U at each volumetric centroid
  Allocate( Uexact(NI,NT) )
  Allocate( am(NI) )
  Allocate( bm(NI) )
  Allocate( cm(NI) )

  Allocate( L1(NI,NT) )
  Allocate( L2(NI,NT) )
  Allocate( Aplus(NI,NT) )
  Allocate( Aminus(NI,NT) )

  Allocate( Fplus(NI,NT) )
  Allocate( Fminus(NI,NT) )

  Allocate( Delta(NI,NT) ) 
  Allocate( RHS(NI,NT) )
  
  Allocate( V(NV) )
 
  Call D1(NI, dx, X )   ! Subroutine to compute the vector of Points.
  Call V1(NV, dx, V)

  !Check Points - Comment this in the full Version
  Write(6,'(AA)') 'Points'

  do i=1,NI
     !write(6,'(2D16.8)') X(i)
     write(*,*) 'I, X(i)', i, X(i)
  end do
  Write(6,*) ''

  !!-------------------------------------

  Write(6,'(AA)') 'Volume Centroids'
  do i=1,NV
     !write(6,'(2D16.8)') X(i)
     write(*,*) 'I, V(i)', i, V(i)
  end do
  Write(6,*) ''
  !Get rid of this in the full Version, obv

!  write(6,*) NI

  !Establish the B.C.'s-------------------------------------------------
  Do i=1,NI
     xU(i,0)=0.0  !specifdy U at edges
  End Do
  Do i=1,NV
     vU(i,0)=0.0  !specify U at centers
  End do


  write(*,*) '----------------Initial Conditions on U at volume centroid--------------'
  write(*,*) '                  0. <= x <= 1.0'
  count1=1
  count2=int(real(1.0/dx))+1
  Do i=count1,count2
     vU(i,0)=1.5
     write(*,*) 'i, V(i,0), U(x,0)', i, V(i), vU(i,0)
  End do
  write(*,*) '                  1.0 < x <= 3.0'
  count1=int(real(1.0/dx))+2
  count2=int(real(3.0/dx))+2
  Do i=count1,count2
     vU(i,0)=0.5
     write(*,*) 'i, V(i,0), vU(x,0)', i, V(i), vU(i,0)
  End do
  write(*,*) ''  

  write(*,*) '----------------Boundary Conditions on U at volume face--------------'
  write(*,*) '                  0. <= x <= 1.0'
  count1=1
  count2=int(real(1.0/dx))+1
  Do i=count1,count2
     xU(i,0)=1.5
     write(*,*) 'i, X(i,0), xU(x,0)', i, X(i), xU(i,0)
  End do
  write(*,*) '                  1.0 < x <= 3.0'
  count1=int(real(1.0/dx))+2
  count2=int(real(3.0/dx))+1
  Do i=count1,count2
     xU(i,0)=0.5
     write(*,*) 'i, X(i,0), xU(x,0)', i, X(i), xU(i,0)
  End do
  write(*,*) ''  

!pause
 
  !End of the IVP B.C. portion------------------------------------------------



!!--------------Simple Upwind Difference (backward space difference) - Working
!   Write (*,*) '------------------Begin Simple Upwind Difference------------------------'
!   Write (*,*) '-------Program Setup to Find Solutions at the Volume Centroids----------'
!   Do t=1,NT
!      vU(1,t)=vU(2,t-1)     ! Note Phantoms now update
!      vU(NV,t)=vU(NV-1,t-1)
!      Write(*,*) 'Begin Time step t =', t
!      do i=2,NV-1   
!         vU(i,t)=vU(i,t-1)-0.5*dt*((vU(i,t-1)**2.)-(vU(i-1,t-1)**2.))/dx  !! Watch your signs
!      write(*,*) 'i, vol center loc,   U(x,t)   ',i, v(i), vU(i,t)
!      end do
!         !!Update your xU's - How do you update the Cell walls??

!      write(*,*) 'End timestep t'
!   end do

! !!--------------Lax-Friedrichs (Central Space Difference)
  
!   Write (*,*) '------------------Begin Lax-Friedrichs Method---------------------------'
!   Do t=1,NT
!      Write(*,*) 'Begin Time step t'
!      vU(1,t)=vU(2,t-1)
!      vU(NV,t)=vU(NV-1,t-1)
!      do i=2,NV-1
!         vU(i,t) =0.5*(vU(i+1,t-1)+vU(i-1,t-1))-0.25*(dt/dx)*((vU(i+1,t-1))**2.-(vU(i-1,t-1))**2.)
!         write(*,*) ' i, X(i)   u(x,t)   ', i, V(i), vU(i,t)
!      end do
!      write(*,*) 'End timestep t'
!   end do


!!--------------Lax-Wendroff (Central Space Difference)   ------MUST DEBUG THIS CODE!!*!
  
!   Write (*,*) '------------------Begin Lax-Wendroff Method---------------------------'
!   Do t=1,NT
!      Write(*,*) 'Begin Time step t'
!      vU(1,t)=vU(2,t-1)
!      vU(NV,t)=vU(NV-1,t-1)
!      do i=2,NV-2
!         L1(i,t-1)=0.25*(dt/dx)*((vU(i+1,t-1))**2.-(vU(i-1,t-1))**2.)
!         !write(*,*) 'L1(i,t-1)', L1(i,t-1)

!         Aplus(i,t-1)=0.5*( vU(i,t-1) + vU(i+1,t-1)   )
!         !write(*,*) 'Aplus', Aplus(i,t-1)

!         Aminus(i,t-1)=0.5*( vU(i,t-1) + vU(i-1,t-1)   )
!         !write(*,*) 'Aminus', Aminus(i,t-1)

!         Fplus(i,t-1)=(vU(i+1,t-1))**2.-(vU(i,t-1))**2.  !* remember to include the 1/2 's *!
!         Fminus(i,t-1)=(vU(i,t-1))**2.-(vU(i-1,t-1))**2.

!         L2(i,t-1)=0.25*((dt/dx)**2.)*((Aplus(i,t-1)*Fplus(i,t-1))-(Aminus(i,t-1)*Fminus(i,t-1)))
!         !write(*,*) 'L2', L2(i,t-1)
!         vU(i,t) =vU(i,t-1) - L1(i,t-1) + L2(i,t-1)
!         write(*,*) ' X(i)   u(x,t)   ', V(i), vU(i,t)
!      end do
!      write(*,*) 'End timestep t'
!     ! pause
!   end do


  Write (*,*) '------------------Begin Implicit Euler Method--------------------------'
  ! NI is the number of elements in the diagonal
  ! am is the subdiagonal
  ! bm is the diagonal
  ! cm is the superdiagonal
  ! U(:,t-1) is the RHS
  ! U(:,t) is the Solution at time t
  ! t is indexed from 0 to 30.  NT=30
  Do i=1,NI-2  !!! Make the Delta matrix which excludes the bc's
     Delta(i,0)=vU(i+1,0)
  End do
  
  Do t=1,NT
  ! Compute the Subdiagonal, a, Diagonal, b, And Superdiagonal, c, Terms
     !Delta(1,t-1)=1.5!U(1,t-2)
     !Delta(NI,t-1)=0.5!U(NI,t-2)

     do i=2,NI-2
        am(i)=(-dt/(2.*dx))*Delta(i-1,t-1)
     End do

     do i=1,NI-3
        cm(i)=(dt/(2.*dx))*Delta(i+1,t-1)
     End do


     do i=1,NI-2

       
       Fplus(i,t-1)=((vU(i+2,t-1))**2.)/2.  !* remember to include the 1/2 's later*!
       Fminus(i,t-1)=((vU(i,t-1))**2.)/2.


       bm(i)=1.

    end do

    !am(1)=0.
    !cm(NI)=0.
    !Fplus(NI,t-1)=0.
    !Fminus(1,t-1)=0.

    !! Solve for the delta U(:,t)

    !! Call my solver
    !Call tlmtri (am(:),bm(:),cm(:),-.5*(dt/dx)*(Fplus(:,t-1)-Fminus(:,t-1)),Delta(:,t),NI-2)

    !! Call Swaff Solver
    Call trid (am,bm,cm,(-.5*(dt/dx)*(Fplus(:,t-1)-Fminus(:,t-1))),Delta(:,t),NI-2)
    
    !!Convert to U(:,t)
    Do i=1,NI-2
       Delta(i,t)=Delta(i,t-1)+Delta(i,t)
       vU(i+1,t)=Delta(i,t)
    End Do
    !Delta(:,t)=Delta(:,t-1)
    vU(1,t)=1.5
    vU(NI,t)=0.5

    write(6,*) 'Timestep', T
    Do i=1,NI
       write(*,*) 'new U(i,t)', vU(i,t) 
    End do    
 End do

!!----------------------------End Euler Scheme----------------------------------------





!!----------------------------Exact Solution----------------------------------------

 !Do t=1,NT
 
 t=40*dt
    Do i=1,NI
       Uexact(i,40)=vU(i,0)
    End Do
 !End Do
 print*
 Write(6,*) 'exact solution at t=40'
 Do i=1,NI
    write(*,*) ' Uexact ', Uexact(i,40)
 End do



!!----------------------------End Exact Solution----------------------------------------


 Write(6,*) 'v, vU (velocity at the FV centroids)'
  ! output Gnuplot file
  OPEN(UNIT=1, FILE='LaxWendroff.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  !IF (ios /= 0) THEN
  !  PRINT *, "Error opening Gnu plot file "
  !  STOP
  !END IF
  Do t=1, 5
     DO i=2, NV-1
        !WRITE(UNIT=1,FMT='(E19.10," ",E19.10," 0.0")') x(i,it1(t)),u(i,it1(t))
        WRITE(UNIT=1,FMT='(E19.10," ",E19.10)') v(i),vU(i,it1(t))
     END DO
     WRITE(UNIT=1,FMT='(" ")')
  End Do

  
 Write(6,*) 'v, U exact'
  Do i=1,NI
     WRITE(UNIT=1,FMT='(E19.10," ",E19.10)') v(i)+dt*Nt*0.5*(1.5+.5),Uexact(i,40) 
  End do
  CLOSE(UNIT=1)




!! Plotting in gnu plot
!! once in the proper directory in terminal, open a second tab
!! type:  gnuplot   to open gnuplot
!! type:  plot "Gnu.dat" w l        to plot






  Write(6,*) 'deallocating'


!   IF (ALLOCATED(Points)) DEALLOCATE(Points)
!   IF (ALLOCATED(S)) DEALLOCATE(S)
!   IF (ALLOCATED(T)) DEALLOCATE(T) 
!   IF (ALLOCATED(To)) DEALLOCATE(To) 
!   IF (ALLOCATED(Treal)) DEALLOCATE(Treal) 
  
  Write(6,*) 'deallocating x'
  IF (ALLOCATED(X)) DEALLOCATE(X)
  Write(6,*) 'deallocating xU'
  IF (ALLOCATED(xU)) DEALLOCATE(xU) 
  Write(6,*) 'deallocating vU'
  IF (ALLOCATED(vU)) DEALLOCATE(vU) 
  Write(6,*) 'deallocating Uexact'
  IF (ALLOCATED(Uexact)) DEALLOCATE(Uexact) 
  Write(6,*) 'deallocating a,b,c'
  IF (ALLOCATED(am)) DEALLOCATE(am) 
  IF (ALLOCATED(bm)) DEALLOCATE(bm) 
  IF (ALLOCATED(cm)) DEALLOCATE(cm)  

  Write(6,*) 'deallocating L1'
  IF (ALLOCATED(L1)) DEALLOCATE(L1) 
  IF (ALLOCATED(L2)) DEALLOCATE(L2) 
  Write(6,*) 'deallocating Aplus'
  IF (ALLOCATED(Aplus)) DEALLOCATE(Aplus)  
  IF (ALLOCATED(Aminus)) DEALLOCATE(Aminus) 

  Write(6,*) 'deallocating Fplus'
  IF (ALLOCATED(Fplus)) DEALLOCATE(Fplus) 
  IF (ALLOCATED(Fminus)) DEALLOCATE(Fminus)  

  
  Write(6,*) 'deallocating Delta'
  IF (ALLOCATED(Delta)) DEALLOCATE(Delta) 
  IF (ALLOCATED(RHS)) DEALLOCATE(RHS) 

  
  Write(6,*) 'deallocating V'
  IF (ALLOCATED(V)) DEALLOCATE(V)

  Write(6,*) 'End Program'
END PROGRAM solver


