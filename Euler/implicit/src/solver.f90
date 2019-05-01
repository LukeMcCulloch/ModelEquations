!!
!! Luke McCulloch
!! solver.f90
!! 11-9-2011
!! Qasi 1D Euler Solver



PROGRAM solver

  use precise, only : defaultp   ! Module to handle precision
  use GEO                        ! Module to write the Points on the surface
  use Swafftridag                ! Swaffs version

  IMPLICIT NONE  

  integer, parameter :: WP=defaultp
  real, parameter :: e=2.718281828459045
  INTEGER, DIMENSION(5), PARAMETER :: it1 = (/ 1, 5, 10, 15, 40 /)   ! GNU output times

  integer::count_0, count1, count2, count_rate, count_max, walltime, tstart, tend

  Integer :: NI, i, t, NT
  Real(WP) :: dx, L, dt, Umx, radius
  REAL(wp), Allocatable, DIMENSION(:) :: X, ratiodtdx  ! X locations
  REAL(wp), Allocatable, DIMENSION(:) :: am, bm, cm ! Diagonal Matrices
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: RHS, Delta, U, Uexact, Aplus, Aminus, L1, L2, fplus, fminus ! Velocity, scheme components

  ! Simple start up message
  WRITE(6,'(A)') ' CFD1 Quasi 1D Euler Solver'
  Write(6,*)     ' 11/9/2011 Implementation By Luke McCulloch'
  write(6,*) ''
  Write(6,*) ''
  write(6,*) '-----------------------------------------------------'
  write(6,*) '-----------------------------------------------------'

  NI=51
  dx=.2
  L=10.
  dt = 0.5 
  NT=40

  Allocate( X(NI) )
  Allocate( U(NI,NT) )
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
 

  !Subroutine to Get the Grid!
  Call D1(NI, dx, X ) 

  !Check Points - Comment this in the full Version
!  Write(6,'(AA)') 'Points'

!!$  do i=1,NI
!!$     !write(6,'(2D16.8)') X(i)
!!$     write(*,*) 'I, X(i)', i, X(i)
!!$  end do
!!$  Write(6,*) ''
  !Get rid of this in the full Version, obv

!  write(6,*) NI

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
     write(*,*) 'i, x(0), u(x,0)', i, X(i), U(i,0)
  End do
  write(*,*) '                  1.0 < x <= 3.0'
  count1=int(real(1.0/dx))+2
  count2=int(real(3.0/dx))+1
  Do i=count1,count2
     u(i,0)=0.5
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




  Write (*,*) '------------------Begin Implicit Method--------------------------'
  ! NI is the number of elements in the diagonal--------------------------------!
  ! am is the subdiagonal
  ! bm is the diagonal
  ! cm is the superdiagonal
  ! U(:,t-1) is the RHS
  ! U(:,t) is the Solution at time t
  ! t is indexed from 0 to 30.  NT=30

  !Initialize Delta
  Do i=1,NI-2 
     Delta(i,0)=U(i+1,0)
  End do
  
  !Begin the Time Loop!
  Do t=1,NT

  ! Compute the Subdiagonal, a, Diagonal, b, And Superdiagonal, c, Terms

     !Compute the sub!
     do i=2,NI-2
        am(i)=(-dt/(2.*dx))*Delta(i-1,t-1)
     End do

     !Compute the Super!
     do i=1,NI-3
        cm(i)=(dt/(2.*dx))*Delta(i+1,t-1)
     End do

     !!Computer the RHS and the Main Diagonal
     do i=1,NI-2
       !* remember to include the 1/2 's later*!
       Fplus(i,t-1)=((U(i+2,t-1))**2.)/2. 
       Fminus(i,t-1)=((U(i,t-1))**2.)/2.
       bm(i)=1.
    end do


    !! Call Swaff Solver
    Call trid (am,bm,cm,(-.5*(dt/dx)*(Fplus(:,t-1)-Fminus(:,t-1))),Delta(:,t),NI-2)
    
    !!Convert to U(:,t)
    Do i=1,NI-2
    Delta(i,t)=Delta(i,t-1)+Delta(i,t)
    U(i+1,t)=Delta(i,t)
    End Do
    !Delta(:,t)=Delta(:,t-1)
    U(1,t)=1.5
    U(NI,t)=0.5

    write(6,*) 'Timestep', T
    Do i=1,NI
       write(*,*) 'new U(i,t)', U(i,t) 
    End do    
 End do

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



!!----------------------------End Exact Solution----------------------------------------


  ! output Gnuplot file
  OPEN(UNIT=1, FILE='Gnu.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  !IF (ios /= 0) THEN
  !  PRINT *, "Error opening Gnu plot file "
  !  STOP
  !END IF
  Do t=1, 5
     DO i=1, NI
        !WRITE(UNIT=1,FMT='(E19.10," ",E19.10," 0.0")') x(i,it1(t)),u(i,it1(t))
        WRITE(UNIT=1,FMT='(E19.10," ",E19.10)') x(i),U(i,it1(t))
     END DO
     WRITE(UNIT=1,FMT='(" ")')
  End Do
  Do i=1,NI
     WRITE(UNIT=1,FMT='(E19.10," ",E19.10)') x(i)+dt*Nt*0.5*(1.5+.5),Uexact(i,40) 
  End do
  CLOSE(UNIT=1)




!! Plotting in gnu plot
!! once in the proper directory in terminal, open a second tab
!! type:  gnuplot   to open gnuplot
!! type:  plot "Gnu.dat" w l        to plot



  IF (ALLOCATED(X)) DEALLOCATE(X)
  IF (ALLOCATED(U)) DEALLOCATE(U) 
  IF (ALLOCATED(Uexact)) DEALLOCATE(Uexact) 
  IF (ALLOCATED(am)) DEALLOCATE(am) 
  IF (ALLOCATED(bm)) DEALLOCATE(bm) 
  IF (ALLOCATED(cm)) DEALLOCATE(cm)  

  IF (ALLOCATED(L1)) DEALLOCATE(L1) 
  IF (ALLOCATED(L2)) DEALLOCATE(L2) 
  IF (ALLOCATED(Aplus)) DEALLOCATE(Aplus)  
  IF (ALLOCATED(Aminus)) DEALLOCATE(Aminus) 

  IF (ALLOCATED(Fplus)) DEALLOCATE(Fplus) 
  IF (ALLOCATED(Fminus)) DEALLOCATE(Fminus)  

  IF (ALLOCATED(Delta)) DEALLOCATE(Delta) 
  IF (ALLOCATED(RHS)) DEALLOCATE(RHS) 

END PROGRAM solver


