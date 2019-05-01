!!
!! Luke McCulloch
!! solver.f90
!! 8-25-2011
!! Serial Solution Various equations - This will be Inviscid Burgers.
!! Euler Solution Contains the Swaff Tri Diag OFFICIAL Correct Indicies!
!!   Variable Listing - Far from Updated hahah




PROGRAM solver

  use precise, only : defaultp   ! Module to handle precision
  use geo                        ! Module to write the Points on the surface
  !!use Source                     ! Module to compute an arbitrary source term.
  use inv                        ! My tri dag inversion module from Spring 2011
  use Swafftridag                ! Swaffs version

  IMPLICIT NONE  

  integer, parameter :: WP=defaultp
  real, parameter :: e=2.718281828459045
  INTEGER, DIMENSION(5), PARAMETER :: it1 = (/ 1, 5, 10, 15, 40 /)   ! GNU output times

  integer::count_0, count1, count2, count_rate, count_max, walltime, tstart, tend

  Integer :: NI, i, t, NT, m, mmax
  Real(WP) :: dx, L, dt, Umx, norm
  Real(WP), DIMENSION(2) :: ux           ! Unit Vectors
  REAL(wp), Allocatable, DIMENSION(:) :: X, ratiodtdx ! X locations
  REAL(wp), Allocatable, DIMENSION(:) :: am, bm, cm ! Diagonal Matrices
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: RHS, Delta, U, Uexact, Aplus, Aminus, L1, L2, fplus, fminus, F, L2norm ! Velocity, scheme components
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: Unm  

  !! FV variables
  Integer :: NV
  REAL(wp), Allocatable, DIMENSION(:) :: V 
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: vU 

  ! Simple start up message
  WRITE(6,'(A)') 'CFD Burgers Schemes Code'
  Write(6,*)     ' 10-09-2011 Implementation By Luke McCulloch'
  write(6,*) ''
  Write(6,*)     ' Implicit Solver FV Inviscid Burgers Equation in 1D'
  write(6,*) '--------------------------------------------------------'


  Umx=-1.0e10
  dx=0.05
  L=3.
  NV=int(real(L/dx))   !! NI and NV have been reversed for convienience.
  NV=NV+1
  NI=NV+1
  dt = 0.03333333 
  NT=40

  Allocate( X(NI) )
  Allocate( U(NI,NT) )
  Allocate( Uexact(NI,NT) )
  Allocate( am(NI-3) )
  Allocate( bm(NI-2) )
  Allocate( cm(NI-3) )

  Allocate( L1(NI,NT) )
  Allocate( L2(NI,NT) )
  Allocate( Aplus(NI,NT) )
  Allocate( Aminus(NI,NT) )

  Allocate( Fplus(NI-2,NT) )
  Allocate( Fminus(NI-2,NT) )
  Allocate( F(NI-2,NT) )

  Allocate( Delta(NI-2,NT) ) 
  Allocate( RHS(NI,NT) )
  Allocate( Unm(NI,NT) )

  Allocate( L2norm(NT,10) )
  !! Allocate the FV variables
  Allocate( vU(NV,NT) )
  Allocate( V(NV) )  

 
  Call D1(NV, dx, V )   ! Subroutine to compute the vector of Points.
  Call V1(NI, dx, X)    ! Now, X's will be at Volume Centroids


  !Check Points - Comment this in the full Version
  Write(6,'(AA)') 'Volume Centroids'
  do i=1,NI
     write(*,*) 'I, X(i)', i, X(i)
  end do
  Write(6,*) ''

  !!-------------------------------------

  Write(6,'(AA)') 'Points'
  do i=1,NV
     write(*,*) 'I, V(i)', i, V(i)
  end do
  Write(6,*) ''
  !Get rid of this in the full Version, obv

  !Establish the B.C.'s-------------------------------------------------
  Do i=1,NI
     U(i,0)=0.0  ! U at the Volume centers now
  End Do
  Do i=1,NV
     vU(i,0)=0.0  !specify U at Volume boundaries now
  End do  

  write(*,*) '----------------Boundary Conditions on the Boundary Faces--------------'
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
  count2=int(real(3.0/dx))+2
  Do i=count1,count2
     u(i,0)=0.5
     Unm(i,0)=0.5
     write(*,*) 'i, x(0), u(x,0)', i, X(i), U(i,0)
  End do
  write(*,*) ''  
 
  !End of the B.C. portion------------------------------------------------


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
     ! Note:  Contains the OFFICIAL Swaff Tridag Correct indicies of the sub, main, and super diagonals!
     !
     !!Do the Newton Iterations from 1-10
     mmax=10
     do m=1,mmax
!        write(*,*)  'Begin m loop'
 
        norm=0.0
        do i=2,NI-2
           !!am(i)=(-dt/(2.*dx))*Delta(i-1,t-1)
           am(i)=-(dt/dx)*Unm(i,t-1)   !! Range = 3,NI-1
        End do
     
        do i=1,NI-3
           !!cm(i)=(dt/(2.*dx))*Delta(i+1,t-1)
           cm(i)=0.0    !! Range = 2,NI-2
        End do



        
        do i=1,NI-2               
           Fplus(i,t-1) =((Unm(i+1,t-1))**2.)  !! Range = 2,NI-1
           Fminus(i,t-1)=((Unm(i,t-1))**2.)    !! Range = 1,NI-2
           F(i,t-1)=Unm(i+1,t-1)-U(i+1,t-1)+0.5*dt*(Fplus(i,t-1)-Fminus(i,t-1))/dx
           !write(*,*) 'i, Umn(i,t-1),  Umn(i+1,t-1)', i, Unm(i,t-1),  Unm(i+1,t-1)
           bm(i)=1.0+(dt/dx)*Unm(i+1,t-1)   !! Range = 2,NI-1
        
        end do

        
        !! Call my solver
        !Call tlmtri (am(:),bm(:),cm(:),-.5*(dt/dx)*(Fplus(:,t-1)-Fminus(:,t-1)),Delta(:,t),NI-2)
        
        !! Call Swaff Solver
        Call trid (am,bm,cm,-F(:,t-1),Delta(:,t),NI-2)
        norm=0.0
        
        Do i=1,NI-2
           Unm(i+1,t-1)=Delta(i,t)+Unm(i+1,t-1)  !!Caareful, Unm range = [1,NI]
           !!Compute the L2 norm whilst you are down here:
           !write(*,*) 'Delta', Delta(i,t)
           norm=Delta(i,t)**2.+norm
        End Do

        L2norm(t,m)=sqrt(norm)/real(NI-2)
  


     Unm(1,t)=1.5
     Unm(NI,t)=0.5  
       
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
     
    ! write(6,*) 'Timestep', T

  End do    !! End Timestep

!!----------------------------End Euler Scheme----------------------------------------

write(*,*) 'L2norm(5,m)'
do m=1,10
   write(*,*) L2norm(5,m)
end do
write(*,*)
write(*,*) 'L2norm(10,m)'
do m=1,10
   write(*,*) L2norm(10,m)
end do
write(*,*)
write(*,*) 'L2norm(15,m)'
do m=1,10
   write(*,*) L2norm(15,m)
end do
write(*,*)
write(*,*) 'L2norm(40,m)'
do m=1,10
   write(*,*) L2norm(40,m)
end do
write(*,*)
!!----------------------------Exact Solution----------------------------------------

 t=40*dt
    Do i=1,NI
       Uexact(i,40)=U(i,0)
    End Do
 print*
 Write(6,*) 'solving exact solution at t=40'


!!----------------------------End Exact Solution----------------------------------------


  ! output Gnuplot file-------------------------Solution---------------------------
  OPEN(UNIT=1, FILE='EulerImplicit_FV.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
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

!!-------------------------------------------------------L2NORM----------------------
  ! output Gnuplot file
  OPEN(UNIT=2, FILE='L2norm.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  Do t=2,5
     DO m=1, 10
        WRITE(UNIT=2,FMT='(I5," ",E19.10)') m,L2norm(it1(t),m)
     END DO
     WRITE(UNIT=2,FMT='(" ")')
  End do
  CLOSE(UNIT=2)



!! Plotting in gnu plot
!! once in the proper directory in terminal, open a second tab
!! type:  gnuplot   to open gnuplot
!! type:  plot "Gnu.dat" w l        to plot


  IF (ALLOCATED(X)) DEALLOCATE(X)
  IF (ALLOCATED(U)) DEALLOCATE(U) 
  IF (ALLOCATED(Unm)) DEALLOCATE(Unm)   
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
  IF (ALLOCATED(F)) DEALLOCATE(F) 

  IF (ALLOCATED(Delta)) DEALLOCATE(Delta) 
  IF (ALLOCATED(RHS)) DEALLOCATE(RHS) 


  IF (ALLOCATED(L2norm)) DEALLOCATE(L2norm)
  !!Deallocate FV variables

  IF (ALLOCATED(vU)) DEALLOCATE(vU) 
  IF (ALLOCATED(V)) DEALLOCATE(V) 

END PROGRAM solver


