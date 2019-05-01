!!
!! Luke McCulloch
!! solver.f90
!! 11-9-2011
!! Qasi 1D Euler Solver


PROGRAM solver

  use precise, only : defaultp   ! Module to handle precision
  use GEO                        ! Module to write the Points on the surface
  !use Swafftridag                ! Swaffs version
  use blocktri

  IMPLICIT NONE  

  integer, parameter :: WP=defaultp
  real(wp), parameter :: e=2.718281828459045
  REAL(wp), PARAMETER :: pi=4.0*ATAN(1.0) 
  INTEGER, DIMENSION(5), PARAMETER :: it1 = (/ 1, 5, 10, 15, 40 /)   ! GNU output times


  integer::count_0, count1, count2, count_rate, count_max, walltime, tstart, tend

  Integer :: NI, i, j, k, t, NT, m, mmax
  Integer :: min, max 
  Real(WP) :: dx, L, dt, norm, freestreamMach, gamma
  real(WP) :: start_time, energy, pressure
  REAL(WP), DIMENSION(3,3) :: Identity
  Real(WP), Allocatable, Dimension(:,:,:) :: Q, Qf, qc
  Real(WP), Allocatable, Dimension(:,:,:) :: dfdq, dhdq
  Real(WP), Allocatable, Dimension(:,:) :: H, F 
  REAL(wp), Allocatable, DIMENSION(:) :: X, r, S, Vol, deltaS, deltaVol, p
  REAL(wp), Allocatable, DIMENSION(:) :: rface, Sface  
  REAL(wp), Allocatable, DIMENSION(:,:,:) :: am, bm, cm
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: L2, L2norm, newtonL2norm
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: deltaQ, RHS



  !! FV variables
  Integer :: NV
  REAL(wp), Allocatable, DIMENSION(:) :: V 
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: vU 

  !*Test Case*!
  REAL(WP), DIMENSION(3,3,4) :: Dtest
  REAL(WP), DIMENSION(3,3,4) :: Utest
  REAL(WP), DIMENSION(3,3,4) :: Ltest
  REAL(WP), DIMENSION(3,4) :: Btest 
  REAL(WP), DIMENSION(3,4) :: Xtest

  ! Simple start up message
  write(6,*) ''
  Write(6,*) ''
  write(6,*) '-----------------------------------------------------'
  write(6,*) '-----------------------------------------------------'
  write(6,*) '                                                     '
  Write(6,*) '                                                     '
  WRITE(6,'(A)') '  CFD1 Quasi 1D Euler Solver                     '
  Write(6,*)     ' 11/9/2011 Implementation By Luke McCulloch      '
  write(6,*) '                                                     '
  Write(6,*) '                                                     '
  write(6,*) '-----------------------------------------------------'
  write(6,*) '-----------------------------------------------------'


  !Umx=-1.0e10

  dx=0.2
  L=10.
  NV=int(real(L/dx))  
  NV=NV+1
  NI=NV+1
  dt = 0.5 
  NT=150




  Allocate( X(NI) )

!!-----------------------------
  Allocate( Q(3,NI,0:nt) )
  Allocate( Qf(3,NV,0:nt) )
  Allocate( qc(3,NI,0:nt) )
  !*Big Chance - there are really only 51 faces*!
  Allocate( F(3,NV) )
  Allocate( H(3,NI) )

  !*Should be NI-2*!
  Allocate( deltaQ(3,NI) )
  Allocate( RHS(3,NI) )

  !*Big Chance - there are really only 51 faces*!
  Allocate( dfdq(3,3,NV) )
  Allocate( dhdq(3,3,NI) )

  !!Allocate( rface(NV) )  
  !!Allocate( Sface(NV) ) 
  !*r now = rface*! 
  Allocate( r(NV) )  
  Allocate( S(NV) )  
  Allocate( Vol(NI) ) 
  Allocate( deltaS(NI) ) 
  Allocate( deltaVol(NI) )   

  Allocate( P(NI) ) 

  !*Block Tri full matrices*!
  Allocate( am(3,3,NI-2) )
  Allocate( bm(3,3,NI-2) )
  Allocate( cm(3,3,NI-2) )

!!-----------------------------
  !Allocate( U(NI,NT) )
  !Allocate( Uexact(NI,NT) )

  Allocate( L2(NI,0:NT) )

  Allocate( L2norm(1,0:NT) )
  Allocate( newtonL2norm(0:NT,15) )

  Allocate( vU(NV,0:NT) )
  Allocate( V(NV) )  




  !* Subroutines to compute grid points  *!
  Call D1(NV, dx, V )   ! Face points
  Call V1(NI, dx, X)    ! Volume Centers



  !Check Points - Comment this in the full Version
  Write(6,'(AA)') 'Volume Centroids'
  do i=1,NI
     write(*,*) 'I, X(i)', i, X(i)
  end do
  Write(6,*) ''


  Write(6,'(AA)') 'Face Points'
  do i=1,NV
     write(*,*) 'I, V(i)', i, V(i)
  end do
  Write(6,*) ''
  !pause


  !* Identity Matrix *!
  Do j=1,3
     Do i=1,3
        Identity(i,j) = 0.0
     End do
  End do
  Identity(1,1) = 1.0
  Identity(2,2) = 1.0
  Identity(3,3) = 1.0
  !*-----------------*!


  ! x<0 @ i=1 so set x=0
     !!r(1)=1.398+0.347*tanh(-4.0)
  !* Compute r(x) *!
  Do i=1,NV
     r(i)=1.398+0.347*tanh(0.8*V(i)-4.0)
  End do

  
  Write(6,'(AA)') 'r Points'
  do i=1,NV
     write(*,*) 'I, r(i)', i, r(i)
  end do
  Write(6,*) ''
  !pause

  ! x<0 @ i=1 so set x=0
  !r(1)=1.398+0.347*tanh(-4.0)
  !* Compute r(x) *!
!!$  Do i=1,NV
!!$     rface(i)=1.398+0.347*tanh(0.8*V(i)-4.0)
!!$     Sface(i)=pi*rface(i)**2
!!$  End do



  !* Set the Cell Face Tube geometry  *!
  Do i=1,NV
     S(i)=pi*r(i)**2
  End do

  Do i=2,NI-1
     !*Frustum Volume*!
!!     if (i==2) then
!!        dx=dx/2.0
!!     end if
     Vol(i)=(pi*dx/3.0)*( (r(i)**2)+(r(i-1)**2)+r(i)*r(i-1) )
     deltaS(i)=S(i)-S(i-1)
     deltaVol(i)=Vol(i)
  End do
  deltaS(1)=0.0!deltaS(2)
  deltaVol(1)=0.0!deltaVol(2)
  deltaS(NI)=deltaS(NI-1)
  deltaVol(NI)=deltaVol(NI-1)

  do i=1,NI
     write(*,*) 'deltaS', deltaS(i)
  end do
  
  Write(6,*) 'setup inflow'
  freestreamMach=1.5
  gamma=1.4
  !* --------------- Compute Q inflow at volume centers and faces-------------- *!
  do i=1,NI
     q(1,i,0)=1.0
     q(2,i,0)=freestreamMach
     q(3,i,0)=(1./(gamma*(gamma-1.0))+0.5*freestreamMach**2)
  End do
  do i=1,NV
     qf(1,i,0)=0.5*(q(1,i,0)+q(1,i+1,0))
     qf(2,i,0)=0.5*(q(2,i,0)+q(2,i+1,0))
     qf(3,i,0)=0.5*(q(3,i,0)+q(3,i+1,0)) 
  end do
  !* End ---------------------------------------------*!
  p(1)=(gamma-1.0)*q(3,1,0)-0.5*(gamma-1.0)*(q(2,1,0)**2)/q(1,1,0)
  !* --------------- Compute F inflow --------------- *!
  f(1,1)=S(i)*qf(2,1,0)
  f(2,1)=S(i)*(((qf(2,1,0)**2)/qf(1,1,0))*(3.-gamma)/2.)+(gamma-1.0)*qf(3,1,0)
  f(3,1)=S(i)*(qf(2,1,0)/qf(1,1,0))*(gamma*qf(3,1,0) - 0.5*(gamma-1.0)*((qf(2,1,0)**2)/qf(1,1,0)))
  do j=1,3
     f(j,1)=S(1)*f(j,1)
  end do
  !* End ---------------------------------------------*!

  !* --------------- Compute H inflow --------------- *!
  h(1,1)=0.0
  h(2,1)=deltaS(1)*p(1)*freestreamMach
  !h(2,1)=freestreamMach  
  h(3,1)=0.0
  !* End -------------------------------------------- *!
  !* End of the B.C. / I.V. portion------------------ *!


  !------------------------------------------------------------
  write(*,*) 'T=', 0
  t=0
  Do i=1,NI
     write(*,*) 'q1,q2,q3', q(1,i,t), q(2,i,t), q(3,i,t)
  End do
  write(*,*) ''
  !------------------------------------------------------------

  Write (*,*) '------------------Begin Implicit Euler Method--------------------------'
  ! max is the number of elements in the diagonal
  ! am is the subdiagonal
  ! bm is the diagonal
  ! cm is the superdiagonal

  !* Begin the time loop  *!
  Do t=0,NT

     !*Set Constant q RHS term for NM loop*!
     Do j=1,3
        do i=1,NI
           qc(j,i,t)=q(j,i,t)
        end do
     end do

     !*Newton's Method Iterative Loop*!
     do m=1,15
     norm=0.0



        !*Write Q on the faces*!
        do i=1,NV
           qf(1,i,t)=0.5*(q(1,i,t)+q(1,i+1,t))
           qf(2,i,t)=0.5*(q(2,i,t)+q(2,i+1,t))
           qf(3,i,t)=0.5*(q(3,i,t)+q(3,i+1,t)) 
        end do
        !* Enforce Boundary Nodes For Newton Loops*!
!!$        Do j=1,3
!!$           q(j,1,t)=q(j,1,0)
!!$           q(j,NI,t)=q(j,NI-1,t)
!!$           !qf(j,1,t)=qf(j,1,0)
!!$           !*The condition below "KILL" my solution*!
!!$           !qf(j,NV,t)=qf(j,NV-1,t)
!!$        End do

        

        
        !* Update F,H and Jacobians throughout the domain *!
        Do i=1,NV    
           !* --------------- Compute F ---------------------------------------------------- *!
           pressure = (gamma-1.0)*(qf(3,i,t)-(0.5*(qf(2,i,t)**2)/qf(1,i,t)))
           energy = ((1.0/(gamma-1.0))*pressure/qf(1,i,t))+0.5*((qf(2,i,t)/qf(1,i,t))**2)
           f(1,i)=S(i)*qf(2,i,t)
           !f(2,i)=S(i)*(  ( ( (qf(2,i,t)**2)/qf(1,i,t) )*(3.0-gamma)/2.0 )+(gamma-1)*qf(3,i,t)  )
           f(2,i)=S(i)*(  ( ( (qf(2,i,t)**2)/qf(1,i,t) ) +pressure))
           !f(3,i)=S(i)*(qf(2,i,t)/qf(1,i,t))*(gamma*qf(3,i,t) - 0.5*(gamma-1)*((qf(2,i,t)**2)/qf(1,i,t)))
           f(3,i)=S(i)*(qf(2,i,t)/qf(1,i,t))*(qf(1,i,t)*energy+pressure)
           !* End ---------------------------------------------------------------------------*!
        End do
        Do i=1,NI
           !* --------------- Compute H ---------------------------------------------------- *!
           h(1,i)=0.0
           h(2,i)=deltaS(i)*((gamma-1.0)*q(3,i,t)-0.5*(gamma-1.0)*((q(2,i,t)**2)/q(1,i,t)))
           h(3,i)=0.0
           !* End ---------------------------------------------------------------------------*!             
        End do
        
        
        !* Compute flux jacobians  *!
        Do i=1,NV
           !*  A=df/dq  *!
           dfdq(1,1,i)=0.0
           dfdq(1,2,i)=1.0
           dfdq(1,3,i)=0.0
           dfdq(2,1,i)=0.5*(gamma-3.0)*(qf(2,i,t)**2)/(qf(1,i,t)**2)
           dfdq(2,2,i)=(3.0-gamma)*qf(2,i,t)/qf(1,i,t)
           dfdq(2,3,i)=(gamma-1.0)
           pressure = (gamma-1.0)*(qf(3,i,t)-(0.5*(qf(2,i,t)**2)/qf(1,i,t)))
           energy = ((1.0/(gamma-1.0))*pressure/qf(1,i,t))+0.5*((qf(2,i,t)/qf(1,i,t))**2)
           !!dfdq(3,1,i)=(gamma-1.0)*(qf(2,i,t)**3/(qf(1,i,t)**3)) - qf(2,i,t)*qf(3,i,t)*(gamma)/(qf(1,i,t)**2)
           !dfdq(3,1,i)=-gamma*(qf(2,i,t)/qf(1,i,t))*energy+(gamma-1.0)*((qf(2,i,t)/qf(1,i,t))**3)

           dfdq(3,1,i)=-gamma*energy*qf(2,i,t)/qf(1,i,t)+(gamma-1.0)*((qf(2,i,t)/qf(1,i,t))**3)

           !!dfdq(3,2,i)=((1.0/qf(1,i,t))*(gamma*qf(3,i,t)-0.5*(gamma-1.0)*(qf(2,i,t)**2)/qf(1,i,t)))
           dfdq(3,2,i)=gamma*energy-(3.0*0.5)*(gamma-1.0)*((qf(2,i,t)/qf(1,i,t))**2)
           dfdq(3,3,i)=(gamma*qf(2,i,t)/qf(1,i,t))
        End do
        Do i=1,NI 
           !*  A=dh/dq  *!
           dhdq(1,1,i)=0.0
           dhdq(1,2,i)=0.0        
           dhdq(1,3,i)=0.0
           dhdq(2,1,i)=deltaS(i)*0.5*(gamma-1.0)*(q(2,i,t)**2)/(q(1,i,t)**2)
           dhdq(2,2,i)=deltaS(i)*(1.0-gamma)*q(2,i,t)/q(1,i,t)
           dhdq(2,3,i)=deltaS(i)*(gamma-1.0)
           dhdq(3,1,i)=0.0
           dhdq(3,2,i)=0.0
           dhdq(3,3,i)=0.0
        End do
        
        !--------------------------------------!
        !!*Set up the Block Tridagonal System*!!
        !   Beware: Indicies set such that     !
        !    Vbls follow convention            !
        !        which includes the            !
        !        out of scope                  !
        !         sub and super diagonals      !
        !         diagonals                    !
                   min=1
                   max=NI-2
        !--------------------------------------!
        
        !* Compute the subdiagonal *!
        do k=1,3
           do j=1,3
              do i=min+1,max  
                 am(j,k,i)=-S(i)*dfdq(j,k,i)   
                 am(j,k,min)=0.0   !Key to indici logic
              end do
           end do
        End do
     
!!$     do i=1,NI
!!$        do k=1,3
!!$           do j=1,3
!!$              write(*,*) 'am components', dt, deltaVol(i+1), S(i), dfdq(j,k,i), Vol(i)
!!$           end do
!!$        end do
!!$     end do
        
        !* Compute the main diagonal  *!
        do k=1,3
           do j=1,3
              do i=min,max
                 bm(j,k,i)= (2.0*deltaVol(i+1)/dt)*Identity(j,k) - &
                      S(i)*dfdq(j,k,i)+S(i+1)*dfdq(j,k,i+1) - &
                      2*dhdq(j,k,i+1)
              End do
           End Do
        End Do
        
        
        !* Compute the superdiagonal *!
        do k=1,3   
           do j=1,3
              do i=min,max-1  
                 cm(j,k,i)=S(i+1)*dfdq(j,k,i+1)   
                 cm(j,k,max)=0.0   !Key to indici logic  
              End do
           End Do
        End Do
        
        
        !* Compute the RHS *!
        do j=1,3
           do i=min,max   
              RHS(j,i)= -2.0*(deltaVol(i+1)*((q(j,i+1,t)-qc(j,i+1,t))/dt)+(F(j,i+1)-F(j,i)-H(j,i+1)) )
              !write(*,*) 'deltaVoli Fi+1/2, Fi-1/2, Hi', deltaVol(i+1),F(j,i+1),F(j,i),H(j,i+1)
           End do
           !write(*,*)''
        end do
        
        !! Call block tri diagonal Solver
!!$        write(*,*) 'Elements into the solver'
!!$            do i=min,max
!!$        do j=1,3  
!!$           do k=1,3
!!$  
!!$                 write(*,*) j,k,i,' am=',am(j,k,i),'bm= ', bm(j,k,i),' cm=', cm(j,k,i)
!!$              End do
!!$           End do
!!$        End do
        Call solve (am,bm,cm,RHS,1,50)
        
        
        !* Sticking with block indices*!
        Do j=1,3
           Do i=min,max
              deltaQ(j,i+1)=RHS(j,i)
           End do
        End do
        
        !*BC, #5 in assignment*!
        Do j=1,3
           deltaQ(j,1)=0.0
           deltaQ(j,NI)=0.0
        end do
        
        norm=0.0
        
        !*-- Update Q --------------------------------------------------------------------*!
        !* Real Indicies *!
        Do j=1,3
           Do i=1,NI
              q(j,i,t)=deltaQ(j,i)+q(j,i,t)
           End Do
        End do


        !* Enforce Boundary Nodes For Newton Loops*!
        Do j=1,3
           q(j,1,t)=q(j,1,0)
           q(j,NI,t)=q(j,NI-1,t)
           qf(j,1,t)=qf(j,1,0)
           !*The condition below "KILL" my solution*!
           qf(j,NV,t)=qf(j,NV-1,t)
        End do

        !*Update Q on the faces*!
!!$        do i=1,NV
!!$           qf(1,i,t)=0.5*(q(1,i,t)+q(1,i+1,t))
!!$           qf(2,i,t)=0.5*(q(2,i,t)+q(2,i+1,t))
!!$           qf(3,i,t)=0.5*(q(3,i,t)+q(3,i+1,t)) 
!!$        end do
!!$        
        !*-- Update norm --------------------------------------------------------------------*!
        !* Real Indicies *!
           Do i=2,NI-1
              norm=deltaQ(1,i)**2.+deltaQ(2,i)**2+deltaQ(3,i)**2+norm
           End Do
        !*-- End -------------------------------------------------------------------------*!
        L2norm(1,t)=sqrt(norm)/real(NI-2)
        newtonL2norm(t,m)=L2norm(1,t)
        
        !------------------------------------------------------------
        write(*,*) 'T=', t, 'm=',m
        Do i=1,NI
           write(*,*) 'r,q(1-3),deltq(1-3)', i, q(1,i,t), q(2,i,t), q(3,i,t), deltaQ(1,i), deltaQ(2,i), deltaQ(3,i)
        End do
        write(*,*) ''
        !------------------------------------------------------------        
          
     !*End Newton's Method iteration*!
     End do
     
     !*Update the new time level cell centers*!
     Do j=1,3
        do i=1,NI
           q(j,i,t+1)=q(j,i,t)
        end do
     end do
!!$        !* Enforce cell center Boundary Nodes for time loops*!
        Do j=1,3
           q(j,1,t+1)=q(j,1,0)
           q(j,NI,t+1)=q(j,NI-1,t)
        End do
        !Update Face at new time level
        do i=1,NV
           qf(1,i,t+1)=0.5*(q(1,i,t+1)+q(1,i+1,t+1))
           qf(2,i,t+1)=0.5*(q(2,i,t+1)+q(2,i+1,t+1))
           qf(3,i,t+1)=0.5*(q(3,i,t+1)+q(3,i+1,t+1)) 
        end do
        !Enforce BC on the Face!
!!$        Do j=1,3
!!$           qf(j,1,t+1)=qf(j,1,0)
!!$           qf(j,NV,t+1)=qf(j,NV-1,t+1)
!!$        End do
 
     
  !*End Timestep*!
  End do

write(6,*) 'End Euler Scheme'
!!----------------------------End Euler Scheme----------------------------------------

write(*,*) 'L2norm(1,t)'
do t=1,NT
   write(*,*) t, L2norm(1,t)
end do
write(*,*)

write(*,*)



write(6,*) 'Replace q terms with velocity and pressure'
!*Replace q's with velocity and pressure*!
do i=1,NI
   !*Velocity*!
   q(2,i,NT)=q(2,i,NT)/q(1,i,NT)
   !*Pressure*!
   q(3,i,NT)=(gamma-1.0)*q(3,i,NT) - 0.5*(gamma-1.0)*(q(2,i,t)**2)/q(1,i,t)
end do




write(*,*) 'L2norm(5,m)'
do m=1,15
   write(*,*) newtonL2norm(0,m)
end do
write(*,*)
write(*,*) 'L2norm(10,m)'
do m=1,15
   write(*,*) newtonL2norm(10,m)
end do
write(*,*)
write(*,*) 'L2norm(15,m)'
do m=1,15
   write(*,*) newtonL2norm(20,m)
end do
write(*,*)
write(*,*) 'L2norm(40,m)'
do m=1,15
   write(*,*) newtonL2norm(50,m)
end do
write(*,*)


write(6,*) 'Write outputs'


  write(6,*) 'Write x,q1'
  ! output Gnuplot file-------------------------Solution---------------------------
  OPEN(UNIT=1, FILE='q1.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  DO i=1, NI
     WRITE(UNIT=1,FMT='(E19.10," ",E19.10)') x(i),q(1,i,NT)
  End do
  WRITE(UNIT=1,FMT='(" ")')
  !Do i=1,NI
  !   WRITE(UNIT=1,FMT='(E19.10," ",E19.10)') x(i)+dt*Nt*0.5*(1.5+.5),Uexact(i,40) 
  !End do
  CLOSE(UNIT=1)

  write(6,*) 'Write x,q2'
  ! output Gnuplot file-------------------------Solution---------------------------
  OPEN(UNIT=3, FILE='q2.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  DO i=1, NI
     WRITE(UNIT=3,FMT='(E19.10," ",E19.10)') x(i),q(2,i,NT)
  End do
  WRITE(UNIT=3,FMT='(" ")')
!  Do i=1,NI
!     WRITE(UNIT=3,FMT='(E19.10," ",E19.10)') x(i)+dt*Nt*0.5*(1.5+.5),Uexact(i,40) 
!  End do
  CLOSE(UNIT=3)

  write(6,*) 'Write x,q3'
  ! output Gnuplot file-------------------------Solution---------------------------
  OPEN(UNIT=4, FILE='q3.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  DO i=1, NI
     WRITE(UNIT=4,FMT='(E19.10," ",E19.10)') x(i),q(3,i,NT)
  End do
  WRITE(UNIT=4,FMT='(" ")')
!  Do i=1,NI
!     WRITE(UNIT=4,FMT='(E19.10," ",E19.10)') x(i)+dt*Nt*0.5*(1.5+.5),Uexact(i,40) 
!  End do
  CLOSE(UNIT=4)

  write(6,*) 'Write t,L2norm'
!!-------------------------------------------------------L2NORM----------------------
  ! output Gnuplot file
  OPEN(UNIT=2, FILE='L2norm.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  Do t=1,nt
     WRITE(UNIT=2,FMT='(I5,"",E19.10)') t, L2norm(1,t)
  End do
  WRITE(UNIT=2,FMT='(" ")')
  CLOSE(UNIT=2)

  write(6,*) 'Write x,r'
  ! output Gnuplot file-------------------------Solution---------------------------
  OPEN(UNIT=5, FILE='r.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
     DO i=1, NV
        WRITE(UNIT=5,FMT='(E19.10," ",E19.10)') x(i),r(i)
     End do
     WRITE(UNIT=5,FMT='(" ")')
!  Do i=1,NI
!     WRITE(UNIT=4,FMT='(E19.10," ",E19.10)') x(i)+dt*Nt*0.5*(1.5+.5),Uexact(i,40) 
!  End do
  CLOSE(UNIT=5)

  write(6,*) 'Write it1(t),newtonL2norm'
!!-------------------------------------------------------L2NORM----------------------
  ! output Gnuplot file
  OPEN(UNIT=6, FILE='newtonL2norm.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  Do t=3,5
     DO m=1,15
        WRITE(UNIT=6,FMT='(I5," ",E19.10)') m,newtonL2norm(it1(t),m)
     END DO
     WRITE(UNIT=6,FMT='(" ")')
  End do
  CLOSE(UNIT=6)

  write(6,*) 'Write it1(t),newtonL2norm, Done'

!! Plotting in gnu plot
!! once in the proper directory in terminal, open a second tab
!! type:  gnuplot   to open gnuplot
!! type:  plot "Gnu.dat" w l        to plot

  write(6,*) 'deallocating'

  
  write(6,*) 'deallocating 1'
  IF (ALLOCATED(X)) DEALLOCATE(X)
  !IF (ALLOCATED(U)) DEALLOCATE(U) 
  !IF (ALLOCATED(Unm)) DEALLOCATE(Unm)   
  !IF (ALLOCATED(Uexact)) DEALLOCATE(Uexact) 
  IF (ALLOCATED(am)) DEALLOCATE(am) 
  IF (ALLOCATED(bm)) DEALLOCATE(bm) 
  IF (ALLOCATED(cm)) DEALLOCATE(cm)  

  
  write(6,*) 'deallocating 2'
  !IF (ALLOCATED(L1)) DEALLOCATE(L1) 
  IF (ALLOCATED(L2)) DEALLOCATE(L2) 
  !IF (ALLOCATED(Aplus)) DEALLOCATE(Aplus)  
  !IF (ALLOCATED(Aminus)) DEALLOCATE(Aminus) 

  !IF (ALLOCATED(Fplus)) DEALLOCATE(Fplus) 
  !IF (ALLOCATED(Fminus)) DEALLOCATE(Fminus)  
  IF (ALLOCATED(RHS)) DEALLOCATE(RHS) 

  !IF (ALLOCATED(Delta)) DEALLOCATE(Delta) 


  write(6,*) 'deallocating 3'
  IF (ALLOCATED(L2norm)) DEALLOCATE(L2norm)
  IF (ALLOCATED(newtonL2norm)) DEALLOCATE(newtonL2norm)

  
  write(6,*) 'deallocating 4'
  IF (ALLOCATED(vU)) DEALLOCATE(vU) 
  IF (ALLOCATED(V)) DEALLOCATE(V) 

  
  write(6,*) 'deallocating 5'
!!$  IF (ALLOCATED(rface)) DEALLOCATE(rface) 
!!$  IF (ALLOCATED(Sface)) DEALLOCATE(Sface)
  IF (ALLOCATED(r)) DEALLOCATE(r) 
  IF (ALLOCATED(S)) DEALLOCATE(S)
  IF (ALLOCATED(Vol)) DEALLOCATE(Vol)
  IF (ALLOCATED(deltaS)) DEALLOCATE(deltaS)
  IF (ALLOCATED(deltaVol)) DEALLOCATE(deltaVol)

  
  write(6,*) 'deallocating qc'
  IF (ALLOCATED(qc)) DEALLOCATE(qc)
  write(6,*) 'deallocating qf'
  IF (ALLOCATED(Qf)) DEALLOCATE(Qf) 
  write(6,*) 'deallocating q'
  IF (ALLOCATED(Q)) DEALLOCATE(Q) 
  write(6,*) 'deallocating f'
  IF (ALLOCATED(f)) DEALLOCATE(f) 
  write(6,*) 'deallocating h'
  IF (ALLOCATED(h)) DEALLOCATE(h) 
  write(6,*) 'deallocating p'
  IF (ALLOCATED(p)) DEALLOCATE(p)   

  
  write(6,*) 'deallocating 7'
  IF (ALLOCATED(deltaQ)) DEALLOCATE(deltaQ) 
  IF (ALLOCATED(dfdq)) DEALLOCATE(dfdq)   
  IF (ALLOCATED(dhdq)) DEALLOCATE(dhdq)   



END PROGRAM solver


