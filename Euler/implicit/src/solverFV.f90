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

  Integer :: NI, i, j, k, t, NT,  mmax
  Integer :: min, max 
  Real(WP) :: dx, L, dt, norm, freestreamMach, gamma
  REAL(WP), DIMENSION(3,3) :: Identity
  Real(WP), Allocatable, Dimension(:,:,:) :: Q
  Real(WP), Allocatable, Dimension(:,:,:) :: dfdq, dhdq
  Real(WP), Allocatable, Dimension(:,:) :: H, F 
  REAL(wp), Allocatable, DIMENSION(:) :: X, r, S, Vol, deltaS, deltaVol, p
  REAL(wp), Allocatable, DIMENSION(:,:,:) :: am, bm, cm
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: L2, L2norm
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: deltaQ, RHS



  !REAL(WP), DIMENSION(3,52,0:100) :: Q
  !REAL(WP), DIMENSION(3,3,52) :: dfdq
  !REAL(WP), DIMENSION(3,3,52) :: dhdq

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
  NT=500




  Allocate( X(NI) )

!!-----------------------------
  Allocate( Q(3,NI,0:nt) )
  Allocate( F(3,NI) )
  Allocate( H(3,NI) )

  Allocate( deltaQ(3,NI) )
  Allocate( RHS(3,NI) )

  Allocate( dfdq(3,3,NI) )
  Allocate( dhdq(3,3,NI) )


  Allocate( r(NI) )  
  Allocate( S(NI) )  
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

  Allocate( L2norm(3,0:NT) )

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


  Write(6,'(AA)') 'Points'
  do i=1,NV
     write(*,*) 'I, X(i)', i, V(i)
  end do
  Write(6,*) ''


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
     r(1)=1.398+0.347*tanh(-4.0)
  !* Compute r(x) *!
  Do i=2,NI
     r(i)=1.398+0.347*tanh(0.8*x(i)-4.0)
  End do




  !* Set the Tube geometry  *!
  Do i=1,NI
     S(i)=pi*r(i)**2
  End do

  Do i=2,NI
     !*Frustum Volume*!
     if (i==2) then
        dx=dx/2.0
     end if
     Vol(i)=(pi*dx/3.0)*( (r(i)**2)+(r(i-1)**2)+r(i)*r(i-1) )
     deltaS(i)=S(i)-S(i-1)
     deltaVol(i)=Vol(i)
  End do
  deltaS(1)=0.0!deltaS(2)
  deltaVol(1)=0.0!deltaVol(2)

  do i=1,NI
     write(*,*) 'deltaS', deltaS(i)
  end do
  
  freestreamMach=1.5
  gamma=1.4
  !* --------------- Compute Q inflow -------------- *!
  do i=1,NI
     q(1,i,0)=1.0
     q(2,i,0)=freestreamMach
     q(3,i,0)=(1./(gamma*(gamma-1.0))+0.5*freestreamMach**2)
  End do
  !* End ---------------------------------------------*!
  p(1)=(gamma-1.0)*q(3,1,0)-0.5*(gamma-1.0)*(q(2,1,0)**2)/q(1,1,0)
  !* --------------- Compute F inflow --------------- *!
  f(1,1)=S(i)*q(2,1,0)
  f(2,1)=S(i)*(((q(2,1,0)**2)/q(1,1,0))*(3.-gamma)/2.)+(gamma-1.0)*q(3,1,0)
  f(3,1)=S(i)*(q(2,1,0)/q(1,1,0))*(gamma*q(3,1,0) - 0.5*(gamma-1.0)*((q(2,1,0)**2)/q(1,1,0)))
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

     !* Enforce Boundary Nodes *!
     Do j=1,3
     !   q(j,1,t)=q(j,1,0)
        q(j,NI,t)=q(j,NI-1,t)
     End do
 
     norm=0.0

     !* Update F,H and Jacobians throughout the domain *!
     Do i=1,NI             
        !* --------------- Compute F ---------------------------------------------------- *!
        f(1,i)=S(i)*q(2,i,t)
        f(2,i)=S(i)*(  ( ( (q(2,i,t)**2)/q(1,i,t) )*(3.0-gamma)/2.0 )+(gamma-1)*q(3,i,t)  )
        f(3,i)=S(i)*(q(2,i,t)/q(1,i,t))*(gamma*q(3,i,t) - 0.5*(gamma-1)*((q(2,i,t)**2)/q(1,i,t)))
        !* End ---------------------------------------------------------------------------*!
        
        !* --------------- Compute H ---------------------------------------------------- *!
        h(1,i)=0.0
        h(2,i)=deltaS(i)*((gamma-1.0)*q(3,i,t)-0.5*(gamma-1.0)*((q(2,i,t)**2)/q(1,i,t)))
        h(3,i)=0.0
        !* End ---------------------------------------------------------------------------*!                              
     End do   


     !* Compute flux jacobians  *!
     Do i=1,NI 
        !*  A=df/dq  *!
        dfdq(1,1,i)=0.0
        dfdq(1,2,i)=1.0
        dfdq(1,3,i)=0.0
        dfdq(2,1,i)=0.5*(gamma-3.0)*(q(2,i,t)**2)/(q(1,i,t)**2)
        dfdq(2,2,i)=(3.0-gamma)*q(2,i,t)/q(1,i,t)
        dfdq(2,3,i)=(gamma-1.0)
        dfdq(3,1,i)=(gamma-1.0)*(q(2,i,t)**3/(q(1,i,t)**3)) - q(2,i,t)*q(3,i,t)*(gamma)/(q(1,i,t)**2)
        dfdq(3,2,i)=((1.0/q(1,i,t))*(gamma*q(3,i,t)-0.5*(gamma-1.0)*(q(2,i,t)**2)/q(1,i,t)))
        dfdq(3,3,i)=(gamma*q(2,i,t)/q(1,i,t))

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
     !        out of scope sub and super    !
              min=1
              max=NI-2
     !--------------------------------------!

     !* Compute the subdiagonal *!
     do k=1,3
        do j=1,3
           do i=min+1,max  
              am(j,k,i)=-(dt/deltaVol(i+1))*S(i+1)*dfdq(j,k,i)   
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
              bm(j,k,i)=(Identity(j,k)) + (dt/deltaVol(i+1))*( S(i+1)*dfdq(j,k,i+1) - dhdq(j,k,i+1) )
           End do
        End Do
     End Do


     !* Compute the superdiagonal *!
     do k=1,3   
        do j=1,3
           do i=min,max-1
              cm(j,k,i)=0.0
              cm(j,k,max)=0.0      
           End do
        End Do
     End Do
     
   
     !* Compute the RHS *!
     do j=1,3
        do i=min,max   
           RHS(j,i)=(-dt/deltaVol(i+1))*(F(j,i+1)-F(j,i)-H(j,i+1))
        End do
     end do
 
!!$     do i=min,5   
!!$        do k=1,3
!!$           do j=1,3
!!$              write(*,*) 'Sub,Main,Sup', am(j,k,i), bm(j,k,i), cm(j,k,i)
!!$           end do
!!$        end do
!!$        write(*,*) ''
!!$     end do
  !------------------------------------------------------------
!!$     Do j=1,3
!!$        do k=1,3
!!$           write(*,*) 'dfdq', dfdq(j,k,1)
!!$        End do
!!$     end do
!!$     write(*,*) ''
!!$     Do j=1,3
!!$        do k=1,3
!!$           write(*,*) 'dhdq', dhdq(j,k,1)
!!$        End do
!!$     end do
!!$     write(*,*) ''
  !------------------------------------------------------------

     !! Call block tri diagonal Solver
     Call solve (am,bm,cm,RHS,1,50)


     !* Sticking with block indices*!
     Do j=1,3
        Do i=min,max
           deltaQ(j,i+1)=RHS(j,i)
        End do
     End do

     Do j=1,3
        deltaQ(j,1)=0.0
        deltaQ(j,NI)=0.0
     end do

     norm=0.0

     !*-- Update Q --------------------------------------------------------------------*!
     !* Real Indicies *!
     Do j=1,3
        Do i=1,NI
           q(j,i,t+1)=deltaQ(j,i)+q(j,i,t)
        End Do
     End do
     !*-- End -------------------------------------------------------------------------*!
     !*-- Update norm --------------------------------------------------------------------*!
     !* Real Indicies *!
     Do j=1,3
        Do i=2,NI-1
           norm=deltaQ(1,i)**2.+norm
        End Do
     End do
     !*-- End -------------------------------------------------------------------------*!

  !------------------------------------------------------------
  write(*,*) 'T=', t
  Do i=1,NI
     write(*,*) 'r,q(1-3),deltq(1-3)', i, q(1,i,t), q(2,i,t), q(3,i,t), deltaQ(1,i), deltaQ(2,i), deltaQ(3,i)
  End do
  write(*,*) ''
  !------------------------------------------------------------


     L2norm(1,t)=sqrt(norm)/real(NI-2)

  End do    !! End Timestep

!!----------------------------End Euler Scheme----------------------------------------

write(*,*) 'L2norm(1,t)'
do t=1,NT
   write(*,*), t, L2norm(1,t)
end do
write(*,*)

write(*,*)


!*Replace q's with velocity and pressure*!
do i=1,NI
   !*Velocity*!
   q(2,i,NT)=q(2,i,NT)/q(1,i,NT)
   !*Pressure*!
   q(3,i,NT)=(gamma-1.0)*q(3,i,NT) - 0.5*(gamma-1.0)*(q(2,i,t)**2)/q(1,i,t)
end do







!!$!* Test matricies *!
!!$!Non - Swafford Form!
!!$ do i=1,3
!!$    do j=1,3
!!$       do k=1,4          
!!$          Dtest(i,j,k)=0.0               
!!$          if(i==j) then
!!$             Dtest(i,j,k)=1.0
!!$          End if 
!!$       End do
!!$    end do
!!$ end do
!!$
!!$ do i=1,3
!!$    do j=1,3
!!$       do k=1,4
!!$          Utest(i,j,k)=1.0   
!!$          Ltest(i,j,k)=1.0  
!!$          Ltest(i,j,1)=0.0  
!!$          Utest(i,j,4)=0.0  
!!$       End do
!!$    end do
!!$ end do
!!$
!!$  
!!$
!!$     Btest(1,1)=4.0
!!$     Btest(2,1)=4.0
!!$     Btest(3,1)=4.0
!!$
!!$     Btest(1,4)=4.0
!!$     Btest(2,4)=4.0
!!$     Btest(3,4)=4.0
!!$
!!$     Btest(1,2)=7.0
!!$     Btest(2,2)=7.0
!!$     Btest(3,2)=7.0
!!$
!!$     Btest(1,3)=7.0
!!$     Btest(2,3)=7.0
!!$     Btest(3,3)=7.0
!!$  
!!$  write(*,*) '1'
!!$  !* Call block tri diagonal Solver *!
!!$  !Call solve (Ltest,Dtest,Utest,Btest,1,4)
!!$
!!$  write(*,*) '2'
!!$
!!$  do i=1,3
!!$     do k=1,4
!!$        Xtest(i,k)=0.0
!!$     end do
!!$  end do
!!$
!!$  do i=1,3
!!$     do k=1,4
!!$        Xtest(i,k)=Btest(i,k)
!!$     end do
!!$  end do
!!$
!!$  do i=1,3
!!$     do k=1,4
!!$        write(*,*) Xtest(i,k)
!!$     end do
!!$  end do
!!$
!!$  write(*,*)




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

!!-------------------------------------------------------L2NORM----------------------
  ! output Gnuplot file
  OPEN(UNIT=2, FILE='L2norm.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  Do t=1,nt
     WRITE(UNIT=2,FMT='(I5,"",E19.10)') t, L2norm(1,t)
  End do
  WRITE(UNIT=2,FMT='(" ")')
  CLOSE(UNIT=2)

  ! output Gnuplot file-------------------------Solution---------------------------
  OPEN(UNIT=5, FILE='r.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
     DO i=1, NI
        WRITE(UNIT=5,FMT='(E19.10," ",E19.10)') x(i),r(i)
     End do
     WRITE(UNIT=5,FMT='(" ")')
!  Do i=1,NI
!     WRITE(UNIT=4,FMT='(E19.10," ",E19.10)') x(i)+dt*Nt*0.5*(1.5+.5),Uexact(i,40) 
!  End do
  CLOSE(UNIT=5)


!! Plotting in gnu plot
!! once in the proper directory in terminal, open a second tab
!! type:  gnuplot   to open gnuplot
!! type:  plot "Gnu.dat" w l        to plot


  IF (ALLOCATED(X)) DEALLOCATE(X)
  !IF (ALLOCATED(U)) DEALLOCATE(U) 
  !IF (ALLOCATED(Unm)) DEALLOCATE(Unm)   
  !IF (ALLOCATED(Uexact)) DEALLOCATE(Uexact) 
  IF (ALLOCATED(am)) DEALLOCATE(am) 
  IF (ALLOCATED(bm)) DEALLOCATE(bm) 
  IF (ALLOCATED(cm)) DEALLOCATE(cm)  

  !IF (ALLOCATED(L1)) DEALLOCATE(L1) 
  IF (ALLOCATED(L2)) DEALLOCATE(L2) 
  !IF (ALLOCATED(Aplus)) DEALLOCATE(Aplus)  
  !IF (ALLOCATED(Aminus)) DEALLOCATE(Aminus) 

  !IF (ALLOCATED(Fplus)) DEALLOCATE(Fplus) 
  !IF (ALLOCATED(Fminus)) DEALLOCATE(Fminus)  
  IF (ALLOCATED(RHS)) DEALLOCATE(RHS) 

  !IF (ALLOCATED(Delta)) DEALLOCATE(Delta) 


  IF (ALLOCATED(L2norm)) DEALLOCATE(L2norm)


  IF (ALLOCATED(vU)) DEALLOCATE(vU) 
  IF (ALLOCATED(V)) DEALLOCATE(V) 

  IF (ALLOCATED(r)) DEALLOCATE(r) 
  IF (ALLOCATED(S)) DEALLOCATE(S)
  IF (ALLOCATED(Vol)) DEALLOCATE(Vol)
  IF (ALLOCATED(deltaS)) DEALLOCATE(deltaS)
  IF (ALLOCATED(deltaVol)) DEALLOCATE(deltaVol)

  IF (ALLOCATED(q)) DEALLOCATE(q) 
  IF (ALLOCATED(f)) DEALLOCATE(f) 
  IF (ALLOCATED(h)) DEALLOCATE(h) 
  IF (ALLOCATED(p)) DEALLOCATE(p)   

  IF (ALLOCATED(deltaQ)) DEALLOCATE(deltaQ) 
  IF (ALLOCATED(dfdq)) DEALLOCATE(dfdq)   
  IF (ALLOCATED(dhdq)) DEALLOCATE(dhdq)   



END PROGRAM solver


