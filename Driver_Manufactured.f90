program Driver_Manufactured
  use Zeitintegration
  implicit none
  REAL(KIND=RP)                                      :: t=0.0_rp,tend=1.0_RP,CFL=0.4_RP,dt,a
  INTEGER,parameter                                  :: n=1,anz=1
  REAL(KIND=RP),DIMENSION(:,:,:,:,:,:,:),allocatable :: u, usolution
  REAL(KIND=RP),DIMENSION(1:n+1,1:n+1)               :: D
  CHARACTER(len=2)                                   :: whichflux='PI'
  REAL(KIND=RP),DIMENSION(1:5,1:anz)                 :: errors
  INTEGER, DIMENSION(1:anz)                          :: nq
  INTEGER                                            :: l,m,o,k,i,j,start=3
  nq=2**(/ (I,I=start,start+anz-1) /)
  DO k=1,anz
    allocate(u(1:Nq(k),1:nq(k),1:nq(k),1:n+1,1:n+1,1:n+1,1:5),usolution(1:Nq(k),1:nq(k),1:nq(k),1:n+1,1:n+1,1:n+1,1:5))
    call Vorbereiten(n,nq(k),D,t)
    call Initialcondition(u,NQ(k),N)
    call lambdaMaxGlobal(u,a,NQ(k),N)
    dt=CFL/(3*a)*(dx/real(N+1,KIND=RP))
    !-ffpe-trap=denormal,invalid,zero,overflow,underflow
    DO while(tend-t>epsilon(dt))
      print*,'t'
      print*,t
      print*,'dt'
      print*,dt
      print*,'sum(energy)'
      print*,sum(U(:,:,:,:,:,:,5))
      call lambdaMaxGlobal(u,a,NQ(k),N)
      dt=CFL/(3*a)*(dx/real(N+1))
      IF(t+dt>tend) dt=tend-t
      call RungeKutta5explizit(u,nq(k),n,5,dt,D,whichflux)
      !print*,u(1,1,1,1,:,:,1)
      t=t+dt
    END DO
    ! Berechne Fehler und Loesung
    call computeSolution(usolution,NQ(k),N) 
    call computeError(u,usolution,NQ(k),N,errors(:,k)) 
    print*,errors(1,:)
    print*,errors(2,:)
    print*,errors(3,:)
    print*,errors(4,:)
    print*,errors(5,:)
    ! Setzte alles wieder auf 0
    deallocate(u,usolution)
    t=0.0_RP
  END DO 
end program Driver_Manufactured
