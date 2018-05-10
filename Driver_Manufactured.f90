program Driver_Manufactured
    use Zeitintegration
    implicit none
    REAL(KIND=RP) :: t=0.0_rp,tend=1.0_RP,CFL=0.4_RP,dt,a
    INTEGER,parameter :: n=0,anz=3
    REAL(KIND=RP),Dimension(:,:,:,:,:,:,:),allocatable :: u, usolution
    REAL(KIND=RP),DIMENSION(1:n+1,1:n+1)    ::D
    real(kind=RP), dimension(1:5,1:anz) ::errors
    INTEGER, dimension(1:anz) :: nq
    INTEGER :: l,m,o,k,i,j
    nq=2**(/ (I,I=1,anz) /)
    do k=1,anz
    allocate(u(1:Nq(k),1:nq(k),1:nq(k),1:n+1,1:n+1,1:n+1,1:5),usolution(1:Nq(k),1:nq(k),1:nq(k),1:n+1,1:n+1,1:n+1,1:5))
    call Vorbereiten(n,nq(k),D,t)
    call Initialcondition(u)
    call lambdaMaxGlobal(u,a)
    dt=CFL/(3*a)*(dx/real(N+1,kind=RP))
!-ffpe-trap=denormal,invalid,zero,overflow,underflow
    do while(tend-t>epsilon(dt))
        print*,'t'
        print*,t
        print*,'dt'
        print*,dt
        print*,'sum(energy)'
        print*,sum(U(:,:,:,:,:,:,5))
        call lambdaMaxGlobal(u,a)
        dt=CFL/(3*a)*(dx/real(N+1))
        if(t+dt>tend) dt=tend-t
        call RungeKutta5explizit(u,nq(k),n,5,dt,D)
        !print*,u(1,1,1,1,:,:,1)
        t=t+dt
    end do
! Berechne Fehler und Loesung
       call computeSolution(usolution) 
       call computeError(u,usolution,errors(:,k)) 
       print*,errors(1,:)
! Setzte alles wieder auf 0
       deallocate(u,usolution)
       t=0.0_RP
       end do 
end program Driver_Manufactured
