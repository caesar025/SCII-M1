program Driver
    use Zeitintegration
    implicit none
    REAL(KIND=RP) :: t=0.0_rp,tend=10.0_RP,CFL=0.4_RP,dt,a
    INTEGER,parameter :: n=0,nq=4
    REAL(KIND=RP),Dimension(:,:,:,:,:,:,:),allocatable :: u, usolution
    REAL(KIND=RP),DIMENSION(1:n+1,1:n+1)    ::D
    real(kind=RP), dimension(1:5) :: EOC,Fehler
    INTEGER :: l,m,o,k,i,j
    allocate(u(1:Nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:5),usolution(1:Nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:5))
    Fehler=0.0_RP
    EOC=0.0_RP
    call Vorbereiten(n,nq,D,t)
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
        call RungeKutta5explizit(u,nq,n,5,dt,D)
        !print*,u(1,1,1,1,:,:,1)
        t=t+dt
    end do
       call computeSolution(usolution) 
        
end program Driver
