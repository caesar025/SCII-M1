program Driver
    use Zeitintegration
    implicit none
    REAL(KIND=RP) :: t=0.0_rp,tend=1.0_RP,CFL=0.2,dt,a
    INTEGER,parameter :: n=4,nq=10
    REAL(KIND=RP),Dimension(:,:,:,:,:,:,:),allocatable :: u,xyz
    REAL(KIND=RP),Dimension(:),allocatable :: xi,xl
    REAL(KIND=RP),Dimension(:,:),allocatable :: xin
    REAL(KIND=RP),DIMENSION(1:n+1,1:n+1)    ::D
    real(kind=RP), dimension(1:5) :: EOC,Fehler
    INTEGER :: l,m,o,k,i,j
    allocate(u(1:Nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:5),xyz(1:Nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:3))
    allocate(xi(1:n+1),xl(1:nq))
    allocate(xin(1:n+1,1:nq))
    Fehler=0.0_RP
    EOC=0.0_RP
    call Vorbereiten(n,nq,D)
    call LegendreGaussLobattoNodesAndWeights(n,xi,w)
   !! Bestimme GL punkte in jeder Zelle
    do k=0,nq-1
    xl(k+1)=(k+1.0_rp/2)*dx
    do i=1,n+1
    xin(i,k+1)=xl(k+1)+dx/2*xi(i)
    end do
    end do
  !! Bestimme alle punkte.  
    do o=1,nq
    do l=1,nq
    do m=1,nq
    do k=1,n+1
    do j=1,n+1
    do i=1,n+1
    xyz(m,l,o,i,j,k,1)=xin(i,m)
    xyz(m,l,o,i,j,k,2)=xin(j,l)
    xyz(m,l,o,i,j,k,3)=xin(k,o)
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    call Initialcondition(xyz,u)
    call lambdamaxglobal(u,a)
    dt=CFL/a*(dx/real(N+1,kind=RP))
    do while(tend-t>epsilon(dt))
        print*,t
        call lambdamaxglobal(u,a)
        dt=CFL*(dx/(a*real(N+1,kind=RP)))
        if(t+dt>tend) dt=tend-t
        call RungeKutta5explizit(u,nq,n,5,dt,D)
        t=t+dt
    end do
                




    
end program Driver
