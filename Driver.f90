program Driver
    use Zeitintegration
    implicit none
    REAL(KIND=RP) :: t=0.0_rp,tend,delx,dely,delz,gam=1.4_rp
    INTEGER,parameter :: n=4,nq=10
    REAL(KIND=RP),Dimension(:,:,:,:,:,:,:),allocatable :: u,xyz
    REAL(KIND=RP),Dimension(:),allocatable :: xi,xl
    REAL(KIND=RP),Dimension(:,:),allocatable :: xin
    INTEGER :: l,m,o,k,i,j



    call LegendreGaussLobattoNodesAndWeights(n,xi,w)

   !! Bestimme GL punkte in jeder Zelle
    do k=0,nq-1
    xl(k+1)=(k+1.0_rp/2)*delx
    do i=1,n+1
    xin(i,k+1)=xl(k+1)+delx/2*xi(i)
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


   u(:,:,:,:,:,:,1)=2.0_rp+sin(pi*(xyz(:,:,:,:,:,:,1)+xyz(:,:,:,:,:,:,2)+xyz(:,:,:,:,:,:,3)))/10.0_rp
   u(:,:,:,:,:,:,2)=2.0_rp+sin(pi*(xyz(:,:,:,:,:,:,1)+xyz(:,:,:,:,:,:,2)+xyz(:,:,:,:,:,:,3)))/10.0_rp
   u(:,:,:,:,:,:,3)=2.0_rp+sin(pi*(xyz(:,:,:,:,:,:,1)+xyz(:,:,:,:,:,:,2)+xyz(:,:,:,:,:,:,3)))/10.0_rp
   u(:,:,:,:,:,:,4)=2.0_rp+sin(pi*(xyz(:,:,:,:,:,:,1)+xyz(:,:,:,:,:,:,2)+xyz(:,:,:,:,:,:,3)))/10.0_rp
   u(:,:,:,:,:,:,5)=u(:,:,:,:,:,:,1)*u(:,:,:,:,:,:,1)
                




    
end program Driver
