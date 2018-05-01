module Zeitintegration
    use Quadraturroutinen
    REAL(KIND=RP),DIMENSION(:),allocatable      ::x,w,xmit,xges
    REAL(KIND=RP)                               ::gk=9.812_RP,dx,gamma=1.4_RP
contains
    subroutine Vorbereiten(n,nq,Dval,Sval)
        implicit none
        INTEGER,INTENT(IN)      ::n,nq
        INTEGER                 ::j,i
        REAL(KIND=RP),DIMENSION(1:n+1,1:n+1),INTENT(out)    :: Dval,Sval
        allocate(x(1:n+1),w(1:n+1),xges(1:nq*(n+1)),xmit(1:nq+1))
        call LegendreGaussLobattoNodesandWeights(N,x,w)
        call DifferentiationsmatrixBerechnen(x,Dval,N+1)
        Sval=0.0_RP

        Sval(1,1)=1.0_RP/w(1)
        Sval(n+1,n+1)=-1.0_RP/w(n+1)
        dx=1.0_RP/real(nq,kind=RP)
        call linspace(dx/2.0_RP,1.0_RP-dx/2.0_RP,NQ,xmit)
        do i=1,NQ
            xges((i-1)*(N+1)+1:i*(N+1))=xmit(i)+dx/2.0_RP*x
        enddo
    end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Operator aus Semidiskreterdarstellung
    function R(u,t,N,NQ,D) result(solution)
        implicit none
        integer, intent(in) :: n,NQ
        real(kind=RP), dimension(1:NQ,1:nq,1:nq,1:n+1,1:(N+1),1:n+1,1:5)                 :: solution
        real(kind=RP), intent(in)                                           :: t
        real(kind=RP), intent(in), dimension(1:nq*nq,1:(n+1),1:n+1,1:3)     :: u
        REAL(KIND=RP),intent(in),dimension(:,:)                             :: D



    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine computeL(u,D,dir,result)
        implicit none
        REAL(KIND=RP),DIMENSION(:,:,:,:,:,:,:),intent(in)           ::u
        INTEGER,intent(in)                                      ::dir
        REAL(KIND=RP),DIMENSION(:,:),intent(in)                 ::D
        !u dimensions:(nummer zelle x,zelle y,zelle z,x,y,z,variable)
        REAL(KIND=RP),dimension(:,:,:,:,:,:,:),allocatable,intent(out)::result
        !local Variables
        INTEGER                                             ::var,k,j,i,o,l,m,n,nq
        REAL(KIND=RP),dimension(:,:),allocatable                        ::Fsharp
        nq=size(u,dim=1)
        n=size(u,dim=4)
        select case(dir)
            case(1)
                do var=1,5
                do k=1,n+1
                    do j=1,n+1
                        do i=1,n+1
                            do o=1,nq
                                do l=1,nq
                                    do m=1,nq
                                        call computeFsharp(u(m,l,o,i,j,k,:),u(m,l,o,:,j,k,:),dir,'ST',Fsharp)
                                        result(m,l,o,i,j,k,var)=2*dot_product(D(i,:),Fsharp(:,var))
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
                enddo

        end select
    end subroutine
    subroutine calculateEuler3DFlux(u,dir,result)
    !Subroutine gets the values of u and return the flux for the spcified direction
    !dir=1,2,3 stands for x,y,z direction
        implicit none
        REAL(KIND=RP),dimension(:,:,:,:),intent(in)     :: u
        INTEGER                         ,intent(in)     :: dir
        REAL(KIND=RP),dimension(:,:,:,:),intent(out),allocatable    :: result
        !local variables beyond here
        INTEGER                                         :: n
        REAL(KIND=RP),dimension(:,:,:),allocatable    ::p
        n=size(u,dim=1)
        ALLOCATE(result(1:n+1,1:n+1,1:n+1,5))
        ALLOCATE(p(1:n+1,1:n+1,1:n+1))
        p=(gamma-1.0_RP)*(u(:,:,:,5)-0.5_RP*(u(:,:,:,2)**2+u(:,:,:,3)**2+u(:,:,:,4)**2)/u(:,:,:,1))
        SELECT CASE (dir)
            CASE(1)
                result(:,:,:,1)=u(:,:,:,2)
                result(:,:,:,2)=u(:,:,:,2)**2/u(:,:,:,1)+p
                result(:,:,:,3)=u(:,:,:,2)*u(:,:,:,3)/u(:,:,:,1)
                result(:,:,:,4)=u(:,:,:,2)*u(:,:,:,4)/u(:,:,:,1)
                result(:,:,:,5)=(u(:,:,:,5)+p)*u(:,:,:,2)/u(:,:,:,1)
            CASE(2)
                result(:,:,:,1)=u(:,:,:,3)
                result(:,:,:,2)=u(:,:,:,2)*u(:,:,:,3)/u(:,:,:,1)
                result(:,:,:,3)=u(:,:,:,3)**2/u(:,:,:,1)+p
                result(:,:,:,4)=u(:,:,:,3)*u(:,:,:,4)/u(:,:,:,1)
                result(:,:,:,5)=(u(:,:,:,5)+p)*u(:,:,:,3)/u(:,:,:,1)
            CASE(3)
                result(:,:,:,1)=u(:,:,:,4)
                result(:,:,:,2)=u(:,:,:,2)*u(:,:,:,4)/u(:,:,:,1)
                result(:,:,:,3)=u(:,:,:,3)*u(:,:,:,4)/u(:,:,:,1)
                result(:,:,:,4)=u(:,:,:,4)**2/u(:,:,:,1)+p
                result(:,:,:,5)=(u(:,:,:,5)+p)*u(:,:,:,4)/u(:,:,:,1)
            CASE DEFAULT
                print*,'Specify 1,2 or 3'
        end SELECT
    end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine computeFsharp(u1,u2,dir,whichflux,result)
    !Subroutine computes the Volume flux Fsharp
        implicit none
        REAL(KIND=RP),intent(in) ,dimension(5)                       ::u1
        REAL(KIND=RP),intent(in),dimension(:,:)           ::u2
        CHARACTER(len=*),intent(in)                     ::whichflux
        INTEGER         ,intent(in)                     ::dir
        REAL(KIND=RP),intent(out),dimension(:,:),allocatable   :: result
        !local variables go beyond here
        INTEGER                                         ::n
        REAL(KIND=RP),dimension(:),allocatable                   ::p2
        REAL(KIND=RP)                               ::p1
        n=size(u2,dim=1)
        allocate(result(n,5))
        p1=(gamma-1.0_RP)*(u1(5)-0.5_RP*(u1(2)**2+u1(3)**2+u1(4)**2)/u1(1))
        p2=(gamma-1.0_RP)*(u2(:,5)-0.5_RP*(u2(:,2)**2+u2(:,3)**2+u2(:,4)**2)/u2(:,1))
        SELECT CASE(whichflux)
            CASE('ST')
                SELECT CASE(dir)
                    CASE(1)
                        result(:,1)=(u1(2)+u2(:,2))*0.5_RP
                        result(:,2)=(u1(2)**2/u1(1)+u2(:,2)**2/u2(:,1))*0.5_RP+(p1+p2)*0.5_RP
                        result(:,3)=(u1(2)*u1(3)/u1(1)+u2(:,2)*u2(:,3)/u2(:,1))*0.5_RP
                        result(:,4)=(u1(2)*u1(4)/u1(1)+u2(:,2)*u2(:,4)/u2(:,1))*0.5_RP
                        result(:,5)=(u1(2)/u1(1)*(u1(5)+p1)+u2(:,2)/u2(:,1)*(u2(:,5)+p2))*0.5_RP
                    CASE(2)
                        result(:,1)=(u1(3)+u2(:,3))*0.5_RP
                        result(:,2)=(u1(2)*u1(3)/u1(1)+u2(:,2)*u2(:,3)/u2(:,1))*0.5_RP
                        result(:,3)=(u1(3)**2/u1(1)+u2(:,3)**2/u2(:,1))*0.5_RP+(p1+p2)*0.5_RP
                        result(:,4)=(u1(3)*u1(4)/u1(1)+u2(:,3)*u2(:,4)/u2(:,1))*0.5_RP
                        result(:,5)=(u1(3)/u1(1)*(u1(5)+p1)+u2(:,3)/u2(:,1)*(u2(:,5)+p2))*0.5_RP
                    CASE(3)
                        result(:,1)=(u1(4)+u2(:,4))*0.5_RP
                        result(:,2)=(u1(2)*u1(4)/u1(1)+u2(:,2)*u2(:,4)/u2(:,1))*0.5_RP
                        result(:,3)=(u1(3)*u1(4)/u1(1)+u2(:,3)*u2(:,4)/u2(:,1))*0.5_RP
                        result(:,4)=(u1(4)**2/u1(1)+u2(:,4)**2/u2(:,1))*0.5_RP+(p1+p2)*0.5_RP
                        result(:,5)=(u1(4)/u1(1)*(u1(5)+p1)+u2(:,3)/u2(:,1)*(u2(:,5)+p2))*0.5_RP
                end SELECT
        END SELECT
    end subroutine
   !#subroutine Initialcondition(x,n,nq,u,x2,y)
        !implicit none

    !end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine computeLocalLaxFriedrich(uL,uR,dir,result)
    ! Compute loacal LaxFriedirch 
    ! Returns Matrix with the localLaxFriedrich at every Interface point
         implicit none
         REAL(KIND=RP), INTENT(IN), DIMENSION(:,:,:),allocatable :: uL, uR
         REAL(KIND=RP), INTENT(out), DIMENSION(:,:,:),allocatable:: result
         INTEGER, INTENT(IN):: dir 
         REAL(KIND=RP), DIMENSION(:,:,:),allocatable :: FL,FR 
         REAL(KIND=RP), DIMENSION(:,:),allocatable :: lamMax
         INTEGER :: n 

        n=size(uL,dim=1)
        ALLOCATE(result(1:n+1,1:n+1,5))
        ALLOCATE(FR(1:n+1,1:n+1,5))
        ALLOCATE(FL(1:n+1,1:n+1,5))
        ALLOCATE(lamMax(1:n+1,1:n+1))

         call calculateEulerRandFlux(uL,dir,FL)
         call calculateEulerRandFlux(uR,dir,FR)
         call lambdaMax(uL,uR,dir,lamMax)

       result(:,:,1)=(FL(:,:,1)+FR(:,:,1)-lamMax*(uR(:,:,1)-uL(:,:,1)))/2 
       result(:,:,2)=(FL(:,:,2)+FR(:,:,2)-lamMax*(uR(:,:,2)-uL(:,:,2)))/2 
       result(:,:,3)=(FL(:,:,3)+FR(:,:,3)-lamMax*(uR(:,:,3)-uL(:,:,3)))/2 
       result(:,:,4)=(FL(:,:,4)+FR(:,:,4)-lamMax*(uR(:,:,4)-uL(:,:,4)))/2 
       result(:,:,5)=(FL(:,:,5)+FR(:,:,5)-lamMax*(uR(:,:,5)-uL(:,:,5)))/2 
    
    
    
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lambdaMax (uL,uR,dir,result)
         !Computes the max of the eigenvalues at every Interface
        
        implicit none 
        REAL(KIND=RP),dimension(:,:,:),intent(in)     :: uR, uL
        INTEGER                         ,intent(in)     :: dir 
        REAL(KIND=RP),dimension(:,:),intent(OUT),allocatable   ::result
        REAL(KIND=RP),dimension(:,:),allocatable    :: pR,hR,cR,pL,hL,cL
        REAL(KIND=RP),dimension(:,:,:),allocatable :: lambda 
        INTEGER                         :: l,k,n

        n=size(uL,dim=1)



        ALLOCATE(result(1:n+1,1:n+1))
        ALLOCATE(lambda(1:n+1,1:n+1,6))
        ALLOCATE(pR(1:n+1,1:n+1))
        ALLOCATE(hR(1:n+1,1:n+1))
        ALLOCATE(cR(1:n+1,1:n+1)) 
        ALLOCATE(pL(1:n+1,1:n+1))
        ALLOCATE(hL(1:n+1,1:n+1))
        ALLOCATE(cL(1:n+1,1:n+1))
        
        pR=(gamma-1.0_RP)*(uR(:,:,5)-0.5_RP*(uR(:,:,2)**2+uR(:,:,3)**2+uR(:,:,4)**2)/uR(:,:,1)) 
        hR=(uR(:,:,5)+pR)/uR(:,:,1)
        cR=sqrt((gamma-1)*(hR-((uR(:,:,2)/uR(:,:,1))**2+(uR(:,:,3)/uR(:,:,1))**2+(uR(:,:,4)/uR(:,:,1))**2)/2))
        lambda(:,:,1)=uR(:,:,dir+1)/uR(:,:,1)
        lambda(:,:,2)=uR(:,:,dir+1)/uR(:,:,1)-cR
        lambda(:,:,3)=uR(:,:,dir+1)/uR(:,:,1)+cR
        pL=(gamma-1.0_RP)*(uL(:,:,5)-0.5_RP*(uL(:,:,2)**2+uL(:,:,3)**2+uL(:,:,4)**2)/uL(:,:,1)) 
        hL=(uL(:,:,5)+pL)/uL(:,:,1)
        cL=sqrt((gamma-1)*(hL-((uL(:,:,2)/uL(:,:,1))**2+(uL(:,:,3)/uL(:,:,1))**2+(uL(:,:,4)/uL(:,:,1))**2)/2))
        lambda(:,:,4)=uL(:,:,dir+1)/uL(:,:,1)
        lambda(:,:,5)=uL(:,:,dir+1)/uL(:,:,1)-cL
        lambda(:,:,6)=uL(:,:,dir+1)/uL(:,:,1)+cL
        lambda=abs(lambda)


        do l=1,n+1
        do k=1,n+1
        result(l,k)=maxval(lambda(l,k,:))
        enddo
        enddo
        
        end subroutine 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calculateEulerRandFlux(u,dir,result)
        IMPLICIT NONE
    !Subroutine gets the values of u and return the flux for the spcified direction
    !dir=1,2,3 stands for x,y,z direction
        REAL(KIND=RP),dimension(:,:,:),intent(in)     :: u
        INTEGER                         ,intent(in)     :: dir
        REAL(KIND=RP),dimension(:,:,:),intent(out),allocatable    :: result
        !local variables beyond here
        INTEGER                                         :: n
        REAL(KIND=RP),dimension(:,:),allocatable    ::p
        n=size(u,dim=1)
        ALLOCATE(result(1:n+1,1:n+1,5))
        ALLOCATE(p(1:n+1,1:n+1))
        p=(gamma-1.0_RP)*(u(:,:,5)-0.5_RP*(u(:,:,2)**2+u(:,:,3)**2+u(:,:,4)**2)/u(:,:,1))
        SELECT CASE (dir)
            CASE(1)
                result(:,:,1)=u(:,:,2)
                result(:,:,2)=u(:,:,2)**2/u(:,:,1)+p
                result(:,:,3)=u(:,:,2)*u(:,:,3)/u(:,:,1)
                result(:,:,4)=u(:,:,2)*u(:,:,4)/u(:,:,1)
                result(:,:,5)=(u(:,:,5)+p)*u(:,:,2)/u(:,:,1)
            CASE(2)
                result(:,:,1)=u(:,:,3)
                result(:,:,2)=u(:,:,2)*u(:,:,3)/u(:,:,1)
                result(:,:,3)=u(:,:,3)**2/u(:,:,1)+p
                result(:,:,4)=u(:,:,3)*u(:,:,4)/u(:,:,1)
                result(:,:,5)=(u(:,:,5)+p)*u(:,:,3)/u(:,:,1)
            CASE(3)
                result(:,:,1)=u(:,:,4)
                result(:,:,2)=u(:,:,2)*u(:,:,4)/u(:,:,1)
                result(:,:,3)=u(:,:,3)*u(:,:,4)/u(:,:,1)
                result(:,:,4)=u(:,:,4)**2/u(:,:,1)+p
                result(:,:,5)=(u(:,:,5)+p)*u(:,:,4)/u(:,:,1)
            CASE DEFAULT
                print*,'Specify 1,2 or 3'
        end SELECT
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine RungeKutta5explizit(ustar,t,nq,n,numvar,dt,Dval,Sval)
        IMPLICIT NONE
        !ustar=bekannte Werte
        INTEGER,INTENT(IN)                                          :: numvar,n,nq
        REAL(KIND=RP),intent(in),dimension(:,:)         :: Dval,Sval
        REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:nq,1:nq,1:nq,1:(n+1),1:n+1,1:n+1,1:numvar)  :: ustar
        REAL(KIND=RP),INTENT(IN)                                    :: t
        !local
        REAL(KIND=RP),DIMENSION(1:nq,1:nq,1:nq,1:(n+1),1:n+1,1:n+1,1:numvar)    :: g
        INTEGER                                                     :: step=1
        REAL(KIND=RP),DIMENSION(5)                                  :: a,b,c
        REAL(KIND=RP),INTENT(IN)                                    :: dt
        g=R(ustar,t,n,nq,Dval)
        a(1)=0.0_RP
        b(1)=0.0_RP
        c(1)=1432997174477.0_RP/9575080441755.0_RP
        a(2)=-567301805773.0_RP/1357537059087.0_RP
        b(2)=1432997174477.0_RP/9575080441755.0_RP
        c(2)=5161836677717.0_RP/13612068292357.0_RP
        a(3)=-2404267990393.0_RP/2016746695238.0_RP
        b(3)=2526269341429.0_RP/6820363962896.0_RP
        c(3)=1720146321549.0_RP/2090206949498.0_RP
        a(4)=-3550918686646.0_RP/2091501179385.0_RP
        b(4)=2006345519317.0_RP/3224310063776.0_RP
        c(4)=3134564353537.0_RP/4481467310338.0_RP
        a(5)=-1275806237668.0_RP/842570457699.0_RP
        b(5)=2802321613138.0_RP/2924317926251.0_RP
        c(5)=2277821191437.0_RP/14882151754819.0_RP
        do step=1,5
            g=a(step)*g+R(ustar,t+b(step)*dt,n,nq,Dval)
            ustar=ustar+c(step)*dt*g
        enddo
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module
