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
        lambda=sqrt(lame/rho)
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
    function R(u,t,N,NQ,D,S) result(solution)
        implicit none
        integer, intent(in) :: n,NQ
        real(kind=RP), dimension(1:NQ*NQ,1:n+1,1:(N+1),1:3)                 :: solution
        real(kind=RP), intent(in)                                           :: t
        real(kind=RP), intent(in), dimension(1:nq*nq,1:(n+1),1:n+1,1:3)     :: u
        REAL(KIND=RP),intent(in),dimension(:,:)                             :: D



    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine computeL(u,dir,result)
        implicit none
        REAL(KIND=RP),DIMENSION(:,:,:,:,:,:,:),intent(in)           ::u
        INTEGER,intent(in)                                      ::dir
        !u dimensions:(nummer zelle x,zelle y,zelle z,x,y,z,variable)
        REAL(KIND=RP),dimension(:,:,:,:,:,:,:),allocatable,intent(out)::result
        !local Variables
        INTEGER
        select case(dir)
        case(1)                                                          ::i,j,k,m
        do m=1,nq**3
            do i=1,n+1
                do j=1,n+1
                    do k=1,n+1



    subroutine calculateEuler3DFlux(u,dir,result)
    !Subroutine gets the values of u and return the flux for the spcified direction
    !dir=1,2,3 stands for x,y,z direction
        implicit none
        REAL(KIND=RP),dimension(:,:,:,5),intent(in)     :: u
        INTEGER                         ,intent(in)     :: dir
        REAL(KIND=RP),dimension(:,:,:,:),intent(out),allocatable    :: result
        !local variables beyond here
        INTEGER                                         :: n
        REAL(KIND=RP),dimension(:,:,:,:),allocatable    ::p
        n=size(u,dim=1)
        ALLOCATE(result(1:n+1,1:n+1,1:n+1,5))
        ALLOCATE(p(1:n+1,1:n+1,1:n+1,5))
        p=(gamma-1.0_RP)*(u(:,:,:,5)-0.5_RP*(u(:,:,:,2)**2+u(:,:,:,3)**2+u(:,:,:,4)**2)/u(:,:,:,1)
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
        REAL(KIND=RP),intent(in),dimension(:,5)           ::u2
        CHARACTER(len=*),intent(in)                     ::whichflux
        INTEGER         ,intent(in)                     ::dir
        REAL(KIND=RP),intent(out),dimension(:),allocatable   :: result
        !local variables go beyond here
        INTEGER                                         ::n
        REAL(KIND=RP),dimension(:),allocatable                   ::p2
        REAL(KIND=RP)                               ::p1
        n=size(u2,dim=1)
        allocate(result(n))
        SELECT CASE(whichflux)
        p1=(gamma-1.0_RP)*(u1(5)-0.5_RP*(u1(2)**2+u1(3)**2+u1(4)**2)/u1(1)
        p2=(gamma-1.0_RP)*(u2(:,5)-0.5_RP*(u2(:,2)**2+u2(:,3)**2+u2(:,4)**2)/u2(:,1)
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
    subroutine Initialcondition(x,n,nq,u,x2,y)
        implicit none

    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine computeLocalLaxFriedrich(uL,uR)
        IMPLICIT NONE
        REAL(KIND=RP)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine RungeKutta5explizit(ustar,t,nq,n,numvar,dt,Dval,Sval)
        IMPLICIT NONE
        !ustar=bekannte Werte
        INTEGER,INTENT(IN)                                          :: numvar,n,nq
        REAL(KIND=RP),intent(in),dimension(:,:)         :: Dval,Sval
        REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:nq*nq,1:numvar,1:(n+1),1:n+1,1:n+1)  :: ustar
        REAL(KIND=RP),INTENT(IN)                                    :: t
        !local
        REAL(KIND=RP),DIMENSION(1:nq*nq,1:numvar,1:(n+1),1:n+1,1:n+1)     :: g
        INTEGER                                                     :: step=1
        REAL(KIND=RP),DIMENSION(5)                                  :: a,b,c
        REAL(KIND=RP),INTENT(IN)                                    :: dt
        g=R(ustar,t,n,nq,Dval,Sval)
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
            g=a(step)*g+R(ustar,t+b(step)*dt,n,nq,Dval,Sval)
            ustar=ustar+c(step)*dt*g
        enddo
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module
