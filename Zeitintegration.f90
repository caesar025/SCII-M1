module Zeitintegration
    use Quadraturroutinen
    REAL(KIND=RP),DIMENSION(:),allocatable      ::x,w,xmit,xges
    REAL(KIND=RP)                               ::gk=9.812_RP,dx,gamma=1.4_RP
contains
    subroutine Vorbereiten(n,nq,Dval)
        implicit none
        INTEGER,INTENT(IN)      ::n,nq
        INTEGER                 ::i
        REAL(KIND=RP),DIMENSION(1:n+1,1:n+1),INTENT(out)    :: Dval
        allocate(x(1:n+1),w(1:n+1),xges(1:nq*(n+1)),xmit(1:nq+1))
        call LegendreGaussLobattoNodesandWeights(N,x,w)
        call DifferentiationsmatrixBerechnen(x,Dval,N+1)
        dx=1.0_RP/real(nq,kind=RP)
        call linspace(dx/2.0_RP,1.0_RP-dx/2.0_RP,NQ,xmit)
        do i=1,NQ
            xges((i-1)*(N+1)+1:i*(N+1))=xmit(i)+dx/2.0_RP*x
        enddo
    end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Operator aus Semidiskreterdarstellung
    function R(u,N,NQ,D) result(solution)
        implicit none
        integer, intent(in) :: n,NQ
        real(kind=RP), dimension(1:NQ,1:nq,1:nq,1:n+1,1:(N+1),1:n+1,1:5)                 :: solution
        real(kind=RP), intent(in), dimension(1:nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:5)     :: u
        REAL(KIND=RP),intent(in),dimension(:,:)                             :: D
        REAL(KIND=RP),dimension(:,:,:,:,:,:,:),allocatable                  ::L1,l2,l3
        call computeL(u,D,1,L1)
        call computeL(u,D,2,L2)
        call computeL(u,D,3,l3)
        solution=8.0_RP/(dx**3)*(-0.25_RP*dx**2*l1-0.25_RP*dx**2*l2-0.25_RP*dx**2*l3)

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
        REAL(KIND=RP),dimension(:,:,:),allocatable                        ::FRand0,FRand1,uR,uL
        nq=size(u,dim=1)
        n=size(u,dim=4)-1
        allocate(result(1:nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:5))
        allocate(ul(1:n+1,1:n+1,1:5),ur(1:n+1,1:n+1,1:5))
        select case(dir)
        case(1)
           do o=1,nq
               do l=1,nq
                    do m=1,nq
                    do k=1,n+1
                        do j=1,n+1
                            do i=1,n+1
					call computeFsharp(u(m,l,o,i,j,k,:),u(m,l,o,:,j,k,:),dir,'ST',Fsharp)
                                    do var=1,5 !! besser
                                        result(m,l,o,i,j,k,var)=2*dot_product(D(i,:),Fsharp(:,var))

                                    enddo
                            enddo
                        enddo
                    enddo
                            !Randbedingungen 
                        if(m==1) then
                              uL=u(nq,l,o,n+1,:,:,:)
                        ELSE
                              uL=u(m-1,l,o,n+1,:,:,:)
                        endif
                        if(m==nq) then
                                   uR=u(1,l,o,1,:,:,:)
                        ELSE
                                   uR=u(m+1,l,o,1,:,:,:)
                        endif

                        call computeLocalLaxFriedrich(uL,u(m,l,o,1,:,:,:),dir,0,FRand0)
                        call computeLocalLaxFriedrich(u(m,l,o,n+1,:,:,:),uR,dir,1,FRand1)


                            result(m,l,o,1,:,:,:)=result(m,l,o,1,:,:,:)-FRand0
                            result(m,l,o,n+1,:,:,:)=result(m,l,o,n+1,:,:,:)+FRand1


                enddo
                enddo
           enddo

       case(2)
         do o=1,nq
             do l=1,nq
                  do m=1,nq
                  do k=1,n+1
                      do j=1,n+1
                          do i=1,n+1
					call computeFsharp(u(m,l,o,i,j,k,:),u(m,l,o,i,:,k,:),dir,'ST',Fsharp)
                                    do var=1,5 !! besser
                                        
                                        result(m,l,o,i,j,k,var)=2*dot_product(D(j,:),Fsharp(:,var))
                                    enddo
                            enddo
                      enddo
                  enddo
                            !Randbedingungen 
                            if(l==1) then 
                                   uL=u(m,nq,o,:,n+1,:,:)
                            ELSE
                                   uL=u(m,l-1,o,:,n+1,:,:)
                            endif 
                            if(l==nq) then
                                   uR=u(m,1,o,:,1,:,:)
                            ELSE
                                   uR=u(m,l+1,o,:,1,:,:)
                            endif 
                            call computeLocalLaxFriedrich(uL,u(m,l,o,:,1,:,:),dir,0,FRand0)
                            call computeLocalLaxFriedrich(u(m,l,o,:,n+1,:,:),uR,dir,1,FRand1)



                            result(m,l,o,:,1,:,:)=result(m,l,o,:,1,:,:)-FRand0
                            result(m,l,o,:,n+1,:,:)=result(m,l,o,:,n+1,:,:)+FRand1


                    enddo
                enddo
                enddo
        case(3)
           do o=1,nq
               do l=1,nq
                    do m=1,nq
                    do k=1,n+1
                        do j=1,n+1
                            do i=1,n+1
					call computeFsharp(u(m,l,o,i,j,k,:),u(m,l,o,i,j,:,:),dir,'ST',Fsharp)
                                    do var=1,5 !! besser
                                        
                                        result(m,l,o,i,j,k,var)=2*dot_product(D(k,:),Fsharp(:,var))
                                    enddo
                            enddo
                        enddo
                    enddo
                            !Randbedingungen 
                            if(o==1) then 
                                   uL=u(m,l,nq,:,:,n+1,:)
                            ELSE
                                   uL=u(m,l,o-1,:,:,n+1,:)
                            endif 
                            if(o==nq) then
                                   uR=u(m,l,1,:,:,1,:)
                            ELSE
                                   uR=u(m,l,o+1,:,:,1,:)
                            endif 
                            call computeLocalLaxFriedrich(uL,u(m,l,o,:,:,1,:),dir,0,FRand0)
                            call computeLocalLaxFriedrich(u(m,l,o,:,:,n+1,:),uR,dir,1,FRand1)

                            
                            result(m,l,o,:,:,1,:)=result(m,l,o,:,:,1,:)-FRand0
                            result(m,l,o,:,:,n+1,:)=result(m,l,o,:,:,n+1,:)+FRand1
                            

                    enddo
                enddo
            enddo
        end select
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        n=size(u,dim=1)-1
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
        REAL(KIND=RP),dimension(:),allocatable                   ::p2,c2,h2
        REAL(KIND=RP)                               ::p1,c1,h1
        n=size(u2,dim=1)-1
        allocate(result(1:n+1,5))
        allocate(p2(1:n+1))
        allocate(c2(1:n+1))
        allocate(h2(1:n+1))
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
            CASE('PI')
                !! Pirozzoli
                c1=sqrt((gamma)*p1/(u1(1)))
                c2=sqrt((gamma)*p2/(u2(:,1)))
                h1=c1+p1/u1(1)
                h2=c2+p2/u2(:,1)
                SELECT CASE(dir)
                    CASE(1)
                        result(:,1)=(u1(1)+u2(:,1))*(u1(2)/u1(1)+u2(:,2)/u2(:,1))*0.25_RP
                        result(:,2)=((u1(1)+u2(:,1))*(u1(2)/u1(1)+u2(:,2)/u2(:,1))**2)*0.125_RP+(u1(1)+u2(:,1))*0.5_RP
                        result(:,3)=(u1(1)+u2(:,1))*(u1(2)/u1(1)+u2(:,2)/u2(:,1))*(u1(3)/u1(1)+u2(:,3)/u2(:,1))*0.125_RP
                        result(:,4)=(u1(1)+u2(:,1))*(u1(2)/u1(1)+u2(:,2)/u2(:,1))*(u1(4)/u1(1)+u2(:,4)/u2(:,1))*0.125_RP
                        result(:,5)=(u1(1)+u2(:,1))*(u1(2)/u1(1)+u2(:,2)/u2(:,1))*(h1+h2)*0.125_RP
                    CASE(2)
                        result(:,1)=(u1(1)+u2(:,1))*(u1(3)/u1(1)+u2(:,3)/u2(:,1))*0.25_RP
                        result(:,2)=(u1(1)+u2(:,1))*(u1(2)/u1(1)+u2(:,2)/u2(:,1))*(u1(3)/u1(1)+u2(:,3)/u2(:,1))*0.125_RP
                        result(:,3)=((u1(1)+u2(:,1))*(u1(3)/u1(1)+u2(:,3)/u2(:,1))**2)*0.125_RP+(u1(1)+u2(:,1))*0.5_RP
                        result(:,4)=(u1(1)+u2(:,1))*(u1(3)/u1(1)+u2(:,3)/u2(:,1))*(u1(4)/u1(1)+u2(:,4)/u2(:,1))*0.125_RP
                        result(:,5)=(u1(1)+u2(:,1))*(u1(3)/u1(1)+u2(:,3)/u2(:,1))*(h1+h2)*0.125_RP
                    CASE(3)
                        result(:,1)=(u1(1)+u2(:,1))*(u1(4)/u1(1)+u2(:,4)/u2(:,1))*0.25_RP
                        result(:,2)=(u1(1)+u2(:,1))*(u1(2)/u1(1)+u2(:,2)/u2(:,1))*(u1(4)/u1(1)+u2(:,4)/u2(:,1))*0.125_RP
                        result(:,3)=(u1(1)+u2(:,1))*(u1(3)/u1(1)+u2(:,3)/u2(:,1))*(u1(4)/u1(1)+u2(:,4)/u2(:,1))*0.125_RP
                        result(:,4)=((u1(1)+u2(:,1))*(u1(4)/u1(1)+u2(:,4)/u2(:,1))**2)*0.125_RP+(u1(1)+u2(:,1))*0.5_RP
                        result(:,5)=(u1(1)+u2(:,1))*(u1(4)/u1(1)+u2(:,4)/u2(:,1))*(h1+h2)*0.125_RP
                end SELECT
        END SELECT
    end subroutine
    subroutine Initialcondition(xyz,u)
        implicit none
        REAL(KIND=RP),DIMENSION(:,:,:,:,:,:,:),INTENT(IN),allocatable::xyz
        REAL(KIND=RP),DIMENSION(:,:,:,:,:,:,:),INTENT(INOUT),allocatable::u
        u(:,:,:,:,:,:,1)=2.0_rp+sin(pi*(xyz(:,:,:,:,:,:,1)+xyz(:,:,:,:,:,:,2)+xyz(:,:,:,:,:,:,3)))/10.0_rp
        u(:,:,:,:,:,:,2)=1.0_RP
        u(:,:,:,:,:,:,3)=1.0_RP
        u(:,:,:,:,:,:,4)=1.0_RP
        u(:,:,:,:,:,:,5)=u(:,:,:,:,:,:,1)*u(:,:,:,:,:,:,1)
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine computeLocalLaxFriedrich(uL,uR,dir,pos,result)
    ! Compute loacal LaxFriedirch 
    ! Returns Matrix with the localLaxFriedrich at every Interface point
         implicit none
         REAL(KIND=RP), INTENT(IN), DIMENSION(:,:,:) :: uL, uR
         REAL(KIND=RP), INTENT(out), DIMENSION(:,:,:),allocatable:: result
         INTEGER, INTENT(IN):: dir,pos 
         REAL(KIND=RP), DIMENSION(:,:,:),allocatable :: FL,FR 
         REAL(KIND=RP), DIMENSION(:,:),allocatable :: lamMax
         INTEGER :: n 

        n=size(uL,dim=1)-1
        ALLOCATE(result(1:n+1,1:n+1,5))
        ALLOCATE(FR(1:n+1,1:n+1,5))
        ALLOCATE(FL(1:n+1,1:n+1,5))
        ALLOCATE(lamMax(1:n+1,1:n+1))

         call calculateEulerRandFlux(uL,dir,FL)
         call calculateEulerRandFlux(uR,dir,FR)
         call lambdaMax(uL,uR,dir,lamMax)
         SELECT CASE(pos)
         CASE(0)
                 result(:,:,1)=(FL(:,:,1)+FR(:,:,1)-lamMax*(uR(:,:,1)-uL(:,:,1)))/2-FR(:,:,1)
                 result(:,:,2)=(FL(:,:,2)+FR(:,:,2)-lamMax*(uR(:,:,2)-uL(:,:,2)))/2-FR(:,:,2)
                 result(:,:,3)=(FL(:,:,3)+FR(:,:,3)-lamMax*(uR(:,:,3)-uL(:,:,3)))/2-FR(:,:,3)
                 result(:,:,4)=(FL(:,:,4)+FR(:,:,4)-lamMax*(uR(:,:,4)-uL(:,:,4)))/2-FR(:,:,4)
                 result(:,:,5)=(FL(:,:,5)+FR(:,:,5)-lamMax*(uR(:,:,5)-uL(:,:,5)))/2-FR(:,:,5)
         CASE(1)

                 result(:,:,1)=(FL(:,:,1)+FR(:,:,1)-lamMax*(uR(:,:,1)-uL(:,:,1)))/2-FL(:,:,1)
                 result(:,:,2)=(FL(:,:,2)+FR(:,:,2)-lamMax*(uR(:,:,2)-uL(:,:,2)))/2-FL(:,:,2)
                 result(:,:,3)=(FL(:,:,3)+FR(:,:,3)-lamMax*(uR(:,:,3)-uL(:,:,3)))/2-FL(:,:,3)
                 result(:,:,4)=(FL(:,:,4)+FR(:,:,4)-lamMax*(uR(:,:,4)-uL(:,:,4)))/2-FL(:,:,4)
                 result(:,:,5)=(FL(:,:,5)+FR(:,:,5)-lamMax*(uR(:,:,5)-uL(:,:,5)))/2-FL(:,:,5)
       END SELECT
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

        n=size(uL,dim=1)-1


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
    subroutine lambdaMaxGlobal(u,lambdamax)
    !computes the max eigenvalue glabally for the timestep
    implicit none
        REAL(KIND=RP),DIMENSION(:,:,:,:,:,:,:),INTENT(IN)           ::u
        REAL(KIND=RP)                         ,INTENT(OUT)          ::lambdamax
        REAL(KIND=RP),DIMENSION(:,:,:,:,:,:),allocatable          ::p,c
        INTEGER                                                     ::n,nq
        !local variables beyond here
        nq=size(u,dim=1)
        n=size(u,dim=4)-1
        allocate(p(1:nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1),c(1:nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1))
        p=(gamma-1.0_RP)*(u(:,:,:,:,:,:,5)-0.5_RP*(u(:,:,:,:,:,:,2)**2+u(:,:,:,:,:,:,3)**2+u(:,:,:,:,:,:,4)**2)/u(:,:,:,:,:,:,1))
        c=sqrt((gamma)*p/(u(:,:,:,:,:,:,1)))

        !! max from abs(eigenvalue)
        lambdamax=max(maxval(abs(u(:,:,:,:,:,:,2)/u(:,:,:,:,:,:,1)+c)),maxval(abs(u(:,:,:,:,:,:,3)/u(:,:,:,:,:,:,1)+c)),&
            maxval(abs(u(:,:,:,:,:,:,4)/u(:,:,:,:,:,:,1)+c)),maxval(abs(u(:,:,:,:,:,:,2)/u(:,:,:,:,:,:,1)-c)),&
            maxval(abs(u(:,:,:,:,:,:,3)/u(:,:,:,:,:,:,1)-c)),maxval(abs(u(:,:,:,:,:,:,4)/u(:,:,:,:,:,:,1)-c)))


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
        n=size(u,dim=1)-1
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
    subroutine RungeKutta5explizit(ustar,nq,n,numvar,dt,Dval)
        IMPLICIT NONE
        !ustar=bekannte Werte
        INTEGER,INTENT(IN)                                          :: numvar,n,nq
        REAL(KIND=RP),intent(in),dimension(:,:)         :: Dval
        REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:nq,1:nq,1:nq,1:(n+1),1:n+1,1:n+1,1:numvar)  :: ustar
        !local
        REAL(KIND=RP),DIMENSION(1:nq,1:nq,1:nq,1:(n+1),1:n+1,1:n+1,1:numvar)    :: g
        INTEGER                                                     :: step=1
        REAL(KIND=RP),DIMENSION(5)                                  :: a,b,c
        REAL(KIND=RP),INTENT(IN)                                    :: dt
        g=R(ustar,n,nq,Dval)
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
            g=a(step)*g+R(ustar,n,nq,Dval)
            ustar=ustar+c(step)*dt*g
        enddo
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module
