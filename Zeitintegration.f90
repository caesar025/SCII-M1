MODULE Zeitintegration
  USE DGtoolbox
  REAL(KIND=RP),DIMENSION(:),ALLOCATABLE             :: x,w,xmit,xges
  REAL(KIND=RP)                                      :: gk=9.812_RP,dx,gamma=1.4_RP,mu=0.001_RP,Pr=0.72_RP,Rkonst=287.058_RP
  REAL(KIND=RP),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: xyz

CONTAINS
  !
  SUBROUTINE Vorbereiten(N,NQ,Dval)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                         :: N,NQ
    REAL(KIND=RP),INTENT(OUT),DIMENSION(1:N+1,1:N+1) :: Dval
    REAL(KIND=RP),DIMENSION(:),ALLOCATABLE           :: xi,xl                 ! Hilfsvariablen
    REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE         :: xin                 ! Hilfsvariablen
    ! local variables
    INTEGER                                          :: l,m,o,k,i,j                           ! Laufvariablen
    !
    allocate(xi(1:n+1),xl(1:nq))
    allocate(xin(1:n+1,1:nq))
    ALLOCATE(x(1:N+1),w(1:N+1),xges(1:NQ*(N+1)),xmit(1:NQ+1))
    allocate(xyz(1:nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:3))
    CALL LegendreGaussLobattoNodesAndWeights(N,x,w)
    Dval= baryzdiffmatrix(x,n)
    dx=1.0_RP/REAL(nq,KIND=RP)
    CALL linspace(dx/2.0_RP,1.0_RP-dx/2.0_RP,NQ,xmit)
    DO i=1,NQ
      xges((i-1)*(N+1)+1:i*(N+1))=xmit(i)+dx/2.0_RP*x
    ENDDO ! i

    CALL LegendreGaussLobattoNodesAndWeights(N,xi,w)
    !! Bestimme GL punkte in jeder Zelle
    DO k=0,nq-1
      xl(k+1)=(k+1.0_rp/2)*dx
      DO i=1,n+1
        xin(i,k+1)=xl(k+1)+dx/2*xi(i)
      END DO
    END DO
    !! Bestimme alle punkte.  
    DO o=1,nq
      DO l=1,nq
        DO m=1,nq
          DO k=1,n+1
            DO j=1,n+1
              DO i=1,n+1
                xyz(m,l,o,i,j,k,1)=xin(i,m)
                xyz(m,l,o,i,j,k,2)=xin(j,l)
                xyz(m,l,o,i,j,k,3)=xin(k,o)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !!! Speicher t intern 
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Operator aus Semidiskreterdarstellung
  FUNCTION R(u,N,NQ,D,whichflux) result(solution)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                                                 :: N,NQ
    REAL(KIND=RP),INTENT(IN),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: u
    REAL(KIND=RP),INTENT(IN),DIMENSION(1:N+1,1:N+1)                          :: D
    REAL(KIND=RP)           ,DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: solution
    CHARACTER(len=2),INTENT(IN)                                              :: whichflux
    ! local
    REAL(KIND=RP)           ,DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) ::L1,L2,L3
    !
    CALL computeL(u,D,1,L1,N,NQ,whichflux)
    CALL computeL(u,D,2,L2,N,NQ,whichflux)
    CALL computeL(u,D,3,L3,N,NQ,whichflux)
    solution=8.0_RP/(dx**3)*((-0.25_RP*dx**2)*L1-(0.25_RP*dx**2)*L2-(0.25_RP*dx**2)*L3)
  END FUNCTION
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION Rmanu(u,N,NQ,D,t,whichflux) result(solution)
    ! Setze Rechteseite zusammen mit manufactuerd solution
    IMPLICIT NONE
    INTEGER, INTENT(IN)                                                        :: n,NQ
    REAL(KIND=RP), DIMENSION(1:NQ,1:nq,1:nq,1:n+1,1:(N+1),1:n+1,1:5)           :: solution    ! Rechte Seite
    REAL(KIND=RP), DIMENSION(1:NQ,1:nq,1:nq,1:n+1,1:(N+1),1:n+1,1:5)           :: res                       ! Quellterme
    REAL(KIND=RP), INTENT(IN), DIMENSION(1:nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:5) :: u             ! U
    REAL(KIND=RP),INTENT(IN),DIMENSION(:,:)                                    :: D                        ! Diff-Matrix
    CHARACTER(len=2),INTENT(IN)                                                :: whichflux
    REAL(KIND=RP),INTENT(IN)                                                   :: t                    ! Startzeit
    REAL(KIND=RP)           ,DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5)   :: L1,L2,L3
    CALL computeL(u,D,1,L1,N,NQ,whichflux)
    CALL computeL(u,D,2,L2,N,NQ,whichflux)
    CALL computeL(u,D,3,L3,N,NQ,whichflux)
    call Residuum(NQ,N,t,res)
    solution=8.0_RP/(dx**3)*(-0.25_RP*dx**2*l1-0.25_RP*dx**2*l2-0.25_RP*dx**2*l3)
    solution=solution+res

  END FUNCTION
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Residuum (NQ,N,t,result)
    ! Berechne Quellterme
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                                        :: n,NQ
    REAL(KIND=RP),INTENT(IN)                                        :: t                    ! Startzeit
    REAL(KIND=RP),DIMENSION(1:NQ,1:nq,1:nq,1:n+1,1:(N+1),1:n+1,1:5) :: result                       ! Quellterme
    REAL(KIND=RP)                                                   :: c1,c2,c3,c4,c5               ! Hilfsvariablen
    INTEGER                                                         :: o,l,m,k,j,i


    c1=pi/10.0_RP
    c2=-1.0_RP/5.0_RP*pi+pi/20.0_rp*(1.0_rp+5.0_RP*gamma)
    c3=pi/100.0_RP*(gamma-1)
    c4=(-16.0_RP*pi+pi*(9.0_RP+15.0_RP*gamma))/20.0_RP
    c5=(3*pi*gamma-2*pi)/100.0_RP

    DO o=1,nq
      DO l=1,nq
        DO m=1,nq
          DO k=1,n+1
            DO j=1,n+1
              DO i=1,n+1
                result(m,l,o,i,j,k,1)=c1*cos(pi*(xyz(m,l,o,i,j,k,1)+xyz(m,l,o,i,j,k,2)+xyz(m,l,o,i,j,k,3)-2.0_RP*t))    
                result(m,l,o,i,j,k,2)=c2*cos(pi*(xyz(m,l,o,i,j,k,1)+xyz(m,l,o,i,j,k,2)+xyz(m,l,o,i,j,k,3)-2.0_RP*t))&
                  +c3*cos(2.0_RP*pi*(xyz(m,l,o,i,j,k,1)+xyz(m,l,o,i,j,k,2)+xyz(m,l,o,i,j,k,3)-2.0_RP*t))    
                result(m,l,o,i,j,k,3)=result(m,l,o,i,j,k,2)
                result(m,l,o,i,j,k,4)=result(m,l,o,i,j,k,2)
                result(m,l,o,i,j,k,5)=c4*cos(pi*(xyz(m,l,o,i,j,k,1)+xyz(m,l,o,i,j,k,2)+xyz(m,l,o,i,j,k,3)-2.0_RP*t))&
                  +c5*cos(2.0_RP*pi*(xyz(m,l,o,i,j,k,1)+xyz(m,l,o,i,j,k,2)+xyz(m,l,o,i,j,k,3)-2.0_RP*t))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO





  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeL(u,D,dir,result,N,NQ,whichflux)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                                                  :: N,NQ
    INTEGER      ,INTENT(IN)                                                  :: dir
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: u
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:N+1,1:N+1)                          :: D
    CHARACTER(len=2),INTENT(IN)                                               :: whichflux
    !u DIMENSIONs:(nummer zelle x,zelle y,zelle z,x,y,z,variable)
    REAL(KIND=RP),INTENT(OUT),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: result
    !local Variables
    INTEGER                                  :: var,k,j,i,o,l,m
    REAL(KIND=RP),DIMENSION(1:N+1,1:5)       :: Fsharp
    REAL(KIND=RP),DIMENSION(1:N+1,1:N+1,1:5) :: FRand0,FRand1,uR,uL
    SELECT CASE(dir)
    CASE(1)
      DO o=1,NQ
        DO l=1,NQ
          DO m=1,NQ
            DO k=1,N+1
              DO j=1,N+1
                DO i=1,N+1
                  CALL computeFsharp(u(m,l,o,i,j,k,:),u(m,l,o,:,j,k,:),dir,whichflux,Fsharp,N)
                  DO var=1,5 !! besser
                    result(m,l,o,i,j,k,var)=2.0_RP*DOt_product(D(i,:),Fsharp(:,var))
                  ENDDO ! var
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
            !Randbedingungen
            IF (m==1) THEN
              uL=u(nq,l,o,N+1,:,:,:)
            ELSE
              uL=u(m-1,l,o,N+1,:,:,:)
            ENDIF
            IF (m==NQ) THEN
              uR=u(1,l,o,1,:,:,:)
            ELSE
              uR=u(m+1,l,o,1,:,:,:)
            ENDIF
            CALL computeLocalLaxFriedrich(uL,u(m,l,o,1,:,:,:),dir,0,whichflux,FRand0,N)
            CALL computeLocalLaxFriedrich(u(m,l,o,N+1,:,:,:),uR,dir,1,whichflux,FRand1,N)
            result(m,l,o,1,:,:,:)=result(m,l,o,1,:,:,:)-FRand0
            result(m,l,o,N+1,:,:,:)=result(m,l,o,N+1,:,:,:)+FRand1
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
    CASE(2)
      DO o=1,NQ
        DO l=1,NQ
          DO m=1,NQ
            DO k=1,N+1
              DO j=1,N+1
                DO i=1,N+1
                  CALL computeFsharp(u(m,l,o,i,j,k,:),u(m,l,o,i,:,k,:),dir,whichflux,Fsharp,N)
                  DO var=1,5 !! besser
                    result(m,l,o,i,j,k,var)=2.0_RP*DOt_product(D(j,:),Fsharp(:,var))
                  ENDDO ! var
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
            !Randbedingungen
            IF (l==1) THEN
              uL=u(m,nq,o,:,N+1,:,:)
            ELSE
              uL=u(m,l-1,o,:,N+1,:,:)
            ENDIF
            IF (l==nq) THEN
              uR=u(m,1,o,:,1,:,:)
            ELSE
              uR=u(m,l+1,o,:,1,:,:)
            ENDIF
            CALL computeLocalLaxFriedrich(uL,u(m,l,o,:,1,:,:),dir,0,whichflux,FRand0,N)
            CALL computeLocalLaxFriedrich(u(m,l,o,:,N+1,:,:),uR,dir,1,whichflux,FRand1,N)
            result(m,l,o,:,1,:,:)=result(m,l,o,:,1,:,:)-FRand0
            result(m,l,o,:,N+1,:,:)=result(m,l,o,:,N+1,:,:)+FRand1
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
    CASE(3)
      DO o=1,NQ
        DO l=1,NQ
          DO m=1,NQ
            DO k=1,N+1
              DO j=1,N+1
                DO i=1,N+1
                  CALL computeFsharp(u(m,l,o,i,j,k,:),u(m,l,o,i,j,:,:),dir,whichflux,Fsharp,N)
                  DO var=1,5 !! besser
                    result(m,l,o,i,j,k,var)=2.0_RP*DOt_product(D(k,:),Fsharp(:,var))
                  ENDDO ! var
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
            !Randbedingungen
            IF (o==1) THEN
              uL=u(m,l,nq,:,:,N+1,:)
            ELSE
              uL=u(m,l,o-1,:,:,N+1,:)
            ENDIF
            IF (o==nq) THEN
              uR=u(m,l,1,:,:,1,:)
            ELSE
              uR=u(m,l,o+1,:,:,1,:)
            ENDIF
            CALL computeLocalLaxFriedrich(uL,u(m,l,o,:,:,1,:),dir,0,whichflux,FRand0,N)
            CALL computeLocalLaxFriedrich(u(m,l,o,:,:,N+1,:),uR,dir,1,whichflux,FRand1,N)
            result(m,l,o,:,:,1,:)=result(m,l,o,:,:,1,:)-FRand0
            result(m,l,o,:,:,N+1,:)=result(m,l,o,:,:,N+1,:)+FRand1
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
    END SELECT
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    SUBROUTINE calculateEuler3DFlux(u,dir,result)
  !    !SUBROUTINE gets the values of u and return the flux for the spcified direction
  !    !dir=1,2,3 stands for x,y,z direction
  !      IMPLICIT NONE
  !      REAL(KIND=RP),DIMENSION(:,:,:,:),INTENT(IN)     :: u
  !      INTEGER                         ,INTENT(IN)     :: dir
  !      REAL(KIND=RP),DIMENSION(:,:,:,:),INTENT(OUT),allocatable    :: result
  !      !local variables beyond here
  !      INTEGER                                         :: n
  !      REAL(KIND=RP),DIMENSION(:,:,:),allocatable    ::p
  !      n=size(u,dim=1)-1
  !      ALLOCATE(result(1:N+1,1:N+1,1:N+1,5))
  !      ALLOCATE(p(1:N+1,1:N+1,1:N+1))
  !      p=(gamma-1.0_RP)*(u(:,:,:,5)-0.5_RP*(u(:,:,:,2)**2+u(:,:,:,3)**2+u(:,:,:,4)**2)/u(:,:,:,1))
  !      SELECT CASE (dir)
  !        CASE(1)
  !          result(:,:,:,1)=u(:,:,:,2)
  !          result(:,:,:,2)=u(:,:,:,2)**2/u(:,:,:,1)+p
  !          result(:,:,:,3)=u(:,:,:,2)*u(:,:,:,3)/u(:,:,:,1)
  !          result(:,:,:,4)=u(:,:,:,2)*u(:,:,:,4)/u(:,:,:,1)
  !          result(:,:,:,5)=(u(:,:,:,5)+p)*u(:,:,:,2)/u(:,:,:,1)
  !        CASE(2)
  !          result(:,:,:,1)=u(:,:,:,3)
  !          result(:,:,:,2)=u(:,:,:,2)*u(:,:,:,3)/u(:,:,:,1)
  !          result(:,:,:,3)=u(:,:,:,3)**2/u(:,:,:,1)+p
  !          result(:,:,:,4)=u(:,:,:,3)*u(:,:,:,4)/u(:,:,:,1)
  !          result(:,:,:,5)=(u(:,:,:,5)+p)*u(:,:,:,3)/u(:,:,:,1)
  !        CASE(3)
  !          result(:,:,:,1)=u(:,:,:,4)
  !          result(:,:,:,2)=u(:,:,:,2)*u(:,:,:,4)/u(:,:,:,1)
  !          result(:,:,:,3)=u(:,:,:,3)*u(:,:,:,4)/u(:,:,:,1)
  !          result(:,:,:,4)=u(:,:,:,4)**2/u(:,:,:,1)+p
  !          result(:,:,:,5)=(u(:,:,:,5)+p)*u(:,:,:,4)/u(:,:,:,1)
  !        CASE DEFAULT
  !          PRINT*,'Specify 1,2 or 3'
  !      END SELECT
  !    END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeviscosFlux(u,dv,dir,n,result)
    !computes the viscos part of the flux analytically
    IMPLICIT NONE
    INTEGER       ,INTENT(IN)                               :: N,dir
    REAL(KIND=RP) ,INTENT(IN) ,DIMENSION(1:n+1,1:n+1,1:5)   :: u, dv
    REAL(KIND=RP) ,INTENT(OUT),DIMENSION(1:n+1,1:n+1,1:5)   :: result
    REAL(KIND=RP)             ,DIMENSION(1:N+1,1:N+1)       ::P
    result(:,:,1)=0.0_RP
    p=(gamma-1.0_RP)*(u(:,:,5)-0.5_RP*(u(:,:,2)*u(:,:,2)+u(:,:,3)*u(:,:,3)+u(:,:,4)*u(:,:,4))/u(:,:,1))
    SELECT CASE(dir)
    CASE(1)
      result(:,:,2)=mu*4.0_RP/3.0_RP*dv(:,:,1)
      result(:,:,3)=mu*dv(:,:,1)
      result(:,:,4)=mu*dv(:,:,1)
      result(:,:,5)=mu*(dv(:,:,1)*(4.0_RP/3.0_RP*u(:,:,2)/u(:,:,1)+u(:,:,3)/u(:,:,1)+u(:,:,4)/u(:,:,1))+&
        p/(u(:,:,1)*Rkonst))
    CASE(2)
      result(:,:,3)=mu*4.0_RP/3.0_RP*dv(:,:,2)
      result(:,:,2)=mu*dv(:,:,2)
      result(:,:,4)=mu*dv(:,:,2)
      result(:,:,5)=mu*(dv(:,:,2)*(4.0_RP/3.0_RP*u(:,:,3)/u(:,:,1)+u(:,:,2)/u(:,:,1)+u(:,:,4)/u(:,:,1))+&
        p/(u(:,:,1)*Rkonst))
    CASE(3)
      result(:,:,4)=mu*4.0_RP/3.0_RP*dv(:,:,3)
      result(:,:,2)=mu*dv(:,:,3)
      result(:,:,3)=mu*dv(:,:,3)
      result(:,:,5)=mu*(dv(:,:,3)*(4.0_RP/3.0_RP*u(:,:,4)/u(:,:,1)+u(:,:,3)/u(:,:,1)+u(:,:,2)/u(:,:,1))+&
        p/(u(:,:,1)*Rkonst))
    END SELECT
  END SUBROUTINE
  SUBROUTINE computeFsharp(u1,u2,dir,whichflux,result,N)
    !SUBROUTINE computes the Volume flux Fsharp
    IMPLICIT NONE
    INTEGER         ,INTENT(IN)                       :: N
    REAL(KIND=RP)   ,INTENT(IN) ,DIMENSION(1:5)       :: u1
    REAL(KIND=RP)   ,INTENT(IN) ,DIMENSION(1:N+1,1:5) :: u2
    CHARACTER(len=*),INTENT(IN)                       :: whichflux
    INTEGER         ,INTENT(IN)                       :: dir
    REAL(KIND=RP)   ,INTENT(OUT),DIMENSION(1:N+1,1:5) :: result
    !local variables go beyond here
    REAL(KIND=RP)               ,DIMENSION(1:N+1)     :: p2,h2
    REAL(KIND=RP)                                     :: p1,h1
    !
    p1=(gamma-1.0_RP)*(u1(5)  -0.5_RP*(u1(2)*u1(2)    +u1(3)*u1(3)    +u1(4)*u1(4))    /u1(1)  )
    p2=(gamma-1.0_RP)*(u2(:,5)-0.5_RP*(u2(:,2)*u2(:,2)+u2(:,3)*u2(:,3)+u2(:,4)*u2(:,4))/u2(:,1))
    SELECT CASE(whichflux)
    CASE('ST')
      !! Standard DG
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
      END SELECT
    CASE('PI')
      !! Pirozzoli
      !! h = e + p/rho = (rho e + p)/rho = (u(5)+p)/u(1)
      h1=(u1(5)  +p1)   /u1(1)
      h2=(u2(:,5)+p2(:))/u2(:,1)
      SELECT CASE(dir)
      CASE(1)
        result(:,1)=(u1(1)+u2(:,1))*(u1(2)/u1(1)+u2(:,2)/u2(:,1))*0.25_RP
        result(:,2)=result(:,1)    *(u1(2)/u1(1)+u2(:,2)/u2(:,1))*0.5_RP  +(p1+p2)*0.5_RP
        result(:,3)=result(:,1)    *(u1(3)/u1(1)+u2(:,3)/u2(:,1))*0.5_RP
        result(:,4)=result(:,1)    *(u1(4)/u1(1)+u2(:,4)/u2(:,1))*0.5_RP
        result(:,5)=result(:,1)    *(h1+h2)*0.5_RP
      CASE(2)
        result(:,1)=(u1(1)+u2(:,1))*(u1(3)/u1(1)+u2(:,3)/u2(:,1))*0.25_RP
        result(:,2)=result(:,1)    *(u1(2)/u1(1)+u2(:,2)/u2(:,1))*0.5_RP
        result(:,3)=result(:,1)    *(u1(3)/u1(1)+u2(:,3)/u2(:,1))*0.5_RP  +(p1+p2)*0.5_RP
        result(:,4)=result(:,1)    *(u1(4)/u1(1)+u2(:,4)/u2(:,1))*0.5_RP
        result(:,5)=result(:,1)    *(h1+h2)*0.5_RP
      CASE(3)
        result(:,1)=(u1(1)+u2(:,1))*(u1(4)/u1(1)+u2(:,4)/u2(:,1))*0.25_RP
        result(:,2)=result(:,1)    *(u1(2)/u1(1)+u2(:,2)/u2(:,1))*0.5_RP
        result(:,3)=result(:,1)    *(u1(3)/u1(1)+u2(:,3)/u2(:,1))*0.5_RP
        result(:,4)=result(:,1)    *(u1(4)/u1(1)+u2(:,4)/u2(:,1))*0.5_RP  +(p1+p2)*0.5_RP
        result(:,5)=result(:,1)    *(h1+h2)*0.5_RP
      END SELECT
    END SELECT
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Initialcondition(u,NQ,N)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                                                    :: NQ,N
    REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: u
    u(:,:,:,:,:,:,1)=2.0_RP+SIN(pi*(xyz(:,:,:,:,:,:,1)+xyz(:,:,:,:,:,:,2)+xyz(:,:,:,:,:,:,3)))/10.0_RP
    u(:,:,:,:,:,:,2)=1.0_RP
    u(:,:,:,:,:,:,3)=1.0_RP
    u(:,:,:,:,:,:,4)=1.0_RP
    u(:,:,:,:,:,:,5)=u(:,:,:,:,:,:,1)*u(:,:,:,:,:,:,1)
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeSolution(u,NQ,N,t)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                                                    :: NQ,N
    REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: u
    REAL(KIND=RP),INTENT(IN)                                                :: t
    u(:,:,:,:,:,:,1)=2.0_rp+sin(pi*(xyz(:,:,:,:,:,:,1)+xyz(:,:,:,:,:,:,2)+xyz(:,:,:,:,:,:,3)-2*t))/10.0_rp
    u(:,:,:,:,:,:,2)=1.0_RP
    u(:,:,:,:,:,:,3)=1.0_RP
    u(:,:,:,:,:,:,4)=1.0_RP
    u(:,:,:,:,:,:,5)=u(:,:,:,:,:,:,1)*u(:,:,:,:,:,:,1)
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeLocalLaxFriedrich(uL,uR,dir,pos,whichflux,result,N)
    ! Compute local LaxFriedrich
    ! Returns Matrix with the localLaxFriedrich at every Interface point
    IMPLICIT NONE
    INTEGER         ,INTENT(IN)                             :: N,dir,pos
    CHARACTER(len=2),INTENT(IN)                             :: whichflux
    REAL(KIND=RP)   ,INTENT(IN) ,DIMENSION(1:N+1,1:N+1,1:5) :: uL,uR
    REAL(KIND=RP)   ,INTENT(OUT),DIMENSION(1:N+1,1:N+1,1:5) :: result
    ! local
    REAL(KIND=RP),DIMENSION(1:N+1,1:N+1,1:5) :: FL,FR,FPI
    REAL(KIND=RP),DIMENSION(1:N+1,1:N+1)     :: lamMax
    CALL calculateEulerRandFlux(uL,dir,FL,N)
    CALL calculateEulerRandFlux(uR,dir,FR,N)
    SELECT CASE(whichflux)
    CASE('ST')

      CALL lambdaMax(uL,uR,dir,lamMax,N)
      SELECT CASE(pos)
      CASE(0)
        result(:,:,1)=(FL(:,:,1)+FR(:,:,1)-lamMax*(uR(:,:,1)-uL(:,:,1)))/2.0_RP-FR(:,:,1)
        result(:,:,2)=(FL(:,:,2)+FR(:,:,2)-lamMax*(uR(:,:,2)-uL(:,:,2)))/2.0_RP-FR(:,:,2)
        result(:,:,3)=(FL(:,:,3)+FR(:,:,3)-lamMax*(uR(:,:,3)-uL(:,:,3)))/2.0_RP-FR(:,:,3)
        result(:,:,4)=(FL(:,:,4)+FR(:,:,4)-lamMax*(uR(:,:,4)-uL(:,:,4)))/2.0_RP-FR(:,:,4)
        result(:,:,5)=(FL(:,:,5)+FR(:,:,5)-lamMax*(uR(:,:,5)-uL(:,:,5)))/2.0_RP-FR(:,:,5)
      CASE(1)
        result(:,:,1)=(FL(:,:,1)+FR(:,:,1)-lamMax*(uR(:,:,1)-uL(:,:,1)))/2.0_RP-FL(:,:,1)
        result(:,:,2)=(FL(:,:,2)+FR(:,:,2)-lamMax*(uR(:,:,2)-uL(:,:,2)))/2.0_RP-FL(:,:,2)
        result(:,:,3)=(FL(:,:,3)+FR(:,:,3)-lamMax*(uR(:,:,3)-uL(:,:,3)))/2.0_RP-FL(:,:,3)
        result(:,:,4)=(FL(:,:,4)+FR(:,:,4)-lamMax*(uR(:,:,4)-uL(:,:,4)))/2.0_RP-FL(:,:,4)
        result(:,:,5)=(FL(:,:,5)+FR(:,:,5)-lamMax*(uR(:,:,5)-uL(:,:,5)))/2.0_RP-FL(:,:,5)
      END SELECT
    CASE('PI')
      CALL calculateEulerRandFluxPI(uL,uR,dir,FPI,N)
      CALL lambdaMax(uL,uR,dir,lamMax,N)
      SELECT CASE(pos)
      CASE(0)
        result(:,:,1)=FPI(:,:,1)-lamMax*(uR(:,:,1)-uL(:,:,1))/2-FR(:,:,1)
        result(:,:,2)=FPI(:,:,2)-lamMax*(uR(:,:,2)-uL(:,:,2))/2-FR(:,:,2)
        result(:,:,3)=FPI(:,:,3)-lamMax*(uR(:,:,3)-uL(:,:,3))/2-FR(:,:,3)
        result(:,:,4)=FPI(:,:,4)-lamMax*(uR(:,:,4)-uL(:,:,4))/2-FR(:,:,4)
        result(:,:,5)=FPI(:,:,5)-lamMax*(uR(:,:,5)-uL(:,:,5))/2-FR(:,:,5)
      CASE(1)
        result(:,:,1)=FPI(:,:,1)-lamMax*(uR(:,:,1)-uL(:,:,1))/2-FL(:,:,1)
        result(:,:,2)=FPI(:,:,2)-lamMax*(uR(:,:,2)-uL(:,:,2))/2-FL(:,:,2)
        result(:,:,3)=FPI(:,:,3)-lamMax*(uR(:,:,3)-uL(:,:,3))/2-FL(:,:,3)
        result(:,:,4)=FPI(:,:,4)-lamMax*(uR(:,:,4)-uL(:,:,4))/2-FL(:,:,4)
        result(:,:,5)=FPI(:,:,5)-lamMax*(uR(:,:,5)-uL(:,:,5))/2-FL(:,:,5)
      END SELECT
    END SELECT
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE  calculateEulerRandFluxPI(u1,u2,dir,result,N)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                             :: N,dir
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:N+1,1:N+1,1:5) :: u1,u2
    REAL(KIND=RP),INTENT(OUT),DIMENSION(1:N+1,1:N+1,1:5) :: result
    ! local
    REAL(KIND=RP),DIMENSION(1:N+1,1:N+1) :: h1,h2,p1,p2
    !
    ! Pirozzoli
    p1=(gamma-1.0_RP)*(u1(:,:,5)  -0.5_RP*(u1(:,:,2)*u1(:,:,2)    +u1(:,:,3)*u1(:,:,3)    +u1(:,:,4)*u1(:,:,4))    /u1(:,:,1)  )
    p2=(gamma-1.0_RP)*(u2(:,:,5)-0.5_RP*(u2(:,:,2)*u2(:,:,2)+u2(:,:,3)*u2(:,:,3)+u2(:,:,4)*u2(:,:,4))/u2(:,:,1))
    ! h = e + p/rho = (rho e + p)/rho = (u(5)+p)/u(1)
    h1=(u1(:,:,5)+p1)/u1(:,:,1)
    h2=(u2(:,:,5)+p2)/u2(:,:,1)
    SELECT CASE(dir)
    CASE(1)
      result(:,:,1)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,2)/u1(:,:,1)+u2(:,:,2)/u2(:,:,1))*0.25_RP
      result(:,:,2)=((u1(:,:,1)+u2(:,:,1))*(u1(:,:,2)/u1(:,:,1)+u2(:,:,2)/u2(:,:,1))**2)*0.125_RP+(u1(:,:,1)&
        +u2(:,:,1))*0.5_RP
      result(:,:,3)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,2)/u1(:,:,1)+u2(:,:,2)/u2(:,:,1))*(u1(:,:,3)/u1(:,:,1)&
        +u2(:,:,3)/u2(:,:,1))*0.125_RP
      result(:,:,4)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,2)/u1(:,:,1)+u2(:,:,2)/u2(:,:,1))*(u1(:,:,4)/u1(:,:,1)&
        +u2(:,:,4)/u2(:,:,1))*0.125_RP
      result(:,:,5)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,2)/u1(:,:,1)+u2(:,:,2)/u2(:,:,1))*(h1+h2)*0.125_RP
    CASE(2)
      result(:,:,1)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,3)/u1(:,:,1)+u2(:,:,3)/u2(:,:,1))*0.25_RP
      result(:,:,2)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,2)/u1(:,:,1)+u2(:,:,2)/u2(:,:,1))*(u1(:,:,3)/u1(:,:,1)&
        +u2(:,:,3)/u2(:,:,1))*0.125_RP
      result(:,:,3)=((u1(:,:,1)+u2(:,:,1))*(u1(:,:,3)/u1(:,:,1)+u2(:,:,3)/u2(:,:,1))**2)*0.125_RP&
        +(u1(:,:,1)+u2(:,:,1))*0.5_RP
      result(:,:,4)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,3)/u1(:,:,1)+u2(:,:,3)/u2(:,:,1))*(u1(:,:,4)/u1(:,:,1)&
        +u2(:,:,4)/u2(:,:,1))*0.125_RP
      result(:,:,5)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,3)/u1(:,:,1)+u2(:,:,3)/u2(:,:,1))*(h1+h2)*0.125_RP
    CASE(3)
      result(:,:,1)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,4)/u1(:,:,1)+u2(:,:,4)/u2(:,:,1))*0.25_RP
      result(:,:,2)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,2)/u1(:,:,1)+u2(:,:,2)/u2(:,:,1))*(u1(:,:,4)/u1(:,:,1)&
        +u2(:,:,4)/u2(:,:,1))*0.125_RP
      result(:,:,3)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,3)/u1(:,:,1)+u2(:,:,3)/u2(:,:,1))*(u1(:,:,4)/u1(:,:,1)&
        +u2(:,:,4)/u2(:,:,1))*0.125_RP
      result(:,:,4)=((u1(:,:,1)+u2(:,:,1))*(u1(:,:,4)/u1(:,:,1)+u2(:,:,4)/u2(:,:,1))**2)*0.125_RP&
        +(u1(:,:,1)+u2(:,:,1))*0.5_RP
      result(:,:,5)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,4)/u1(:,:,1)+u2(:,:,4)/u2(:,:,1))*(h1+h2)*0.125_RP
    END SELECT
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE lambdaMax(uL,uR,dir,result,N)
    !Computes the max of the eigenvalues at every Interface
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                             :: dir,N
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:N+1,1:N+1,1:5) :: uR,uL
    REAL(KIND=RP),INTENT(OUT),DIMENSION(1:N+1,1:N+1)     :: result
    ! local variables
    INTEGER                                  :: l,k
    REAL(KIND=RP),DIMENSION(1:N+1,1:N+1)     :: pR,hR,cR,pL,hL,cL
    REAL(KIND=RP),DIMENSION(1:N+1,1:N+1,1:6) :: lambda
    !
    pR=(gamma-1.0_RP)*(uR(:,:,5)-0.5_RP*(uR(:,:,2)**2+uR(:,:,3)**2+uR(:,:,4)**2)/uR(:,:,1))
    ! h = e + p/rho = (rho e + p)/rho = (u(5)+p)/u(1)
    hR=(uR(:,:,5)+pR)/uR(:,:,1)
    ! c = sqrt(gamma p / rho) = sqrt(gamma p / u(1))
    cR=SQRT(gamma*pR/uR(:,:,1))
    !        cR=sqrt((gamma-1)*(hR-((uR(:,:,2)/uR(:,:,1))**2+(uR(:,:,3)/uR(:,:,1))**2+(uR(:,:,4)/uR(:,:,1))**2)/2))
    lambda(:,:,1)=uR(:,:,dir+1)/uR(:,:,1)
    lambda(:,:,2)=uR(:,:,dir+1)/uR(:,:,1)-cR
    lambda(:,:,3)=uR(:,:,dir+1)/uR(:,:,1)+cR
    pL=(gamma-1.0_RP)*(uL(:,:,5)-0.5_RP*(uL(:,:,2)**2+uL(:,:,3)**2+uL(:,:,4)**2)/uL(:,:,1))
    ! h = e + p/rho = (rho e + p)/rho = (u(5)+p)/u(1)
    hL=(uL(:,:,5)+pL)/uL(:,:,1)
    ! c = sqrt(gamma p / rho) = sqrt(gamma p / u(1))
    cL=SQRT(gamma*pL/uL(:,:,1))
    !        cL=sqrt((gamma-1)*(hL-((uL(:,:,2)/uL(:,:,1))**2+(uL(:,:,3)/uL(:,:,1))**2+(uL(:,:,4)/uL(:,:,1))**2)/2))
    lambda(:,:,4)=uL(:,:,dir+1)/uL(:,:,1)
    lambda(:,:,5)=uL(:,:,dir+1)/uL(:,:,1)-cL
    lambda(:,:,6)=uL(:,:,dir+1)/uL(:,:,1)+cL
    !
    lambda=ABS(lambda)
    DO l=1,N+1
      DO k=1,N+1
        result(l,k)=MAXVAL(lambda(l,k,:))
      ENDDO ! k
    ENDDO ! l
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE lambdaMaxGlobal(u,lambdamax,NQ,N)
    !computes the max eigenvalue glabally for the timestep
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                                                 :: NQ,N
    REAL(KIND=RP),INTENT(IN),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: u
    REAL(KIND=RP),INTENT(OUT)                                                :: lambdamax
    ! local variables
    REAL(KIND=RP),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1) :: p,c
    !
    p=(gamma-1.0_RP)*(u(:,:,:,:,:,:,5)-0.5_RP*(u(:,:,:,:,:,:,2)**2+u(:,:,:,:,:,:,3)**2+u(:,:,:,:,:,:,4)**2)/u(:,:,:,:,:,:,1))
    c=sqrt(gamma*p/u(:,:,:,:,:,:,1))
    !! max from abs(eigenvalue)
    lambdamax=max(MAXVAL(ABS(u(:,:,:,:,:,:,2)/u(:,:,:,:,:,:,1)+c)),MAXVAL(ABS(u(:,:,:,:,:,:,3)/u(:,:,:,:,:,:,1)+c)),&
      MAXVAL(ABS(u(:,:,:,:,:,:,4)/u(:,:,:,:,:,:,1)+c)),MAXVAL(ABS(u(:,:,:,:,:,:,2)/u(:,:,:,:,:,:,1)-c)),&
      MAXVAL(abs(u(:,:,:,:,:,:,3)/u(:,:,:,:,:,:,1)-c)),MAXVAL(ABS(u(:,:,:,:,:,:,4)/u(:,:,:,:,:,:,1)-c)))
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE calculateEulerRandFlux(u,dir,result,N)
    IMPLICIT NONE
    !SUBROUTINE gets the values of u and return the flux for the spcified direction
    !dir=1,2,3 stands for x,y,z direction
    INTEGER      ,INTENT(IN)                             :: N,dir
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:N+1,1:N+1,1:5) :: u
    REAL(KIND=RP),INTENT(OUT),DIMENSION(1:N+1,1:N+1,1:5) :: result
    !local variables beyond here
    REAL(KIND=RP),DIMENSION(1:N+1,1:N+1) :: p
    !
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
      PRINT*,'Specify 1,2 or 3'
    END SELECT
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE RungeKutta5explizit(ustar,nq,n,numvar,dt,Dval,t,whichflux)
    IMPLICIT NONE
    !ustar=bekannte Werte
    INTEGER      ,INTENT(IN)                                                         :: numvar,n,nq
    CHARACTER(len=*),INTENT(IN)                                                      :: whichflux
    REAL(KIND=RP),INTENT(IN)   ,DIMENSION(1:N+1,1:N+1)                               :: Dval
    REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:nq,1:nq,1:nq,1:N+1,1:N+1,1:N+1,1:numvar) :: ustar
    !local
    REAL(KIND=RP),DIMENSION(1:nq,1:nq,1:nq,1:(N+1),1:N+1,1:N+1,1:numvar)    :: g
    INTEGER                                                                 :: step
    REAL(KIND=RP),DIMENSION(5)                                              :: a,b,c
    REAL(KIND=RP),INTENT(IN)                                                :: dt,t
    g=Rmanu(ustar,n,nq,Dval,t,whichflux)
    !a(1)=0.0_RP
    !b(1)=0.0_RP
    !c(1)=1432997174477.0_RP/9575080441755.0_RP
    !a(2)=-567301805773.0_RP/1357537059087.0_RP
    !b(2)=1432997174477.0_RP/9575080441755.0_RP
    !c(2)=5161836677717.0_RP/13612068292357.0_RP
    !a(3)=-2404267990393.0_RP/2016746695238.0_RP
    !b(3)=2526269341429.0_RP/6820363962896.0_RP
    !c(3)=1720146321549.0_RP/2090206949498.0_RP
    !a(4)=-3550918686646.0_RP/2091501179385.0_RP
    !b(4)=2006345519317.0_RP/3224310063776.0_RP
    !c(4)=3134564353537.0_RP/4481467310338.0_RP
    !a(5)=-1275806237668.0_RP/842570457699.0_RP
    !b(5)=2802321613138.0_RP/2924317926251.0_RP
    !c(5)=2277821191437.0_RP/14882151754819.0_RP




    a=(/0.0_rp, -567301805773.0_rp/1357537059087.0_rp,&
      -2404267990393.0_rp/2016746695238.0_rp, -3550918686646.0_rp/2091501179385.0_rp,&
      -1275806237668.0_rp/842570457699.0_rp/)
    b=(/0.0_rp, 1432997174477.0_rp/9575080441755.0_rp,&
      2526269341429.0_rp/6820363962896.0_rp, 2006345519317.0_rp/3224310063776.0_rp,&
      2802321613138.0_rp/2924317926251.0_rp /)
    c=(/1432997174477.0_rp/9575080441755.0_rp, 5161836677717.0_rp/13612068292357.0_rp,&
      1720146321549.0_rp/2090206949498.0_rp, 3134564353537.0_rp/4481467310338.0_rp,&
      2277821191437.0_rp/14882151754819.0_rp /)



    DO step=1,5
      g=a(step)*g+Rmanu(ustar,n,nq,Dval,t*b(step)*dt,whichflux)
      ustar=ustar+c(step)*dt*g
    ENDDO ! step
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeError (u,usolution,NQ,N,result)
    IMPLICIT NONE 
    INTEGER, INTENT(IN)                                                      :: NQ,N
    REAL(KIND=RP),INTENT(IN),DIMENSION(1:nq,1:nq,1:nq,1:N+1,1:N+1,1:N+1,1:5) :: u,usolution
    REAL(KIND=RP),INTENT(OUT),DIMENSION(5) :: result
    result(1)=maxval(abs(u(:,:,:,:,:,:,1)-usolution(:,:,:,:,:,:,1)))
    result(2)=maxval(abs(u(:,:,:,:,:,:,2)-usolution(:,:,:,:,:,:,2)))
    result(3)=maxval(abs(u(:,:,:,:,:,:,3)-usolution(:,:,:,:,:,:,3)))
    result(4)=maxval(abs(u(:,:,:,:,:,:,4)-usolution(:,:,:,:,:,:,4)))
    result(5)=maxval(abs(u(:,:,:,:,:,:,5)-usolution(:,:,:,:,:,:,5)))

    deallocate(xyz,x,w,xmit,xges)
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeEOC (errors,n,nq,result)
    IMPLICIT NONE
    REAL(KIND=RP),INTENT(IN),DIMENSION(:,:)  :: errors
    REAL(KIND=RP),INTENT(IN),DIMENSION(:)  :: n,nq
    REAL(KIND=RP),INTENT(OUT),DIMENSION(:,:),ALLOCATABLE  :: result
    INTEGER::k,anz
    anz=size(errors,dim=2)
    allocate(result(1:5,1:anz))
    result(:,:)=0.0_RP
    DO k=1,anz
      result(:,k+1)=log(errors(:,k+1)/errors(:,k))/log(real(nq(k)*(n+1),rp)/((n+1)*nq(k+1)))
    END DO

  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE
