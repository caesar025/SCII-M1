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
    ALLOCATE(xi(1:n+1),xl(1:nq))
    ALLOCATE(xin(1:n+1,1:nq))
    ALLOCATE(x(1:N+1),w(1:N+1),xges(1:NQ*(N+1)),xmit(1:NQ+1))
    ALLOCATE(xyz(1:nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:3))
    CALL LegendreGaussLobattoNodesAndWeights(N,x,w)
    Dval= baryzdiffmatrix(x,n)
    dx=2.0_RP/REAL(nq,KIND=RP)
    !CALL linspace(dx/2.0_RP,1.0_RP-dx/2.0_RP,NQ,xmit)
    !DO i=1,NQ
    !  xges((i-1)*(N+1)+1:i*(N+1))=xmit(i)+dx/2.0_RP*x
    !ENDDO ! i

    CALL LegendreGaussLobattoNodesAndWeights(N,xi,w)
    !! Bestimme GL punkte in jeder Zelle
    DO k=0,nq-1
      xl(k+1)=-1+(k+1.0_rp/2.0_RP)*dx
      DO i=1,n+1
        xin(i,k+1)=xl(k+1)+dx/2.0_RP*xi(i)
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
    !print*, xyz(1,2,1,1,:,1,2)
    !stop
    !!! Speicher t intern
  END SUBROUTINE Vorbereiten
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
  END FUNCTION R
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION Rmanu(u,N,NQ,D,t,whichflux,vis,dt) result(solution)
    ! Setze Rechteseite zusammen mit manufactuerd solution
    IMPLICIT NONE
    INTEGER, INTENT(IN)                                                        :: n,NQ
    REAL(KIND=RP), DIMENSION(1:NQ,1:nq,1:nq,1:n+1,1:(N+1),1:n+1,1:5)           :: solution    ! Rechte Seite
    REAL(KIND=RP), DIMENSION(1:NQ,1:nq,1:nq,1:n+1,1:(N+1),1:n+1,1:5)           :: res                       ! Quellterme
    REAL(KIND=RP), INTENT(IN), DIMENSION(1:nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:5) :: u             ! U
    REAL(KIND=RP),INTENT(IN),DIMENSION(1:N+1,1:N+1)                                    :: D                        ! Diff-Matrix
    CHARACTER(len=2),INTENT(IN)                                                :: whichflux,vis
    REAL(KIND=RP),INTENT(IN)                                                   :: t,dt                    ! Startzeit,zeitschritt
    REAL(KIND=RP)           ,DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5)   :: L1,L2,L3,dux,duy,duz,L1vis,L2vis,L3vis
    CALL computeL(u,D,1,L1,N,NQ,whichflux)
    CALL computeL(u,D,2,L2,N,NQ,whichflux)
    CALL computeL(u,D,3,L3,N,NQ,whichflux)
    SELECT CASE(vis)
    CASE('AD')
    solution=8.0_RP/(dx**3)*(-0.25_RP*dx*dx*l1-0.25_RP*dx*dx*l2-0.25_RP*dx*dx*l3)
    call Residuum(NQ,N,t,vis,res)
    solution=solution!+res !dt*res noch mal ueberpruefen !!!!
    CASE('VI')
    solution=8.0_RP/(dx**3)*(-0.25_RP*dx*dx*l1-0.25_RP*dx*dx*l2-0.25_RP*dx*dx*l3)
    call Residuum(NQ,N,t,vis,res)
    call computeGradient(u,n,nq,D,dux,duy,duz)
    call computeLviscous(u,dux,duy,duz,D,1,N,NQ,L1vis)
    call computeLviscous(u,dux,duy,duz,D,2,N,NQ,L2vis)
    call computeLviscous(u,dux,duy,duz,D,3,N,NQ,L3vis)
    solution=solution+res
    solution=solution-8.0_RP/(dx**3)*(0.25_RP*dx*dx*l1vis+0.25_RP*dx*dx*l2vis+0.25_RP*dx*dx*l3vis)
    END SELECT
  END FUNCTION Rmanu
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Residuum (NQ,N,t,vis,result)
    ! Berechne Quellterme
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                                        :: n,NQ
    REAL(KIND=RP),INTENT(IN)                                        :: t                    ! Startzeit
    REAL(KIND=RP),DIMENSION(1:NQ,1:nq,1:nq,1:n+1,1:(N+1),1:n+1,1:5) :: result                       ! Quellterme
    CHARACTER(len=2),INTENT(IN)                                     :: vis                !Schalter ob viskos oder advektiv
    REAL(KIND=RP)                                                   :: c1,c2,c3,c4,c5,ro,rox,Px,p               ! Hilfsvariablen
    INTEGER                                                         :: o,l,m,k,j,i

    c1=pi/10.0_RP
    c2=-1.0_RP/5.0_RP*pi+pi/20.0_rp*(1.0_rp+5.0_RP*gamma)
    c3=pi/100.0_RP*(gamma-1)
    c4=(-16.0_RP*pi+pi*(9.0_RP+15.0_RP*gamma))/20.0_RP
    c5=(3*pi*gamma-2.0*pi)/100.0_RP

    DO o=1,nq
      DO l=1,nq
        DO m=1,nq
          DO k=1,n+1
            DO j=1,n+1
              DO i=1,n+1
                ro=2.0_RP+1.0_RP/10.0_RP*sin(pi*(xyz(m,l,o,i,j,k,1)+xyz(m,l,o,i,j,k,2)+xyz(m,l,o,i,j,k,3)-2.0_RP*t))
                rox=cos(pi*(xyz(m,l,o,i,j,k,1)+xyz(m,l,o,i,j,k,2)+xyz(m,l,o,i,j,k,3)-2.0_RP*t))*pi/10.0_RP
                Px=(gamma-1.0_RP)*((2.0_RP*ro-3.0_RP/2.0_RP)*rox)
                P=(gamma-1.0_RP)*(ro**2-3.0_RP/2.0_RP*ro)
                result(m,l,o,i,j,k,1)=rox
                result(m,l,o,i,j,k,2)=Px+rox
                result(m,l,o,i,j,k,3)=result(m,l,o,i,j,k,2)
                result(m,l,o,i,j,k,4)=result(m,l,o,i,j,k,2)
                result(m,l,o,i,j,k,5)=2.0_RP*ro*rox+3.0_RP*Px
               ! result(m,l,o,i,j,k,1)=c1*cos(pi*(xyz(m,l,o,i,j,k,1)+xyz(m,l,o,i,j,k,2)+xyz(m,l,o,i,j,k,3)-2.0_RP*t))
               ! result(m,l,o,i,j,k,2)=c2*cos(pi*(xyz(m,l,o,i,j,k,1)+xyz(m,l,o,i,j,k,2)+xyz(m,l,o,i,j,k,3)-2.0_RP*t))&
               !   +c3*cos(2.0_RP*pi*(xyz(m,l,o,i,j,k,1)+xyz(m,l,o,i,j,k,2)+xyz(m,l,o,i,j,k,3)-2.0_RP*t))
               ! result(m,l,o,i,j,k,3)=result(m,l,o,i,j,k,2)
               ! result(m,l,o,i,j,k,4)=result(m,l,o,i,j,k,2)
               ! result(m,l,o,i,j,k,5)=c4*cos(pi*(xyz(m,l,o,i,j,k,1)+xyz(m,l,o,i,j,k,2)+xyz(m,l,o,i,j,k,3)-2.0_RP*t))&
               !   +c5*cos(2.0_RP*pi*(xyz(m,l,o,i,j,k,1)+xyz(m,l,o,i,j,k,2)+xyz(m,l,o,i,j,k,3)-2.0_RP*t))
                if (vis=='VI') then
                  result(m,l,o,i,j,k,5)=result(m,l,o,i,j,k,5)+3.0_RP*mu/(Pr*Rkonst)*(Px*ro-rox*p)/(ro**2)

               end if
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE Residuum
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
                    result(m,l,o,i,j,k,var)=2.0_RP*dot_product(D(i,:),Fsharp(:,var))
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
                    result(m,l,o,i,j,k,var)=2.0_RP*dot_product(D(j,:),Fsharp(:,var))
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
                    result(m,l,o,i,j,k,var)=2.0_RP*dot_product(D(k,:),Fsharp(:,var))
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
  END SUBROUTINE computeL
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
  SUBROUTINE computeGradient(u,n,nq,D,dux,duy,duz)
    IMPLICIT NONE
    !Gleichung 57 im skript
    INTEGER      ,INTENT(IN)                                                        :: n,nq
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:5)       :: u
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:n+1,1:n+1)                                :: D
    REAL(KIND=RP),INTENT(OUT),DIMENSION(1:nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:5)       :: dux,duy,duz
    !local variables
    INTEGER                                                                         :: i,j,k,l,m,o,var
    if (.not.allocated(w)) then
        print*,'w not allocated'
        stop
    endif
    do m=1,nq
        do l=1,nq
            do o=1,nq

                do i=1,n+1
                    do j=1,n+1
                        do k=1,n+1
                            do var=1,5
! TODO (dorn_ni#1#): maybe put dux duy duz together but than rank >7. check for compiler version on cheops

                                    dux(m,l,o,i,j,k,var)=-dot_product(D(i,:),u(m,l,o,:,j,k,var)*w)*1.0_RP/w(i)*2.0_RP/dx!+surface term
                                    duy(m,l,o,i,j,k,var)=-dot_product(D(j,:),u(m,l,o,i,:,k,var)*w)*1.0_RP/w(j)*2.0_RP/dx!+surface term
                                    duz(m,l,o,i,j,k,var)=-dot_product(D(k,:),u(m,l,o,i,j,:,var)*w)*1.0_RP/w(k)*2.0_RP/dx!+surface term
                            end do !var
                        end do !k
                    end do !j
                end do !i
                !surfaceterms and boundary conditions
                !x-direction
                if(m==1) then
                dux(m,l,o,1,:,:,:)=(u(m,l,o,N+1,:,:,:)+u(m+1,l,o,1,:,:,:))*0.5_RP+dux(m,l,o,1,:,:,:)
                dux(m,l,o,N+1,:,:,:)=-(u(nq,l,o,N+1,:,:,:)+u(m,l,o,1,:,:,:))*0.5_RP+dux(m,l,o,N+1,:,:,:)
                elseif(m==nq) then
                dux(m,l,o,1,:,:,:)=(u(m,l,o,N+1,:,:,:)+u(1,l,o,1,:,:,:))*0.5_RP+dux(m,l,o,1,:,:,:)
                dux(m,l,o,N+1,:,:,:)=-(u(m-1,l,o,N+1,:,:,:)+u(m,l,o,1,:,:,:))*0.5_RP+dux(m,l,o,N+1,:,:,:)
                else
                dux(m,l,o,1,:,:,:)=(u(m,l,o,N+1,:,:,:)+u(m+1,l,o,1,:,:,:))*0.5_RP+dux(m,l,o,1,:,:,:)
                dux(m,l,o,N+1,:,:,:)=-(u(m-1,l,o,N+1,:,:,:)+u(m,l,o,1,:,:,:))*0.5_RP+dux(m,l,o,N+1,:,:,:)
                end if
                !y-direction
                if(l==1) then
                duy(m,l,o,:,1,:,:)=(u(m,l,o,:,N+1,:,:)+u(m,l+1,o,:,1,:,:))*0.5_RP+duy(m,l,o,:,1,:,:)
                duy(m,l,o,:,N+1,:,:)=-(u(m,nq,o,:,N+1,:,:)+u(m,l,o,:,1,:,:))*0.5_RP+duy(m,l,o,:,N+1,:,:)
                elseif(l==nq) then
                duy(m,l,o,:,1,:,:)=(u(m,l,o,:,N+1,:,:)+u(m,1,o,:,1,:,:))*0.5_RP+duy(m,l,o,:,1,:,:)
                duy(m,l,o,:,N+1,:,:)=-(u(m,l-1,o,:,N+1,:,:)+u(m,l,o,:,1,:,:))*0.5_RP+duy(m,l,o,:,N+1,:,:)
                else
                duy(m,l,o,:,1,:,:)=(u(m,l,o,:,N+1,:,:)+u(m,l+1,o,:,1,:,:))*0.5_RP+duy(m,l,o,:,1,:,:)
                duy(m,l,o,:,N+1,:,:)=-(u(m,l-1,o,:,N+1,:,:)+u(m,l,o,:,1,:,:))*0.5_RP+duy(m,l,o,:,N+1,:,:)
                end if
                !z-direction
                if(o==1) then
                duz(m,l,o,:,:,1,:)=(u(m,l,o,:,:,N+1,:)+u(m,l,o+1,:,:,1,:))*0.5_RP+duz(m,l,o,:,:,1,:)
                duz(m,l,o,:,:,N+1,:)=-(u(m,l,nq,:,:,N+1,:)+u(m,l,o,:,:,1,:))*0.5_RP+duz(m,l,o,:,:,N+1,:)
                elseif(o==nq) then
                duz(m,l,o,:,:,1,:)=(u(m,l,o,:,:,N+1,:)+u(m,l,1,:,:,1,:))*0.5_RP+duz(m,l,o,:,:,1,:)
                duz(m,l,o,:,:,N+1,:)=-(u(m,l,o-1,:,:,N+1,:)+u(m,l,o,:,:,1,:))*0.5_RP+duz(m,l,o,:,:,N+1,:)
                else
                duz(m,l,o,:,:,1,:)=(u(m,l,o,:,:,N+1,:)+u(m,l,o+1,:,:,1,:))*0.5_RP+duz(m,l,o,:,:,1,:)
                duz(m,l,o,:,:,N+1,:)=-(u(m,l,o-1,:,:,N+1,:)+u(m,l,o,:,:,1,:))*0.5_RP+duz(m,l,o,:,:,N+1,:)
                end if
            end do !o
        enddo!l
    end do!m
  END SUBROUTINE
  SUBROUTINE computeLviscous(u,dux,duy,duz,D,dir,N,NQ,result)
    !puts all of the viscous components together
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                                                  :: N,NQ
    INTEGER      ,INTENT(IN)                                                  :: dir
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: u,dux,duy,duz
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:N+1,1:N+1)                          :: D
    !u DIMENSIONs:(nummer zelle x,zelle y,zelle z,x,y,z,variable)
    REAL(KIND=RP),INTENT(OUT),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: result
    !local Variables
    INTEGER                                  :: var,k,j,i,o,l,m
    REAL(KIND=RP),DIMENSION(1:N+1,1:5)       :: Fviscous
    REAL(KIND=RP),DIMENSION(1:N+1,1:N+1,1:5) :: uR,uL,fvisl1,fvisl2,fvisr1,fvisr2
    REAL(KIND=RP),DIMENSION(1:n+1,1:5,1:3)   :: du
    REAL(KIND=RP),DIMENSION(1:n+1,1:n+1,1:5,1:3)   :: duRandl,duRandr

    SELECT CASE(dir)
    CASE(1)
      DO o=1,NQ
        DO l=1,NQ
          DO m=1,NQ
            DO k=1,N+1
              DO j=1,N+1
                DO i=1,N+1
                  du(:,:,1)=dux(m,l,o,:,j,k,:)
                  du(:,:,2)=duy(m,l,o,:,j,k,:)
                  du(:,:,3)=duz(m,l,o,:,j,k,:)
                  CALL computeviscousFlux(u(m,l,o,:,j,k,:),du,dir,n,Fviscous)
                  DO var=1,5 !! besser
                    result(m,l,o,i,j,k,var)=dot_product(D(i,:),Fviscous(:,var)*w)*1.0_RP/w(i)
                  ENDDO ! var
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
            !Randbedingungen
            IF (m==1) THEN
              uL=u(nq,l,o,N+1,:,:,:)
              duRandl(:,:,:,1)=dux(nq,l,o,N+1,:,:,:)
              duRandl(:,:,:,2)=duy(nq,l,o,N+1,:,:,:)
              duRandl(:,:,:,3)=duz(nq,l,o,N+1,:,:,:)
            ELSE
              uL=u(m-1,l,o,N+1,:,:,:)
              duRandl(:,:,:,1)=dux(m-1,l,o,N+1,:,:,:)
              duRandl(:,:,:,2)=duy(m-1,l,o,N+1,:,:,:)
              duRandl(:,:,:,3)=duz(m-1,l,o,N+1,:,:,:)
            ENDIF
            IF (m==NQ) THEN
              uR=u(1,l,o,1,:,:,:)
              duRandr(:,:,:,1)=dux(1,l,o,1,:,:,:)
              duRandr(:,:,:,2)=duy(1,l,o,1,:,:,:)
              duRandr(:,:,:,3)=duz(1,l,o,1,:,:,:)
            ELSE
              uR=u(m+1,l,o,1,:,:,:)
              duRandr(:,:,:,1)=dux(m+1,l,o,1,:,:,:)
              duRandr(:,:,:,2)=duy(m+1,l,o,1,:,:,:)
              duRandr(:,:,:,3)=duz(m+1,l,o,1,:,:,:)
            ENDIF
            CALL computeviscousFluxRand(uL,duRandl,dir,n,Fvisl1)
            CALL computeviscousFluxRand(uR,duRandr,dir,n,Fvisr2)
            duRandl(:,:,:,1)=dux(m,l,o,N+1,:,:,:)
            duRandl(:,:,:,2)=duy(m,l,o,N+1,:,:,:)
            duRandl(:,:,:,3)=duz(m,l,o,N+1,:,:,:)
            duRandr(:,:,:,1)=dux(m,l,o,1,:,:,:)
            duRandr(:,:,:,2)=duy(m,l,o,1,:,:,:)
            duRandr(:,:,:,3)=duz(m,l,o,1,:,:,:)
            CALL computeviscousFluxRand(u(m,l,o,N+1,:,:,:),duRandl,dir,n,Fvisl2)
            CALL computeviscousFluxRand(u(m,l,o,1,:,:,:),duRandr,dir,n,Fvisr1)
            result(m,l,o,1,:,:,:)=result(m,l,o,1,:,:,:)+(Fvisl1+Fvisr1)*0.5_RP
            result(m,l,o,N+1,:,:,:)=result(m,l,o,N+1,:,:,:)-(Fvisl2+Fvisr2)*0.5_RP
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
                  du(:,:,1)=dux(m,l,o,i,:,k,:)
                  du(:,:,2)=duy(m,l,o,i,:,k,:)
                  du(:,:,3)=duz(m,l,o,i,:,k,:)
                  CALL computeviscousFlux(u(m,l,o,i,:,k,:),du,dir,n,Fviscous)
                  DO var=1,5 !! besser
                    result(m,l,o,i,j,k,var)=dot_product(D(j,:),Fviscous(:,var)*w)*1.0_RP/w(j)
                  ENDDO ! var
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
            !Randbedingungen
            IF (l==1) THEN
              uL=u(m,nq,o,:,N+1,:,:)
              duRandl(:,:,:,1)=dux(m,nq,o,:,N+1,:,:)
              duRandl(:,:,:,2)=duy(m,nq,o,:,N+1,:,:)
              duRandl(:,:,:,3)=duz(m,nq,o,:,N+1,:,:)
            ELSE
              uL=u(m,l-1,o,:,N+1,:,:)
              duRandl(:,:,:,1)=dux(m,l-1,o,:,N+1,:,:)
              duRandl(:,:,:,2)=duy(m,l-1,o,:,N+1,:,:)
              duRandl(:,:,:,3)=duz(m,l-1,o,:,N+1,:,:)
            ENDIF
            IF (l==NQ) THEN
              uR=u(m,1,o,:,1,:,:)
              duRandr(:,:,:,1)=dux(m,1,o,:,1,:,:)
              duRandr(:,:,:,2)=duy(m,1,o,:,1,:,:)
              duRandr(:,:,:,3)=duz(m,1,o,:,1,:,:)
            ELSE
              uR=u(m,l+1,o,:,1,:,:)
              duRandr(:,:,:,1)=dux(m,l+1,o,:,1,:,:)
              duRandr(:,:,:,2)=duy(m,l+1,o,:,1,:,:)
              duRandr(:,:,:,3)=duz(m,l+1,o,:,1,:,:)
            ENDIF
            CALL computeviscousFluxRand(uL,duRandl,dir,n,Fvisl1)
            CALL computeviscousFluxRand(uR,duRandr,dir,n,Fvisr2)
            duRandl(:,:,:,1)=dux(m,l,o,:,N+1,:,:)
            duRandl(:,:,:,2)=duy(m,l,o,:,N+1,:,:)
            duRandl(:,:,:,3)=duz(m,l,o,:,N+1,:,:)
            duRandr(:,:,:,1)=dux(m,l,o,:,1,:,:)
            duRandr(:,:,:,2)=duy(m,l,o,:,1,:,:)
            duRandr(:,:,:,3)=duz(m,l,o,:,1,:,:)
            CALL computeviscousFluxRand(u(m,l,o,:,N+1,:,:),duRandl,dir,n,Fvisl2)
            CALL computeviscousFluxRand(u(m,l,o,:,1,:,:),duRandr,dir,n,Fvisr1)
            result(m,l,o,:,1,:,:)=result(m,l,o,:,1,:,:)+(Fvisl1+Fvisr1)*0.5_RP
            result(m,l,o,:,N+1,:,:)=result(m,l,o,:,N+1,:,:)-(Fvisl2+Fvisr2)*0.5_RP
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
                  du(:,:,1)=dux(m,l,o,i,j,:,:)
                  du(:,:,2)=duy(m,l,o,i,j,:,:)
                  du(:,:,3)=duz(m,l,o,i,j,:,:)
                  CALL computeviscousFlux(u(m,l,o,i,j,:,:),du,dir,n,Fviscous)
                  DO var=1,5 !! besser
                    result(m,l,o,i,j,k,var)=dot_product(D(k,:),Fviscous(:,var)*w)*1.0_RP/w(k)
                  ENDDO ! var
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
            !Randbedingungen
            IF (o==1) THEN
              uL=u(m,l,nq,:,:,N+1,:)
              duRandl(:,:,:,1)=dux(m,l,nq,:,:,N+1,:)
              duRandl(:,:,:,2)=duy(m,l,nq,:,:,N+1,:)
              duRandl(:,:,:,3)=duz(m,l,nq,:,:,N+1,:)
            ELSE
              uL=u(m,l,o-1,:,:,N+1,:)
              duRandl(:,:,:,1)=dux(m,l,o-1,:,:,N+1,:)
              duRandl(:,:,:,2)=duy(m,l,o-1,:,:,N+1,:)
              duRandl(:,:,:,3)=duz(m,l,o-1,:,:,N+1,:)
            ENDIF
            IF (o==NQ) THEN
              uR=u(m,l,1,:,:,1,:)
              duRandr(:,:,:,1)=dux(m,l,1,:,:,1,:)
              duRandr(:,:,:,2)=duy(m,l,1,:,:,1,:)
              duRandr(:,:,:,3)=duz(m,l,1,:,:,1,:)
            ELSE
              uR=u(m,l,o+1,:,:,1,:)
              duRandr(:,:,:,1)=dux(m,l,o+1,:,:,1,:)
              duRandr(:,:,:,2)=duy(m,l,o+1,:,:,1,:)
              duRandr(:,:,:,3)=duz(m,l,o+1,:,:,1,:)
            ENDIF
            CALL computeviscousFluxRand(uL,duRandl,dir,n,Fvisl1)
            CALL computeviscousFluxRand(uR,duRandr,dir,n,Fvisr2)
            duRandl(:,:,:,1)=dux(m,l,o,:,:,N+1,:)
            duRandl(:,:,:,2)=duy(m,l,o,:,:,N+1,:)
            duRandl(:,:,:,3)=duz(m,l,o,:,:,N+1,:)
            duRandr(:,:,:,1)=dux(m,l,o,:,:,1,:)
            duRandr(:,:,:,2)=duy(m,l,o,:,:,1,:)
            duRandr(:,:,:,3)=duz(m,l,o,:,:,1,:)
            CALL computeviscousFluxRand(u(m,l,o,:,:,N+1,:),duRandl,dir,n,Fvisl2)
            CALL computeviscousFluxRand(u(m,l,o,:,:,1,:),duRandr,dir,n,Fvisr1)
            result(m,l,o,:,:,1,:)=result(m,l,o,:,:,1,:)+(Fvisl1+Fvisr1)*0.5_RP
            result(m,l,o,:,:,N+1,:)=result(m,l,o,:,:,N+1,:)-(Fvisl2+Fvisr2)*0.5_RP
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
    END SELECT
  END SUBROUTINE computeLviscous
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeviscousFluxRand(u,du,dir,n,result)
    !computes the viscous part of the flux analytically
    IMPLICIT NONE
    INTEGER       ,INTENT(IN)                               :: N,dir
    REAL(KIND=RP) ,INTENT(IN) ,DIMENSION(1:n+1,1:n+1,1:5)   :: u
    REAL(KIND=RP) ,INTENT(IN) ,DIMENSION(1:n+1,1:n+1,1:5,1:3)   :: du !(Punkte x Punkte x Variable x Ableitungsrichtung)
    REAL(KIND=RP) ,INTENT(OUT),DIMENSION(1:n+1,1:n+1,1:5)   :: result
    REAL(KIND=RP)             ,DIMENSION(1:N+1,1:N+1)       ::P
    REAL(KIND=RP)             ,DIMENSION(1:n+1,1:n+1,1:3)   ::dTemp,dv1,dv2,dv3,dP
    result(:,:,1)=0.0_RP
    p=(gamma-1.0_RP)*(u(:,:,5)-0.5_RP*(u(:,:,2)*u(:,:,2)+u(:,:,3)*u(:,:,3)+u(:,:,4)*u(:,:,4))/u(:,:,1))
    dP(:,:,1)=(gamma-1.0_RP)*du(:,:,1,1)*(2.0_RP*u(:,:,1)-3.0_RP/2.0_RP)
    dP(:,:,2)=(gamma-1.0_RP)*du(:,:,1,2)*(2.0_RP*u(:,:,1)-3.0_RP/2.0_RP)
    dP(:,:,3)=(gamma-1.0_RP)*du(:,:,1,3)*(2.0_RP*u(:,:,1)-3.0_RP/2.0_RP)
    dTemp(:,:,1)=(dp(:,:,1)*u(:,:,1)-du(:,:,1,1)*p)/u(:,:,1)**2
    dTemp(:,:,2)=(dp(:,:,2)*u(:,:,1)-du(:,:,1,2)*p)/u(:,:,1)**2
    dTemp(:,:,3)=(dp(:,:,3)*u(:,:,1)-du(:,:,1,3)*p)/u(:,:,1)**2
    dv1(:,:,1)=(du(:,:,1,1)*u(:,:,2)/u(:,:,1))/u(:,:,1)
    dv1(:,:,2)=(du(:,:,1,2)*u(:,:,2)/u(:,:,1))/u(:,:,1)
    dv1(:,:,3)=(du(:,:,1,3)*u(:,:,2)/u(:,:,1))/u(:,:,1)
    dv2(:,:,1)=(du(:,:,1,1)*u(:,:,3)/u(:,:,1))/u(:,:,1)
    dv2(:,:,2)=(du(:,:,1,2)*u(:,:,3)/u(:,:,1))/u(:,:,1)
    dv2(:,:,3)=(du(:,:,1,3)*u(:,:,3)/u(:,:,1))/u(:,:,1)
    dv3(:,:,1)=(du(:,:,1,1)*u(:,:,4)/u(:,:,1))/u(:,:,1)
    dv3(:,:,2)=(du(:,:,1,2)*u(:,:,4)/u(:,:,1))/u(:,:,1)
    dv3(:,:,3)=(du(:,:,1,3)*u(:,:,4)/u(:,:,1))/u(:,:,1)
    SELECT CASE(dir)
    CASE(1)
      result(:,:,2)=mu*(2*du(:,:,1,1)-2.0_RP/3.0_RP*(dv1(:,:,1)+dv2(:,:,2)+dv3(:,:,3)))
      result(:,:,3)=mu*(dv1(:,:,2)+dv2(:,:,1))
      result(:,:,4)=mu*(dv1(:,:,3)+dv3(:,:,1))
      result(:,:,5)=mu*(u(:,:,2)/u(:,:,1)*(2.0_rp*dv2(:,:,2)-2.0_RP/3.0_RP*(dv1(:,:,1)+dv2(:,:,2)+dv3(:,:,3)))+&
                    u(:,:,3)/u(:,:,1)*(dv2(:,:,1)+dv1(:,:,2))+u(:,:,4)/u(:,:,1)*(dv3(:,:,1)+dv1(:,:,3)))+mu/(Pr*Rkonst)*dTemp(:,:,1)
    CASE(2)
      result(:,:,2)=mu*(dv1(:,:,2)+dv2(:,:,1))
      result(:,:,3)=mu*(2.0_RP*dv2(:,:,2)-2.0_RP/3.0_RP*(dv1(:,:,1)+dv2(:,:,2)+dv3(:,:,3)))
      result(:,:,4)=mu*(dv2(:,:,3)+dv3(:,:,2))
      result(:,:,5)=mu*(u(:,:,3)/u(:,:,1)*(2.0_RP*dv2(:,:,2)-2.0_RP/3.0_RP*(dv1(:,:,1)+dv2(:,:,2)+dv3(:,:,3)))+&
                    u(:,:,2)/u(:,:,1)*(dv2(:,:,1)+dv1(:,:,2))+u(:,:,4)/u(:,:,1)*(dv3(:,:,2)+dv2(:,:,3)))+mu/(Pr*Rkonst)*dTemp(:,:,2)
    CASE(3)
      result(:,:,2)=mu*(dv1(:,:,3)+dv3(:,:,1))
      result(:,:,3)=mu*(dv2(:,:,3)+dv3(:,:,2))
      result(:,:,4)=mu*(2.0_RP*dv3(:,:,3)-2.0_RP/3.0_RP*(dv1(:,:,1)+dv2(:,:,2)+dv3(:,:,3)))
      result(:,:,5)=mu*(u(:,:,4)/u(:,:,1)*(2.0_RP*dv3(:,:,3)-2.0_RP/3.0_RP*(dv1(:,:,1)+dv2(:,:,2)+dv3(:,:,3)))+&
                    u(:,:,2)/u(:,:,1)*(dv1(:,:,3)+dv3(:,:,1))+u(:,:,3)/u(:,:,1)*(dv2(:,:,3)+dv3(:,:,2)))+mu/(Pr*Rkonst)*dTemp(:,:,3)
    END SELECT
  END SUBROUTINE computeviscousFluxRand
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeviscousFlux(u,du,dir,n,result)
    !computes the viscous part of the flux analytically
    IMPLICIT NONE
    INTEGER       ,INTENT(IN)                               :: N,dir
    REAL(KIND=RP) ,INTENT(IN) ,DIMENSION(1:n+1,1:5)   :: u
    REAL(KIND=RP) ,INTENT(IN) ,DIMENSION(1:n+1,1:5,1:3)   :: du !(Punkte x Variable x Ableitungsrichtung)
    REAL(KIND=RP) ,INTENT(OUT),DIMENSION(1:n+1,1:5)   :: result
    REAL(KIND=RP)             ,DIMENSION(1:N+1)       ::P
    REAL(KIND=RP)             ,DIMENSION(1:n+1,1:3)   ::dTemp,dv1,dv2,dv3,dp
    result(:,1)=0.0_RP
    p=(gamma-1.0_RP)*(u(:,5)-0.5_RP*(u(:,2)*u(:,2)+u(:,3)*u(:,3)+u(:,4)*u(:,4))/u(:,1))
    dP(:,1)=(gamma-1.0_RP)*du(:,1,1)*(2.0_RP*u(:,1)-3.0_RP/2.0_RP)
    dP(:,2)=(gamma-1.0_RP)*du(:,1,2)*(2.0_RP*u(:,1)-3.0_RP/2.0_RP)
    dP(:,3)=(gamma-1.0_RP)*du(:,1,3)*(2.0_RP*u(:,1)-3.0_RP/2.0_RP)
    dTemp(:,1)=(dp(:,1)*u(:,1)-du(:,1,1)*p)/u(:,1)**2
    dTemp(:,2)=(dp(:,2)*u(:,1)-du(:,1,2)*p)/u(:,1)**2
    dTemp(:,3)=(dp(:,3)*u(:,1)-du(:,1,3)*p)/u(:,1)**2
    dv1(:,1)=(du(:,1,1)*u(:,2)/u(:,1))/u(:,1)
    dv1(:,2)=(du(:,1,2)*u(:,2)/u(:,1))/u(:,1)
    dv1(:,3)=(du(:,1,3)*u(:,2)/u(:,1))/u(:,1)
    dv2(:,1)=(du(:,1,1)*u(:,3)/u(:,1))/u(:,1)
    dv2(:,2)=(du(:,1,2)*u(:,3)/u(:,1))/u(:,1)
    dv2(:,3)=(du(:,1,3)*u(:,3)/u(:,1))/u(:,1)
    dv3(:,1)=(du(:,1,1)*u(:,4)/u(:,1))/u(:,1)
    dv3(:,2)=(du(:,1,2)*u(:,4)/u(:,1))/u(:,1)
    dv3(:,3)=(du(:,1,3)*u(:,4)/u(:,1))/u(:,1)
    SELECT CASE(dir)
    CASE(1)
      result(:,2)=mu*(2.0_RP*dv1(:,1)-2.0_RP/3.0_RP*(dv1(:,1)+dv2(:,2)+dv3(:,3)))
      result(:,3)=mu*(dv1(:,2)+dv2(:,1))
      result(:,4)=mu*(dv1(:,3)+dv3(:,1))
      result(:,5)=mu*(u(:,2)/u(:,1)*(2.0_rp*dv2(:,2)-2.0_RP/3.0_RP*(dv1(:,1)+dv2(:,2)+dv3(:,3)))+&
                    u(:,3)/u(:,1)*(dv2(:,1)+dv1(:,2))+u(:,4)/u(:,1)*(dv3(:,1)+dv1(:,3)))+mu/(Pr*Rkonst)*dTemp(:,1)
    CASE(2)
      result(:,2)=mu*(dv1(:,2)+dv2(:,1))
      result(:,3)=mu*(2.0_RP*dv2(:,2)-2.0_RP/3.0_RP*(dv1(:,1)+dv2(:,2)+dv3(:,3)))
      result(:,4)=mu*(dv2(:,3)+dv3(:,2))
      result(:,5)=mu*(u(:,3)/u(:,1)*(2.0_RP*dv2(:,2)-2.0_RP/3.0_RP*(dv1(:,1)+dv2(:,2)+dv3(:,3)))+&
                    u(:,2)/u(:,1)*(dv2(:,1)+dv1(:,2))+u(:,4)/u(:,1)*(dv3(:,2)+dv2(:,3)))+mu/(Pr*Rkonst)*dTemp(:,2)
    CASE(3)
      result(:,2)=mu*(dv1(:,3)+dv3(:,1))
      result(:,3)=mu*(dv2(:,3)+dv3(:,2))
      result(:,4)=mu*(2.0_RP*dv3(:,3)-2.0_RP/3.0_RP*(dv1(:,1)+dv2(:,2)+dv3(:,3)))
      result(:,5)=mu*(u(:,4)/u(:,1)*(2.0_RP*dv3(:,3)-2.0_RP/3.0_RP*(dv1(:,1)+dv2(:,2)+dv3(:,3)))+&
                    u(:,2)/u(:,1)*(dv1(:,3)+dv3(:,1))+u(:,3)/u(:,1)*(dv2(:,3)+dv3(:,2)))+mu/(Pr*Rkonst)*dTemp(:,3)
    END SELECT
  END SUBROUTINE computeviscousFlux
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        result(:,5)=(u1(4)/u1(1)*(u1(5)+p1)+u2(:,4)/u2(:,1)*(u2(:,5)+p2))*0.5_RP!fehler
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
  END SUBROUTINE computeFsharp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Initialcondition(u,NQ,N)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                                                    :: NQ,N
    REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: u
    u(:,:,:,:,:,:,1)=2.0_RP+SIN(pi*(xyz(:,:,:,:,:,:,1)+xyz(:,:,:,:,:,:,2)+xyz(:,:,:,:,:,:,3)))/10.0_RP
    u(:,:,:,:,:,:,2)=u(:,:,:,:,:,:,1)
    u(:,:,:,:,:,:,3)=u(:,:,:,:,:,:,1)
    u(:,:,:,:,:,:,4)=u(:,:,:,:,:,:,1)
    u(:,:,:,:,:,:,5)=u(:,:,:,:,:,:,1)*u(:,:,:,:,:,:,1)
  END SUBROUTINE Initialcondition
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeSolution(u,NQ,N,t)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                                                    :: NQ,N
    REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: u
    REAL(KIND=RP),INTENT(IN)                                                    :: t
    u(:,:,:,:,:,:,1)=2.0_rp+sin(pi*(xyz(:,:,:,:,:,:,1)+xyz(:,:,:,:,:,:,2)+xyz(:,:,:,:,:,:,3)-2*t))/10.0_rp
    u(:,:,:,:,:,:,2)=1.0_RP*u(:,:,:,:,:,:,1)
    u(:,:,:,:,:,:,3)=1.0_RP*u(:,:,:,:,:,:,1)
    u(:,:,:,:,:,:,4)=1.0_RP*u(:,:,:,:,:,:,1)
    u(:,:,:,:,:,:,5)=u(:,:,:,:,:,:,1)*u(:,:,:,:,:,:,1)
  END SUBROUTINE computeSolution
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
        result(:,:,1)=FPI(:,:,1)-lamMax*(uR(:,:,1)-uL(:,:,1))/2.0_RP-FR(:,:,1)
        result(:,:,2)=FPI(:,:,2)-lamMax*(uR(:,:,2)-uL(:,:,2))/2.0_RP-FR(:,:,2)
        result(:,:,3)=FPI(:,:,3)-lamMax*(uR(:,:,3)-uL(:,:,3))/2.0_RP-FR(:,:,3)
        result(:,:,4)=FPI(:,:,4)-lamMax*(uR(:,:,4)-uL(:,:,4))/2.0_RP-FR(:,:,4)
        result(:,:,5)=FPI(:,:,5)-lamMax*(uR(:,:,5)-uL(:,:,5))/2.0_RP-FR(:,:,5)
      CASE(1)
        result(:,:,1)=FPI(:,:,1)-lamMax*(uR(:,:,1)-uL(:,:,1))/2.0_RP-FL(:,:,1)
        result(:,:,2)=FPI(:,:,2)-lamMax*(uR(:,:,2)-uL(:,:,2))/2.0_RP-FL(:,:,2)
        result(:,:,3)=FPI(:,:,3)-lamMax*(uR(:,:,3)-uL(:,:,3))/2.0_RP-FL(:,:,3)
        result(:,:,4)=FPI(:,:,4)-lamMax*(uR(:,:,4)-uL(:,:,4))/2.0_RP-FL(:,:,4)
        result(:,:,5)=FPI(:,:,5)-lamMax*(uR(:,:,5)-uL(:,:,5))/2.0_RP-FL(:,:,5)
      END SELECT
    END SELECT
  END SUBROUTINE computeLocalLaxFriedrich
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
    h1=(u1(:,:,5)  +p1)   /u1(:,:,1)
      h2=(u2(:,:,5)+p2)/u2(:,:,1)
    SELECT CASE(dir)
      CASE(1)
        result(:,:,1)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,2)/u1(:,:,1)+u2(:,:,2)/u2(:,:,1))*0.25_RP
        result(:,:,2)=result(:,:,1)    *(u1(:,:,2)/u1(:,:,1)+u2(:,:,2)/u2(:,:,1))*0.5_RP  +(p1+p2)*0.5_RP
        result(:,:,3)=result(:,:,1)    *(u1(:,:,3)/u1(:,:,1)+u2(:,:,3)/u2(:,:,1))*0.5_RP
        result(:,:,4)=result(:,:,1)    *(u1(:,:,4)/u1(:,:,1)+u2(:,:,4)/u2(:,:,1))*0.5_RP
        result(:,:,5)=result(:,:,1)    *(h1+h2)*0.5_RP
      CASE(2)
        result(:,:,1)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,3)/u1(:,:,1)+u2(:,:,3)/u2(:,:,1))*0.25_RP
        result(:,:,2)=result(:,:,1)    *(u1(:,:,2)/u1(:,:,1)+u2(:,:,2)/u2(:,:,1))*0.5_RP
        result(:,:,3)=result(:,:,1)    *(u1(:,:,3)/u1(:,:,1)+u2(:,:,3)/u2(:,:,1))*0.5_RP  +(p1+p2)*0.5_RP
        result(:,:,4)=result(:,:,1)    *(u1(:,:,4)/u1(:,:,1)+u2(:,:,4)/u2(:,:,1))*0.5_RP
        result(:,:,5)=result(:,:,1)    *(h1+h2)*0.5_RP
      CASE(3)
        result(:,:,1)=(u1(:,:,1)+u2(:,:,1))*(u1(:,:,4)/u1(:,:,1)+u2(:,:,4)/u2(:,:,1))*0.25_RP
        result(:,:,2)=result(:,:,1)    *(u1(:,:,2)/u1(:,:,1)+u2(:,:,2)/u2(:,:,1))*0.5_RP
        result(:,:,3)=result(:,:,1)    *(u1(:,:,3)/u1(:,:,1)+u2(:,:,3)/u2(:,:,1))*0.5_RP
        result(:,:,4)=result(:,:,1)    *(u1(:,:,4)/u1(:,:,1)+u2(:,:,4)/u2(:,:,1))*0.5_RP  +(p1+p2)*0.5_RP
        result(:,:,5)=result(:,:,1)    *(h1+h2)*0.5_RP
      END SELECT
  END SUBROUTINE calculateEulerRandFluxPI
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
    !result(:,:)=maxval(lambda(:,:,:))
  END SUBROUTINE lambdaMax
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
    if (any(gamma*p/u(:,:,:,:,:,:,1)<0)) then
    print*,1
    endif
    c=sqrt(gamma*p/u(:,:,:,:,:,:,1))
    !! max from abs(eigenvalue)
    lambdamax=max(MAXVAL(ABS(u(:,:,:,:,:,:,2)/u(:,:,:,:,:,:,1)+c)),MAXVAL(ABS(u(:,:,:,:,:,:,3)/u(:,:,:,:,:,:,1)+c)),&
      MAXVAL(ABS(u(:,:,:,:,:,:,4)/u(:,:,:,:,:,:,1)+c)),MAXVAL(ABS(u(:,:,:,:,:,:,2)/u(:,:,:,:,:,:,1)-c)),&
      MAXVAL(abs(u(:,:,:,:,:,:,3)/u(:,:,:,:,:,:,1)-c)),MAXVAL(ABS(u(:,:,:,:,:,:,4)/u(:,:,:,:,:,:,1)-c)),&
      MAXVAL(ABS(u(:,:,:,:,:,:,2)/u(:,:,:,:,:,:,1))),MAXVAL(ABS(u(:,:,:,:,:,:,3)/u(:,:,:,:,:,:,1))),&
      MAXVAL(ABS(u(:,:,:,:,:,:,4)/u(:,:,:,:,:,:,1))))
  END SUBROUTINE lambdaMaxGlobal
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
    p=(gamma-1.0_RP)*(u(:,:,5)  -0.5_RP*(u(:,:,2)*u(:,:,2)    +u(:,:,3)*u(:,:,3)    +u(:,:,4)*u(:,:,4))    /u(:,:,1)  )
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
  END SUBROUTINE calculateEulerRandFlux
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE RungeKutta5explizit(ustar,nq,n,numvar,dt,Dval,t,whichflux,vis)
    IMPLICIT NONE
    !ustar=bekannte Werte
    INTEGER      ,INTENT(IN)                                                         :: numvar,n,nq
    CHARACTER(len=2),INTENT(IN)                                                      :: whichflux,vis
    REAL(KIND=RP),INTENT(IN)   ,DIMENSION(1:N+1,1:N+1)                               :: Dval
    REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:nq,1:nq,1:nq,1:N+1,1:N+1,1:N+1,1:numvar) :: ustar
    !local
    REAL(KIND=RP),DIMENSION(1:nq,1:nq,1:nq,1:(N+1),1:N+1,1:N+1,1:numvar)    :: g
    INTEGER                                                                 :: step
    REAL(KIND=RP),DIMENSION(5)                                              :: a,b,c
    REAL(KIND=RP),INTENT(IN)                                                :: dt,t
    g=Rmanu(ustar,n,nq,Dval,t,whichflux,vis,dt)
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
      g=a(step)*g+Rmanu(ustar,n,nq,Dval,t+b(step)*dt,whichflux,vis,b(step)*dt)
      ustar=ustar+c(step)*dt*g
    ENDDO ! step
  END SUBROUTINE RungeKutta5explizit
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeError (u,usolution,NQ,N,result)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                                                      :: NQ,N
    REAL(KIND=RP),INTENT(IN),DIMENSION(1:nq,1:nq,1:nq,1:N+1,1:N+1,1:N+1,1:5) :: u,usolution
    REAL(KIND=RP),INTENT(OUT),DIMENSION(1:5)                                 :: result
    result(1)=maxval(abs(u(:,:,:,:,:,:,1)-usolution(:,:,:,:,:,:,1)))
    result(2)=maxval(abs(u(:,:,:,:,:,:,2)-usolution(:,:,:,:,:,:,2)))
    result(3)=maxval(abs(u(:,:,:,:,:,:,3)-usolution(:,:,:,:,:,:,3)))
    result(4)=maxval(abs(u(:,:,:,:,:,:,4)-usolution(:,:,:,:,:,:,4)))
    result(5)=maxval(abs(u(:,:,:,:,:,:,5)-usolution(:,:,:,:,:,:,5)))
    print*, result
    DEALLOCATE(xyz,x,w,xmit,xges)
  END SUBROUTINE computeError
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeEOC (errors,n,nq,anz,result)
    IMPLICIT NONE
    INTEGER,INTENT(IN)                             :: anz,n
    REAL(KIND=RP),INTENT(IN),DIMENSION(1:5,1:anz)  :: errors
    INTEGER,INTENT(IN),DIMENSION(1:anz)            :: nq
    REAL(KIND=RP),INTENT(OUT),DIMENSION(1:5,1:anz) :: result
    INTEGER                                        :: k
    result=0.0_RP
    DO k=1,anz-1
      result(:,k+1)=log(errors(:,k+1)/errors(:,k))/log(real(nq(k),rp)/(nq(k+1)))
    END DO

  END SUBROUTINE computeEOC
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Zeitintegration
