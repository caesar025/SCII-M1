MODULE Zeitintegration
  USE DGtoolbox
  USE mpi
  REAL(KIND=RP),DIMENSION(:),ALLOCATABLE             :: x,w,xmit,xges
  REAL(KIND=RP)                                      :: gk=9.812_RP,dx,gamma=1.4_RP,mu=0.001_RP,Pr=0.72_RP,Rkonst=287.058_RP
  INTEGER                                            :: NQX,NQY,NQZ,id,num_procs,ierr
  INTEGER,DIMENSION(MPI_STATUS_SIZE)                 :: stat
  INTEGER,DIMENSION(:,:),ALLOCATABLE                 :: part_map
  INTEGER,DIMENSION(:,:,:),ALLOCATABLE              :: procs_map
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
  !FUNCTION R(u,N,NQ,D,whichflux) result(solution)
  !  IMPLICIT NONE
  !  INTEGER      ,INTENT(IN)                                                 :: N,NQ
  !  REAL(KIND=RP),INTENT(IN),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: u
  !  REAL(KIND=RP),INTENT(IN),DIMENSION(1:N+1,1:N+1)                          :: D
  !  REAL(KIND=RP)           ,DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: solution
  !  CHARACTER(len=2),INTENT(IN)                                              :: whichflux
  !  ! local
  !  REAL(KIND=RP)           ,DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) ::L1,L2,L3
  !  !
  !  CALL computeL(u,D,1,L1,dux,N,NQ,whichflux)
  !  CALL computeL(u,D,2,L2,dux,N,NQ,whichflux)
  !  CALL computeL(u,D,3,L3,dux,N,NQ,whichflux)
  !  solution=8.0_RP/(dx**3)*((-0.25_RP*dx**2)*L1-(0.25_RP*dx**2)*L2-(0.25_RP*dx**2)*L3)
  !END FUNCTION R
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION Rmanu(usub,N,NQ,D,t,whichflux,vis,dt,xyzsub) result(solution)
    ! Setze Rechteseite zusammen mit manufactuerd solution
    IMPLICIT NONE
    INTEGER, INTENT(IN)                                                           :: n,NQ
    REAL(KIND=RP), DIMENSION(1:NQX,1:NQY,1:NQZ,1:n+1,1:(N+1),1:n+1,1:5)           :: solution    ! Rechte Seite
    REAL(KIND=RP), DIMENSION(1:NQX,1:NQY,1:NQZ,1:n+1,1:(N+1),1:n+1,1:5)           :: res                       ! Quellterme
    REAL(KIND=RP), INTENT(IN), DIMENSION(1:NQX,1:NQY,1:NQZ,1:n+1,1:n+1,1:n+1,1:5) :: usub
    REAL(KIND=RP), INTENT(IN), DIMENSION(1:NQX,1:NQY,1:NQZ,1:n+1,1:n+1,1:n+1,1:3) :: xyzsub
    REAL(KIND=RP),INTENT(IN),DIMENSION(1:N+1,1:N+1)                               :: D                        ! Diff-Matrix
    CHARACTER(len=2),INTENT(IN)                                                   :: whichflux,vis
    REAL(KIND=RP),INTENT(IN)                                                      :: t,dt                    ! Startzeit,zeitschritt
    REAL(KIND=RP)           ,DIMENSION(1:NQX,1:NQY,1:NQZ,1:N+1,1:N+1,1:N+1,1:5)   :: L1,L2,L3,dux,duy,duz,L1vis,L2vis,L3vis
    INTEGER                                                                       :: i !num_procs=#Processoren ierr=Errorflag
    REAL(KIND=RP)                                                                 :: per_proc
    !we assume nq**3 is divisible by num_procs
    IF(mod(nq**3,num_procs)/=0) THEN
      call MPI_Finalize(ierr)
      print*,'Number of elements(',nq**3,') is not divided by ',num_procs
    END IF
    CALL computeL(usub,D,1,L1,dux,N,NQ,whichflux)
    CALL computeL(usub,D,2,L2,duy,N,NQ,whichflux)
    CALL computeL(usub,D,3,L3,duZ,N,NQ,whichflux)
    SELECT CASE(vis)
    CASE('AD')
      solution=8.0_RP/(dx**3)*(-0.25_RP*dx*dx*l1-0.25_RP*dx*dx*l2-0.25_RP*dx*dx*l3)
      call Residuum(NQ,N,t,vis,res,xyzsub) 
      solution=solution+res
    CASE('VI')
      solution=8.0_RP/(dx**3)*(-0.25_RP*dx*dx*l1-0.25_RP*dx*dx*l2-0.25_RP*dx*dx*l3)
      call Residuum(NQ,N,t,vis,res,xyzsub)
      call computeLviscous(usub,dux,duy,duz,D,1,N,NQ,L1vis)
      call computeLviscous(usub,dux,duy,duz,D,2,N,NQ,L2vis)
      call computeLviscous(usub,dux,duy,duz,D,3,N,NQ,L3vis)
      solution=solution+res
      solution=solution+2.0_RP/(dx)*(l1vis+l2vis+l3vis)
    END SELECT
  END FUNCTION Rmanu
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Residuum (NQ,N,t,vis,result,xyzsub)
    ! Berechne Quellterme
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                                        :: n,NQ
    REAL(KIND=RP),INTENT(IN)                                        :: t                    ! Startzeit
    REAL(KIND=RP),DIMENSION(1:NQx,1:nqy,1:nqZ,1:n+1,1:(N+1),1:n+1,1:5) :: result                       ! Quellterme
    REAL(KIND=RP),DIMENSION(1:NQx,1:nqy,1:nqZ,1:n+1,1:(N+1),1:n+1,1:3) :: xyzsub                      
    CHARACTER(len=2),INTENT(IN)                                     :: vis                !Schalter ob viskos oder advektiv
    REAL(KIND=RP)                                                   :: c1,c2,c3,c4,c5,ro,rox,Px,p               ! Hilfsvariablen
    INTEGER                                                         :: o,l,m,k,j,i

    c1=pi/10.0_RP
    c2=-1.0_RP/5.0_RP*pi+pi/20.0_rp*(1.0_rp+5.0_RP*gamma)
    c3=pi/100.0_RP*(gamma-1)
    c4=(-16.0_RP*pi+pi*(9.0_RP+15.0_RP*gamma))/20.0_RP
    c5=(3*pi*gamma-2.0*pi)/100.0_RP

    DO o=1,NQZ
      DO l=1,NQY
        DO m=1,NQX
          DO k=1,n+1
            DO j=1,n+1
              DO i=1,n+1
                ro=2.0_RP+1.0_RP/10.0_RP*sin(pi*(xyzsub(m,l,o,i,j,k,1)+xyzsub(m,l,o,i,j,k,2)+xyzsub(m,l,o,i,j,k,3)-2.0_RP*t))
                rox=cos(pi*(xyzsub(m,l,o,i,j,k,1)+xyzsub(m,l,o,i,j,k,2)+xyzsub(m,l,o,i,j,k,3)-2.0_RP*t))*pi/10.0_RP
                Px=(gamma-1.0_RP)*((2.0_RP*ro-3.0_RP/2.0_RP)*rox)
                P=(gamma-1.0_RP)*(ro**2-3.0_RP/2.0_RP*ro)
                result(m,l,o,i,j,k,1)=rox
                result(m,l,o,i,j,k,2)=Px+rox
                result(m,l,o,i,j,k,3)=result(m,l,o,i,j,k,2)
                result(m,l,o,i,j,k,4)=result(m,l,o,i,j,k,2)
                result(m,l,o,i,j,k,5)=2.0_RP*ro*rox+3.0_RP*Px
                IF (vis=='VI') THEN
                  result(m,l,o,i,j,k,5)=result(m,l,o,i,j,k,5)+3.0_RP*mu/(Pr*Rkonst)*(Px*ro-rox*p)/(ro**2)
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE Residuum
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeL(usub,D,dir,result,du,N,NQ,whichflux)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                                                     :: N,NQ
    INTEGER      ,INTENT(IN)                                                     :: dir
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:NQX,1:NQY,1:NQZ,1:N+1,1:N+1,1:N+1,1:5) :: usub
    REAL(KIND=RP),DIMENSION(:,:,:,:,:,:),allocatable                             :: urand
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:N+1,1:N+1)                             :: D
    CHARACTER(len=2),INTENT(IN)                                                  :: whichflux
    !u DIMENSIONs:(nummer zelle x,zelle y,zelle z,x,y,z,variable)
    REAL(KIND=RP),INTENT(OUT),DIMENSION(1:NQX,1:NQY,1:NQZ,1:N+1,1:N+1,1:N+1,1:5) :: result,du
    !local Variables
    INTEGER                                                                      :: var,k,j,i,o,l,m,left,right
    REAL(KIND=RP),DIMENSION(1:N+1,1:5)                                           :: Fsharp
    REAL(KIND=RP),DIMENSION(1:N+1,1:N+1,1:5)                                     :: FRand0,FRand1,uR,uL
    SELECT CASE(dir)
    CASE(1)
      allocate(urand(1:2,1:NQY,1:NQZ,1:N+1,1:N+1,1:5))
      DO o=1,NQZ
        DO l=1,NQY
          DO m=1,NQX
            DO k=1,N+1
              DO j=1,N+1
                DO i=1,N+1
                  CALL computeFsharp(usub(m,l,o,i,j,k,:),usub(m,l,o,:,j,k,:),dir,whichflux,Fsharp,N)
                  DO var=1,5 !! besser
                    result(m,l,o,i,j,k,var)=2.0_RP*dot_product(D(i,:),Fsharp(:,var))
                  ENDDO ! var
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
      
      if(nq/=NQX) then
        call comm_rand_bed(usub,urand,1,N) 
      endif
      DO o=1,NQZ
        DO l=1,NQY
          DO m=1,NQX
            IF (m==1) THEN
              if(nq==nqx) then
                uL=usub(NQX,l,o,N+1,:,:,:)
              else
                uL=urand(1,l,o,:,:,:)
              !  print*, UL
              endif
            ELSE
              uL=usub(m-1,l,o,N+1,:,:,:)
            ENDIF
            IF (m==NQX) THEN
              if(nq==nqx) then
                uR=usub(1,l,o,1,:,:,:)
              else
                uR=urand(2,l,o,:,:,:)
              end if
            ELSE
              uR=usub(m+1,l,o,1,:,:,:)
            ENDIF
            !Randbedingungen
            CALL computeLocalLaxFriedrich(uL,usub(m,l,o,1,:,:,:),dir,0,whichflux,FRand0,N)
            CALL computeLocalLaxFriedrich(usub(m,l,o,N+1,:,:,:),uR,dir,1,whichflux,FRand1,N)
            result(m,l,o,1,:,:,:)=result(m,l,o,1,:,:,:)-FRand0
            result(m,l,o,N+1,:,:,:)=result(m,l,o,N+1,:,:,:)+FRand1
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o 
  call computeGradientEinzeln(usub,urand,n,nq,D,du,dir)
    CASE(2)
      allocate(urand(1:2,1:NQX,1:NQZ,1:N+1,1:N+1,1:5))
      DO o=1,NQZ
        DO l=1,NQY
          DO m=1,NQX
            DO k=1,N+1
              DO j=1,N+1
                DO i=1,N+1
                  CALL computeFsharp(usub(m,l,o,i,j,k,:),usub(m,l,o,i,:,k,:),dir,whichflux,Fsharp,N)
                  DO var=1,5 !! besser
                    result(m,l,o,i,j,k,var)=2.0_RP*dot_product(D(j,:),Fsharp(:,var))
                  ENDDO ! var
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
      if(nq/=NQY) then
        call comm_rand_bed(usub,urand,2,N) 
      endif
      DO o=1,NQZ
        DO l=1,NQY
          DO m=1,NQX
            !Randbedingungen
            IF (l==1) THEN
            if(nq==nqy) then
              uL=usub(m,nqy,o,:,N+1,:,:)
            else
              uL=urand(1,m,o,:,:,:)
            endif
            ELSE
              uL=usub(m,l-1,o,:,N+1,:,:)
            ENDIF
            IF (l==nqy) THEN
            if(nq==nqy) then
              uR=usub(m,1,o,:,1,:,:)
            else
              uR=urand(2,m,o,:,:,:)
            end if
            ELSE
              uR=usub(m,l+1,o,:,1,:,:)
            ENDIF
            !Randbedingungen
            CALL computeLocalLaxFriedrich(uL,usub(m,l,o,:,1,:,:),dir,0,whichflux,FRand0,N)
            CALL computeLocalLaxFriedrich(usub(m,l,o,:,N+1,:,:),uR,dir,1,whichflux,FRand1,N)
            result(m,l,o,:,1,:,:)=result(m,l,o,:,1,:,:)-FRand0
            result(m,l,o,:,N+1,:,:)=result(m,l,o,:,N+1,:,:)+FRand1
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
  call computeGradientEinzeln(usub,urand,n,nq,D,du,dir)
    CASE(3)
      allocate(urand(1:2,1:NQX,1:NQY,1:N+1,1:N+1,1:5))
      DO o=1,NQZ
        DO l=1,NQY
          DO m=1,NQX
            DO k=1,N+1
              DO j=1,N+1
                DO i=1,N+1
                  CALL computeFsharp(usub(m,l,o,i,j,k,:),usub(m,l,o,i,j,:,:),dir,whichflux,Fsharp,N)
                  DO var=1,5 !! besser
                    result(m,l,o,i,j,k,var)=2.0_RP*dot_product(D(k,:),Fsharp(:,var))
                  ENDDO ! var
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
      if(nq/=NQZ) then
        call comm_rand_bed(usub,urand,3,N) 
      endif
      DO o=1,NQZ
        DO l=1,NQY
          DO m=1,NQX
            !Randbedingungen
            IF (o==1) THEN
              if(nq==nqz) then
                uL=usub(m,l,nqz,:,:,N+1,:)
              else
                uL=urand(1,m,l,:,:,:)
              endif
            ELSE
              uL=usub(m,l,o-1,:,:,N+1,:)
            ENDIF
            IF (o==nqz) THEN
              if(nq==nqZ) then
                uR=usub(m,l,1,:,:,1,:)
              else
                uR=urand(2,m,l,:,:,:)
              end if
            ELSE
              uR=usub(m,l,o+1,:,:,1,:)
            ENDIF
            !Randbedingungen
            CALL computeLocalLaxFriedrich(uL,usub(m,l,o,:,:,1,:),dir,0,whichflux,FRand0,N)
            CALL computeLocalLaxFriedrich(usub(m,l,o,:,:,N+1,:),uR,dir,1,whichflux,FRand1,N)
            result(m,l,o,:,:,1,:)=result(m,l,o,:,:,1,:)-FRand0
            result(m,l,o,:,:,N+1,:)=result(m,l,o,:,:,N+1,:)+FRand1
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
  call computeGradientEinzeln(usub,urand,n,nq,D,du,dir)
    END SELECT

    DEALLOCATE(urand)
  END SUBROUTINE computeL
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine comm_rand_bed (u,urand,richtung,N)
    implicit none
    INTEGER, INTENT(IN)                                :: richtung
    REAL(KIND=RP),DIMENSION(:,:,:,:,:,:),INTENT(INOUT) :: urand
    REAL(KIND=RP),DIMENSION(:,:,:,:,:,:,:),INTENT(IN)  :: u
    INTEGER                                            :: links,rechts,i,l
    INTEGER,INTENT(in)                                 :: N
    ! urand(links=1/rechts=2,zelle,zelle,matrix am Rand, matrix am Rand, variablen)
    SELECT CASE(richtung)
    case(1)
     do i=0,num_procs-1
       if(id==i) then 
         links=part_map(id,0)    
         rechts=part_map(id,1)    
         call MPI_SSEND(u(1,:,:,1,:,:,:),NQY*NQZ*(N+1)**2*5,MPI_DOUBLE_PRECISION,links,1,MPI_COMM_WORLD,ierr)
         call MPI_SSEND(u(NQX,:,:,n+1,:,:,:),NQY*NQZ*(N+1)**2*5,MPI_DOUBLE_PRECISION,rechts,2,MPI_COMM_WORLD,ierr)
       ELSEIF(id==part_map(i,0).AND.id==part_map(i,1)) then
         call MPI_RECV(urand(2,:,:,:,:,:),NQY*NQZ*(N+1)**2*5,MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,stat,ierr)
         call MPI_RECV(urand(1,:,:,:,:,:),NQY*NQZ*(N+1)**2*5,MPI_DOUBLE_PRECISION,i,2,MPI_COMM_WORLD,stat,ierr)
         
       elseif(id==part_map(i,0)) then
         call MPI_RECV(urand(2,:,:,:,:,:),NQY*NQZ*(N+1)**2*5,MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,stat,ierr)
       ELSEIF(id==part_map(i,1)) then 
         call MPI_RECV(urand(1,:,:,:,:,:),NQY*NQZ*(N+1)**2*5,MPI_DOUBLE_PRECISION,i,2,MPI_COMM_WORLD,stat,ierr)
       end if
     end do
    case(2)
      do i=0,num_procs-1
        if(id==i) then 
          links=part_map(id,2)    
          rechts=part_map(id,3)    
          call MPI_SEND(u(:,1,:,:,1,:,:),NQX*NQZ*(N+1)**2*5,MPI_DOUBLE_PRECISION,links,1,MPI_COMM_WORLD,ierr)
          call MPI_SEND(u(:,NQY,:,:,n+1,:,:),NQX*NQZ*(N+1)**2*5,MPI_DOUBLE_PRECISION,rechts,2,MPI_COMM_WORLD,ierr)
        ELSEIF(id==part_map(i,2).AND.id==part_map(i,3)) then
          call MPI_RECV(urand(2,:,:,:,:,:),NQX*NQZ*(N+1)**2*5,MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,stat,ierr)
          call MPI_RECV(urand(1,:,:,:,:,:),NQX*NQZ*(N+1)**2*5,MPI_DOUBLE_PRECISION,i,2,MPI_COMM_WORLD,stat,ierr)
        ELSEIF(id==part_map(i,2)) then
          call MPI_RECV(urand(2,:,:,:,:,:),NQX*NQZ*(N+1)**2*5,MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,stat,ierr)
        ELSEIF(id==part_map(i,3)) then 
          call MPI_RECV(urand(1,:,:,:,:,:),NQX*NQZ*(N+1)**2*5,MPI_DOUBLE_PRECISION,i,2,MPI_COMM_WORLD,stat,ierr)
        end if
    end do
    case(3)
    do i=0,num_procs-1
        if(id==i) then 
          links=part_map(id,4)    
          rechts=part_map(id,5)    
          call MPI_SEND(u(:,:,1,:,:,1,:),NQX*NQY*(N+1)**2*5,MPI_DOUBLE_PRECISION,links,1,MPI_COMM_WORLD,ierr)
          call MPI_SEND(u(:,:,NQZ,:,:,n+1,:),NQX*NQY*(N+1)**2*5,MPI_DOUBLE_PRECISION,rechts,2,MPI_COMM_WORLD,ierr)
        ELSEIF(id==part_map(i,4).AND.id==part_map(i,5)) then
          call MPI_RECV(urand(2,:,:,:,:,:),NQX*NQY*(N+1)**2*5,MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,stat,ierr)
          call MPI_RECV(urand(1,:,:,:,:,:),NQX*NQY*(N+1)**2*5,MPI_DOUBLE_PRECISION,i,2,MPI_COMM_WORLD,stat,ierr)
        ELSEIF(id==part_map(i,4)) then
          call MPI_RECV(urand(2,:,:,:,:,:),NQX*NQY*(N+1)**2*5,MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,stat,ierr)
      ELSEIF(id==part_map(i,5)) then 
          call MPI_RECV(urand(1,:,:,:,:,:),NQX*NQY*(N+1)**2*5,MPI_DOUBLE_PRECISION,i,2,MPI_COMM_WORLD,stat,ierr)
        end if
    end do
  END SELECT
  end subroutine
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
  SUBROUTINE computeGradientEinzeln(u,urand,n,nq,D,du,richtung)
    IMPLICIT NONE
    !Gleichung 57 im skript
    INTEGER      ,INTENT(IN)                                                     :: n,nq
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:nqx,1:nqy,1:nqz,1:n+1,1:n+1,1:n+1,1:5) :: u
    REAL(KIND=RP),DIMENSION(:,:,:,:,:,:),INTENT(IN)                              :: urand
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:n+1,1:n+1)                             :: D
    REAL(KIND=RP),INTENT(OUT),DIMENSION(1:nqx,1:nqy,1:nqz,1:n+1,1:n+1,1:n+1,1:5) :: du
    REAL(KIND=RP),DIMENSION(1:N+1,1:N+1,1:5)                                     :: uR,uL
    !local variables
    INTEGER                                                                      :: i,j,k,l,m,o,var,richtung

    IF (.not.allocated(w)) THEN
      print*,'w not allocated'
      stop
    ENDif

    DO m=1,nqx
      DO l=1,nqY
        DO o=1,nqZ
          DO i=1,n+1
            DO j=1,n+1
              DO k=1,n+1
                DO var=1,5
                  SELECT CASE (richtung)
                  Case(1)
                  du(m,l,o,i,j,k,var)=-dot_product(D(:,i),u(m,l,o,:,j,k,var)*w)!+surface term 
                  Case(2)
                  du(m,l,o,i,j,k,var)=-dot_product(D(:,j),u(m,l,o,i,:,k,var)*w)!+surface term
                  Case(3)
                  du(m,l,o,i,j,k,var)=-dot_product(D(:,k),u(m,l,o,i,j,:,var)*w)!+surface term
                end SELECT
                END DO !var
              END DO !k
            END DO !j
          END DO !i

          !surfaceterms and boundary conditions
          !x-direction
          SELECT CASE (richtung)
          CASE(1)
          IF(m==1) THEN
              if(nq==nqx) then
                uL=u(NQX,l,o,N+1,:,:,:)
              else
                uL=urand(1,l,o,:,:,:)
              endif
            du(m,l,o,1,:,:,:)=-(u(m,l,o,1,:,:,:)+uL)*1.0_RP/2.0_RP+du(m,l,o,1,:,:,:)
            du(m,l,o,N+1,:,:,:)=+(u(m,l,o,N+1,:,:,:)+u(m+1,l,o,1,:,:,:))*1.0_RP/2.0_RP+du(m,l,o,N+1,:,:,:)
          ELSEif(m==nqx) THEN
              if(nq==nqx) then
                uR=u(1,l,o,1,:,:,:)
              else
                uR=urand(2,l,o,:,:,:)
              end if
            du(m,l,o,1,:,:,:)=-(u(m-1,l,o,N+1,:,:,:)+u(m,l,o,1,:,:,:))*1.0_RP/2.0_RP+du(m,l,o,1,:,:,:)
            du(m,l,o,N+1,:,:,:)=+(u(m,l,o,N+1,:,:,:)+uR)*1.0_RP/2.0_RP+du(m,l,o,N+1,:,:,:)
          ELSE
            du(m,l,o,1,:,:,:)=-(u(m-1,l,o,N+1,:,:,:)+u(m,l,o,1,:,:,:))*1.0_RP/2.0_RP+du(m,l,o,1,:,:,:)
            du(m,l,o,N+1,:,:,:)=+(u(m,l,o,N+1,:,:,:)+u(m+1,l,o,1,:,:,:))*1.0_RP/2.0_RP+du(m,l,o,N+1,:,:,:)
          END IF
          CASE(2)
          !y-direction
          IF(l==1) THEN
            if(nq==nqy) then
              uL=u(m,nqy,o,:,N+1,:,:)
            else
              uL=urand(1,m,o,:,:,:)
            endif
            du(m,l,o,:,1,:,:)=-(u(m,l,o,:,1,:,:)+uL)*1.0_RP/2.0_RP+du(m,l,o,:,1,:,:)
            du(m,l,o,:,N+1,:,:)=+(u(m,l,o,:,N+1,:,:)+u(m,l+1,o,:,1,:,:))*1.0_RP/2.0_RP+du(m,l,o,:,N+1,:,:)
          ELSEif(l==nqy) THEN
            if(nq==nqy) then
              uR=u(m,1,o,:,1,:,:)
            else
              uR=urand(2,m,o,:,:,:)
            end if
            du(m,l,o,:,1,:,:)=-(u(m,l,o,:,1,:,:)+u(m,l-1,o,:,N+1,:,:))*1.0_RP/2.0_RP+du(m,l,o,:,1,:,:)
            du(m,l,o,:,N+1,:,:)=+(u(m,l,o,:,N+1,:,:)+uR)*1.0_RP/2.0_RP+du(m,l,o,:,N+1,:,:)
          ELSE
            du(m,l,o,:,1,:,:)=-(u(m,l-1,o,:,N+1,:,:)+u(m,l,o,:,1,:,:))*1.0_RP/2.0_RP+du(m,l,o,:,1,:,:)
            du(m,l,o,:,N+1,:,:)=+(u(m,l,o,:,N+1,:,:)+u(m,l+1,o,:,1,:,:))*1.0_RP/2.0_RP+du(m,l,o,:,N+1,:,:)
          END IF
          CASE(3)
          !z-direction
          IF(o==1) THEN
              if(nq==nqz) then
                uL=u(m,l,nqz,:,:,N+1,:)
              else
                uL=urand(1,m,l,:,:,:)
              endif
            du(m,l,o,:,:,1,:)=-(uL+u(m,l,o,:,:,1,:))*1.0_RP/2.0_RP+du(m,l,o,:,:,1,:)
            du(m,l,o,:,:,N+1,:)=(u(m,l,o,:,:,N+1,:)+u(m,l,o+1,:,:,1,:))*1.0_RP/2.0_RP+du(m,l,o,:,:,N+1,:)
          ELSEif(o==nqz) THEN
              if(nq==nqZ) then
                uR=u(m,l,1,:,:,1,:)
              else
                uR=urand(2,m,l,:,:,:)
              end if
            du(m,l,o,:,:,1,:)=-(u(m,l,o-1,:,:,N+1,:)+u(m,l,o,:,:,1,:))*1.0_RP/2.0_RP+du(m,l,o,:,:,1,:)
            du(m,l,o,:,:,N+1,:)=(u(m,l,o,:,:,N+1,:)+uR)*1.0_RP/2.0_RP+du(m,l,o,:,:,N+1,:)
          ELSE
            du(m,l,o,:,:,1,:)=-(u(m,l,o-1,:,:,N+1,:)+u(m,l,o,:,:,1,:))*1.0_RP/2.0_RP+du(m,l,o,:,:,1,:)
            du(m,l,o,:,:,N+1,:)=(u(m,l,o,:,:,N+1,:)+u(m,l,o+1,:,:,1,:))*1.0_RP/2.0_RP+du(m,l,o,:,:,N+1,:)
          END IF
        end SELECT
        END DO !o
      ENDdo!l
    END DO!m

    DO m=1,nqx
      DO l=1,nqy
        DO o=1,nqz
          DO i=1,n+1
            DO j=1,n+1
              DO k=1,n+1
                SELECT CASE(richtung)
                Case(1)
                du(m,l,o,i,j,k,:)=du(m,l,o,i,j,k,:)*2.0_RP/(w(i)*dx)
                Case(2)
                du(m,l,o,i,j,k,:)=du(m,l,o,i,j,k,:)*2.0_RP/(w(j)*dx)
                Case(3)
                du(m,l,o,i,j,k,:)=du(m,l,o,i,j,k,:)*2.0_RP/(w(k)*dx)
              END SELECT
              END DO !k
            END DO !j
          END DO !i
        END DO !o
      ENDdo!l
    END DO!m
    !print*,du(1,:,:,:,1,1,1)
    !stop
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeGradient(u,n,nq,D,dux,duy,duz)
    IMPLICIT NONE
    !Gleichung 57 im skript
    INTEGER      ,INTENT(IN)                                                        :: n,nq
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:5)       :: u
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:n+1,1:n+1)                                :: D
    REAL(KIND=RP),INTENT(OUT),DIMENSION(1:nq,1:nq,1:nq,1:n+1,1:n+1,1:n+1,1:5)       :: dux,duy,duz
    !local variables
    INTEGER                                                                         :: i,j,k,l,m,o,var

    IF (.not.allocated(w)) THEN
      print*,'w not allocated'
      stop
    ENDif
    DO m=1,nq
      DO l=1,nq
        DO o=1,nq
          DO i=1,n+1
            DO j=1,n+1
              DO k=1,n+1
                DO var=1,5
                  ! TODO (dorn_ni#1#): maybe put dux duy duz together but than rank >7. check for compiler version on cheops
                  !dux(m,l,o,i,j,k,var)=-dot_product(w,matmul(D,u(m,l,o,:,j,k,var)))*(w(j)*w(k))!+surface term
                  !duy(m,l,o,i,j,k,var)=-dot_product(w,matmul(D,u(m,l,o,i,:,k,var)))*(w(i)*w(k))!+surface term
                  !duz(m,l,o,i,j,k,var)=-dot_product(w,matmul(D,u(m,l,o,i,j,:,var)))*(w(j)*w(i))!+surface term

                  dux(m,l,o,i,j,k,var)=-dot_product(D(:,i),u(m,l,o,:,j,k,var)*w)!+surface term
                  duy(m,l,o,i,j,k,var)=-dot_product(D(:,j),u(m,l,o,i,:,k,var)*w)!+surface term
                  duz(m,l,o,i,j,k,var)=-dot_product(D(:,k),u(m,l,o,i,j,:,var)*w)!+surface term
                END DO !var
              END DO !k
            END DO !j
          END DO !i

          !surfaceterms and boundary conditions
          !x-direction

          IF(m==1) THEN
            dux(m,l,o,1,:,:,:)=-(u(m,l,o,1,:,:,:)+u(nq,l,o,N+1,:,:,:))*1.0_RP/2.0_RP+dux(m,l,o,1,:,:,:)
            dux(m,l,o,N+1,:,:,:)=+(u(m,l,o,N+1,:,:,:)+u(m+1,l,o,1,:,:,:))*1.0_RP/2.0_RP+dux(m,l,o,N+1,:,:,:)
          ELSEif(m==nq) THEN
            dux(m,l,o,1,:,:,:)=-(u(m-1,l,o,N+1,:,:,:)+u(m,l,o,1,:,:,:))*1.0_RP/2.0_RP+dux(m,l,o,1,:,:,:)
            dux(m,l,o,N+1,:,:,:)=+(u(m,l,o,N+1,:,:,:)+u(1,l,o,1,:,:,:))*1.0_RP/2.0_RP+dux(m,l,o,N+1,:,:,:)
          ELSE
            dux(m,l,o,1,:,:,:)=-(u(m-1,l,o,N+1,:,:,:)+u(m,l,o,1,:,:,:))*1.0_RP/2.0_RP+dux(m,l,o,1,:,:,:)
            dux(m,l,o,N+1,:,:,:)=+(u(m,l,o,N+1,:,:,:)+u(m+1,l,o,1,:,:,:))*1.0_RP/2.0_RP+dux(m,l,o,N+1,:,:,:)
          END IF
          !y-direction
          IF(l==1) THEN
            duy(m,l,o,:,1,:,:)=-(u(m,l,o,:,1,:,:)+u(m,nq,o,:,N+1,:,:))*1.0_RP/2.0_RP+duy(m,l,o,:,1,:,:)
            duy(m,l,o,:,N+1,:,:)=+(u(m,l,o,:,N+1,:,:)+u(m,l+1,o,:,1,:,:))*1.0_RP/2.0_RP+duy(m,l,o,:,N+1,:,:)
          ELSEif(l==nq) THEN
            duy(m,l,o,:,1,:,:)=-(u(m,l,o,:,1,:,:)+u(m,l-1,o,:,N+1,:,:))*1.0_RP/2.0_RP+duy(m,l,o,:,1,:,:)
            duy(m,l,o,:,N+1,:,:)=+(u(m,l,o,:,N+1,:,:)+u(m,1,o,:,1,:,:))*1.0_RP/2.0_RP+duy(m,l,o,:,N+1,:,:)
          ELSE
            duy(m,l,o,:,1,:,:)=-(u(m,l-1,o,:,N+1,:,:)+u(m,l,o,:,1,:,:))*1.0_RP/2.0_RP+duy(m,l,o,:,1,:,:)
            duy(m,l,o,:,N+1,:,:)=+(u(m,l,o,:,N+1,:,:)+u(m,l+1,o,:,1,:,:))*1.0_RP/2.0_RP+duy(m,l,o,:,N+1,:,:)
          END IF
          !!!!!!!!!!!!!!!duz is wrong!!!!!!!!!!!!!!!!!
          ! TODO (dorn_ni#1#): z direction

          !z-direction
          IF(o==1) THEN
            duz(m,l,o,:,:,1,:)=-(u(m,l,nq,:,:,N+1,:)+u(m,l,o,:,:,1,:))*1.0_RP/2.0_RP+duz(m,l,o,:,:,1,:)
            duz(m,l,o,:,:,N+1,:)=(u(m,l,o,:,:,N+1,:)+u(m,l,o+1,:,:,1,:))*1.0_RP/2.0_RP+duz(m,l,o,:,:,N+1,:)
          ELSEif(o==nq) THEN
            duz(m,l,o,:,:,1,:)=-(u(m,l,o-1,:,:,N+1,:)+u(m,l,o,:,:,1,:))*1.0_RP/2.0_RP+duz(m,l,o,:,:,1,:)
            duz(m,l,o,:,:,N+1,:)=(u(m,l,o,:,:,N+1,:)+u(m,l,1,:,:,1,:))*1.0_RP/2.0_RP+duz(m,l,o,:,:,N+1,:)
          ELSE
            duz(m,l,o,:,:,1,:)=-(u(m,l,o-1,:,:,N+1,:)+u(m,l,o,:,:,1,:))*1.0_RP/2.0_RP+duz(m,l,o,:,:,1,:)
            duz(m,l,o,:,:,N+1,:)=(u(m,l,o,:,:,N+1,:)+u(m,l,o+1,:,:,1,:))*1.0_RP/2.0_RP+duz(m,l,o,:,:,N+1,:)
          END IF
        END DO !o
      ENDdo!l
    END DO!m

    DO m=1,nq
      DO l=1,nq
        DO o=1,nq
          DO i=1,n+1
            DO j=1,n+1
              DO k=1,n+1
                dux(m,l,o,i,j,k,:)=dux(m,l,o,i,j,k,:)*2.0_RP/(w(i)*dx)
                duy(m,l,o,i,j,k,:)=duy(m,l,o,i,j,k,:)*2.0_RP/(w(j)*dx)
                duz(m,l,o,i,j,k,:)=duz(m,l,o,i,j,k,:)*2.0_RP/(w(k)*dx)
              END DO !k
            END DO !j
          END DO !i
        END DO !o
      ENDdo!l
    END DO!m
    !print*,dux(1,1,1,:,1,1,1)
    !print*, xyz(1,1,1,:,1,1,1)*xyz(1,1,1,:,1,1,1)*xyz(1,1,1,:,1,1,1)*xyz(1,1,1,:,1,1,1)*5.0_rp
    !stop
  END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeLviscous(u,dux,duy,duz,D,dir,N,NQ,result)
    !puts all of the viscous components together
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)                                                     :: N,NQ
    INTEGER      ,INTENT(IN)                                                     :: dir
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:NQX,1:NQY,1:NQZ,1:N+1,1:N+1,1:N+1,1:5) :: u,dux,duy,duz
    REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:N+1,1:N+1)                             :: D
    !u DIMENSIONs:(nummer zelle x,zelle y,zelle z,x,y,z,variable)
    REAL(KIND=RP),INTENT(OUT),DIMENSION(1:NQX,1:NQY,1:NQZ,1:N+1,1:N+1,1:N+1,1:5) :: result
    !local Variables
    REAL(KIND=RP),DIMENSION(:,:,:,:,:,:),allocatable :: urand,duxrand,duyrand,duzrand
    INTEGER                                          :: var,k,j,i,o,l,m
    REAL(KIND=RP),DIMENSION(1:N+1,1:5)               :: Fviscous
    REAL(KIND=RP),DIMENSION(1:N+1,1:N+1,1:5)         :: uR,uL,fvisl1,fvisl2,fvisr1,fvisr2
    REAL(KIND=RP),DIMENSION(1:n+1,1:5,1:3)           :: du
    REAL(KIND=RP),DIMENSION(1:n+1,1:n+1,1:5,1:3)     :: duRandl,duRandr

    SELECT CASE(dir)
    CASE(1)
      allocate(urand(1:2,1:NQY,1:NQZ,1:N+1,1:N+1,1:5))
      allocate(duxrand(1:2,1:NQY,1:NQZ,1:N+1,1:N+1,1:5))
      allocate(duyrand(1:2,1:NQY,1:NQZ,1:N+1,1:N+1,1:5))
      allocate(duzrand(1:2,1:NQY,1:NQZ,1:N+1,1:N+1,1:5))
      DO o=1,NQZ
        DO l=1,NQY
          DO m=1,NQX
            DO k=1,N+1
              DO j=1,N+1
                CALL computeviscousFlux(u(m,l,o,:,j,k,:),du,dir,n,Fviscous)
                DO i=1,N+1
                  du(:,:,1)=dux(m,l,o,:,j,k,:)
                  du(:,:,2)=duy(m,l,o,:,j,k,:)
                  du(:,:,3)=duz(m,l,o,:,j,k,:)
                  DO var=1,5
                    result(m,l,o,i,j,k,var)=-dot_product(D(:,i),Fviscous(:,var)*w)
                  ENDDO ! var
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
      if(nq/=NQX) then
        call comm_rand_bed(u,urand,1,N) 
        call comm_rand_bed(dux,duxrand,1,N) 
        call comm_rand_bed(duy,duyrand,1,N) 
        call comm_rand_bed(duz,duzrand,1,N) 
      endif

      DO o=1,NQZ
        DO l=1,NQY
          DO m=1,NQX
            !Randbedingungen
            IF (m==1) THEN
              if(nq==nqx) then
                uL=u(NQX,l,o,N+1,:,:,:)
              duRandl(:,:,:,1)=dux(nq,l,o,N+1,:,:,:)
              duRandl(:,:,:,2)=duy(nq,l,o,N+1,:,:,:)
              duRandl(:,:,:,3)=duz(nq,l,o,N+1,:,:,:)
              else
                uL=urand(1,l,o,:,:,:)
              duRandl(:,:,:,1)=duxrand(1,l,o,:,:,:)
              duRandl(:,:,:,2)=duyrand(1,l,o,:,:,:)
              duRandl(:,:,:,3)=duzrand(1,l,o,:,:,:)
              endif
            ELSE
              uL=u(m-1,l,o,N+1,:,:,:)
              duRandl(:,:,:,1)=dux(m-1,l,o,N+1,:,:,:)
              duRandl(:,:,:,2)=duy(m-1,l,o,N+1,:,:,:)
              duRandl(:,:,:,3)=duz(m-1,l,o,N+1,:,:,:)
            ENDIF
            IF (m==NQX) THEN
              if(nq==nqx) then
              uR=u(1,l,o,1,:,:,:)
              duRandr(:,:,:,1)=dux(1,l,o,1,:,:,:)
              duRandr(:,:,:,2)=duy(1,l,o,1,:,:,:)
              duRandr(:,:,:,3)=duz(1,l,o,1,:,:,:)
              else
              uR=urand(2,l,o,:,:,:)
              duRandr(:,:,:,1)=duxrand(2,l,o,:,:,:)
              duRandr(:,:,:,2)=duyrand(2,l,o,:,:,:)
              duRandr(:,:,:,3)=duzrand(2,l,o,:,:,:)
              end if
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
            result(m,l,o,1,:,:,:)=result(m,l,o,1,:,:,:)-(Fvisl1+Fvisr1)*0.5_RP
            result(m,l,o,N+1,:,:,:)=result(m,l,o,N+1,:,:,:)+(Fvisl2+Fvisr2)*0.5_RP
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
      DO m=1,nqX
        DO l=1,nqY
          DO o=1,nqZ
            DO i=1,n+1
              DO j=1,n+1
                DO k=1,n+1
                  result(m,l,o,i,j,k,:)=result(m,l,o,i,j,k,:)/(w(i))
                END DO !k
              END DO !j
            END DO !i
          END DO !o
        ENDdo!l
      END DO!m
    CASE(2)
      allocate(urand(1:2,1:NQX,1:NQZ,1:N+1,1:N+1,1:5))
      allocate(duxrand(1:2,1:NQX,1:NQZ,1:N+1,1:N+1,1:5))
      allocate(duyrand(1:2,1:NQX,1:NQZ,1:N+1,1:N+1,1:5))
      allocate(duzrand(1:2,1:NQX,1:NQZ,1:N+1,1:N+1,1:5))
      DO o=1,NQZ
        DO l=1,NQY
          DO m=1,NQX
            DO k=1,N+1
              DO j=1,N+1
                DO i=1,N+1
                  du(:,:,1)=dux(m,l,o,i,:,k,:)
                  du(:,:,2)=duy(m,l,o,i,:,k,:)
                  du(:,:,3)=duz(m,l,o,i,:,k,:)
                  CALL computeviscousFlux(u(m,l,o,i,:,k,:),du,dir,n,Fviscous)
                  DO var=1,5 !! besser
                    result(m,l,o,i,j,k,var)=-dot_product(D(:,j),Fviscous(:,var)*w)
                  ENDDO ! var
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
      if(nq/=NQY) then
        call comm_rand_bed(u,urand,2,N) 
        call comm_rand_bed(dux,duxrand,2,N) 
        call comm_rand_bed(duy,duyrand,2,N) 
        call comm_rand_bed(duz,duzrand,2,N) 
      endif
            !Randbedingungen
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
      DO o=1,NQZ
        DO l=1,NQY
          DO m=1,NQX
            IF (l==1) THEN
              if(nq==nqy) then
              uL=u(m,nqy,o,:,N+1,:,:)
              duRandl(:,:,:,1)=dux(m,nq,o,:,N+1,:,:)
              duRandl(:,:,:,2)=duy(m,nq,o,:,N+1,:,:)
              duRandl(:,:,:,3)=duz(m,nq,o,:,N+1,:,:)
              else
                uL=urand(1,m,o,:,:,:)
              duRandl(:,:,:,1)=duxrand(1,m,o,:,:,:)
              duRandl(:,:,:,2)=duyrand(1,m,o,:,:,:)
              duRandl(:,:,:,3)=duzrand(1,m,o,:,:,:)
              endif
            ELSE
              uL=u(m,l-1,o,:,N+1,:,:)
              duRandl(:,:,:,1)=dux(m,l-1,o,:,N+1,:,:)
              duRandl(:,:,:,2)=duy(m,l-1,o,:,N+1,:,:)
              duRandl(:,:,:,3)=duz(m,l-1,o,:,N+1,:,:)
            ENDIF
            IF (l==NQY) THEN
              if(nq==nqy) then
              uR=u(m,1,o,:,1,:,:)
              duRandr(:,:,:,1)=dux(m,1,o,:,1,:,:)
              duRandr(:,:,:,2)=duy(m,1,o,:,1,:,:)
              duRandr(:,:,:,3)=duz(m,1,o,:,1,:,:)
              else
              uR=urand(2,m,o,:,:,:)
              duRandr(:,:,:,1)=duxrand(2,m,o,:,:,:)
              duRandr(:,:,:,2)=duyrand(2,m,o,:,:,:)
              duRandr(:,:,:,3)=duzrand(2,m,o,:,:,:)
              end if
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
            result(m,l,o,:,1,:,:)=result(m,l,o,:,1,:,:)-(Fvisl1+Fvisr1)*0.5_RP
            result(m,l,o,:,N+1,:,:)=result(m,l,o,:,N+1,:,:)+(Fvisl2+Fvisr2)*0.5_RP
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
      DO m=1,nqX
        DO l=1,nqY
          DO o=1,nqZ
            DO i=1,n+1
              DO j=1,n+1
                DO k=1,n+1
                  result(m,l,o,i,j,k,:)=result(m,l,o,i,j,k,:)/(w(j))
                END DO !k
              END DO !j
            END DO !i
          END DO !o
        ENDdo!l
      END DO!m
    CASE(3)
      allocate(urand(1:2,1:NQX,1:NQY,1:N+1,1:N+1,1:5))
      allocate(duxrand(1:2,1:NQX,1:NQY,1:N+1,1:N+1,1:5))
      allocate(duyrand(1:2,1:NQX,1:NQY,1:N+1,1:N+1,1:5))
      allocate(duzrand(1:2,1:NQX,1:NQY,1:N+1,1:N+1,1:5))
      DO o=1,NQZ
        DO l=1,NQY
          DO m=1,NQX
            DO k=1,N+1
              DO j=1,N+1
                DO i=1,N+1
                  du(:,:,1)=dux(m,l,o,i,j,:,:)
                  du(:,:,2)=duy(m,l,o,i,j,:,:)
                  du(:,:,3)=duz(m,l,o,i,j,:,:)
                  CALL computeviscousFlux(u(m,l,o,i,j,:,:),du,dir,n,Fviscous)
                  DO var=1,5 !! besser
                    result(m,l,o,i,j,k,var)=-dot_product(D(:,k),Fviscous(:,var)*w)
                  ENDDO ! var
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
      if(nq/=NQZ) then
        call comm_rand_bed(u,urand,3,N) 
        call comm_rand_bed(dux,duxrand,3,N) 
        call comm_rand_bed(duy,duyrand,3,N) 
        call comm_rand_bed(duz,duzrand,3,N) 
      endif
            !Randbedingungen
      DO o=1,NQZ
        DO l=1,NQY
          DO m=1,NQX
            IF (o==1) THEN
              if(nq==nqZ) then
              uL=u(m,l,nq,:,:,N+1,:)
              duRandl(:,:,:,1)=dux(m,l,nq,:,:,N+1,:)
              duRandl(:,:,:,2)=duy(m,l,nq,:,:,N+1,:)
              duRandl(:,:,:,3)=duz(m,l,nq,:,:,N+1,:)
              else
                uL=urand(1,m,l,:,:,:)
              duRandl(:,:,:,1)=duxrand(1,m,l,:,:,:)
              duRandl(:,:,:,2)=duyrand(1,m,l,:,:,:)
              duRandl(:,:,:,3)=duzrand(1,m,l,:,:,:)
              endif
            ELSE
              uL=u(m,l,o-1,:,:,N+1,:)
              duRandl(:,:,:,1)=dux(m,l,o-1,:,:,N+1,:)
              duRandl(:,:,:,2)=duy(m,l,o-1,:,:,N+1,:)
              duRandl(:,:,:,3)=duz(m,l,o-1,:,:,N+1,:)
            ENDIF
            IF (o==NQZ) THEN
              if(nq==nqZ) then
              uR=u(m,l,1,:,:,1,:)
              duRandr(:,:,:,1)=dux(m,l,1,:,:,1,:)
              duRandr(:,:,:,2)=duy(m,l,1,:,:,1,:)
              duRandr(:,:,:,3)=duz(m,l,1,:,:,1,:)
              else
              uR=urand(2,m,l,:,:,:)
              duRandr(:,:,:,1)=duxrand(2,m,l,:,:,:)
              duRandr(:,:,:,2)=duyrand(2,m,l,:,:,:)
              duRandr(:,:,:,3)=duzrand(2,m,l,:,:,:)
              end if
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
            result(m,l,o,:,:,1,:)=result(m,l,o,:,:,1,:)-(Fvisl1+Fvisr1)*0.5_RP
            result(m,l,o,:,:,N+1,:)=result(m,l,o,:,:,N+1,:)+(Fvisl2+Fvisr2)*0.5_RP
          ENDDO ! m
        ENDDO ! l
      ENDDO ! o
      DO m=1,nqx
        DO l=1,nqy
          DO o=1,nqZ
            DO i=1,n+1
              DO j=1,n+1
                DO k=1,n+1
                  result(m,l,o,i,j,k,:)=result(m,l,o,i,j,k,:)/(w(k))
                END DO !k
              END DO !j
            END DO !i
          END DO !o
        ENDdo!l
      END DO!m
    END SELECT
  END SUBROUTINE computeLviscous
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computeviscousFluxRand(u,du,dir,n,result)
    !computes the viscous part of the flux analytically
    IMPLICIT NONE
    INTEGER       ,INTENT(IN)                                 :: N,dir
    REAL(KIND=RP) ,INTENT(IN) ,DIMENSION(1:n+1,1:n+1,1:5)     :: u
    REAL(KIND=RP) ,INTENT(IN) ,DIMENSION(1:n+1,1:n+1,1:5,1:3) :: du !(Punkte x Punkte x Variable x Ableitungsrichtung)
    REAL(KIND=RP) ,INTENT(OUT),DIMENSION(1:n+1,1:n+1,1:5)     :: result
    REAL(KIND=RP)             ,DIMENSION(1:N+1,1:N+1)         :: P
    REAL(KIND=RP)             ,DIMENSION(1:n+1,1:n+1,1:3)     :: dTemp,dv1,dv2,dv3,dP,drv1,drv2,drv3
    result(:,:,1)=0.0_RP
    p=(gamma-1.0_RP)*(u(:,:,5)-0.5_RP*(u(:,:,2)*u(:,:,2)+u(:,:,3)*u(:,:,3)+u(:,:,4)*u(:,:,4))/u(:,:,1))
    !dP(:,:,1)=(gamma-1.0_RP)*du(:,:,1,1)*(2.0_RP*u(:,:,1)-3.0_RP/2.0_RP)
    !dP(:,:,2)=(gamma-1.0_RP)*du(:,:,1,2)*(2.0_RP*u(:,:,1)-3.0_RP/2.0_RP)
    !dP(:,:,3)=(gamma-1.0_RP)*du(:,:,1,3)*(2.0_RP*u(:,:,1)-3.0_RP/2.0_RP)
    !dv1(:,:,1)=(du(:,:,1,1)*u(:,:,2)/u(:,:,1))/u(:,:,1)
    !dv1(:,:,2)=(du(:,:,1,2)*u(:,:,2)/u(:,:,1))/u(:,:,1)
    !dv1(:,:,3)=(du(:,:,1,3)*u(:,:,2)/u(:,:,1))/u(:,:,1)
    !dv2(:,:,1)=(du(:,:,1,1)*u(:,:,3)/u(:,:,1))/u(:,:,1)
    !dv2(:,:,2)=(du(:,:,1,2)*u(:,:,3)/u(:,:,1))/u(:,:,1)
    !dv2(:,:,3)=(du(:,:,1,3)*u(:,:,3)/u(:,:,1))/u(:,:,1)
    !dv3(:,:,1)=(du(:,:,1,1)*u(:,:,4)/u(:,:,1))/u(:,:,1)
    !dv3(:,:,2)=(du(:,:,1,2)*u(:,:,4)/u(:,:,1))/u(:,:,1)
    !dv3(:,:,3)=(du(:,:,1,3)*u(:,:,4)/u(:,:,1))/u(:,:,1)
    dv1(:,:,1)=(du(:,:,2,1)-u(:,:,2)*du(:,:,1,1))/u(:,:,1)
    dv1(:,:,2)=(du(:,:,2,2)-u(:,:,2)*du(:,:,1,2))/u(:,:,1)
    dv1(:,:,3)=(du(:,:,2,3)-u(:,:,2)*du(:,:,1,3))/u(:,:,1)
    dv2(:,:,1)=(du(:,:,3,1)-u(:,:,3)*du(:,:,1,1))/u(:,:,1)
    dv2(:,:,2)=(du(:,:,3,2)-u(:,:,3)*du(:,:,1,2))/u(:,:,1)
    dv2(:,:,3)=(du(:,:,3,3)-u(:,:,3)*du(:,:,1,3))/u(:,:,1)
    dv3(:,:,1)=(du(:,:,4,1)-u(:,:,4)*du(:,:,1,1))/u(:,:,1)
    dv3(:,:,2)=(du(:,:,4,2)-u(:,:,4)*du(:,:,1,2))/u(:,:,1)
    dv3(:,:,3)=(du(:,:,4,3)-u(:,:,4)*du(:,:,1,3))/u(:,:,1)
    drv1(:,:,1)=du(:,:,1,1)*(u(:,:,2)/u(:,:,1))**2+2.0_RP*u(:,:,2)*dv1(:,:,1)
    drv1(:,:,2)=du(:,:,1,2)*(u(:,:,2)/u(:,:,1))**2+2.0_RP*u(:,:,2)*dv1(:,:,2)
    drv1(:,:,3)=du(:,:,1,3)*(u(:,:,2)/u(:,:,1))**2+2.0_RP*u(:,:,2)*dv1(:,:,3)
    drv2(:,:,1)=du(:,:,1,1)*(u(:,:,3)/u(:,:,1))**2+2.0_RP*u(:,:,3)*dv2(:,:,1)
    drv2(:,:,2)=du(:,:,1,2)*(u(:,:,3)/u(:,:,1))**2+2.0_RP*u(:,:,3)*dv2(:,:,2)
    drv2(:,:,3)=du(:,:,1,3)*(u(:,:,3)/u(:,:,1))**2+2.0_RP*u(:,:,3)*dv2(:,:,3)
    drv3(:,:,1)=du(:,:,1,1)*(u(:,:,4)/u(:,:,1))**2+2.0_RP*u(:,:,4)*dv3(:,:,1)
    drv3(:,:,2)=du(:,:,1,2)*(u(:,:,4)/u(:,:,1))**2+2.0_RP*u(:,:,4)*dv3(:,:,2)
    drv3(:,:,3)=du(:,:,1,3)*(u(:,:,4)/u(:,:,1))**2+2.0_RP*u(:,:,4)*dv3(:,:,3)
    dP(:,:,1)=(gamma-1.0_RP)*(du(:,:,5,1)-0.5_RP*(drv1(:,:,1)+drv2(:,:,1)+drv3(:,:,1)))
    dP(:,:,2)=(gamma-1.0_RP)*(du(:,:,5,1)-0.5_RP*(drv1(:,:,2)+drv2(:,:,2)+drv3(:,:,2)))
    dP(:,:,3)=(gamma-1.0_RP)*(du(:,:,5,1)-0.5_RP*(drv1(:,:,3)+drv2(:,:,3)+drv3(:,:,3)))
    dTemp(:,:,1)=(dp(:,:,1)*u(:,:,1)-du(:,:,1,1)*p)/u(:,:,1)**2
    dTemp(:,:,2)=(dp(:,:,2)*u(:,:,1)-du(:,:,1,2)*p)/u(:,:,1)**2
    dTemp(:,:,3)=(dp(:,:,3)*u(:,:,1)-du(:,:,1,3)*p)/u(:,:,1)**2
    SELECT CASE(dir)
    CASE(1)
      result(:,:,2)=mu*(2*dv1(:,:,1)-2.0_RP/3.0_RP*(dv1(:,:,1)+dv2(:,:,2)+dv3(:,:,3)))
      result(:,:,3)=mu*(dv1(:,:,2)+dv2(:,:,1))
      result(:,:,4)=mu*(dv1(:,:,3)+dv3(:,:,1))
      result(:,:,5)=mu*(u(:,:,2)/u(:,:,1)*(2.0_rp*dv1(:,:,1)-2.0_RP/3.0_RP*(dv1(:,:,1)+dv2(:,:,2)+dv3(:,:,3)))+&
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
    INTEGER       ,INTENT(IN)                           :: N,dir
    REAL(KIND=RP) ,INTENT(IN) ,DIMENSION(1:n+1,1:5)     :: u
    REAL(KIND=RP) ,INTENT(IN) ,DIMENSION(1:n+1,1:5,1:3) :: du !(Punkte x Variable x Ableitungsrichtung)
    REAL(KIND=RP) ,INTENT(OUT),DIMENSION(1:n+1,1:5)     :: result
    REAL(KIND=RP)             ,DIMENSION(1:N+1)         :: P
    REAL(KIND=RP)             ,DIMENSION(1:n+1,1:3)     :: dTemp,dv1,dv2,dv3,dp,drv1,drv2,drv3
    result(:,1)=0.0_RP
    p=(gamma-1.0_RP)*(u(:,5)-0.5_RP*(u(:,2)*u(:,2)+u(:,3)*u(:,3)+u(:,4)*u(:,4))/u(:,1))
    !dP(:,1)=(gamma-1.0_RP)*du(:,1,1)*(2.0_RP*u(:,1)-3.0_RP/2.0_RP)
    !dP(:,2)=(gamma-1.0_RP)*du(:,1,2)*(2.0_RP*u(:,1)-3.0_RP/2.0_RP)
    !dP(:,3)=(gamma-1.0_RP)*du(:,1,3)*(2.0_RP*u(:,1)-3.0_RP/2.0_RP)
    !dTemp(:,1)=(dp(:,1)*u(:,1)-du(:,1,1)*p)/u(:,1)**2
    !dTemp(:,2)=(dp(:,2)*u(:,1)-du(:,1,2)*p)/u(:,1)**2
    !dTemp(:,3)=(dp(:,3)*u(:,1)-du(:,1,3)*p)/u(:,1)**2
    ! Fehler glaube ich
    ! dv1(:,1)=(du(:,1,1)*u(:,2)/u(:,1))/u(:,1)
    ! dv1(:,2)=(du(:,1,2)*u(:,2)/u(:,1))/u(:,1)
    ! dv1(:,3)=(du(:,1,3)*u(:,2)/u(:,1))/u(:,1)
    ! dv2(:,1)=(du(:,1,1)*u(:,3)/u(:,1))/u(:,1)
    ! dv2(:,2)=(du(:,1,2)*u(:,3)/u(:,1))/u(:,1)
    ! dv2(:,3)=(du(:,1,3)*u(:,3)/u(:,1))/u(:,1)
    ! dv3(:,1)=(du(:,1,1)*u(:,4)/u(:,1))/u(:,1)
    ! dv3(:,2)=(du(:,1,2)*u(:,4)/u(:,1))/u(:,1)
    ! dv3(:,3)=(du(:,1,3)*u(:,4)/u(:,1))/u(:,1)
    dv1(:,1)=(du(:,2,1)-u(:,2)*du(:,1,1))/u(:,1)
    dv1(:,2)=(du(:,2,2)-u(:,2)*du(:,1,2))/u(:,1)
    dv1(:,3)=(du(:,2,3)-u(:,2)*du(:,1,3))/u(:,1)
    dv2(:,1)=(du(:,3,1)-u(:,3)*du(:,1,1))/u(:,1)
    dv2(:,2)=(du(:,3,2)-u(:,3)*du(:,1,2))/u(:,1)
    dv2(:,3)=(du(:,3,3)-u(:,3)*du(:,1,3))/u(:,1)
    dv3(:,1)=(du(:,4,1)-u(:,4)*du(:,1,1))/u(:,1)
    dv3(:,2)=(du(:,4,2)-u(:,4)*du(:,1,2))/u(:,1)
    dv3(:,3)=(du(:,4,3)-u(:,4)*du(:,1,3))/u(:,1)
    drv1(:,1)=du(:,1,1)*(u(:,2)/u(:,1))**2+2.0_RP*u(:,2)*dv1(:,1)
    drv1(:,2)=du(:,1,2)*(u(:,2)/u(:,1))**2+2.0_RP*u(:,2)*dv1(:,2)
    drv1(:,3)=du(:,1,3)*(u(:,2)/u(:,1))**2+2.0_RP*u(:,2)*dv1(:,3)
    drv2(:,1)=du(:,1,1)*(u(:,3)/u(:,1))**2+2.0_RP*u(:,3)*dv2(:,1)
    drv2(:,2)=du(:,1,2)*(u(:,3)/u(:,1))**2+2.0_RP*u(:,3)*dv2(:,2)
    drv2(:,3)=du(:,1,3)*(u(:,3)/u(:,1))**2+2.0_RP*u(:,3)*dv2(:,3)
    drv3(:,1)=du(:,1,1)*(u(:,4)/u(:,1))**2+2.0_RP*u(:,4)*dv3(:,1)
    drv3(:,2)=du(:,1,2)*(u(:,4)/u(:,1))**2+2.0_RP*u(:,4)*dv3(:,2)
    drv3(:,3)=du(:,1,3)*(u(:,4)/u(:,1))**2+2.0_RP*u(:,4)*dv3(:,3)
    dP(:,1)=(gamma-1.0_RP)*(du(:,5,1)-0.5_RP*(drv1(:,1)+drv2(:,1)+drv3(:,1)))
    dP(:,2)=(gamma-1.0_RP)*(du(:,5,1)-0.5_RP*(drv1(:,2)+drv2(:,2)+drv3(:,2)))
    dP(:,3)=(gamma-1.0_RP)*(du(:,5,1)-0.5_RP*(drv1(:,3)+drv2(:,3)+drv3(:,3)))
    dTemp(:,1)=(dp(:,1)*u(:,1)-du(:,1,1)*p)/u(:,1)**2
    dTemp(:,2)=(dp(:,2)*u(:,1)-du(:,1,2)*p)/u(:,1)**2
    dTemp(:,3)=(dp(:,3)*u(:,1)-du(:,1,3)*p)/u(:,1)**2
    SELECT CASE(dir)
    CASE(1)
      result(:,2)=mu*(2.0_RP*dv1(:,1)-2.0_RP/3.0_RP*(dv1(:,1)+dv2(:,2)+dv3(:,3)))
      result(:,3)=mu*(dv1(:,2)+dv2(:,1))
      result(:,4)=mu*(dv1(:,3)+dv3(:,1))
      result(:,5)=mu*(u(:,2)/u(:,1)*(2.0_rp*dv1(:,1)-2.0_RP/3.0_RP*(dv1(:,1)+dv2(:,2)+dv3(:,3)))+&
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
    ! u(:,:,:,:,:,:,1)=xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)
    ! u(:,:,:,:,:,:,2)=xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)
    ! u(:,:,:,:,:,:,3)=xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)
    ! u(:,:,:,:,:,:,4)=xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)
    ! u(:,:,:,:,:,:,5)=xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*xyz(:,:,:,:,:,:,1)*3.0_RP
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
    DEALLOCATE(xyz)
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
    IF (any(gamma*p/u(:,:,:,:,:,:,1)<0)) THEN
      print*,'kacke'
    ENDif
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
 ! SUBROUTINE RungeKutta5explizit(u,nq,n,numvar,dt,Dval,t,whichflux,vis)
 !   IMPLICIT NONE
 !   !u=bekannte Werte
 !   INTEGER      ,INTENT(IN)                                                         :: numvar,n,nq
 !   CHARACTER(len=2),INTENT(IN)                                                      :: whichflux,vis
 !   REAL(KIND=RP),INTENT(IN)   ,DIMENSION(1:N+1,1:N+1)                               :: Dval
 !   REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:nq,1:nq,1:nq,1:N+1,1:N+1,1:N+1,1:numvar) :: u
 !   !local
 !   REAL(KIND=RP),DIMENSION(1:nq,1:nq,1:nq,1:(N+1),1:N+1,1:N+1,1:numvar)    :: g
 !   INTEGER                                                                 :: step
 !   REAL(KIND=RP),DIMENSION(5)                                              :: a,b,c
 !   REAL(KIND=RP),INTENT(IN)                                                :: dt,t
 !   g=Rmanu(u,n,nq,Dval,t,whichflux,vis,dt)
 !   a=(/0.0_rp, -567301805773.0_rp/1357537059087.0_rp,&
 !     -2404267990393.0_rp/2016746695238.0_rp, -3550918686646.0_rp/2091501179385.0_rp,&
 !     -1275806237668.0_rp/842570457699.0_rp/)
 !   b=(/0.0_rp, 1432997174477.0_rp/9575080441755.0_rp,&
 !     2526269341429.0_rp/6820363962896.0_rp, 2006345519317.0_rp/3224310063776.0_rp,&
 !     2802321613138.0_rp/2924317926251.0_rp /)
 !   c=(/1432997174477.0_rp/9575080441755.0_rp, 5161836677717.0_rp/13612068292357.0_rp,&
 !     1720146321549.0_rp/2090206949498.0_rp, 3134564353537.0_rp/4481467310338.0_rp,&
 !     2277821191437.0_rp/14882151754819.0_rp /)



 !   DO step=1,5
 !     g=a(step)*g+Rmanu(u,n,nq,Dval,t+b(step)*dt,whichflux,vis,b(step)*dt)
 !     u=u+c(step)*dt*g
 !   ENDDO ! step
 ! END SUBROUTINE RungeKutta5explizit
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE RungeKutta5explizitParallel(u,nq,n,numvar,dt,Dval,t,whichflux,vis,xyzsub)
    IMPLICIT NONE
    !u=bekannte Werte
    INTEGER      ,INTENT(IN)                                                            :: numvar,n,nq
    CHARACTER(len=2),INTENT(IN)                                                         :: whichflux,vis
    REAL(KIND=RP),INTENT(IN)   ,DIMENSION(1:N+1,1:N+1)                                  :: Dval
    REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:nqx,1:nqy,1:nqZ,1:N+1,1:N+1,1:N+1,1:numvar) :: u
    REAL(KIND=RP),DIMENSION(1:nqx,1:nqy,1:nqZ,1:N+1,1:N+1,1:N+1,1:numvar)               :: g
    REAL(KIND=RP),DIMENSION(1:nqx,1:nqy,1:nqZ,1:N+1,1:N+1,1:N+1,1:3)                    :: xyzsub
    REAL(KIND=RP),INTENT(IN)                                                            :: dt,t
    !local
    INTEGER                                                                             :: step
    INTEGER                                                                             :: i
    REAL(KIND=RP),DIMENSION(5)                                                          :: a,b,c

     
    
    g=Rmanu(u,n,nq,Dval,t,whichflux,vis,dt,xyzsub)
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
      g=a(step)*g+Rmanu(u,n,nq,Dval,t+b(step)*dt,whichflux,vis,b(step)*dt,xyzsub)
      u=u+c(step)*dt*g
    ENDDO ! step

  END SUBROUTINE RungeKutta5explizitParallel
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Driver_zeit (u,nq,n,numvar,dt,D,t,tend,whichflux,vis,CFL,xnum,ynum,znum)
    implicit none
    INTEGER      ,INTENT(IN)                                                         :: numvar,n,nq,xnum,ynum,znum
    CHARACTER(len=2),INTENT(IN)                                                      :: whichflux,vis
    REAL(KIND=RP),INTENT(IN)   ,DIMENSION(1:N+1,1:N+1)                               :: D
    REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:nq,1:nq,1:nq,1:N+1,1:N+1,1:N+1,1:numvar) :: u
    REAL(KIND=RP),DIMENSION(:,:,:,:,:,:,:),allocatable :: usub,xyzsub
    REAL(KIND=RP)                                                         :: dt,t,tend,dtz
    REAL(KIND=RP) :: a,CFL
    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr) !getting the number of processors
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierr)
    !we assume nq**3 is divisible by num_procs

    IF(mod(nq,num_procs)/=0)then
      print*, 'nicht teilbar'
    ENDif 
    IF(mod(nq**3,num_procs)/=0) THEN
      call MPI_Finalize(ierr)
      print*,'Number of elements(',nq**3,') is not divided by ',num_procs
    END IF
    NQX=NQ/xnum
    NQY=NQ/ynum
    NQZ=NQ/znum
    ALLOCATE(part_map(0:num_procs-1,0:5))
    ALLOCATE(procs_map(0:xnum-1,0:ynum-1,0:znum-1))
    allocate(usub(1:nqx,1:nqy,1:nqz,1:n+1,1:n+1,1:n+1,1:5))
    allocate(xyzsub(1:nqx,1:nqy,1:nqz,1:n+1,1:n+1,1:n+1,1:3))
    call split_Gebiet(u,usub,xyzsub,xnum,ynum,znum,n,nq)

    DO while(tend-t>epsilon(dt))

     if(id==0) then 
      print*,'t'
      print*,t
      print*,'dt'
      print*,dt
      print*,'sum(energy)'
      print*,sum(U(:,:,:,:,:,:,5))
      call lambdaMaxGlobal(u,a,NQ,N)
      dt=CFL/(3.0_RP*a)*(dx/real(2*N+1))**2
    ! print*, 'id'
    ! print*, id
    ! print*, 'Dt'
    ! print*, dt
      call MPI_SEND(dt,1,MPI_DOUBLE_PRECISION,1,23,MPI_COMM_WORLD,ierr)
    else
      !!! andern für mehr als 2
      call MPI_RECV(dtz,1,MPI_DOUBLE_PRECISION,0,23,MPI_COMM_WORLD,stat,ierr)
      dt=dtz
    end if

      IF(t+dt>tend) dt=tend-t
      call RungeKutta5explizitParallel(usub,nq,n,5,dt,D,t,whichflux,vis,xyzsub)
      t=t+dt
      call collect_solution(u,usub,xnum,ynum,znum,n,nq)
    END DO
   DEALLOCATE(xyzsub,usub,part_map) 
    DEALLOCATE(x,w,xmit,xges)
    DEALLOCATE(procs_map)
    end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collect_solution(u,usub,xnum,ynum,znum,n,nq)
  implicit none
  REAL(KIND=RP),INTENT(INout),DIMENSION(1:nq,1:nq,1:nq,1:N+1,1:N+1,1:N+1,1:5)       :: u
  REAL(KIND=RP),INTENT(IN),DIMENSION(1:nqx,1:nqy,1:nqZ,1:N+1,1:N+1,1:N+1,1:5) :: usub
  INTEGER                                                                        :: m,l,o,i,xnum,ynum,znum
  INTEGER                                                                        :: n,nq

  !procs_map speichert die Position der processors im gebiet
  if(id/=0) then
      call MPI_SSEND(usub,NQX*NQY*NQZ*(N+1)**3*5,MPI_DOUBLE_PRECISION,0,id,MPI_COMM_WORLD,ierr)
  endif
  if(id==0) then
  i=0
  do m=0,xnum-1
    do l=0,ynum-1
      do o=0,znum-1
        if(i==0) then 
         i=i+1
          continue
        else
      call MPI_RECV(u(m*NQX+1:(m+1)*NQX,l*NQY+1:(l+1)*NQY,o*NQZ+1:(o+1)*NQZ,:,:,:,:)&
        ,NQX*NQY*NQZ*(N+1)**3*5,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,stat,ierr)
         i=i+1
       end if
      enddo 
    enddo 
  enddo 
  
  u(0*NQX+1:(0+1)*NQX,0*NQY+1:(0+1)*NQY,0*NQZ+1:(0+1)*NQZ,:,:,:,:)=usub
end if
  
  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine split_Gebiet(u,usub,xyzsub,xnum,ynum,znum,n,nq)
  implicit none
  REAL(KIND=RP),INTENT(IN),DIMENSION(1:nq,1:nq,1:nq,1:N+1,1:N+1,1:N+1,1:5)       :: u
  REAL(KIND=RP),INTENT(INout),DIMENSION(1:nqx,1:nqy,1:nqZ,1:N+1,1:N+1,1:N+1,1:5) :: usub
  REAL(KIND=RP),INTENT(INout),DIMENSION(1:nqx,1:nqy,1:nqZ,1:N+1,1:N+1,1:N+1,1:3) :: xyzsub
  INTEGER                                                                        :: m,l,o,i,xnum,ynum,znum,mmin1,lmin1,omin1
  INTEGER                                                                        :: mplus1,lplus1,oplus1,n,nq

  !procs_map speichert die Position der processors im gebiet
  i=0
  do m=0,xnum-1
    do l=0,ynum-1
      do o=0,znum-1
        procs_map(m,l,o)=i
       if(id==i) then
       usub=u(m*NQX+1:(m+1)*NQX,l*NQY+1:(l+1)*NQY,o*NQZ+1:(o+1)*NQZ,:,:,:,:)
       xyzsub=xyz(m*NQX+1:(m+1)*NQX,l*NQY+1:(l+1)*NQY,o*NQZ+1:(o+1)*NQZ,:,:,:,:)
       endif 
        i=i+1
      enddo 
    enddo 
  enddo 
  
  i=0
  do m=0,xnum-1
    do l=0,ynum-1
      do o=0,znum-1 
        if(m==0) then
          mmin1=xnum-1
        else
          mmin1=m-1 
        endif
        if(l==0) then
          lmin1=ynum-1
        else
          lmin1=l-1 
        endif
        if(o==0) then
          omin1=znum-1
        else
          omin1=o-1 
        endif
        if(m==xnum-1) then
          mplus1=0
        else
          mplus1=m+1 
        endif
        if(l==ynum-1) then
          lplus1=0
        else
          lplus1=l+1 
        endif
        if(o==znum-1) then
          oplus1=0
        else
          oplus1=o+1 
        endif
          part_map(i,0)=procs_map(mmin1,l,o)
          part_map(i,1)=procs_map(mplus1,l,o)
          part_map(i,2)=procs_map(m,lmin1,o)
          part_map(i,3)=procs_map(m,lplus1,o)
          part_map(i,4)=procs_map(m,l,omin1)
          part_map(i,5)=procs_map(m,l,oplus1)
          
        i=i+1
      enddo 
    enddo 
  enddo 
  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  subroutine beginne_para ()
    implicit none
    
    call MPI_Init(ierr) !starting MPI
    end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ende_para ()
    implicit none
    
    call MPI_Finalize(ierr)
    end subroutine
END MODULE Zeitintegration
