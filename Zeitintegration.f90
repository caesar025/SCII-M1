!!
!  Overall there were several bugs in this portion of the code that I changed when compared to the "master"
!  branch. Also, I include a few questions about the implementation:
!
!  1) Overall, I now have the code compiling and it runs without producing NaN in the energy when using the 
!     Pirozolli flux. However, I did not debug the code at all. I just changed how the memory is managed as
!     well as some source code formatting issues. You still need to run comprehensive convergence tests to
!     make sure the code is really debugged.
!
!  2) Chief amoung them was a mismanagement of memory and how you pass arrays (both the solution and the
!     derivative matrix). Passing an arguement like DIMENSION(:,:) can be unreliable as you rely on the 
!     compiler to manage the memory and make sure that the "right" information is passed. This should be
!     avoided, always be explicit in the array dimensions that you pass to subroutines. For example, the 
!     derivative matrix would be passed as (1:N+1,1:N+1). Also, NEVER pass somthing as an allocatable! The
!     original code often did this and would allocate arrays (like result) and pass them back. But this 
!     produces problems with "scope", in particular you allocate memory and never release it which can
!     cause seg faults because eventually the progra will run out of memory from the heap. I did my best to 
!     alter the code and handle memory correctly, where subroutines are always passed the correct dimension
!     but there may still be issues
!
!  3) The calculation of the average pressure was wrong in the PI flux
!
!  4) The enthalpy h was incorrect in the PI flux routines. It doesn't depend on the sound speed c. Oddly,
!     the computation of h was correctly done in the lambdaMax routine. So this seems to be a copy/paste 
!     error
!
!  5) You can simplify the PI computation and reuse the value of result(:,1), just a hint
!
!  6) Does the standard DG split form converge correctly? If so, this would rule out a bug in the
!     computeFsharp routine and point to an issue in the PI routine.
!
!  7) As a comment, don't play so fast and loose with ALLOCATABLE arrays. They should be used carefully. As
!     such, my suggested modifications only declares allocatables in the driver and deallocates at the end.
!     Otherwise, subroutines and functions are always passed the array size as an arguement. This ensures
!     that you are always aware and in control of how memory is declared and managed.
!
!  8) Your implementation of the split forms is interesting and a tactic I never considered. That is, passing
!     slices of the fluxes and the computing the volume contribution. I like it! Is convergence working??
!
!  9) It appears that when you use the PI flux you never compute FL and FR to create the dissipation term.
!
! 10) Again, there we many crimes committed when managing memory, especially for the array "result" (in several routines). It
!     was often allocated and the passed back from the routine as INTENT(OUT). However, this can be unreliable because of the
!     scope with respect to memory. It is better to always pass an array with a known size, especially if it's INTENT(IN). This
!     is because if FORTRAN sees an INTENT(IN) as an ALLOCATABLE it might just fill that array with "garbage" values rather 
!     than what was intended, e.g. the derivative matrix, because it thinks that this is a "new" array rather than 
!     something that is precomputed and stored.
!
!!
MODULE Zeitintegration
    USE Quadraturroutinen
    REAL(KIND=RP),DIMENSION(:),ALLOCATABLE :: x,w,xmit,xges
    REAL(KIND=RP)                          :: gk=9.812_RP,dx,gamma=1.4_RP
!
CONTAINS
!
    SUBROUTINE Vorbereiten(N,NQ,Dval)
      IMPLICIT NONE
      INTEGER      ,INTENT(IN)                         :: N,NQ
      REAL(KIND=RP),INTENT(OUT),DIMENSION(1:N+1,1:N+1) :: Dval
! local variables
      INTEGER :: i
!
      ALLOCATE(x(1:N+1),w(1:N+1),xges(1:NQ*(N+1)),xmit(1:NQ+1))
      CALL LegendreGaussLobattoNodesandWeights(N,x,w)
      CALL DifferentiationsmatrixBerechnen(x,Dval,N+1)
      dx=1.0_RP/REAL(nq,kind=RP)
      CALL linspace(dx/2.0_RP,1.0_RP-dx/2.0_RP,NQ,xmit)
      DO i=1,NQ
        xges((i-1)*(N+1)+1:i*(N+1))=xmit(i)+dx/2.0_RP*x
      ENDDO ! i
    END SUBROUTINE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Operator aus Semidiskreterdarstellung
    FUNCTION R(u,N,NQ,D) result(solution)
      IMPLICIT NONE
      INTEGER      ,INTENT(IN)                                                 :: N,NQ
      REAL(KIND=RP),INTENT(IN),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: u
      REAL(KIND=RP),INTENT(IN),DIMENSION(1:N+1,1:N+1)                          :: D
      REAL(KIND=RP)           ,DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: solution
! local
      REAL(KIND=RP)           ,DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) ::L1,L2,L3
!
      CALL computeL(u,D,1,L1,N,NQ)
      CALL computeL(u,D,2,L2,N,NQ)
      CALL computeL(u,D,3,L3,N,NQ)
      solution=8.0_RP/(dx**3)*((-0.25_RP*dx**2)*L1-(0.25_RP*dx**2)*L2-(0.25_RP*dx**2)*L3)
    END FUNCTION
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE computeL(u,D,dir,result,N,NQ)
      IMPLICIT NONE
      INTEGER      ,INTENT(IN)                                                  :: N,NQ
      INTEGER      ,INTENT(IN)                                                  :: dir
      REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: u
      REAL(KIND=RP),INTENT(IN) ,DIMENSION(1:N+1,1:N+1)                          :: D
      !u dimensions:(nummer zelle x,zelle y,zelle z,x,y,z,variable)
      REAL(KIND=RP),INTENT(OUT),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: result
      !local Variables
      INTEGER                                  :: var,k,j,i,o,l,m
      REAL(KIND=RP),DIMENSION(1:N+1,1:5)       :: Fsharp
      REAL(KIND=RP),DIMENSION(1:N+1,1:N+1,1:5) :: FRand0,FRand1,uR,uL
      CHARACTER(len=2)                         :: whichflux='PI'
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
    END SUBROUTINE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    SUBROUTINE calculateEuler3DFlux(u,dir,result)
!    !SUBROUTINE gets the values of u and return the flux for the spcified direction
!    !dir=1,2,3 stands for x,y,z direction
!      IMPLICIT NONE
!      REAL(KIND=RP),dimension(:,:,:,:),INTENT(IN)     :: u
!      INTEGER                         ,INTENT(IN)     :: dir
!      REAL(KIND=RP),dimension(:,:,:,:),INTENT(OUT),allocatable    :: result
!      !local variables beyond here
!      INTEGER                                         :: n
!      REAL(KIND=RP),dimension(:,:,:),allocatable    ::p
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
    SUBROUTINE Initialcondition(xyz,u,NQ,N)
      IMPLICIT NONE
      INTEGER      ,INTENT(IN)                                                    :: NQ,N
      REAL(KIND=RP),INTENT(IN)   ,DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: xyz
      REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: u
      u(:,:,:,:,:,:,1)=2.0_RP+SIN(pi*(xyz(:,:,:,:,:,:,1)+xyz(:,:,:,:,:,:,2)+xyz(:,:,:,:,:,:,3)))/10.0_RP
      u(:,:,:,:,:,:,2)=1.0_RP
      u(:,:,:,:,:,:,3)=1.0_RP
      u(:,:,:,:,:,:,4)=1.0_RP
      u(:,:,:,:,:,:,5)=u(:,:,:,:,:,:,1)*u(:,:,:,:,:,:,1)
    END SUBROUTINE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      REAL(KIND=RP),dimension(1:N+1,1:N+1) :: p
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
    SUBROUTINE RungeKutta5explizit(ustar,nq,n,numvar,dt,Dval)
      IMPLICIT NONE
      !ustar=bekannte Werte
      INTEGER      ,INTENT(IN)                                                         :: numvar,n,nq
      REAL(KIND=RP),INTENT(IN)   ,DIMENSION(1:N+1,1:N+1)                               :: Dval
      REAL(KIND=RP),INTENT(INOUT),DIMENSION(1:nq,1:nq,1:nq,1:N+1,1:N+1,1:N+1,1:numvar) :: ustar
      !local
      REAL(KIND=RP),DIMENSION(1:nq,1:nq,1:nq,1:(N+1),1:N+1,1:N+1,1:numvar)    :: g
      INTEGER                                                                 :: step
      REAL(KIND=RP),DIMENSION(5)                                              :: a,b,c
      REAL(KIND=RP),INTENT(IN)                                                :: dt
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
      DO step=1,5
        g=a(step)*g+R(ustar,n,nq,Dval)
        ustar=ustar+c(step)*dt*g
      ENDDO ! step
    END SUBROUTINE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE
