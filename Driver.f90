!!
!  This origianlly uses allocatables but doesn't need too. You fix N and NQ in the beginning so you could directly
!  declare the arrays to be of the proper size. Alternatively, you could use an input file to read in the
!  polynomial order and number of elements in each direction such that you can change these values without
!  having to recompile the source code.
!!
PROGRAM Driver
  USE Zeitintegration
  IMPLICIT NONE
  INTEGER,PARAMETER                                             :: N=2,NQ=10
  REAL(KIND=RP)                                                 :: t=0.0_RP,tend=1.0_RP,CFL=0.1_RP,dt,a
  REAL(KIND=RP),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:5) :: u
  REAL(KIND=RP),DIMENSION(1:NQ,1:NQ,1:NQ,1:N+1,1:N+1,1:N+1,1:3) :: xyz
  REAL(KIND=RP),DIMENSION(1:N+1)                                :: xi
  REAL(KIND=RP),DIMENSION(1:NQ)                                 :: xl
  REAL(KIND=RP),DIMENSION(1:N+1,1:NQ)                           :: xin
  REAL(KIND=RP),DIMENSION(1:N+1,1:N+1)                          :: D
  REAL(KIND=RP),DIMENSION(1:5)                                  :: EOC,Fehler
  INTEGER                                                       :: l,m,o,k,i,j
!
  Fehler=0.0_RP
  EOC   =0.0_RP
  CALL Vorbereiten(N,NQ,D)
  CALL LegendreGaussLobattoNodesandWeights(N,xi,w)
  !! Bestimme GL punkte in jeder Zelle
  DO k=0,NQ-1
    xl(k+1)=(k+1.0_rp/2)*dx
    DO i=1,N+1
      xin(i,k+1)=xl(k+1)+dx/2*xi(i)
    ENDDO ! i
  ENDDO ! k
  !! Bestimme alle punkte.  
  DO o=1,NQ
    DO l=1,NQ
      DO m=1,NQ
        DO k=1,N+1
          DO j=1,N+1
            DO i=1,N+1
              xyz(m,l,o,i,j,k,1)=xin(i,m)
              xyz(m,l,o,i,j,k,2)=xin(j,l)
              xyz(m,l,o,i,j,k,3)=xin(k,o)
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k
      ENDDO ! m
    ENDDO ! l
  ENDDO ! o
  CALL Initialcondition(xyz,u,NQ,N)
  CALL lambdaMaxGlobal(u,a,NQ,N)
  dt=CFL/(3*a)*(dx/real(N+1,kind=RP))
!-ffpe-trap=denormal,invalid,zero,overflow,underflow
  DO WHILE (tend-t>epsilon(dt))
    PRINT*,'t'
    PRINT*,t
    PRINT*,'dt'
    PRINT*,dt
    PRINT*,'sum(energy)'
    PRINT*,sum(U(:,:,:,:,:,:,5))
    CALL lambdaMaxGlobal(u,a,NQ,N)
    dt=CFL/(3*a)*(dx/real(N+1))
    IF (t+dt>tend) dt=tend-t
    CALL RungeKutta5explizit(u,NQ,N,5,dt,D)
    !PRINT*,u(1,1,1,1,:,:,1)
    t=t+dt
  ENDDO ! tend
END PROGRAM Driver
