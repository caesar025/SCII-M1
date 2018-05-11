MODULE DGtoolbox
  USE Diverses
  IMPLICIT NONE
contains


  SUBROUTINE LegendrePolynomialAndDerivative(N,x,l,lstrich)
    REAL(KIND=RP),INTENT(IN)  :: x
    INTEGER,INTENT(IN)        :: N
    REAL(KIND=RP),INTENT(OUT) :: l,lstrich
    REAL(KIND=RP)             :: lnmin1,lstrichnmin1,lnmin2,lstrichnmin2
    INTEGER                   :: k
    IF(N==0) THEN

      l=1.0_RP
      lstrich=0.0_RP
    ELSE IF(N==1) THEN

      l=x
      lstrich=1.0_RP
    ELSE
      lnmin2       = 1.0_RP
      lnmin1       = x
      lstrichnmin2 = 0.0_RP
      lstrichnmin1 = 1.0_RP


      DO k=2,N
        l            = (2*REAL(k,rp)-1)/REAL(k,rp)*x*lnmin1-(REAL(k,rp)-1)/REAL(k,rp)*lnmin2
        lstrich      = lstrichnmin2+(2*REAL(k,rp)-1)*lnmin1
        lnmin2       = lnmin1
        lnmin1       = l
        lstrichnmin2 = lstrichnmin1
        lstrichnmin1 = lstrich
      END DO
    END IF
  END SUBROUTINE LegendrePolynomialAndDerivative

  !!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE LegendreGaussNodesAndWeights(N,x,w)
    INTEGER,INTENT(IN)                       :: N
    REAL(KIND=RP),INTENT(OUT),DIMENSION(0:N) :: x,w
    REAL(KIND=RP)                            :: tol=4*1e-16,l,lstrich,delta
    INTEGER                                  :: k,j,nit=4
    IF(N==0) THEN
      x(0)=0
      w(0)=2
    ELSE IF(N==1) THEN
      x(0)=-sqrt(1.0_rp/3.0_rp)
      w(0)=1
      x(1)=-x(0)
      w(1)=w(0)
    ELSE
      DO j=0,floor((REAL(N,rp)+1.0_rp)/2.0_rp)-1
        x(j)=-cos((2.0_rp*REAL(j,rp)+1.0_rp)/(2.0_rp*REAL(N,rp)+2.0_rp)*pi)

        DO k=0,nit
          call LegendrePolynomialAndDerivative(N+1,x(j),l,lstrich)

          delta=-l/lstrich
          x(j)=x(j)+delta

          IF(abs(delta)<=tol*abs(x(j))) exit
        END DO

        call LegendrePolynomialAndDerivative(N+1,x(j),l,lstrich)
        x(n-j)=-x(j)
        w(j)=2/((1-x(j)*x(j))*(lstrich*lstrich) )
        w(n-j)=w(j)
      END DO
    END IF
    IF(mod(N,2)==0) THEN
      call LegendrePolynomialAndDerivative(N+1,0.0_RP,l,lstrich)
      x(n/2)=0
      w(n/2)=2/(lstrich*lstrich)
    END IF


  END SUBROUTINE LegendreGaussNodesAndWeights
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE qAndLEvaluation(n,x,q,qstrich,l)
    INTEGER,INTENT(IN)        :: N
    REAL(KIND=RP),INTENT(OUT) :: q,qstrich,l
    REAL(KIND=RP),INTENT(IN)  :: x
    REAL(KIND=RP)             :: lnmin1,lstrichnmin1,lnmin2,lstrichnmin2,ln1,lstrich1,lstrich
    INTEGER                   :: k
    lstrichnmin1=1.0_RP
    lstrichnmin2=0.0_RP
    lnmin1=x
    lnmin2=1.0_RP
    DO k=2,n
      l=(2*REAL(k,rp)-1)/REAL(k,rp)*x*lnmin1-(REAL(k,rp)-1)/REAL(k,rp)*lnmin2
      lstrich=lstrichnmin2+(2*REAL(k,rp)-1)*lnmin1
      lnmin2=lnmin1
      lnmin1=l
      lstrichnmin2=lstrichnmin1
      lstrichnmin1=lstrich
    END DO
    k=n+1
    ln1=(2*REAL(k,rp)-1)/REAL(k,rp)*x*lnmin1-(REAL(k,rp)-1)/REAL(k,rp)*lnmin2
    lstrich1=lstrichnmin2+(2*REAL(k,rp)-1)*lnmin1
    q=ln1-lnmin2
    qstrich=lstrich1-lstrichnmin2
  END SUBROUTINE qAndLEvaluation

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE LegendreGaussLobattoNodesAndWeights(N,x,w)
    INTEGER,INTENT(IN)                       :: N
    REAL(KIND=RP),INTENT(OUT),DIMENSION(0:N) :: x,w
    REAL(KIND=RP)                            :: l,lstrich,delta,q,qstrich,tol=4*1e-16
    INTEGER                                  :: k,j,nit=4
    IF(N==1) THEN
      x(0)=-1
      w(0)=1
      x(1)=1
      w(1)=w(0)
    ELSE
      x(0)=-1
      w(0)=2/(n*(n+1))
      x(n)=1
      w(n)=w(0)
      DO j=0,floor((REAL(N,rp)+1.0_rp)/2.0_rp)-1
        x(j)=-cos(((REAL(j,rp)+1.0_rp/4_rp)*pi/N-3/(8.0_rp*n*pi*(REAL(j,rp)+1.0_rp/4_rp))))
        DO k=0,nit
          call qANdLEvaluation(n,x(j),q,qstrich,l)
          delta=-q/qstrich
          x(j)=x(j)+delta
          IF(abs(delta)<=tol*abs(x(j))) exit
        END DO
        call qANdLEvaluation(n,x(j),q,qstrich,l)
        x(n-j)=-x(j)
        w(j)=2.0_rp/((REAL(n,rp)*(REAL(n,rp)+1)*l*l))
        w(n-j)=w(j)
      END DO
    END IF
    IF(mod(N,2)==0) THEN
      call qANdLEvaluation(n,0.0_rp,q,qstrich,l)
      x(n/2)=0
      w(n/2)=2/((REAL(n,rp)*(REAL(n,rp)+1)*l*l))
    END IF


  END SUBROUTINE LegendreGaussLobattoNodesAndWeights
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION f1 (x,n,i)
    INTEGER,INTENT(IN):: N,i
    REAL(KIND=RP),INTENT(IN) :: x
    REAL(KIND=RP) :: f1
    IF(i==1) THEN
      f1=cos(x)
    ELSE IF(i==2) THEN
      f1=1/(1+x*x)

    ELSE IF(i==3) THEN
      f1=x**(2*n-2)
    ELSE IF(i==4) THEN
      f1=x**(2*n)
    ELSE IF(i==5) THEN
      f1=x**(2*n+2)
    END IF

  END FUNCTION f1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION exakt1 (n,i)
    INTEGER,INTENT(IN) :: N,i
    REAL(KIND=RP)      :: exakt1


    IF(i==1) THEN
      exakt1=sin(1.0_rp)-sin(-1.0_rp)
    ELSE IF(i==2) THEN
      exakt1=atan(1.0_rp)-atan(-1.0_rp)
    ELSE IF(i==3) THEN
      exakt1=1/(2*REAL(n,rp)-1)*(1.0_rp)**(2*REAL(n,rp)-1)-1/(2*REAL(n,rp)-1)*(-1.0_rp)**(2*REAL(n,rp)-1)
    ELSE IF(i==4) THEN
      exakt1=1/(2*REAL(n,rp)+1)*(1.0_rp)**(2*REAL(n,rp)+1)-1/(2*REAL(n,rp)+1)*(-1.0_rp)**(2*REAL(n,rp)+1)
    ELSE IF(i==5) THEN
      exakt1=1/(2*REAL(n,rp)+3)*(1.0_rp)**(2*REAL(n,rp)+3)-1/(2*REAL(n,rp)+3)*(-1.0_rp)**(2*REAL(n,rp)+3)
    END IF
  END FUNCTION exakt1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION baryzweight(x,n)

    IMPLICIT NONE
    INTEGER                                 :: i,j
    INTEGER, INTENT(IN)                     :: n
    REAL(KIND=RP),DIMENSION(0:n),INTENT(IN) :: x
    REAL(KIND=RP),DIMENSION(0:n)            :: baryzweight
    REAL(KIND=RP)                           :: pro

    !
    DO i=0,n
      pro = 1.0_RP
      DO j=0,n
        IF (i==j) CYCLE
        pro=pro*(x(i)-x(j))
      END DO
      baryzweight(i)=1/pro

    END DO


  END FUNCTION baryzweight

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION baryzvalue(f,x,n,point)
    IMPLICIT NONE
    INTEGER                      :: i
    INTEGER, INTENT(IN)          :: n
    REAL(KIND=RP),INTENT(IN)     :: point
    REAL(KIND=RP),DIMENSION(0:n) :: f,x,baryzweights,summevec
    REAL(KIND=RP)                :: baryzvalue
    DO i=0,n
      IF(point==x(i)) THEN
        baryzvalue=f(i)
        return
      END IF
    END DO

    baryzweights=baryzweight(x,n)
    DO i=0,n
      summevec(i)=baryzweights(i)/(point-x(i))
    END DO
    baryzvalue= sum(f*summevec)/sum(summevec)
  END FUNCTION baryzvalue
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION baryzdiffmatrix(x,n)
    IMPLICIT NONE
    INTEGER                                  :: i,j
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=RP),DIMENSION(0:n), INTENT(IN) :: x
    REAL(KIND=RP),DIMENSION(0:n)             :: summevec,w
    REAL(KIND=RP),DIMENSION(0:n,0:n)         :: baryzdiffmatrix

    w=baryzweight(x,n)
    DO i=0,n
      summevec(i)=0
      DO j=0,n
        IF(i==j) cycle
        summevec(i)=summevec(i)+w(j)/(x(i)-x(j))
      END DO
    END DO
    summevec=summevec*(-1/w)


    DO i=0,n
      DO j=0,n
        IF(i==j) THEN
          baryzdiffmatrix(i,j)=summevec(i)
        ELSE
          baryzdiffmatrix(i,j)=w(j)/(w(i)*(x(i)-x(j)))
        END IF

      END DO
    END DO

  END FUNCTION baryzdiffmatrix
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION lagrangevalue(x,j,n,point)
    IMPLICIT NONE
    INTEGER                                  :: i
    INTEGER, INTENT(IN)                      :: n,j
    REAL(KIND=RP),INTENT(IN)                 :: point
    REAL(KIND=RP),DIMENSION(0:n), INTENT(IN) :: x
    REAL(KIND=RP)                            :: lagrangevalue
    lagrangevalue=1.0_rp

    DO i=0,n
      IF(i==j) cycle
      lagrangevalue=lagrangevalue*(point-x(i))/(x(j)-x(i))
    END DO



  END FUNCTION lagrangevalue
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE plotdatasave(x,f,n,fileid)
    IMPLICIT NONE
    INTEGER,INTENT(IN)                      :: n,fileid
    REAL(KIND=RP),DIMENSION(1:n),INTENT(IN) :: x,f
    INTEGER                                 :: j
    DO j=1,n
      write(fileid,*) x(j),f(j)
    END DO
  END SUBROUTINE plotdatasave
END MODULE DGtoolbox
