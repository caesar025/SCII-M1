!!
!  These routines appear to be working as you probably USEd them in the past semesters of the DG course. But
!  some of the loop structure is very weird and I have never seen them before in fortran. Often DO loops
!  are created without an index to run over. I DOn't know how this works, perhaps some KIND of implicit
!  loop structure. I wouldn't recommend this as relying on the compiler to "be smart" on how to manage memory
!  is very unreliable. Be very explicit on your loop bound and how you declare memory!! This will avoid weird
!  bugs and possible segmentation faults when you move to the parallel version of the code. Sometimes when
!  I compile the code it runs without producing NaN, but IF I change the spacing in the loops and try to
!  specify the loop bounds the code will NaN again. This points to memory not being handled properly
!  I expect IF you were to compile with an optimization flag, e..g -02, the code would break and NaN again
!  I would suggest to rewrite these routines to explicitly state what vairable the DO loops USE!
!!
MODULE Quadraturroutinen
  USE Diverses
contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE LegendreGaussNodesAndWeights(N,x,w)
    IMPLICIT NONE
    INTEGER,INTENT(IN)          :: N
    REAL(KIND=RP),INTENT(INOUT) :: x(0:N),w(0:N)
    INTEGER                     :: Boundary,j=0,nit=4,k=0
    REAL(KIND=RP)               :: LN,LN1
    REAL(KIND=RP)               :: delta,TOL=4*1e-16
    !
    Boundary=ceiling((real(N,KIND=RP)+1.0_RP)/2.0_RP)-1
    IF(N==0) THEN
      x(0)=0.0_RP
      w(0)=2.0_RP
    ELSEif(N==1) THEN
      x(0)=-sqrt((1.0_RP/3.0_RP))
      w(0)=1.0_RP
      x(1)=-x(0)
      w(1)=w(0)
    ELSE
      DO

        IF (j>Boundary) exit
        x(j)=-cos(((2*real(j,KIND=rp)+1.0_RP)/(2*real(N,KIND=RP)+2.0_RP))*pi)
        DO
          IF (k>nit) exit
          call LegendrePolynomialandDerivative(N+1,x(j),LN,LN1)
          delta=-LN/LN1
          x(j)=x(j)+delta
          IF (abs(delta)<=TOL*abs(x(j))) exit
          k=k+1
        ENDdo
        k=0;
        call LegendrePolynomialAndDerivative(N+1,x(j),LN,LN1)
        x(N-j)=-x(j)
        w(j)=2.0_RP/((1.0_RP-x(j)**2)*(LN1**2))
        w(N-j)=w(j)
        j=j+1
      ENDdo
      j=0
    ENDif
    IF (mod(N,2)==0) THEN
      call LegendrePolynomialAndDerivative(N+1,0.0_RP,LN,LN1)
      x(N/2)=0.0_RP
      w(N/2)=2.0_RP/((LN1)**2)
    ENDif
  END SUBROUTINE LegendreGaussNodesAndWeights
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  SUBROUTINE LegendrePolynomialAndDerivative(N,x,LN,LN1)
    IMPLICIT NONE
    INTEGER,INTENT(IN)          :: N
    REAL(KIND=RP),INTENT(IN)    :: x
    REAL(KIND=RP),INTENT(OUT)   :: LN,LN1
    REAL(KIND=RP),DIMENSION(3)  :: L,L1
    REAL(KIND=RP)               :: k
    IF (N==0) THEN
      LN=1
      LN1=0
    ELSEif (N==1) THEN
      LN=x
      LN1=1
    ELSE
      L(1)=1
      L(2)=x
      L1(1)=0
      L1(2)=1
      k=2;
      DO
        IF (k>N) exit
        L(3)=(2*k-1)/(k)*x*L(2)-(k-1)/(k)*L(1)
        L1(3)=L1(1)+(2*k-1)*L(2)
        L(1)=L(2)
        L(2)=L(3)
        L1(1)=L1(2)
        L1(2)=L1(3)
        k=k+1
      ENDdo
      LN=L(3)
      LN1=L1(3)
    ENDif
  END SUBROUTINE LegendrePolynomialAndDerivative
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  SUBROUTINE qAndLEvaluation(N,x,q,q1,LN)
    IMPLICIT NONE
    INTEGER,INTENT(IN)          :: N
    REAL(KIND=RP),INTENT(IN)    :: x
    REAL(KIND=RP),INTENT(OUT)   :: q,q1,LN
    REAL(KIND=RP),DIMENSION(3)  :: L,L1
    INTEGER                     :: k=2
    IF (N<2) return
    L(1)=1.0_RP
    L(2)=x
    L1(1)=0.0_RP
    L1(2)=1.0_RP
    DO k=1,N
      L(3)=(2.0_RP*real(k,rp)-1.0_RP)/real(k,rp)*x*L(2)-(real(k,rp)-1.0_RP)/real(k,rp)*L(1)
      L1(3)=L1(1)+(2.0_RP*real(k,rp)-1)*L(2)
      !if (k==N) exit
      L(1)=L(2)
      L(2)=L(3)
      L1(1)=L1(2)
      L1(2)=L1(3)
    ENDdo
    k=N+1
    LN=(2.0_RP*real(k,rp)-1.0_RP)/real(k,rp)*x*L(2)-(real(k,rp)-1.0_RP)/real(k,rp)*L(1)
    L1(3)=L1(1)+(2.0_RP*real(k,rp)-1.0_RP)*L(2)
    q=LN-L(1)
    q1=L1(3)-L1(1)
    k=2
  END SUBROUTINE
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  SUBROUTINE LegendreGaussLobattoNodesandWeights(N,x,w)
    IMPLICIT NONE
    INTEGER,INTENT(IN)                                      :: N
    REAL(KIND=RP),INTENT(INOUT)                             :: x(0:N),w(0:N)
    REAL(KIND=RP)                                           :: q,q1,LN
    INTEGER                                                 :: Boundary,k=0,nit=4,j=1
    REAL(KIND=RP)                                           :: delta,TOL=4*1e-16
    Boundary=ceiling((real(N,KIND=RP)+1)/2)-1
    IF (N==1) THEN
      x(0)=-1.0_RP
      w(0)=1.0_RP
      x(1)=1.0_RP
      w(1)=w(0)
    ELSEif (N==2) THEN
      x(0)=-1.0_RP
      x(1)=0.0_RP
      x(2)=1.0_RP
      w(0)=1.0_RP/3.0_RP
      w(2)=w(0)
      w(1)=4.0_RP/3.0_RP
      return
    ELSE
      x(0)=-1.0_RP
      w(0)=2.0_RP/(real(N,KIND=RP)*(real(N,KIND=RP)+1.0_RP))
      x(N)=1.0_RP
      w(N)=w(0)
      DO
        IF (j>Boundary) exit
        x(j)=-cos(((real(j,KIND=RP)+1.0_RP/4.0_RP)*pi)/real(N,KIND=RP)- &
          3.0_RP/(8.0_RP*real(N,KIND=RP)*pi)*1.0_RP/(real(j,KIND=RP)+1.0_RP/4.0_RP))
        DO
          IF (k>nit) exit
          call qAndLEvaluation(N,x(j),q,q1,LN)
          delta=-q/q1
          x(j)=x(j)+delta
          IF (abs(delta)<=TOL*abs(x(j))) exit
          k=k+1
        ENDdo
        call qAndLEvaluation(N-1,x(j),q,q1,LN)
        x(N-j)=-x(j)
        w(j)=2.0_RP/(real(N,KIND=RP)*(real(N,KIND=RP)+1.0_RP)*((LN)**2.0_RP))
        w(N-j)=w(j)
        k=0
        j=j+1
      ENDdo
      j=1
    ENDif
    IF (mod(N,2)==0.and.n/=2) THEN
      call qAndLEvaluation(N-1,0.0_RP,q,q1,LN)
      x(N/2)=0.0_RP
      w(N/2)=2.0_RP/(real(N,KIND=RP)*(real(N,KIND=RP)+1.0_RP)*((LN)**2.0_RP))
    ENDif
  END SUBROUTINE
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  SUBROUTINE BaryzGewichte(x,w)
    IMPLICIT NONE
    REAL(KIND=RP),DIMENSION(:),            INTENT(IN)   :: x
    REAL(KIND=RP),DIMENSION(:),allocatable,INTENT(OUT)  :: w
    INTEGER                                             :: N,j,start
    N=ubound(x,1)
    start=lbound(x,1)
    allocate(w(start:N))
    DO j=start,N
      w(j)=1.0_RP/product(x(j)-x,MASK=x/=x(j))
    ENDdo
  END SUBROUTINE
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  SUBROUTINE LagrangeIntBaryz(xst,res,x,b,f)
    IMPLICIT NONE
    REAL(KIND=RP),DIMENSION(:),            INTENT(IN)   :: xst,f
    REAL(KIND=RP)                         ,INTENT(IN)   :: x
    REAL(KIND=RP)                         ,INTENT(OUT)  :: res
    REAL(KIND=RP),DIMENSION(:),allocatable              :: w
    INTEGER                                             :: N,j,start
    REAL(KIND=RP),DIMENSION(:),allocatable,INTENT(OUT)  :: b
    !
    N=ubound(xst,1)
    start=lbound(xst,1)
    allocate(w(start:N))
    allocate(b(start:N))
    call BaryzGewichte(xst,w)
    DO j=start,N
      IF (x==xst(j)) THEN
        res=f(j)
        return
      ELSE
        b(j)=w(j)/(x-xst(j))
      ENDif
    ENDdo
    res=dot_product(f,b)/sum(b)
  END SUBROUTINE
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  SUBROUTINE DifferentiationsmatrixBerechnen(x,D,N)
    IMPLICIT NONE
    REAL(KIND=RP),DIMENSION(:),INTENT(IN)        :: x
    INTEGER,INTENT(IN)                           :: N
    REAL(KIND=RP),DIMENSION(1:n,1:n),INTENT(OUT) :: D
    INTEGER                                      :: i,j,start
    REAL(KIND=RP),DIMENSION(:),allocatable       :: w
    REAL(KIND=RP),DIMENSION(:,:),allocatable     :: b
    start=lBound(x,1)
    allocate(w(1:N))
    allocate(b(1:N,1:N))
    call BaryzGewichte(x,w)
    DO i=1,N
      DO j=1,N
        IF (i==j) cycle
        D(i,j)=w(j)/(w(i)*(x(i)-x(j)))
        b(i,j)=w(j)/(x(i)-x(j))
      ENDdo
      D(i,i)=-(1.0_RP/w(i))*(sum(b(i,1:i-1))+sum(b(i,i+1:N)))
    ENDdo
  END SUBROUTINE
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  SUBROUTINE LagrangePolynomAuswerten(x,j,l,xst)
    IMPLICIT NONE
    REAL(KIND=RP),DIMENSION(:),INTENT(IN)                   :: x,xst
    INTEGER,INTENT(IN)                                      :: j
    REAL(KIND=RP),DIMENSION(:),allocatable,INTENT(OUT)      :: l
    INTEGER                                                 :: N,start,i,m
    start=lbound(x,1)
    N=ubound(x,1)
    allocate(l(start-1:n-1))
    DO i=start,N
      l(i-1)=1.0_RP
      DO m=1,uBound(xst,1)
        IF (m==j+1) cycle
        l(i-1)=l(i-1)*((x(i)-xst(m))/(xst(j+1)-xst(m)))
      ENDdo
    ENDdo
  END SUBROUTINE
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    SUBROUTINE LagrangePolynomAuswertenEinzel(x,j,l,xst)
      IMPLICIT NONE
      REAL(KIND=RP),INTENT(IN)                                :: x
      REAL(KIND=RP),DIMENSION(:),INTENT(IN)                   :: xst
      INTEGER,INTENT(IN)                                      :: j
      REAL(KIND=RP),INTENT(OUT)                               :: l
      INTEGER                                                 :: m
      l=1.0_RP
      DO m=1,uBound(xst,1)
        IF (m==j+1) cycle
        l=l*((x-xst(m))/(xst(j+1)-xst(m)))
      ENDdo
    END SUBROUTINE
  END module
