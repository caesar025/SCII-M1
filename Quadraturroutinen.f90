module Quadraturroutinen
use Diverses
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LegendreGaussNodesAndWeights(N,x,w)
    implicit none
        INTEGER,INTENT(IN)                                :: N
    REAL(KIND=RP),INTENT(INOUT)                           :: x(0:N),w(0:N)
    INTEGER                                               :: Boundary,j=0,nit=4,k=0
    REAL(KIND=RP)                                         :: LN,LN1
    REAL(KIND=RP)                                         :: delta,TOL=4*1e-16
    Boundary=ceiling((real(N,KIND=RP)+1.0_RP)/2.0_RP)-1
    if(N==0) then
    x(0)=0.0_RP
    w(0)=2.0_RP
    elseif(N==1) then
    x(0)=-sqrt((1.0_RP/3.0_RP))
    w(0)=1.0_RP
    x(1)=-x(0)
    w(1)=w(0)
    else
    do

        if (j>Boundary) exit
        x(j)=-cos(((2*real(j,kind=rp)+1.0_RP)/(2*real(N,KIND=RP)+2.0_RP))*pi)
        do
            if (k>nit) exit
            call LegendrePolynomialandDerivative(N+1,x(j),LN,LN1)
            delta=-LN/LN1
            x(j)=x(j)+delta
            if (abs(delta)<=TOL*abs(x(j))) exit
            k=k+1
        enddo
        k=0;
        call LegendrePolynomialAndDerivative(N+1,x(j),LN,LN1)
        x(N-j)=-x(j)
        w(j)=2.0_RP/((1.0_RP-x(j)**2)*(LN1**2))
        w(N-j)=w(j)
        j=j+1
    enddo
    j=0
    endif
    if (mod(N,2)==0) then
    call LegendrePolynomialAndDerivative(N+1,0.0_RP,LN,LN1)
    x(N/2)=0.0_RP
    w(N/2)=2.0_RP/((LN1)**2)
    endif
end subroutine LegendreGaussNodesAndWeights
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine LegendrePolynomialAndDerivative(N,x,LN,LN1)
    implicit none
    INTEGER,INTENT(IN)          :: N
    REAL(KIND=RP),INTENT(IN)    :: x
    REAL(KIND=RP),INTENT(OUT)   :: LN,LN1
    REAL(KIND=RP),DIMENSION(3)  :: L,L1
    REAL(KIND=RP)               :: k
    if (N==0) then
    LN=1
    LN1=0
    elseif (N==1) then
    LN=x
    LN1=1
    else
    L(1)=1
    L(2)=x
    L1(1)=0
    L1(2)=1
    k=2;
    do
    if (k>N) exit
    L(3)=(2*k-1)/(k)*x*L(2)-(k-1)/(k)*L(1)
    L1(3)=L1(1)+(2*k-1)*L(2)
    L(1)=L(2)
    L(2)=L(3)
    L1(1)=L1(2)
    L1(2)=L1(3)
    k=k+1
    enddo
    LN=L(3)
    LN1=L1(3)
    endif
end subroutine LegendrePolynomialAndDerivative
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine qAndLEvaluation(N,x,q,q1,LN)
    implicit none
    INTEGER,INTENT(IN)          :: N
    REAL(KIND=RP),INTENT(IN)    :: x
    REAL(KIND=RP),INTENT(OUT)   :: q,q1,LN
    REAL(KIND=RP),DIMENSION(3)  :: L,L1
    REAL(KIND=RP)               :: k=2.0_RP
    if (N<2) return
    L(1)=1.0_RP
    L(2)=x
    L1(1)=0.0_RP
    L1(2)=1.0_RP
    do
    if (k>N) exit
    L(3)=(2*k-1)/(k)*x*L(2)-(k-1)/(k)*L(1)
    L1(3)=L1(1)+(2*k-1)*L(2)
    !if (k==N) exit
    L(1)=L(2)
    L(2)=L(3)
    L1(1)=L1(2)
    L1(2)=L1(3)
    k=k+1.0_RP
    enddo
    k=real(N,KIND=RP)+1.0_RP
    LN=(2*k-1)/(k)*x*L(2)-(k-1)/(k)*L(1)
    L1(3)=L1(1)+(2*k-1)*L(2)
    q=LN-L(1)
    q1=L1(3)-L1(1)
    k=2.0_RP
end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine LegendreGaussLobattoNodesandWeights(N,x,w)
    implicit none
    INTEGER,INTENT(IN)                                      :: N
    REAL(KIND=RP),INTENT(INOUT)                             :: x(0:N),w(0:N)
    REAL(KIND=RP)                                           :: q,q1,LN
    INTEGER                                                 :: Boundary,k=0,nit=4,j=1
    REAL(KIND=RP)                                           :: delta,TOL=4*1e-16
    Boundary=ceiling((real(N,KIND=RP)+1)/2)-1
    if (N==1) then
        x(0)=-1.0_RP
        w(0)=1.0_RP
        x(1)=1.0_RP
        w(1)=w(0)
    elseif (N==2) then
        x(0)=-1.0_RP
        x(1)=0.0_RP
        x(2)=1.0_RP
        w(0)=1.0_RP/3.0_RP
        w(2)=w(0)
        w(1)=4.0_RP/3.0_RP
        return
    else
        x(0)=-1.0_RP
        w(0)=2.0_RP/(real(N,KIND=RP)*(real(N,KIND=RP)+1.0_RP))
        x(N)=1.0_RP
        w(N)=w(0)
        do
            if (j>Boundary) exit
            x(j)=-cos(((real(j,KIND=RP)+1.0_RP/4.0_RP)*pi)/real(N,KIND=RP)- &
            3.0_RP/(8.0_RP*real(N,KIND=RP)*pi)*1.0_RP/(real(j,KIND=RP)+1.0_RP/4.0_RP))
            do
                if (k>nit) exit
                call qAndLEvaluation(N,x(j),q,q1,LN)
                delta=-q/q1
                x(j)=x(j)+delta
                if (abs(delta)<=TOL*abs(x(j))) exit
                k=k+1
            enddo
            call qAndLEvaluation(N-1,x(j),q,q1,LN)
            x(N-j)=-x(j)
            w(j)=2.0_RP/(real(N,KIND=RP)*(real(N,KIND=RP)+1.0_RP)*((LN)**2.0_RP))
            w(N-j)=w(j)
            k=0
            j=j+1
        enddo
        j=1
    endif
    if (mod(N,2)==0.and.n/=2) then
        call qAndLEvaluation(N-1,0.0_RP,q,q1,LN)
        x(N/2)=0.0_RP
        w(N/2)=2.0_RP/(real(N,KIND=RP)*(real(N,KIND=RP)+1.0_RP)*((LN)**2.0_RP))
    endif
end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine BaryzGewichte(x,w)
    Implicit none
    REAL(KIND=RP),Dimension(:),            INTENT(IN)   :: x
    REAL(KIND=RP),Dimension(:),allocatable,INTENT(OUT)  :: w
    INTEGER                                             :: N,j,start
    N=ubound(x,1)
    start=lbound(x,1)
    allocate(w(start:N))
    do j=start,N
        w(j)=1.0_RP/product(x(j)-x,MASK=x/=x(j))
    enddo
end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine LagrangeIntBaryz(xst,res,x,b,f)
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
    do j=start,N
    if (x==xst(j)) then
    res=f(j)
    return
    else
    b(j)=w(j)/(x-xst(j))
    endif
    enddo
    res=dot_product(f,b)/sum(b)
end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine DifferentiationsmatrixBerechnen(x,D,N)
    Implicit none
    REAL(KIND=RP),Dimension(:),INTENT(IN)                   :: x
    INTEGER,intent(in)                                      :: N
    REAL(KIND=RP),DIMENSION(1:n,1:n),intent(out)        :: D
    INTEGER                                                 :: i,j,start
    REAL(KIND=RP),Dimension(:),allocatable                  :: w
    REAL(KIND=RP),DIMENSION(:,:),allocatable                :: b
    start=lBound(x,1)
    allocate(w(1:N))
    allocate(b(1:N,1:N))
    call BaryzGewichte(x,w)
    do i=1,N
        do j=1,N
            if (i==j) cycle
                D(i,j)=w(j)/(w(i)*(x(i)-x(j)))
                b(i,j)=w(j)/(x(i)-x(j))
        enddo
        D(i,i)=-(1.0_RP/w(i))*(sum(b(i,1:i-1))+sum(b(i,i+1:N)))
    enddo
end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine LagrangePolynomAuswerten(x,j,l,xst)
    Implicit none
    REAL(KIND=RP),DIMENSION(:),INTENT(IN)                   :: x,xst
    INTEGER,INTENT(IN)                                      :: j
    REAL(KIND=RP),DIMENSION(:),allocatable,INTENT(OUT)      :: l
    INTEGER                                                 :: N,start,i,m
    start=lbound(x,1)
    N=ubound(x,1)
    allocate(l(start-1:n-1))
    do i=start,N
    l(i-1)=1.0_RP
        do m=1,uBound(xst,1)
        if (m==j+1) cycle
            l(i-1)=l(i-1)*((x(i)-xst(m))/(xst(j+1)-xst(m)))
        enddo
    enddo
end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine LagrangePolynomAuswertenEinzel(x,j,l,xst)
    Implicit none
    REAL(KIND=RP),INTENT(IN)                                :: x
    REAL(KIND=RP),DIMENSion(:),intent(in)                   :: xst
    INTEGER,INTENT(IN)                                      :: j
    REAL(KIND=RP),INTENT(OUT)                               :: l
    INTEGER                                                 :: N,start,i,m
    l=1.0_RP
        do m=1,uBound(xst,1)
        if (m==j+1) cycle
            l=l*((x-xst(m))/(xst(j+1)-xst(m)))
        enddo
end subroutine
end module
