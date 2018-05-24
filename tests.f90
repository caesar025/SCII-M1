program tests
  use DGtoolbox
  implicit none

  !Testprogram for SCII-M1
  integer,parameter                                             :: n=2,nq=3,dir=3,pos=1
  REAL(KIND=RP), DIMENSION(1:n+1,1:n+1)                             :: D
  integer                                                         ::i
  real(KIND=RP), DIMENSION(1:n+1)                                 :: x,w,f
  call LegendreGaussLobattoNodesAndWeights(N,x,w)
  D=baryzdiffmatrix(x,n)
 f=x**2
  print*,matmul(D,f)-2*x
  end program
