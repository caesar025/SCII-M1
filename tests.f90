program tests
  use Zeitintegration
  implicit none

  !Testprogram for SCII-M1
  integer,parameter                                             :: n=3,nq=3,dir=3,pos=1
  REAL(KIND=RP),DIMENSION(1:n+1,1:5) ::u2,result
  REAL(KIND=RP),dimension(1:5)        ::u1
  u1=1.0_RP
  u1(2)=2.0_RP
  u1(5)=2.0_RP
  u2=1.0_RP
  u2(:,5)=2.0_RP
  call computeFsharp(u1,u2,dir,'ST',result,N)
  print*,result(1,:)


end program
