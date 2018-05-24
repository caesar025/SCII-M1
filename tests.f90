program tests
  use Zeitintegration
  implicit none

  !Testprogram for SCII-M1
  integer,parameter                                             :: n=3,nq=3,dir=3,pos=1
  REAL(KIND=RP),DIMENSION(1:n+1,1:5) ::u2,result
  REAL(KIND=RP),dimension(1:5)        ::u1
  u1(1)=2.0_RP
  u1(2)=3.0_RP
  u1(3)=4.0_RP
  u1(4)=5.0_RP
  u1(5)=15.0_RP
  u2(:,1)=2.0_RP
  u2(:,2)=3.0_RP
  u2(:,3)=4.0_RP
  u2(:,4)=5.0_RP
  u2(:,5)=15.0_RP
  call computeFsharp(u1,u2,dir,'PI',result,N)
  print*,result(1,:)


end program
