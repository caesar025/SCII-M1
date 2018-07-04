program Driver_Manufactured
  use Zeitintegration
  implicit none
  REAL(KIND=RP)                                      :: t=0.0_rp,tend=0.1_RP,CFL=0.4_RP,dt,a
  INTEGER,parameter                                  :: n=3,anz=3
  REAL(KIND=RP),DIMENSION(:,:,:,:,:,:,:),allocatable :: u, usolution
  REAL(KIND=RP),DIMENSION(:,:,:,:),allocatable       :: uplot,xplot,yplot,zplot
  REAL(KIND=RP),DIMENSION(1:n+1,1:n+1)               :: D
  CHARACTER(len=2)                                   :: whichflux='PI',vis='AD' !whichflux: IF pirozzoli or standard fluxes; vis: viskos or just advective
  CHARACTER(Len=3)                                    ::numChar
  CHARACTER(LEN=17) :: fName  = "Movies/UXXX.tec"
  REAL(KIND=RP),DIMENSION(1:5,1:anz)                 :: errors,EOC
  INTEGER, DIMENSION(1:anz)                          :: nq
  INTEGER                                            :: k,i,start=1,m=0,l,o,li=1,j=0,xnum=1,ynum=1,znum=2
  nq=2**(/ (I,I=start,start+anz-1) /)
  call beginne_para() 
  DO k=1,anz
    allocate(u(1:Nq(k),1:nq(k),1:nq(k),1:n+1,1:n+1,1:n+1,1:5),usolution(1:Nq(k),1:nq(k),1:nq(k),1:n+1,1:n+1,1:n+1,1:5))
    allocate(uplot(1:nq(k)**3,1:n+1,1:n+1,1:n+1),xplot(1:nq(k)**3,1:n+1,1:n+1,1:n+1),yplot(1:nq(k)**3,1:n+1,1:n+1,1:n+1)&
      ,zplot(1:nq(k)**3,1:n+1,1:n+1,1:n+1))
    call Vorbereiten(n,nq(k),D)
    call Initialcondition(u,NQ(k),N)
    !-ffpe-trap=denormal,invalid,zero,overflow,underflow

    call Driver_zeit (u,nq(k),n,5,dt,D,t,tend,whichflux,vis,CFL,xnum,ynum,znum)
    ! Ueberpruefen ob Dichte/Druck negativ werden
    ! IF (ANY(u(:,:,:,:,:,:,1) < 0)) print*, 'Druck/Dichte sind negativ!'
    !print*,u(1,1,1,1,:,:,1)

    ! IF (MODULO(j,10).EQ.0 .and. nq(k)==4) THEN
    !   !!
    !   !  Print solution every 10 timesteps for movies
    !   !!
    !   DO li=1,nq(k)
    !     DO l=1,nq(k)
    !       DO o=1,nq(k)

    !         uplot(o+nq(k)*(l-1)+nq(k)**2*(li-1),:,:,:)=u(o,l,li,:,:,:,1)
    !         xplot(o+nq(k)*(l-1)+nq(k)**2*(li-1),:,:,:)=xyz(o,l,li,:,:,:,1)
    !         yplot(o+nq(k)*(l-1)+nq(k)**2*(li-1),:,:,:)=xyz(o,l,li,:,:,:,2)
    !         zplot(o+nq(k)*(l-1)+nq(k)**2*(li-1),:,:,:)=xyz(o,l,li,:,:,:,3)

    !       ENDdo
    !     ENDdo
    !   ENDdo


    !   m = m + 1
    !   WRITE(numChar,'(i3)')
    !   IF (m.GE.100) THEN
    !     fName(9:11) = numChar
    !   ELSEif(m.GE.10) then
    !     fName(9:9)    = "0"
    !     fName(10:11)  = numChar(2:3)
    !   ELSE
    !     fName(9:10) = "00"
    !     fName(11:11)=numChar(3:3)
    !   END IF
    !   open(unit=15,file=fName)
    !   call ExportToTecplot_3D(xplot,yplot,zplot,uplot,N,NQ(k)**3,15,'rho')
    !   close(15)

    ! ENDif
    ! j=j+1

    ! Berechne Fehler und Loesung
    call computeSolution(usolution,NQ(k),N,tend)
    call computeError(u,usolution,NQ(k),N,errors(:,k))
    if(id==0) then
      print*, 'FEHLER'
      print*,errors(1,:)
      print*,errors(2,:)
      print*,errors(3,:)
      print*,errors(4,:)
      print*,errors(5,:)
    end if
    ! Setzte alles wieder auf 0
    deallocate(u,usolution)
    t=0.0_RP
    deallocate(uplot,xplot,yplot,zplot)
  END DO

  call ende_para() 

  ! call computeEOC(errors,n,nq,anz,EOC)
  ! print*, "EOC"
  !print*, EOC(1,:)
  !print*, EOC(2,:)
  !print*, EOC(3,:)
  !print*, EOC(4,:)
  !print*, EOC(5,:)

  !       DO l=1,nq(1)
  !           DO o=1,nq(1)
  !               uplot(o+o*(l-1),:,:)=u(o,l,m,:,:,1,1)
  !               xplot(o+o*(l-1),:,:)=xyz(o,l,m,:,:,1,1)
  !               yplot(o+o*(l-1),:,:)=xyz(o,l,m,:,:,1,2)
  !           ENDdo
  !       ENDdo
  !   open(unit=15,file='rho.tec')
  !   call ExportToTecplot_2D(xplot,yplot,uplot,N,64,15,'rho')

  !   close(15)
  if(id==0) then
    call computeEOC(errors,n,nq,anz,EOC)
    print*, "EOC"
    print*, EOC(1,:)
    print*, EOC(2,:)
    print*, EOC(3,:)
    print*, EOC(4,:)
    print*, EOC(5,:)
  ENDIF
end program Driver_Manufactured
