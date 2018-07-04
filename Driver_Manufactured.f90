program Driver_Manufactured
    use Zeitintegration
    use mpi
    implicit none
    REAL(KIND=RP)                                      :: t=0.0_rp,tend=0.3_RP,CFL=0.1_RP,dt,a
    INTEGER,parameter                                  :: n=3,anz=2
    REAL(KIND=RP),DIMENSION(:,:,:,:,:,:,:),allocatable :: u, usolution
    REAL(KIND=RP),DIMENSION(:,:,:,:),allocatable       :: uplot,xplot,yplot,zplot
    REAL(KIND=RP),DIMENSION(1:n+1,1:n+1)               :: D
    CHARACTER(len=2)                                   :: whichflux='PI',vis='VI' !whichflux: if pirozzoli or standard fluxes; vis: viskos or just advective
    CHARACTER(Len=3)                                    ::numChar
    CHARACTER(LEN=17) :: fName  = "Movies/UXXX.tec"
    REAL(KIND=RP),DIMENSION(1:5,1:anz)                 :: errors,EOC
    INTEGER, DIMENSION(1:anz)                          :: nq
    INTEGER                                            :: k,i,start=6,m=0,l,o,li=1,j=0,ierralloc
    call MPI_Init(ierr) !starting MPI
    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr) !getting the number of processors
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierr)
    nq=2*(/ (I,I=start,start+anz-1) /)
    DO k=1,anz
        if (id==0) then
            allocate(u(1:Nq(k),1:nq(k),1:nq(k),1:n+1,1:n+1,1:n+1,1:5),STAT=ierralloc)
            allocate(usolution(1:Nq(k),1:nq(k),1:nq(k),1:n+1,1:n+1,1:n+1,1:5))
            allocate(uplot(1:nq(k)**3,1:n+1,1:n+1,1:n+1),xplot(1:nq(k)**3,1:n+1,1:n+1,1:n+1),yplot(1:nq(k)**3,1:n+1,1:n+1,1:n+1)&
                ,zplot(1:nq(k)**3,1:n+1,1:n+1,1:n+1))

            call Vorbereiten(n,nq(k),D)
            call Initialcondition(u,NQ(k),N)
            call lambdaMaxGlobal(u,a,NQ(k),N)
            dt=CFL/(3.0_RP*a)*(dx/real(2*N+1,KIND=RP))**2
        endif
        !-ffpe-trap=denormal,invalid,zero,overflow,underflow
        DO while(tend-t>epsilon(dt))
            if (id==0) then
                print*,'t'
                print*,t
                print*,'dt'
                print*,dt
                print*,'sum(energy)'
                print*,sum(U(:,:,:,:,:,:,5))
                call lambdaMaxGlobal(u,a,NQ(k),N)
                dt=CFL/(3.0_RP*a)*(dx/real(2*N+1))**2
                IF(t+dt>tend) dt=tend-t
            else
                if(.not.allocated(u)) allocate(u(1:nq(k),1:nq(k),1:nq(k),1:n+1,1:n+1,1:n+1,1:5))
            end if
						
            call MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	    			
            call RungeKutta5explizit(u,nq(k),n,5,dt,D,t,whichflux,vis)
            ! Ueberpruefen ob Dichte/Druck negativ werden
            IF (ANY(u(:,:,:,:,:,:,1) < 0)) print*, 'Druck/Dichte sind negativ!'
            !print*,u(1,1,1,1,:,:,1)
            t=t+dt

            IF (MODULO(j,10).EQ.0 .and. nq(k)==4.and.id==0) THEN
                !!
                !  Print solution every 10 timesteps for movies
                !!
                do li=1,nq(k)
                    do l=1,nq(k)
                        do o=1,nq(k)

                            uplot(o+nq(k)*(l-1)+nq(k)**2*(li-1),:,:,:)=u(o,l,li,:,:,:,1)
                            xplot(o+nq(k)*(l-1)+nq(k)**2*(li-1),:,:,:)=xyz(o,l,li,:,:,:,1)
                            yplot(o+nq(k)*(l-1)+nq(k)**2*(li-1),:,:,:)=xyz(o,l,li,:,:,:,2)
                            zplot(o+nq(k)*(l-1)+nq(k)**2*(li-1),:,:,:)=xyz(o,l,li,:,:,:,3)

                        enddo
                    enddo
                enddo


                m = m + 1
                WRITE(numChar,'(i3)')m
                IF (m.GE.100) THEN
                    fName(9:11) = numChar
                ELSEif(m.GE.10) then
                    fName(9:9)    = "0"
                    fName(10:11)  = numChar(2:3)
                ELSE
                    fName(9:10) = "00"
                    fName(11:11)=numChar(3:3)
                END IF
                open(unit=15,file=fName)
                call ExportToTecplot_3D(xplot,yplot,zplot,uplot,N,NQ(k)**3,15,'rho')
                close(15)

            endif
            j=j+1

        END DO
        ! Berechne Fehler und Loesung
        call computeSolution(usolution,NQ(k),N,t)
        call computeError(u,usolution,NQ(k),N,errors(:,k))
        if (id==0) then
            print*, 'FEHLER'
            print*,errors(1,:)
            print*,errors(2,:)
            print*,errors(3,:)
            print*,errors(4,:)
            print*,errors(5,:)
            ! Setzte alles wieder auf 0
            deallocate(usolution)
            deallocate(uplot,xplot,yplot,zplot)
        else
            deallocate (u)
        endif
        t=0.0_RP

    END DO


    ! call computeEOC(errors,n,nq,anz,EOC)
    ! print*, "EOC"
     !print*, EOC(1,:)
     !print*, EOC(2,:)
     !print*, EOC(3,:)
     !print*, EOC(4,:)
     !print*, EOC(5,:)

    !       do l=1,nq(1)
    !           do o=1,nq(1)
    !               uplot(o+o*(l-1),:,:)=u(o,l,m,:,:,1,1)
    !               xplot(o+o*(l-1),:,:)=xyz(o,l,m,:,:,1,1)
    !               yplot(o+o*(l-1),:,:)=xyz(o,l,m,:,:,1,2)
    !           enddo
    !       enddo
    !   open(unit=15,file='rho.tec')
    !   call ExportToTecplot_2D(xplot,yplot,uplot,N,64,15,'rho')

    !   close(15)
    call computeEOC(errors,n,nq,anz,EOC)
    print*, "EOC"
    print*, EOC(1,:)
    print*, EOC(2,:)
    print*, EOC(3,:)
    print*, EOC(4,:)
    print*, EOC(5,:)
    call MPI_Finalize(ierr)
end program Driver_Manufactured
