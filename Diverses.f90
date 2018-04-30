module Diverses
    implicit none
    INTEGER,PARAMETER                                     :: RP = SELECTED_REAL_KIND(15)
    REAL(KIND=RP),PARAMETER                                     :: pi = 4.0_RP*ATAN(1.0_RP)
contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE ExportToTecplot_1D(x,u,N,fUnit,name)
        ! File writer for one dimensional solution vector
        IMPLICIT NONE
        INTEGER ,INTENT(IN) :: N,fUnit
        REAL(KIND=RP),DIMENSION(0:N),INTENT(IN) :: x,u
        CHARACTER(len=*),INTENT(IN) ::name
        ! Local variables
        INTEGER :: j
        WRITE(FUnit,*)name
        DO j = 0,N
            WRITE(fUnit,*)x(j),u(j)
        END DO
    END SUBROUTINE ExportToTecplot_1D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE ExportToTecplot_2D(x,y,u,N,K,fUnit,solutionFile)
        IMPLICIT NONE
        INTEGER ,INTENT(IN) :: N,K,fUnit
        CHARACTER(LEN=*) ,INTENT(IN) :: solutionFile
        REAL(KIND=RP),DIMENSION(K,0:N,0:N),INTENT(IN) :: x,y,u
        ! Local variables
        INTEGER :: i,j,l
        !
        WRITE(fUnit,*)'TITLE = "',solutionFile,'solution.tec"'
        WRITE(fUnit,*)'VARIABLES = "x","y","',solutionFile,'"'
        !
        DO l = 1,K
            WRITE(fUnit,*)"ZONE I =",N+1,",J=",N+1,",F=POINT"
            DO j = 0,N
                DO i = 0,N
                    WRITE(fUnit,*)x(l,i,j),y(l,i,j),u(l,i,j)
                END DO
            END DO
        END DO
        !
        RETURN
    END SUBROUTINE ExportToTecplot_2D
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE MatrixPlot(M)
        Implicit none
        REAL(KIND=RP),Dimension(:,:),INTENT(IN)         ::m
        INTEGER                                         ::i,N,start
        start=lbound(m,1)
        N=uBOund(m,1)
        do i=1,N
            print*,m(i,:)
        enddo
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine linspace(d1,d2,n,y)
        Implicit none
        INTEGER,INTENT(IN)                          ::n
        REAL(KIND=RP),INTENT(IN)                    ::d1,d2
        REAL(KIND=RP),DIMENSION(1:N),INTENT(OUT)    ::y
        INTEGER                                     ::i
        REAL(KIND=RP)                               ::c
        c=d2-d1
        do i=1,n
            y(i)=d1+(i-1)*c/real(n-1,KIND=RP)
        enddo
    end subroutine
end module Diverses
