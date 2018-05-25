MODULE Diverses
    IMPLICIT NONE
    INTEGER      ,PARAMETER :: RP = SELECTED_REAL_KIND(15)
    REAL(KIND=RP),PARAMETER :: pi = 4.0_RP*ATAN(1.0_RP)
CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE ExportToTecplot_1D(x,u,N,fUnit,name)
        ! File writer for one dimensional solution vector
        IMPLICIT NONE
        INTEGER         ,INTENT(IN)                :: N,fUnit
        CHARACTER(len=*),INTENT(IN)                :: name
        REAL(KIND=RP)   ,INTENT(IN),DIMENSION(0:N) :: x,u
        ! Local variables
        INTEGER :: j
        !
        WRITE(FUnit,*)name
        DO j = 0,N
            WRITE(fUnit,*)x(j),u(j)
        END DO ! j
    END SUBROUTINE ExportToTecplot_1D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE ExportToTecplot_2D(x,y,u,N,K,fUnit,solutionFile)
        IMPLICIT NONE
        INTEGER         ,INTENT(IN)                      :: N,K,fUnit
        CHARACTER(LEN=*),INTENT(IN)                      :: solutionFile
        REAL(KIND=RP)   ,INTENT(IN),DIMENSION(1:K,1:N+1,1:N+1) :: x,y,u
        ! array size issue? your computation uses 1:N+1 but this uses 0:N structure, might be a problem
        ! also you'll need a 3D version of this plotting routine, available in the FORTRAN notes
        ! Local variables
        INTEGER :: i,j,l
        !
        WRITE(fUnit,*)'TITLE = "',solutionFile,'solution.tec"'
        WRITE(fUnit,*)'VARIABLES = "x","y","',solutionFile,'"'
        !
        DO l = 1,K
            WRITE(fUnit,*)"ZONE I =",N+1,",J=",N+1,",F=POINT"
            DO j = 1,N+1
                DO i = 1,N+1
                    WRITE(fUnit,*)x(l,i,j),y(l,i,j),u(l,i,j)
                END DO ! i
            END DO ! j
        END DO ! l
        !
        RETURN
    END SUBROUTINE ExportToTecplot_2D
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE MatrixPlot(M,N,l)
        IMPLICIT NONE
        INTEGER      ,INTENT(IN)                    :: N,l
        REAL(KIND=RP),INTENT(IN),DIMENSION(1:n,1:l) :: m
        ! local variables
        INTEGER :: i
        !
        DO i=1,N
            PRINT*,m(i,:)
        ENDDO ! i
    END SUBROUTINE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE linspace(d1,d2,n,y)
        IMPLICIT NONE
        INTEGER      ,INTENT(IN)                 :: n
        REAL(KIND=RP),INTENT(IN)                 :: d1,d2
        REAL(KIND=RP),INTENT(OUT),DIMENSION(1:N) :: y
        ! local variables
        INTEGER                                  :: i
        REAL(KIND=RP)                            :: c
        !
        c=d2-d1
        DO i=1,n
            y(i)=d1+(i-1)*c/real(n-1,KIND=RP)
        ENDDO ! i
    END SUBROUTINE
    SUBROUTINE ExportToTecplot_3D(x,y,z,u,N,K,fUnit,solutionFile)
        IMPLICIT NONE
        INTEGER                                 ,INTENT(IN) :: N,K,fUnit
        CHARACTER(LEN=*)                        ,INTENT(IN) :: solutionFile
        REAL(KIND=RP),DIMENSION(1:K,0:N,0:N,0:N),INTENT(IN) :: x,y,z,u
        !  Local variables
        INTEGER           :: i,j,h,p
        CHARACTER(LEN=32) :: valuesFMT
        !  Format the output of real numbers
        valuesFMT = "(4E13.5)"
        !
        WRITE(fUnit,*)'TITLE = "',solutionFile,'solution.tec"'
        WRITE(fUnit,*)'VARIABLES = "x","y","z","',solutionFile,'"'
        !
        DO h = 1,K
          WRITE(fUnit,*)"ZONE I =",N+1,",J=",N+1,",K=",N+1,",F=POINT"
          DO i = 0,N
            Do j = 0,N
              DO p = 0,N
                WRITE(fUnit,valuesFMT)x(h,i,j,p),y(h,i,j,p),z(h,i,j,p),u(h,i,j,p)
              END DO
              END DO
            END DO
        END DO
        !
        RETURN
    END SUBROUTINE ExportToTecplot_3D
END MODULE Diverses
