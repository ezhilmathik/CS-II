
! Part of TsunAWI
! Copyright (c) 2014
! Alfred Wegener Institut,
! Helmholtz Centre for Polar and Marine Research
! Contact: rakowsky@awi.de


! MODULE: ELEMENTS -- contains the local element definitions, basis functions

MODULE ELEMENTS

    USE PARAMETERS
    IMPLICIT NONE
    SAVE 
 
    REAL(KIND=dp)                                :: Vol2D

    REAL(KIND=dp), DIMENSION(2,3)                :: derivative_zeta_stdbf,    &
                                                   derivative_vel_stdbf

    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:)   :: bafu_xy
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:)     :: voltriangle
  
contains



SUBROUTINE standard_element_definition

    !
    !  DEFINITION OF :
    !    - 2D STANDARD ELEMENTS
    !    - BASISFUNCTIONS ON 2D STANDARD ELEMENTS (stdbafu)
    !     Zeta    : stdbafu(1)=1-x-y   stdbafu(2)=x       stdbafu(3)=y
    !     Velocity: stdbafu(1)=1-2y    stdbafu(2)=2x+2y-1 stdbafu(3)=1-2x   
    !    - SCALARPRODUCTS
    !    - DERIVATIVES OF STD.BASISFUNCTIONS
    !
  
    USE MESH
    IMPLICIT NONE
  
    INTEGER :: i, j
    Vol2D = 1./2.
    !
    derivative_zeta_stdbf      =  0.
    derivative_zeta_stdbf(:,1) = -1.
    derivative_zeta_stdbf(1,2) =  1. 
    derivative_zeta_stdbf(2,3) =  1.

END SUBROUTINE standard_element_definition


SUBROUTINE basisfunctions

    USE MESH
    IMPLICIT NONE
  
    REAL(KIND=dp)                           :: DET2D
    REAL(KIND=dp), DIMENSION(2, 3)          :: derivative_loczeta,      &
                                              derivative_locvel
    REAL(KIND=dp), DIMENSION(2,2)           :: jacobian2D, jacobian2D_inv
    INTEGER                                :: elem, i
  
    ALLOCATE(bafu_xy(6,elem2d))
    ALLOCATE(voltriangle(elem2d))
  
!$OMP PARALLEL
!$OMP DO 
    do elem=1,elem2D
       bafu_xy(:,elem)   = 0.
       voltriangle(elem) = 0.
    enddo
!$OMP END DO
!$OMP END PARALLEL

! bafu_xy = bafuzeta_xy
! Employ: bafuvel_xy(1,el) = -2.*bafuzeta_xy(3,el)
!         bafuvel_xy(2,el) = -2.*bafuzeta_xy(1,el)
!         bafuvel_xy(3,el) = -2.*bafuzeta_xy(2,el)

    DO elem = 1, elem2d
        CALL local_element_def(elem,        & 
            jacobian2D, jacobian2D_inv, DET2D,         &
		    derivative_loczeta)
        DO i = 1, 3
            bafu_xy(i, elem)   = derivative_loczeta(1, i)
            bafu_xy(i+3, elem) = derivative_loczeta(2, i)
        
        ENDDO
        voltriangle(elem) = abs(DET2D) * Vol2D
    ENDDO

  
END SUBROUTINE basisfunctions

SUBROUTINE local_element_def(element,    &
     jacobian2D, jacobian2D_inv, DET, derivative_loczeta)
    !--------------------------------------------------------------------------
    !                     2D LOCAL ELEMENT DEFINITIONS
    !     USES: 
    !       coord_nod2D(2,nod2D)
    !       elem2D_nodes(3, elem2D)
    !       derivative_stabafu_x_2D(2, 3)
    !
    !     OUTPUT:
    !       jacobian2D(2,2), jacobian2D_inv(2,2), DET
    !       derivative_locbafu_x_2D(2, 3)
    !--------------------------------------------------------------------------
    USE MESH

    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)                         :: element
    REAL(KIND=dp), DIMENSION(2,2),  INTENT(OUT)  :: jacobian2D
    REAL(KIND=dp), DIMENSION(2,2),  INTENT(OUT)  :: jacobian2D_inv
    REAL(KIND=dp), INTENT(OUT)                   :: DET
    REAL(KIND=dp), DIMENSION(2, 3), INTENT(OUT)  :: derivative_loczeta
    !
    REAL(KIND=dp), DIMENSION(2, 3)               :: local_cart
    REAL(KIND=dp), DIMENSION(3, 2)               :: der_transp
    INTEGER                                     :: i
    INTEGER                                     :: node
    REAL(kind=dp), parameter                    :: one_third = 1._dp/3._dp      


    DO i = 1, 3

       ! cartesian coordinates
        local_cart(1:2, i) = scaling * coord_nod2D(1:2, elem2D_nodes(i, element))

    END DO

    !  T R A N S F O R M A T I O N - M A T R I X   "jacobian"

    jacobian2D(:, 1) = local_cart(:, 2) - local_cart(:, 1)
    jacobian2D(:, 2) = local_cart(:, 3) - local_cart(:, 1)

    if ( coordinate_type == 1 ) then
       !correct for global grids for jumps across -180/180
       if (abs(jacobian2D(1, 1)) > .5*pi*r_earth) then
          jacobian2D(1, 1) = -sign(1._dp,jacobian2D(1, 1))*(2*pi*r_earth-abs(jacobian2D(1, 1)))
       end if
       if (abs(jacobian2D(1, 2)) > .5*pi*r_earth) then
          jacobian2D(1, 2) = -sign(1._dp,jacobian2D(1, 2))*(2*pi*r_earth-abs(jacobian2D(1, 2)))
       end if

       jacobian2D(1, :) = jacobian2D(1, :) * one_third*sum(cos_nod( elem2D_nodes(1:3,element)))
    END IF	  

    !  I N V E R S E   O F   jacobian

    CALL matrix_inverse_2x2(jacobian2D, jacobian2D_inv, DET)

    der_transp = MATMUL(TRANSPOSE(derivative_zeta_stdbf), jacobian2D_inv)
    derivative_loczeta = TRANSPOSE(der_transp)

END SUBROUTINE local_element_def



SUBROUTINE  matrix_inverse_2x2 (A, AINV, DET)
    !
    !     * * * * * * * * * * * * * * * * * * * * * * * * * *
    !     CALCULATE THE DETERMINATE AND INVERSE OF A(2,2)
    !     * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !        A    = ORIGINAL MATRIX
    !        AINV = INVERSE OF MATRIX A
    !        DET  = DETERMINANT OF A
    !
    IMPLICIT NONE
    !
    REAL(KIND=dp), DIMENSION(2,2), INTENT(IN)  :: A
    REAL(KIND=dp), DIMENSION(2,2), INTENT(OUT) :: AINV
    REAL(KIND=dp), INTENT(OUT)                 :: DET
    !
    INTEGER                                   :: i,j
    !
    DET  = A(1,1) * A(2,2) - A(1,2) * A(2,1)
    IF ( DET .eq. 0.0 )  THEN
        DO j = 1,2
            WRITE(*,*) (A(i,j), i = 1,2)
        END DO
        print *,'SINGULAR 2X2 MATRIX'
        STOP 1
    ELSE
        AINV(1,1) =  A(2,2) / DET
        AINV(1,2) = -A(1,2) / DET
        AINV(2,1) = -A(2,1) / DET
        AINV(2,2) =  A(1,1) / DET
    ENDIF
END SUBROUTINE matrix_inverse_2x2

  
END MODULE ELEMENTS
