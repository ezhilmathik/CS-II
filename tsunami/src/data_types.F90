
! Part of TsunAWI
! Copyright (c) 2014
! Alfred Wegener Institut,
! Helmholtz Centre for Polar and Marine Research
! Contact: rakowsky@awi.de

! MODULE: DATA_TYPES -- contains the declaration of the specified data types used by other modules and subroutines

MODULE DATA_TYPES
    IMPLICIT NONE
    SAVE

  
    TYPE addresstype
        INTEGER                            :: nmb
        INTEGER, DIMENSION(:), allocatable :: addr
    END TYPE addresstype
  
END MODULE DATA_TYPES

