! Part of TsunAWI
! Copyright (c) 2014
! Alfred Wegener Institut,
! Helmholtz Centre for Polar and Marine Research
! Contact: rakowsky@awi.de

MODULE INITIAL_CONDITIONS
  !
  PRIVATE ! all functions except for the ones below are private
  PUBLIC :: start_point_init_file, start_point_init

  ! This module supplies routines which are related to the initial
  ! conditions.

CONTAINS

  SUBROUTINE start_point_init_file
    !
    ! !DESCRIPTION:
    ! This subroutine reads initialization information from file !{start\_point\_init\_file}.

    USE PARAMETERS
    USE MESH
    USE ELEMENTS
    USE SWE

    IMPLICIT NONE

    INTEGER                     :: nd, n, i, nmb
    INTEGER                     :: infile_inicon
    REAL(kind=wp)               :: max_ssh0

    OPEN (newunit=infile_inicon,file=trim(IniCondPath)//trim(init_filename), status='old')

    READ(infile_inicon,*) ssh0(1:nod2D)

    CLOSE(infile_inicon)
    
    max_ssh0=0.
    
!$OMP target teams & !!!distribute PARALLEL DEFAULT(SHARED) &
!$OMP&          PRIVATE(i,n,nmb) &
!$OMP&          REDUCTION(max:max_ssh0)    
    DO i=1,nmb_iter_smooth_inicond
!$OMP distribute parallel DO simd
       DO n=1,nod2D
          nmb=nghbr_nod2D(n)%nmb
          if (nmb > 0) then
             ssh_init(n)=SUM(ssh0(nghbr_nod2D(n)%addr(:)))/real(nmb,wp)
          else
             ssh_init(n) = ssh0(n)
          endif
       END DO
!$OMP END distribute parallel DO simd
!$OMP distribute parallel DO simd
       DO n=1,nod2D
          ssh0(n) = ssh_init(n)
       END DO
!$OMP END distribute parallel DO simd
    END DO
    
    ! Adjust topography
!$OMP distribute parallel DO simd
    DO n=1,nod2D
       max_ssh0 = max(max_ssh0,abs(ssh0(n)))
       nodhn(n) = nodhn(n)-ssh0(n)
    END DO
!$OMP END distribute parallel DO simd
!$OMP END target teams
    
    IF (max_ssh0 == 0.) then 
       print *, 'Initial Condition is all zero'
       STOP 1
    end IF


  END SUBROUTINE start_point_init_file

  !================================================================================

  SUBROUTINE start_point_init

    USE MESH
    USE ELEMENTS
    USE PARAMETERS
    USE SWE

    IMPLICIT NONE

    CHARACTER(LEN=12), parameter :: ideal_scaling_default= 'strasser2010'
    CHARACTER(LEN=6), parameter  :: ideal_shape_default  = 'cosine'
    REAL(kind=wp), parameter     :: mu_default = 3.0e11 ! Default value for rigidity

    INTEGER                            :: n, i,icase,ixy
    REAL(kind=wp)                       :: x, y, r, xp, yp

    REAL(kind=wp) :: x1, x2, y1, y2, d1, d2, dst1, dst2, amp, r_max
    REAL(kind=wp) :: dx, dy, rr, elup, eldn

    REAL(kind=wp) :: ideal_mw, ideal_epi_lon, ideal_epi_lat, ideal_rigidity
    CHARACTER(LEN=100) :: ideal_scaling, ideal_shape
    REAL(kind=wp) :: fault_area, mu, moment_dynecm, slip
    REAL(kind=wp) :: ssh0amp, ttl_volume, lc, dist

    INTEGER :: nmb, n_nonzero
    INTEGER :: infile_nml, i_nml_present
    INTEGER :: outfile

    NAMELIST /model_init_idealised/ ideal_mw, ideal_epi_lon, ideal_epi_lat, &
         ideal_shape, ideal_scaling, ideal_rigidity
     
     ideal_rigidity = -1.
     ideal_scaling  = ''
     ideal_shape    = ''
     ideal_epi_lon  = -9999.
     ideal_epi_lat  = -9999.
     ideal_mw       = -9999.
        
    
     open(newunit=infile_nml, file='namelist.tsunami')
     read(infile_nml, NML=model_init_idealised,iostat=i_nml_present)
     close(infile_nml)

        ! Check namelist input and set to defaults for missing values
     if (ideal_rigidity < 0.) then
        mu = mu_default
        print *,'/model_init_idealised/ Setting missing value for rigidity: mu=',mu_default
     else
        mu = ideal_rigidity
     endif
     if (len_trim(ideal_scaling) == 0) then
        print *,'/model_init_idealised/ Setting missing value ideal_scaling=',ideal_scaling_default
        ideal_scaling = ideal_scaling_default
        
     elseif (trim(ideal_scaling) /= ideal_scaling_default) then
        print *,'/model_init_idealised/ ideal_scaling = ',ideal_scaling
        print *,' is not implemented.'
        STOP 1
     endif
     
     if (len(trim(ideal_shape)) == 0) then
        print *,'/model_init_idealised/ Setting missing value ideal_shape=',ideal_shape_default
        ideal_shape = ideal_shape_default
     endif
     initial_shape = ideal_shape
       
     if (.not. cli_eq_params) then
        if (i_nml_present == 0  .and. ideal_epi_lon > -9999.  &
             .and. ideal_mw > 0..and. ideal_epi_lat > -9999.) then
           eq_epi_x = ideal_epi_lon
           eq_epi_y = ideal_epi_lat
           eq_mag   = ideal_mw
        else
           print *,  '/model_init_idealised/ EQ parameters mag, lon, lat not specified.'
           print *,  '/model_init_idealised/ Plese add them to either namelist or command line.'
           STOP 1
        endif
        
     else 
        print *,  'Earthquake parameters mag, lon, lat from command line will be used'
     endif
       
     
    if (LEN_TRIM(initial_shape)==6 .AND. trim(initial_shape)=='cosine') then
 
       if (coordinate_type /= 1) then
          print *,'Cosine bell only implemented for lon-lat-coordinates' 
          stop 1
       endif
 

       fault_area = 10._wp**(-3.476_wp + 0.952_wp*eq_mag) * 1.d6  ! in [m^2] 

       rr = sqrt(fault_area/pi) ! radius in [m]

       moment_dynecm = 10._wp**((eq_mag + 10.73_wp)*1.5_wp)
       slip = moment_dynecm / (mu * fault_area * 1.e4) *.01_wp  ! in [m]

       ttl_volume = fault_area*slip  ! in [m^3]

       ssh0amp= ttl_volume*pi /(4._wp*rr*rr *(pi-2._wp))

       lc = cos(eq_epi_y*rad)**2
       n_nonzero = 0

!$OMP target teams distribute PARALLEL DO  REDUCTION(+:n_nonzero) PRIVATE(n,dist)
       do n=1,nod2D
          dist = r_earth*r_earth *( lc*(coord_nod2D(1,n)-eq_epi_x*rad)**2 + &
                                       (coord_nod2D(2,n)-eq_epi_y*rad)**2 )
          if (dist <= rr*rr) then
             ssh0(n)   = ssh0amp * cos(sqrt(dist) * pi/(2._wp*rr))
             n_nonzero = n_nonzero+1
          else
             ssh0(n)=0._wp
          end if

       end do
!$OMP END target teams distribute PARALLEL DO

       print *,'====================================='
       print *,'Cosine shaped initial condition'
       print *,'Slip value: ',slip,' m'
       print *,'Amplitude: ',ssh0amp,' m'
       print *,'Radius',rr*0.001,' km'
       print *,'Resolved with',n_nonzero,' grid nodes'
       print *,'====================================='

    else

       print *,'/model_init_idealised/, ideal_shape = ', initial_shape
       print *,' is not implemented.'
       STOP 1
    endif

    IF (ALL(ABS(ssh0)==0._wp)) then
       print *, 'Initial Condition is all zero'
       STOP 1
    end IF
    ! Smoothing of the initial condition

!$OMP target teams PRIVATE(n,nmb) !!!PARALLEL DEFAULT(SHARED) PRIVATE(n,nmb)
    DO i=1,nmb_iter_smooth_inicond
!$OMP distribute parallel DO simd
       DO n=1,nod2D
          nmb=nghbr_nod2D(n)%nmb
          if (nmb > 0) then
             ssh_init(n) = SUM(ssh0(nghbr_nod2D(n)%addr(:)))/real(nmb,wp)
          else
             ssh_init(n) = ssh0(n)
          endif
       END DO
!$OMP END distribute parallel DO simd
!$OMP distribute parallel DO simd
       DO n=1,nod2D
          ssh0(n) = ssh_init(n)
       END DO
!$OMP END distribute parallel DO simd
    END DO

    ! Adjust topography
!$OMP distribute parallel DO simd
    DO n=1,nod2D
       nodhn(n)    = nodhn(n)-ssh0(n)
       ssh_init(n) = ssh0(n)
    END DO
!$OMP END distribute parallel DO simd
!$OMP END target teams !!! PARALLEL 



    ! OPEN(newunit=outfile, file='boundary_cond.out')
    !  DO i=1, nod2D
    !     IF (ssh0(i)==0.) cycle
    !     IF (coord_nod2D(1, i)/pi*180<75. .or. coord_nod2D(1, i)/pi*180>110.) cycle
    !     WRITE(outfile, *) coord_nod2D(1, i)/pi*180, coord_nod2D(2, i)/pi*180., ssh0(i)
    !  END DO
    ! CLOSE(outfile)
  END SUBROUTINE start_point_init

END MODULE INITIAL_CONDITIONS

