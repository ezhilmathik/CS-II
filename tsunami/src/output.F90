! Part of TsunAWI
! Copyright (c) 2014
! Alfred Wegener Institut,
! Helmholtz Centre for Polar and Marine Research
! Contact: rakowsky@awi.de


!---------------------------------------------------------------------
SUBROUTINE ascii_out_init

    USE PARAMETERS
    USE MESH
    USE ELEMENTS
    USE SWE
  
!---- local declarations

    IMPLICIT NONE
    CHARACTER(len=100)    :: ssh_file, vel_file, ttt_file
    INTEGER               :: i

    ssh_file   = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_ssh_000000.out'
    vel_file   = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_vel_000000.out'

!---- open file & write sea surface height  
  
    OPEN(41, file = ssh_file)
    DO i = 1, nod2D
        WRITE(41,*) real(ssh2(i),sp)
    END DO
    CLOSE(41)

!----- open file & write velocity

    if (write_velocity) then
       OPEN(41, file = vel_file)
       DO i = 1, nod2D
          WRITE(41,*)  real(uv_node(1,i),sp), real(uv_node(2,i),sp)
       END DO
       CLOSE(41)
    endif

    if (write_ttt) then
       ttt_file   = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_ttt.out'
       OPEN(41, file = ttt_file)
       DO i = 1, nod2D
          WRITE(41, *) real(ttt(i),sp)
       END DO
       CLOSE(41)
    endif

  END SUBROUTINE ascii_out_init
!=====================================
SUBROUTINE ascii_out

    USE PARAMETERS
    USE MESH
    USE ELEMENTS
    USE SWE
  
!---- local declarations

    IMPLICIT NONE
    CHARACTER(len=100)    :: ssh_file, vel_file
    INTEGER               :: i
    INTEGER, save         :: n_out=0
    CHARACTER(len=6)      :: c_out

    n_out = n_out+1
    write(c_out,'(i6.6)') n_out

    ssh_file   = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_ssh_'//c_out//'.out'
    vel_file   = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_vel_'//c_out//'.out'

!---- open file & write sea surface height  
  
    OPEN(41, file = ssh_file)
    DO i = 1, nod2D
        WRITE(41, *) real(ssh2(i),sp)
    END DO
    CLOSE(41)

!----- open file & write velocity

    if (write_velocity) then
       OPEN(41, file = vel_file)
       DO i = 1, nod2D
          WRITE(41,*)  real(uv_node(1,i),sp), real(uv_node(2,i),sp)
       END DO
       CLOSE(41)
    endif

  END SUBROUTINE ascii_out

!================================================================

SUBROUTINE ascii_out_finalize    

  USE PARAMETERS
  USE MESH
  USE ELEMENTS
  USE SWE
  
!---- local declarations

    IMPLICIT NONE
    CHARACTER(len=100)    :: mwh_file, arr_file
    INTEGER               :: i

    arr_file   = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_arrivaltime.out'
    mwh_file   = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_mwh.out'

!---- open file & write arrival time

    OPEN(41, file = arr_file)
    DO i = 1, nod2D
        WRITE(41, '(2e11.4)') arrival_time(i)
    END DO
    CLOSE(41)

!---- open file & write maximum wave height

    OPEN(41, file = mwh_file)
    DO i = 1, nod2D
        WRITE(41, '(e10.3)') mwh(i)
    END DO
    CLOSE(41)


  END SUBROUTINE ascii_out_finalize



#ifndef NO_NETCDF
SUBROUTINE nc_init

!----- uses

    USE PARAMETERS
    USE MESH
    USE ELEMENTS
    USE SWE
    USE OKADA_FAULT_PARAMETERS

!-----local declarations

    IMPLICIT NONE

    include 'netcdf.inc'

    CHARACTER(len=100)      :: ncfile
    CHARACTER(len=200)      :: attstr
    INTEGER                 :: attValInt
    REAL(KIND=dp)           :: attVal, notyetknown
    CHARACTER(len=2)        :: par_num
    CHARACTER(len=3)        :: StrMeshID
    INTEGER                 :: i
    INTEGER                 :: fileid 			     ! IDs netCDF file
    INTEGER                 :: DimId_iter		     ! dimension: iteration
    INTEGER                 :: DimId_n2D		     ! dimension: nodes
    INTEGER                 :: DimId_el2D		     ! dimension: elements
    INTEGER                 :: Dim1, Dim2, Dim3, Dim4	     ! dimensions: 
    INTEGER                 :: VarId_time, VarId_iter	     ! variables: time, iteration
    INTEGER                 :: VarId_epot, VarId_ekin        ! variables: kinetic, potential energy
    INTEGER                 :: VarId_trian		     ! variable: triangle
    INTEGER                 :: VarId_loc2D, VarId_lon, VarId_lat ! variables: 2D location, longitude, latitude
    INTEGER                 :: VarId_ssh		     ! variable: sea surface height
    INTEGER                 :: VarId_tri_area		     ! variable: area
    INTEGER                 :: VarId_topo, VarId_index	     ! variables: topography, index
    INTEGER                 :: VarId_topo_ini                ! variables: initial topography
    INTEGER                 :: VarId_arrival                 ! variable: arrival time
    INTEGER                 :: VarId_ttt                     ! variable: estimated arrival time ttt
    INTEGER                 :: VarId_mwh  		     ! variable: maximum wave height
    INTEGER                 :: VarId_flux, VarId_vel         ! variables: flux and mean absolute velocity
    INTEGER                 :: VarId_uplift                  ! variable: topo/bathy uplift
    INTEGER                 :: VarId_nodal_vel               ! variable: velocity in nodes

    INTEGER                 :: VarId_nodal_cfl               ! variable: vel / edgelength 
    INTEGER                 :: s       			     ! auxiliari: status counter
    INTEGER, DIMENSION(3)   :: dimarray			     ! auxiliari: array dimension
    INTEGER, DIMENSION(500) :: stat    			     ! auxiliari: status array
        
    REAL, DIMENSION(2,nod2D):: AUX_loc
    REAL(KIND=dp)           ::   X0, Y0, D, L, W, TH, DL, HH, RD

    REAL(KIND=dp)           :: att_mw
    integer                 ::   att_epi_i,  att_epi_j

    ! Set nc output filename

    netcdf_file  = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'.nc'

    ncfile = netcdf_file
    


!----- open file and write global attributes

    s = 1

    stat(s) = NF_CREATE(TRIM(ncfile), 0, fileid) 
    s = s + 1
    attstr  = 'Tsunami Simulation'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'title',       LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1
    attstr  = 'AWI'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'institution', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'Conventions', LEN_TRIM('CF-1.0'),  &
                              'CF-1.0') 

    IF (enable_ruptgen_scenario) THEN

       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'TSID',     LEN_TRIM(TSID), TSID) 
       
       s = s + 1

       attstr  = 'RuptGen V2.0 - Rupture Generator provided by GFZ Potsdam, Germany'
       stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'sourceDescription', LEN_TRIM(attstr),    &
            TRIM(attstr)) 
       s = s + 1

        stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Input_Magnitude',         NF_DOUBLE, 1, eq_mag) 
        s = s + 1
        stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Input_Epicenter_X',       NF_DOUBLE, 1, eq_epi_x) 
        s = s + 1
        stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Input_Epicenter_Y',       NF_DOUBLE, 1, eq_epi_y) 
       ! s = s + 1
       !stat(s) = NF_PUT_ATT_INT(fileid,    NF_GLOBAL, 'Tsunami_Input_Direction_Rupture', NF_FLOAT,  1, eq_dir_rupt) 
        s = s + 1
        attValInt=Epi_i
        stat(s) = NF_PUT_ATT_INT(fileid,    NF_GLOBAL,  'sourceParameterRuptGen_i',          NF_INT,    1, attValInt)    
        s = s + 1

        attValInt=Epi_j
        stat(s) = NF_PUT_ATT_INT(fileid,    NF_GLOBAL,  'sourceParameterRuptGen_j',          NF_INT,    1, attValInt)    
        s = s + 1

        attval=mw
        stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'sourceParameterRuptGen_mw',         NF_DOUBLE, 1, attval) 
        s = s + 1

        attstr  = 'Topography data: SRTM X- and C-Band data Interpolated to X-Band resolution, adjusted in Kuta/Bali by local data'
        stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'meshTopographyDataDescription', LEN_TRIM(attstr),    &
             TRIM(attstr)) 
        s = s + 1
        
        attstr  = 'Bathymetry data: GEBCO 30sec, SONNE Cruise SO 186, SCOTT Cruise, Seamap Data; Padang Data by Franzius Institut'
        stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'meshBathymetryDataDescription', LEN_TRIM(attstr),    &
             TRIM(attstr))

     ELSE IF (enable_okada_scenario) THEN
 
        s = s + 1
       attstr  = 'Okada parameters'
       stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'sourceDescription', LEN_TRIM(attstr),    &
            TRIM(attstr)) 

        stat(s) = NF_PUT_ATT_INT(fileid,    NF_GLOBAL,  'number_of_fault_planes',          NF_INT,    1, nfault) 
        s = s + 1
     
       IF (nfault == 1) THEN
               X0 = AllFaults(1, 1)
               Y0 = AllFaults(2, 1)
               D  = AllFaults(3, 1)
               L  = AllFaults(4, 1)
               W  = AllFaults(5, 1)
               TH = AllFaults(6, 1)
               DL = AllFaults(7, 1)
               HH = AllFaults(8, 1)
               RD = AllFaults(9, 1)

           stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Plane_X0',              NF_DOUBLE, 1, X0) 
           s = s + 1
           stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Plane_Y0',              NF_DOUBLE, 1, Y0) 
           s = s + 1
           stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Slip',                  NF_DOUBLE, 1,  D) 
           s = s + 1
           stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Fault_Length',          NF_DOUBLE, 1,  L) 
           s = s + 1
           stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Fault_Width',           NF_DOUBLE, 1,  W) 
           s = s + 1
           stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Strike_Angle',          NF_DOUBLE, 1, TH) 
           s = s + 1
           stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Dip_Angle',             NF_DOUBLE, 1, DL) 
           s = s + 1
           stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Fault_Depth',           NF_DOUBLE, 1, HH) 
           s = s + 1
           stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Direction_Dislocation', NF_DOUBLE, 1, RD) 
       ELSE
         if (nfault <= 20) then
           DO i = 1, nfault

               WRITE(par_num, '(I2.2)') i

               X0 = AllFaults(1, faults(i))
               Y0 = AllFaults(2, faults(i))
               D  = AllFaults(3, faults(i))
               L  = AllFaults(4, faults(i))
               W  = AllFaults(5, faults(i))
               TH = AllFaults(6, faults(i))
               DL = AllFaults(7, faults(i))
               HH = AllFaults(8, faults(i))
               RD = AllFaults(9, faults(i))

               s = s + 1
               stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Plane_X0_'//par_num,              &
                   NF_DOUBLE, 1, X0) 
               s = s + 1
               stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Plane_Y0_'//par_num,              &
                   NF_DOUBLE, 1, Y0) 
               s = s + 1
               stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Slip_'//par_num,                  &
                   NF_DOUBLE, 1, D) 
               s = s + 1
               stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Fault_Length_'//par_num,          &
                   NF_DOUBLE, 1, L) 
               s = s + 1
               stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Fault_Width_'//par_num,           &
                   NF_DOUBLE, 1, W) 
               s = s + 1
               stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Strike_Angle_'//par_num,          &
                   NF_DOUBLE, 1, TH) 
               s = s + 1
               stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Dip_Angle_'//par_num,             &
                   NF_DOUBLE, 1, DL)
               s = s + 1
               stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Fault_Depth_'//par_num,           &
                   NF_DOUBLE, 1, HH) 
               s = s + 1
               stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'Tsunami_Okada_Direction_Dislocation_'//par_num, &
                   NF_DOUBLE, 1, RD) 
            END DO
         end if
      END IF

   elseif (enable_aifdr_scenario) then

       s = s + 1
       attstr  = 'Source provided by AIFDR'
       stat(s) = NF_PUT_ATT_TEXT(fileId, NF_GLOBAL, 'sourceDescription', len_trim(attstr),    &
                                 TRIM(attstr))
       s = s + 1

       stat(s) = NF_PUT_ATT_TEXT(fileId, NF_GLOBAL, 'source_zone', len_trim(source_zone), TRIM(source_zone))
       s = s + 1

       stat(s) = NF_PUT_ATT_TEXT(fileId, NF_GLOBAL, 'source_zone_id', len_trim(sz), sz)
       s = s + 1

       stat(s) = NF_PUT_ATT_INT(fileId,    NF_GLOBAL,  'sourceParameterEpi_i',   NF_INT,    1, epi_i)    
       s = s + 1

       stat(s) = NF_PUT_ATT_INT(fileId,    NF_GLOBAL,  'sourceParameter_Epi_j',  NF_INT,    1, epi_j)    
       s = s + 1

       stat(s) = NF_PUT_ATT_DOUBLE(fileId, NF_GLOBAL, 'sourceParameter_mw',   NF_DOUBLE, 1, mw)

    elseif (enable_idealised) then

       s = s + 1
       attstr  = 'Source: Idealised, cosine bell'
       stat(s) = NF_PUT_ATT_TEXT(fileId, NF_GLOBAL, 'sourceDescription', len_trim(attstr),    &
                                 TRIM(attstr))

       stat(s) = NF_PUT_ATT_DOUBLE(fileId,    NF_GLOBAL,  'sourceParameterEpi_lon',   NF_DOUBLE,    1, eq_epi_x)    
       s = s + 1

       stat(s) = NF_PUT_ATT_DOUBLE(fileId,    NF_GLOBAL,  'sourceParameter_Epi_lat',  NF_DOUBLE,    1, eq_epi_y)    
       s = s + 1

       stat(s) = NF_PUT_ATT_DOUBLE(fileId, NF_GLOBAL, 'sourceParameter_mag',   NF_DOUBLE, 1, eq_mag)

       

    END IF

    WRITE(StrMeshID, '(I3.3)') MeshID
    attstr  = 'TsunAWI Revision '//revision//' - Mesh: '//StrMeshID
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'source', LEN_TRIM(attstr),     &
                              TRIM(attstr)) 
    s = s + 1
    attstr  = 'no entry'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'history', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1
    attstr  = 'none'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'references', LEN_TRIM(attstr), &
                              TRIM(attstr)) 
    s = s + 1
    attstr  = 'none'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'comment', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1

    attstr  = 'AWI'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'provider', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1

    attstr  = 'TsunAWI'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'model', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1

    attstr  = revision
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'modelVersion', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1


    attstr  = 'regional'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'scope', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1

    notyetknown=-99999.

    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'globalMinSSH',         NF_DOUBLE, 1, notyetknown) 
    s = s + 1
    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'globalMaxSSH',         NF_DOUBLE, 1, notyetknown) 
    s = s + 1

    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'boundingBox_UpperRight_Lon', NF_DOUBLE, 1, real(BoundingBox_xmax,dp)) 
    s = s + 1
    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'boundingBox_UpperRight_Lat', NF_DOUBLE, 1,  real(BoundingBox_ymax,dp)) 
    s = s + 1
    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'boundingBox_LowerLeft_Lon',  NF_DOUBLE, 1,  real(BoundingBox_xmin,dp)) 
    s = s + 1
    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'boundingBox_LowerLeft_Lat',  NF_DOUBLE, 1,  real(BoundingBox_ymin,dp)) 

    
    s = s + 1
    attval=T_out
    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'outputTimestepSize',         NF_DOUBLE, 1, attval) 
    s = s + 1
    
    attval=T_end
    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'integrationTime',         NF_DOUBLE, 1, attval)
    s = s + 1
    
    attval=longestEdge
    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'spatialResolution_Coarse',         NF_DOUBLE, 1, attval)
    s = s + 1
    
    attval=shortestEdge
    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'spatialResolution_Fine',         NF_DOUBLE, 1, attval)
    s = s + 1

    attValInt=nod2D
    stat(s) = NF_PUT_ATT_INT(fileid,    NF_GLOBAL,  'degreeOfFreedom',          NF_INT,    1, attValInt)

    if (enable_benchmark .and. benchmark_ident<=2) then
       s = s + 1
       if  (benchmark_ident==1) attstr  = 'Bathymetry / topography data: Sloping beach experiment'
       if  (benchmark_ident==1) attstr  = 'Bathymetry / topography data: Monai beach channel experiment'

       stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'meshDataDescription', LEN_TRIM(attstr),    &
            TRIM(attstr)) 
    endif
    



!----- DEFINE DIMENSIONS ------------------------------

    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'nodes_2D',     nod2D, DimId_n2D)             
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'elements_2D', elem2D, DimId_el2D)        
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'one',   1, dim1)                           
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'two',   2, dim2)                           
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'three', 3, dim3)                         
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'four',  4, dim4)                          
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'iteration', NF_UNLIMITED, DimId_iter)

    DO i = 1,  s
       IF (stat(i) /= NF_NOERR) then
          PRINT *, 'nc_init, NetCDF error',stat(i),' in dimension definitions, no.', i
          print *, nf_strerror(stat(i))
       endif
    END DO


!----- DEFINE VARIABLES ---------------------------------

    s = 1

!----- coordinates

    dimarray(1) = dim2
    dimarray(2) = DimId_n2D
    
    stat(s) = NF_DEF_VAR(fileid, 'surface_locations', NF_FLOAT, 2, dimarray(1:2), VarId_loc2D) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'longitude',         NF_FLOAT, 1, DimId_n2D,     VarId_lon) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'latitude',          NF_FLOAT, 1, DimId_n2D,     VarId_lat) 
    s = s + 1

!----- grid

    dimarray(1) = dim3
    dimarray(2) = DimId_el2D

    stat(s) = NF_DEF_VAR(fileid, 'triangles',  NF_INT, 2, dimarray(1:2), VarId_trian) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'tri_area', NF_FLOAT, 1, DimId_el2D,    VarId_tri_area) 
    s = s + 1

!----- node indices
    stat(s) = NF_DEF_VAR(fileid, 'node_index', NF_INT, 1, DimId_n2D,     VarId_index) 
    s = s + 1

!----- numbers 
    stat(s) = NF_DEF_VAR(fileid, 'iteration',  NF_INT, 1, DimId_iter,    VarId_iter) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'time',     NF_FLOAT, 1, DimId_iter,    VarId_time) 
    s = s + 1

!----- energy
    if (write_energy) then
       stat(s) = NF_DEF_VAR(fileid, 'epot',  NF_FLOAT, 1, DimId_iter,    VarId_epot) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'ekin',  NF_FLOAT, 1, DimId_iter,    VarId_ekin) 
       s = s + 1
    endif
       
!----- topography
    stat(s) = NF_DEF_VAR(fileid,'initial_topography',NF_FLOAT, 1, DimId_n2D, VarId_topo_ini) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid,'topography',NF_FLOAT, 1, DimId_n2D,    VarId_topo) 
    s = s + 1

!----- F I E L D S 

!-----  scalar variables

    dimarray(1) = DimId_n2D
    dimarray(2) = DimId_iter

    stat(s) = NF_DEF_VAR(fileid, 'ssh', NF_FLOAT, 2, dimarray(1:2), VarId_ssh); 
    s = s + 1

!----- arrival times, max. wave height, max. abs. velocity, max. flux, initial uplift

    stat(s) = NF_DEF_VAR(fileid, 'first_arrival',     NF_FLOAT, 1, DimId_n2D, VarId_arrival)
 
    if (write_ttt) then
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'first_arrival_ttt', NF_FLOAT, 1, DimId_n2D, VarId_ttt) 
    endif
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'max_wave_height',   NF_FLOAT, 1, DimId_n2D, VarId_mwh) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'max_abs_vel',       NF_FLOAT, 1, DimId_n2D, VarId_vel) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'max_flux',          NF_FLOAT, 1, DimId_n2D, VarId_flux) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'initial_uplift',    NF_FLOAT, 1, DimId_n2D, VarId_uplift) 

!----- vector variables

!-----  velocity in nodes
    IF (write_velocity) THEN
       dimarray(1) = dim2
       dimarray(2) = DimId_n2D
       dimarray(3) = DimId_iter

       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'velocity_in_nodes', NF_FLOAT, 3, dimarray(1:3), VarId_nodal_vel)
    END IF

    IF (write_cfl) THEN
       dimarray(1) = DimId_n2D
       dimarray(2) = DimId_iter

       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'cfl', NF_FLOAT, 2, dimarray(1:2), VarId_nodal_cfl)
       
       s = s + 1
       attstr = 'cfl in nodes, velocity vs egdlen/dt'
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_nodal_cfl, 'long_name', LEN_TRIM(attstr), TRIM(attstr))
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_nodal_cfl, 'field',       19, 'cfl, scalar, series') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_nodal_cfl, 'connections', 20, 'triangles, triangles') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_nodal_cfl, 'positions',   17, 'surface_locations')

       ! no units, because the CFL is dimensionless
    END IF
 
    DO i = 1,  s
        IF (stat(i) /= NF_NOERR) then
            PRINT *, 'nc_init, NetCDF error',stat(i),' in variable definition, no.', i
            print *, nf_strerror(stat(i))
         endif
    END DO

!----- DEFINE ATTRIBUTES -----------------------------------------

    s = 1

    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_loc2D,    'long_name',   19, 'surface_coordinates') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_loc2D,    'units',       18, 'degrees east/north') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_lon,      'long_name',    9, 'longitude') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_lon,      'units',       12, 'degrees_east') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_lat,      'long_name',    8, 'latitude') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_lat,      'units',       13, 'degrees_north') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_index,    'long_name',   16, 'nodal attributes') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tri_area, 'long_name',   13, 'triangle area') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssh,      'long_name',   21, 'sea surface elevation') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssh,      'units',        5, 'meter') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssh,      'field',       19, 'ssh, scalar, series') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssh,      'connections', 20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssh,      'positions',   17, 'surface_locations') 
    s = s + 1
    
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo_ini,     'long_name',   18, 'initial_topography') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo_ini,     'field',       18, 'initial_topography') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo_ini,     'connections', 20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo_ini,     'positions',   17, 'surface_locations') 
    s = s + 1
    attstr  = 'meters (before shock bathymetry and topography)'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo_ini,     'units',       LEN_TRIM(attstr), attstr) 
    s = s + 1
    
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'long_name',   10, 'topography') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'field',       10, 'topography') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'connections', 20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'positions',   17, 'surface_locations') 
    s = s + 1
    
    attstr  = 'meters (after shock bathymetry and topography)'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'units',       LEN_TRIM(attstr), attstr) 
    s = s + 1
    
    attstr='initial uplift due to earthquake'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uplift,   'long_name',   LEN_TRIM(attstr), attstr) 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uplift,   'field',       14, 'initial_uplift') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uplift,   'connections', 20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uplift,   'positions',   17, 'surface_locations') 
    s = s + 1
    
    attstr  = 'first wave arrival (sea surface elevation larger than +0.01m)'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_arrival,   'long_name',   LEN_TRIM(attstr), attstr) 
    s = s + 1    
    attstr = 'eta_thd_orig'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_arrival,   'eta_type', LEN_TRIM(attstr),attstr)
    s = s + 1
    attstr  = 'seconds after rupture'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_arrival,   'units',       LEN_TRIM(attstr), attstr) 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_arrival,   'field',       13, 'first_arrival') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_arrival,   'connections', 20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_arrival,   'positions',   17, 'surface_locations') 

    if (write_ttt) then
       s = s + 1
       attstr  = 'first wave arrival (estimated arrival time - ttt)'
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ttt      , 'long_name',   LEN_TRIM(attstr), attstr) 
       s = s + 1
       attstr = 'ttt_est_orig'
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ttt,   'eta_type', LEN_TRIM(attstr),attstr)
       s = s + 1
 
       attstr  = 'seconds after rupture'
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ttt,       'units',       LEN_TRIM(attstr), attstr) 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ttt,       'field',       17, 'first_arrival_ttt') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ttt,       'connections', 20, 'triangles, triangles') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ttt,         'positions',   17, 'surface_locations') 
    endif

    s = s + 1
    attstr  = 'maximum wave height\n (above mean sea level on sea \n and above after shock topography on land)'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'long_name',   LEN_TRIM(attstr), attstr) 
    s = s + 1
    
    attstr  = 'meters'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'units',       LEN_TRIM(attstr), attstr) 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'field',       15, 'max_wave_height') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'connections', 20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'positions',   17, 'surface_locations') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'long_name',   36, 'maximum nodal velocity approximation') 
    s = s + 1
    
    attstr  = 'm/sec'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'units',       LEN_TRIM(attstr), attstr) 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'field',       11, 'max_abs_vel'); s=s+1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'connections', 20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'positions',   17, 'surface_locations') 
    s = s + 1

    attstr  = 'maximum flux in nodes\n (maximum value of the product of velocity and wave height)'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,      'long_name',   LEN_TRIM(attstr), attstr) 
    s = s + 1
    
    attstr  = 'm^2/sec'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,     'units',        LEN_TRIM(attstr), attstr) 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,     'field',         8, 'max_flux') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,     'connections',  20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,     'positions',    17, 'surface_locations') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'long_name',     4, 'time'); 
    s = s + 1
    
    attstr  = 'seconds since rupture'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'units',        LEN_TRIM(attstr), attstr) 
    s = s + 1
    
    attstr  = '0'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'time_origin',  LEN_TRIM(attstr), attstr) 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'field',        12, 'time, series') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_iter,     'long_name',     9, 'iteration') 


    if (write_energy) then
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_epot,   'units',        6, 'Joules') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ekin,   'units',        6, 'Joules') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_epot,   'long_name',   22, 'total potential energy')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ekin,   'long_name',   20, 'total kinetic energy')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_epot,   'field',       14, 'energy, series')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ekin,   'field',       14, 'energy, series')
    endif

    
    DO i = 1, s
        IF (stat(i) /= NF_NOERR) then
            PRINT *, 'nc_init, NetCDF error',stat(i),' in attribute assignments, no.', i
            print *, nf_strerror(stat(i))
         endif
    END DO

    s = 1

    stat(s) = NF_ENDDEF(fileid) 
    

!----- write constant variables

!----- coordinates

    IF (coordinate_type == 1) THEN
        AUX_loc(1, :) = coord_nod2D(1, :) / rad
        AUX_loc(2, :) = coord_nod2D(2, :) / rad
    ELSE
        AUX_loc(1, :) = coord_nod2D(1, :)
        AUX_loc(2, :) = coord_nod2D(2, :)
    END IF

    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_loc2D,   AUX_loc(1:2, 1:nod2D)) 
    s = s + 1

    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_lon,     AUX_loc(1, 1:nod2D) ) 
    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_lat,     AUX_loc(2, 1:nod2D) ) 
    s = s + 1

!----- triangles

    stat(s) = NF_PUT_VAR_INT(fileid,  VarId_trian,    elem2D_nodes(:, :) - 1) 
    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_tri_area, REAL(voltriangle(:), sp)) 
    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_topo_ini, nodhn_init(:)) 
    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_topo,     REAL(nodhn(:), sp)) 
    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_uplift,   REAL(ssh_init(:), sp)) 
    s = s + 1
    stat(s) =NF_PUT_VAR_INT(fileid,   VarId_index,    index_nod2D(:)) 

    if (write_ttt) then
       s = s + 1
       stat(s) = NF_PUT_VAR_REAL(fileid, VarId_ttt,   ttt(:)) 
    endif

    DO i = 1, s
        IF (stat(i) /= NF_NOERR) then
            PRINT *,'nc_init, NetCDF error',stat(i),' in writing variables, no.', i
            print *, nf_strerror(stat(i))
         endif
    END DO

    stat(1) = NF_CLOSE(fileid)

    IF (stat(1) /= NF_NOERR) THEN
        PRINT *, 'nc_init, NetCDF error',stat(1),' in closing NetCDF file'
        print *, nf_strerror(stat(1))
    END IF

END SUBROUTINE nc_init
#endif




#ifndef NO_NETCDF
SUBROUTINE nc_out(iteration)

!----- uses

    USE PARAMETERS
    USE MESH
    USE ELEMENTS
    USE SWE

!----- local declarations

    IMPLICIT NONE

    include 'netcdf.inc'
    
    CHARACTER(len=100) 	         :: ncfile
    INTEGER, INTENT(in) 	 :: iteration
    INTEGER, SAVE 	         :: iter=1
    INTEGER 		         :: FileId
    INTEGER 		         :: VarId_iter, VarId_time ! numbers
    INTEGER 		         :: VarId_topo             ! topography (after rupture and smoothing)
    INTEGER 		         :: VarId_topo_ini         ! inital topography
    INTEGER 		         :: VarId_ssh              ! fields: ssh
    INTEGER 		         :: VarId_arrival          ! fields: Arrival times
    INTEGER 		         :: VarId_mwh              ! fields: maximum wave height
    INTEGER 		         :: VarId_flux, VarId_vel  !fields: flux and mean absolute velocity
    INTEGER		         :: VarId_nodal_vel 	   !fields: vector
    INTEGER		         :: VarId_nodal_cfl 	   !fields: cfl value
    INTEGER		         :: VarId_epot, VarId_ekin !fields: energies

    INTEGER 		         :: i, n, iA, iZ, edg, j
    INTEGER 		         :: pos1, nmb
    INTEGER, DIMENSION(2)      :: pos1vec, nmbvec
    INTEGER, DIMENSION(3)      :: pos1vec3, nmbvec3
 
    INTEGER,    DIMENSION(50)    :: stat                     ! auxiliaries
    INTEGER                      :: s                        ! auxiliaries
    INTEGER                      :: restart_iter, dimid_iter
    REAL,       DIMENSION(nod2D) :: AUX_n2D
    REAL,     DIMENSION(2,nod2D) :: AUX_n2D_2
    REAL                         :: acc, inv_acc
   
    ncfile = netcdf_file

    stat(1) = NF_OPEN(TRIM(ncfile), NF_WRITE, FileId)

    IF (stat(1) /= NF_NOERR) then
        print *,'nc_out, nc-file error:',stat(1)
        print *, nf_strerror(stat(1))
        STOP 1
     endif
    IF (read_restart .AND. iter == 1) THEN
        stat(2) = NF_INQ_DIMID(fileid,  'iteration', dimid_iter)
        stat(2) = NF_INQ_DIMLEN(fileid, dimid_iter,  restart_iter)
        iter    = restart_iter + 1
    END IF

!----- INQUIRE VARIABLE IDs

    s = 1

    stat(s) = NF_INQ_VARID(fileid, "iteration",         VarId_iter) 
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, "time",              VarId_time) 
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, "ssh",                  VarId_ssh) 
    if (write_diag_in_nc_snapshots) then
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "first_arrival",     VarId_arrival) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "max_wave_height",   VarId_mwh) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "max_abs_vel",       VarId_vel) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "max_flux",          VarId_flux) 
    endif

    IF (write_velocity) THEN
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "velocity_in_nodes", VarId_nodal_vel) 
    END IF
    IF (write_cfl) THEN
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "cfl",               VarId_nodal_cfl) 
    END IF
    IF (write_energy) THEN
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "epot",              VarId_epot) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "ekin",              VarId_ekin) 
    END IF
       
    

    DO i = 1, s 
       IF (stat(i) /= NF_NOERR) then
          PRINT *, 'nc_out, NetCDF error',stat(i),' inquiring variable IDs, no.', i
          print *, nf_strerror(stat(i))
       endif
    END DO

!----- WRITE VARIABLES
    s = 1

    stat(s) = NF_PUT_VARA_INT(fileid,  VarId_iter, iter, 1, iteration) 
    s = s + 1
    stat(s) = NF_PUT_VARA_REAL(fileid, VarId_time, iter, 1, REAL(time, sp)) 
    s = s + 1

    if (write_energy) then
       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_epot, iter, 1, REAL(epot, sp)) 
       s = s + 1
       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_ekin, iter, 1, REAL(ekin, sp)) 
       s = s + 1
    endif
       
!----- 2D FIELDS

    pos1vec = (/ 1, iter  /)
    nmbvec  = (/ nod2D, 1 /)

!----- sea surface height
!----- subtract terrain height on land
  
    AUX_n2D = 0.
    if (nc_snapshot_accuracy_ssh <= 0.) then
!$OMP PARALLEL DO 
       DO n = 1, nod2D
          IF (nodhn(n) < 0.) THEN
             AUX_n2D(n) = max(0., ssh1(n) + nodhn(n))
          ELSE
             AUX_n2D(n) = ssh1(n)
          END IF
       END DO
!$OMP END PARALLEL DO
    else
       acc = real(nc_snapshot_accuracy_ssh ,sp)
       inv_acc = 1._sp/nc_snapshot_accuracy_ssh 
!$OMP PARALLEL DO 
       DO n = 1, nod2D
          IF (nodhn(n) < 0.) THEN
             AUX_n2D(n) =  acc* INT(inv_acc*max(0., real(ssh1(n) + nodhn(n),sp)))
          ELSE
             AUX_n2D(n) =  acc* INT(inv_acc*real(ssh1(n),sp))
          END IF
       END DO
!$OMP END PARALLEL DO
    end if
       

    stat(s) = NF_PUT_VARA_REAL(fileid, VarId_ssh, pos1vec, nmbvec, AUX_n2D(:)) 

    if (write_diag_in_nc_snapshots) then
       s = s + 1
       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_arrival, 1, nod2D, arrival_time(:)) 
       s = s + 1
       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_mwh,     1, nod2D, mwh(:))
       s = s + 1
       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_vel,     1, nod2D, max_abs_vel(:))
       s = s + 1
       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_flux,    1, nod2D, max_flux(:))
    endif

!----- VECTOR FIELDS

!----- nodal velocity

    IF (write_velocity) THEN
       AUX_n2D_2(1,:)=uv_node(1,1:nod2D)
       AUX_n2D_2(2,:)=uv_node(2,1:nod2D)
       pos1vec3 = (/ 1, 1, iter /)
       nmbvec3  = (/ 2, nod2D, 1 /)
       s = s + 1
       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_nodal_vel, pos1vec3, nmbvec3, AUX_n2D_2)
    END IF

    IF (write_cfl) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP&            PRIVATE(n,iA,iZ,j,edg)     
       do n = 1,nod2D
          AUX_n2D(n)=0.
          if (depth(n) > small) then

             iA = nodedge_ptr(n)
             iZ = nodedge_ptr(n+1)-1
             ! The CFL is calculated on the edges and interpolated to the nodes:
             ! |vel|*dt / |edge|
             ! edg_length contains 1./edg_length here.
             do j=iA, iZ
                edg = nodedge_addr(j)
                AUX_n2D(n) = AUX_n2D(n) + nodedge_wght(j)* &
                     sqrt(uv2(1,edg)*uv2(1,edg) + uv2(2,edg)*uv2(2,edg)) &
                     * edg_length(edg)                      
             end do
             AUX_n2D(n) = AUX_n2D(n) *dt
             
          end if
       end do
!$OMP END PARALLEL DO    
       s = s + 1
       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_nodal_cfl, pos1vec, nmbvec, AUX_n2D)
    END IF

    DO i = 1, s 
       IF (stat(i) /= NF_NOERR) then
          PRINT *, 'nc_out, NetCDF error',stat(i),' in writing variables, no.', i
          print *, nf_strerror(stat(i))
       end IF
    END DO


!----- CLOSE THE FILE

    stat(1) = NF_CLOSE(fileid)

    IF (stat(1) /= NF_NOERR) THEN
        PRINT *, 'nc_out, NetCDF error',stat(1),' in closing NetCDF file'
        print *, nf_strerror(stat(1))
    END IF

    iter = iter + 1

END SUBROUTINE nc_out
#endif


#ifndef NO_NETCDF
subroutine nc_finalize

!----- uses

    USE PARAMETERS
    USE MESH
    USE ELEMENTS
    USE SWE

!----- local declarations

    IMPLICIT NONE

    include 'netcdf.inc'
    
    CHARACTER(len=100) 	    :: ncfile
    INTEGER 		    :: FileId
    INTEGER                 :: stat(32), s, i           ! auxiliaries
    REAL(KIND=dp)            :: attVal
    INTEGER 		    :: VarId_arrival            ! fields: Arrival times
    INTEGER 		    :: VarId_mwh                ! fields: maximum wave height
    INTEGER 		    :: VarId_flux, VarId_vel    !fields: flux and mean absolute velocity

    ncfile=netcdf_file

    stat(1) = NF_OPEN(TRIM(ncfile), NF_WRITE, FileId)

    if (stat(1) /= NF_NOERR) then
       print *,'nc-file error', stat(1)
       print *, nf_strerror(stat(1))
       stop 1
    endif
    attVal=-99999.
    s=0
    s=s+1; stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'globalMinSSH',         NF_DOUBLE, 1, attVal) 

    attVal=maxval(mwh(:))
    s=s+1; stat(1) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'globalMaxSSH',         NF_DOUBLE, 1, attVal) 

    do i=1,s
       if (stat(i) /= NF_NOERR) then
          PRINT *, 'nc_finalize, NetCDF error',stat(i),' in rewriting global attribute No.', i
          print *, nf_strerror(stat(i))
       endif
    enddo

    if (.not. write_diag_in_nc_snapshots) then
       
       if (verbosity >= 1) print *,'Writing MWH, ETA, Max Vel, Max Flux'

       s = 1
       stat(s) = NF_INQ_VARID(fileid, "first_arrival",     VarId_arrival) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "max_wave_height",   VarId_mwh) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "max_abs_vel",       VarId_vel) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "max_flux",          VarId_flux) 

       DO i = 1, s 
          IF (stat(i) /= NF_NOERR) then
             PRINT *, 'nc_finalize, NetCDF error',stat(i),', inquiring variable IDs, no.', i
             print *, nf_strerror(stat(i))
          endif
       END DO

       s = 1
       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_arrival, 1, nod2D, arrival_time(:)) 
       s = s + 1
       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_mwh,     1, nod2D, mwh(:))
       s = s + 1
       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_vel,     1, nod2D, max_abs_vel(:))
       s = s + 1
       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_flux,    1, nod2D, max_flux(:))
      
       DO i = 1, s 
          IF (stat(i) /= NF_NOERR) then
             PRINT *, 'nc_finalize, NetCDF error,',stat(i),' writing variables, no.', i
             print *, nf_strerror(stat(i))
          end IF
       END DO
    endif


    stat(1) = NF_CLOSE(fileid)
    if (stat(1) /= NF_NOERR) then
       print *, 'nc_finalize, NetCDF error',stat(1),' in closing NetCDF file'
       print *, nf_strerror(stat(1))
    endif

end subroutine nc_finalize
#endif

!=====================================================================================================

!=======================================================================================

SUBROUTINE bin_restart_in

  USE PARAMETERS
  USE MESH
  USE ELEMENTS
  USE SWE
  
  IMPLICIT NONE
  
  INTEGER :: i
  CHARACTER(LEN=100) :: fname
  INTEGER :: restartfile_ssh
  INTEGER :: restartfile_uv
  INTEGER :: restartfile_grad
  INTEGER :: restartfile_diag
  INTEGER :: restartfile_info
  INTEGER :: restartfile_nodhn

       if (enable_nc_out) netcdf_file  = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'.nc'

       fname=TRIM(restart_file)//'_nodhn'
       OPEN(newunit=restartfile_nodhn,file=fname,form='unformatted')
       read(restartfile_nodhn) nodhn(:)
       CLOSE(restartfile_nodhn)

       fname=TRIM(restart_file)//'_ssh'

       OPEN(newunit=restartfile_ssh,file=fname,form='unformatted')
       READ(restartfile_ssh) ssh0(:)
       READ(restartfile_ssh) ssh1(:)
       READ(restartfile_ssh) ssh2(:)
       CLOSE(restartfile_ssh)

       fname=TRIM(restart_file)//'_uv'

       OPEN(newunit=restartfile_uv,file=fname,form='unformatted')
       READ(restartfile_uv) uv0(:,:)
       READ(restartfile_uv) uv1(:,:)
       READ(restartfile_uv) uv2(:,:)
       READ(restartfile_uv) uv_node(:,:)
       CLOSE(restartfile_uv)

       fname=TRIM(restart_file)//'_grad'

       OPEN(newunit=restartfile_grad,file=fname,form='unformatted')
       READ(restartfile_grad) grad(:,:)
       CLOSE(restartfile_grad)


       fname=TRIM(restart_file)//'_diag'

       OPEN(newunit=restartfile_diag, file=fname,form='unformatted')
       READ(restartfile_diag) arrival_time(:)
       READ(restartfile_diag) mwh(:)
       READ(restartfile_diag) max_flux(:)
       READ(restartfile_diag) max_abs_vel(:)
       CLOSE(restartfile_diag)
       
       fname=TRIM(restart_file)//'_info'

       OPEN(newunit=restartfile_info,file=fname,form='unformatted')
       READ(restartfile_info) time0
       CLOSE(restartfile_info)

        time=time0

END SUBROUTINE bin_restart_in

!========================================

SUBROUTINE bin_restart_out

  USE PARAMETERS
  USE MESH
  USE ELEMENTS
  USE SWE
  
  IMPLICIT NONE
  
  INTEGER :: i
  CHARACTER(LEN=100) :: fname
  INTEGER :: restartfile_ssh
  INTEGER :: restartfile_uv
  INTEGER :: restartfile_grad
  INTEGER :: restartfile_diag
  INTEGER :: restartfile_info
  INTEGER :: restartfile_nodhn
  INTEGER :: outfile_status

  fname=TRIM(restart_file)//'_ssh'
  
  OPEN(newunit=restartfile_ssh,file=fname,form='unformatted')
  WRITE(restartfile_ssh) ssh0(:)
  WRITE(restartfile_ssh) ssh1(:)
  WRITE(restartfile_ssh) ssh2(:)
  CLOSE(restartfile_ssh)

  fname=TRIM(restart_file)//'_uv'

  OPEN(newunit=restartfile_uv,file=fname,form='unformatted')
  WRITE(restartfile_uv) uv0(:,:)
  WRITE(restartfile_uv) uv1(:,:)
  WRITE(restartfile_uv) uv2(:,:)
  WRITE(restartfile_uv) uv_node(:,:)
  CLOSE(restartfile_uv)

  fname=TRIM(restart_file)//'_grad'

  OPEN(newunit=restartfile_grad,file=fname,form='unformatted')
  WRITE(restartfile_grad) grad(:,:)
  CLOSE(restartfile_grad)

  fname=TRIM(restart_file)//'_diag'

  OPEN(newunit=restartfile_diag,file=fname,form='unformatted')
  WRITE(restartfile_diag) arrival_time(:)
  WRITE(restartfile_diag) mwh(:)
  WRITE(restartfile_diag) max_flux(:)
  WRITE(restartfile_diag) max_abs_vel(:)
  CLOSE(restartfile_diag)
  
  fname=TRIM(restart_file)//'_info'

  OPEN(newunit=restartfile_info,file=fname,form='unformatted')
  WRITE(restartfile_info) time
  CLOSE(restartfile_info)

  fname=TRIM(restart_file)//'_nodhn'

  OPEN(newunit=restartfile_nodhn,file=fname,form='unformatted')
  WRITE(restartfile_nodhn) nodhn(:)
  CLOSE(restartfile_nodhn)
  
  OPEN(newunit=outfile_status,file='status.txt')
  write(outfile_status,*) 'go on'
  CLOSE(outfile_status)
  
END SUBROUTINE bin_restart_out



SUBROUTINE fill_sat_aux_arrays

  USE MESH

  IMPLICIT NONE

  INTEGER :: i

  OPEN(25,file='/work/ab0042/a270038/GRIDS/AUX/TrBF.IndianOcean_hi')
  OPEN(26,file='/work/ab0042/a270038/GRIDS/AUX/TrElems.IndianOcean_hi')
  do i=1,551
    read(25,*) TrBF_J1(i,:)
    read(26,*) TrElems_J1(i)
  end do
  CLOSE(25)
  CLOSE(26)

  OPEN(25,file='/work/ab0042/a270038/GRIDS/AUX/TrBF.TP.IndianOcean_hi')
  OPEN(26,file='/work/ab0042/a270038/GRIDS/AUX/TrElems.TP.IndianOcean_hi')
  do i=1,378
    read(25,*) TrBF_TP(i,:)
    read(26,*) TrElems_TP(i)
  end do
  CLOSE(25)
  CLOSE(26)
    
END SUBROUTINE fill_sat_aux_arrays

!==================================================================================


SUBROUTINE ascii_out_tidegauge_init

  USE PARAMETERS
  USE MESH
  USE SWE

  IMPLICIT NONE

! The mareogram for each TG is extracted from the nearest grid node
! with a water depth of mindepth. If such a node is more than distmax
! away from the TG location, the TG is skipped.   
  real(kind=dp), parameter        :: distmax=10000.  ! in meters (here: 10km)
  real(kind=dp), parameter        :: mindepth=1.  ! in meters (here: 1m)

  real(kind=dp)                   :: mindepth_gg
  integer                        :: n_all
  integer                        :: tg, n, count, io_stat, pos1, pos2
  real(kind=dp), allocatable      :: lon(:), lat(:), dist(:)
  real(kind=dp)                   :: maxlon, minlon, maxlat, minlat
  character(len=20)              :: tg_id
  character(len=1)               :: char
  character(len=100)             ::  dummy 
  character(len=100), allocatable:: tg_urn(:)
  integer, allocatable           :: aux(:)
  real(kind=dp)                   :: r
  character(len=4)                :: land_id
  
  ! Define identifier for gauges on land (for inundation studies)
  land_id='LAND';

  open(55,file=TRIM(TideGaugeFile))
  io_stat = 0
  n_all = -1
  DO WHILE (io_stat==0)
     read(55,*,iostat=io_stat) char
     n_all = n_all+1
  END DO
  close(55)
  print *,'Total number of tide gauges: ',n_all

  
  allocate(lon(n_all), lat(n_all))
  allocate(aux(n_all), dist(n_all))
  allocate(tg_urn(n_all))

  open(55,file=TRIM(TideGaugeFile))
  DO tg = 1, n_all
     read(55,*) tg_urn(tg), lon(tg), lat(tg)
  END DO
  close(55)

! Find nearest grid node in water. If the distance between
! node and TG location is larger than distmax, reject tidegauge 

  aux(:)  = 0
  dist(:) = 4.*pi*pi
  count = 0

! computing grid bounding box
  maxlon = (maxval(coord_nod2D(1,:)) + distmax/r_earth) / rad
  minlon = (minval(coord_nod2D(1,:)) - distmax/r_earth) / rad
  maxlat = (maxval(coord_nod2D(2,:)) + distmax/r_earth) / rad
  minlat = (minval(coord_nod2D(2,:)) - distmax/r_earth) / rad

  DO tg = 1, n_all
     if (lat(tg) > minlat .and. lat(tg) < maxlat .and. &
         lon(tg) > minlon .and. lon(tg) < maxlon) then

        mindepth_gg=mindepth
        if (tg_urn(tg)(1:4)==land_id) then
           print *,'Land gauge'
           mindepth_gg=-100.0
        end if

        DO n=1,nod2D
           if (nodhn(n) >= mindepth_gg) then
              ! Remind to convert lon, lat to radian
              r = (rad*lon(tg) - coord_nod2D(1,n))**2 + (rad*lat(tg) - coord_nod2D(2,n))**2 
              if ( r < dist(tg)) then
                 aux(tg) = n
                 dist(tg) = r
              endif
           end if
        ENDDO
     
        dist(tg) = sqrt(dist(tg)) *r_earth
        
        if (dist(tg) <= distmax) then
           count = count+1
           if (verbosity >= 5) &
                print *,trim(tg_urn(tg)),': nearest sea node ',aux(tg),', distance ',dist(tg),'m'
        elseif (verbosity >= 5) then
                print *,'Tidegauge ',trim(tg_urn(tg)),' not considered: distance [m] to nearest sea node',dist(tg)
        endif

     elseif (verbosity >= 5) then
        dist(tg) = sqrt(dist(tg)) *r_earth
        print *,'Tidegauge ',trim(tg_urn(tg)),' not considered: outside bounding box of the computational domain'
     endif
     
  ENDDO

  n_tidegauges = count

 if (verbosity >= 1)  print *,'Number of tide gauges taken into account: ',n_tidegauges
  allocate(tidegauge_node_index(n_tidegauges))
  allocate(tidegauge_ascii_filename(n_tidegauges))

  n = 0
  DO tg = 1,n_all
     if (dist(tg) <= distmax) then
        n=n+1

        dummy = tg_urn(tg)
        if (LEN_TRIM(dummy)>8) then
          pos1 = index(tg_urn(tg),':',back=.true.) +1
          pos2 = index(tg_urn(tg),' ',back=.false.) -1
          tg_id = dummy(pos1:pos2)
        else
          tg_id = trim(dummy)
       end if

        tidegauge_node_index(n) = aux(tg)
        tidegauge_ascii_filename(n) = 'TG_'//trim(tg_id)//'_'//trim(TSID)//'.out'

        if (verbosity >= 5) &
             print *, trim(tidegauge_ascii_filename(n))

        if (.not. read_restart) then
           ! write the information in the preambel
           open(66,file=tidegauge_ascii_filename(n))
           write(66,'(a2$)') "% "
           write(66,*) trim(tg_urn(tg))
           write(66,'("% TG longitude, latitude      ",2f12.6)') lon(tg), lat(tg)
           if (tg_urn(tg)(1:4)==land_id) then
              write(66,'("% Index of nearest land grid node ",i12)') tidegauge_node_index(n)
           else 
              write(66,'("% Index of nearest grid node with depth>1m",i12)') tidegauge_node_index(n)
           end if
           write(66,'("% Grid node longitude latitude",2f12.6)') coord_nod2D(1:2,tidegauge_node_index(n))/rad
           write(66,'("% Distance [m] to grid node",f12.4)') dist(tg)
           if (tg_urn(tg)(1:4)==land_id) then
              write(66,'("% Land height [m] at grid node before rupture",f12.4)') &
                   -nodhn_init(tidegauge_node_index(n))
              write(66,'("% Land height [m] at grid node after rupture ",f12.4)') &
                   -nodhn(tidegauge_node_index(n))
           else
              write(66,'("% Water depth [m] at grid node before rupture",f12.4)') &
                   nodhn_init(tidegauge_node_index(n))
              write(66,'("% Water depth [m] at grid node after rupture ",f12.4)') &
                   nodhn(tidegauge_node_index(n))
           end if
           write(66,'(a52)') "% time since rupture [min]   sea surface elevation [m]"
           write(66,'(f15.2,f13.5)') time/60., ssh2(tidegauge_node_index(n))        
           close(66)
        endif
     endif
  ENDDO


END SUBROUTINE ascii_out_tidegauge_init


!==================================================================================


SUBROUTINE ascii_out_tidegauge

  USE PARAMETERS
  USE MESH
  USE SWE

  IMPLICIT NONE

  integer :: tg

  DO tg = 1, n_tidegauges
     open(66,file=tidegauge_ascii_filename(tg),status='old',position='append',action='write')
     write(66,'(f15.2,f13.5)') time/60., ssh2(tidegauge_node_index(tg))
     close(66)
  ENDDO


END SUBROUTINE ascii_out_tidegauge


!====================================================================================

SUBROUTINE write_checksums

  USE PARAMETERS
  USE MESH
  USE ELEMENTS
  USE SWE
  
!---- local declarations

  IMPLICIT NONE
  CHARACTER(len=100)    :: checksum_file
  integer               :: u_checksum
  
  real(kind=wp)          :: ssh2max=-1000._wp, ssh2min=1000._wp, ssh2sum=0._wp
  real(kind=wp)          :: u2max  =-1000._wp, u2min  =1000._wp,   u2sum=0._wp
  real(kind=wp)          :: v2max  =-1000._wp, v2min  =1000._wp,   v2sum=0._wp
  integer               :: n, ed

!$OMP PARALLEL DEFAULT(shared) PRIVATE(n,ed) &
!$OMP          REDUCTION(+:  ssh2sum, u2sum, v2sum) &
!$OMP          REDUCTION(max:ssh2max, u2max, v2max) &
!$OMP          REDUCTION(min:ssh2min, u2min, v2min) 

!$OMP DO 
  do n=1,nod2D
     ssh2max = max(ssh2max,ssh2(n))
     ssh2min = min(ssh2min,ssh2(n))
     ssh2sum = ssh2sum +   ssh2(n)
  enddo
!$OMP END DO NOWAIT
!$OMP DO 
  do ed=1,edg2D
     u2max = max(u2max,uv2(1,ed));  v2max = max(v2max,uv2(2,ed)); 
     u2min = min(u2min,uv2(1,ed));  v2min = min(v2min,uv2(2,ed)); 
     u2sum = u2sum +   uv2(1,ed);   v2sum = v2sum +   uv2(2,ed);     
  enddo
!$OMP END DO
!$OMP END PARALLEL

  checksum_file   = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_checksums.out'
  
  open(newunit=u_checksum, file=checksum_file, status='replace', action='write')
  
  write(u_checksum,*) "max_ssh2=",ssh2max
  write(u_checksum,*) "min_ssh2=",ssh2min
  write(u_checksum,*) "sum_ssh2=",ssh2sum
  write(u_checksum,*) "max_u2=  ",u2max
  write(u_checksum,*) "min_u2=  ",u2min
  write(u_checksum,*) "sum_u2=  ",u2sum
  write(u_checksum,*) "max_v2=  ",v2max
  write(u_checksum,*) "min_v2=  ",v2min
  write(u_checksum,*) "sum_v2=  ",v2sum

  close(u_checksum)

END SUBROUTINE write_checksums


#ifndef NO_NETCDF
SUBROUTINE nc_inundation

!----- uses

    USE PARAMETERS
    USE MESH
    USE ELEMENTS
    USE SWE
    USE OKADA_FAULT_PARAMETERS

!-----local declarations

    IMPLICIT NONE

    include 'netcdf.inc'

! For reduced mesh (inundated area, mwh > threshold on land)    
    REAL(KIND=sp), parameter    :: mwh_threshold=0.01
    INTEGER                     :: inun_nod2D, inun_elem2D
    INTEGER                     :: n, k, el
    INTEGER                     :: inun_nod(nod2D), inun_nod_inv(nod2D)
    INTEGER                     :: inun_elem_nod(3,elem2D), inun_elem(elem2D)
    REAL(KIND=dp)               :: inun_bbox_xmax, inun_bbox_xmin
    REAL(KIND=dp)               :: inun_bbox_ymax, inun_bbox_ymin
    
    CHARACTER(len=100)          :: ncfile
    CHARACTER(len=200)          :: attstr
    INTEGER                     :: attValInt
    REAL(KIND=dp)               :: attVal, notyetknown
    CHARACTER(len=2)            :: par_num
    CHARACTER(len=3)            :: StrMeshID
    INTEGER                     :: i
    INTEGER                     :: fileid 			     ! IDs netCDF file
    INTEGER                     :: DimId_n2D			     ! dimension: nodes
    INTEGER                     :: DimId_el2D			     ! dimension: elements
    INTEGER                     :: Dim1, Dim2, Dim3, Dim4	     ! dimensions: 
    INTEGER                     :: VarId_trian			     ! variable: triangle
    INTEGER                     :: VarId_loc2D, VarId_lon, VarId_lat ! variables: 2D location, longitude, latitude
    INTEGER                     :: VarId_tri_area		     ! variable: area
    INTEGER                     :: VarId_topo, VarId_index	     ! variables: topography, index
    INTEGER                     :: VarId_topo_ini                    ! variables: initial topography
    INTEGER                     :: VarId_mwh  			     ! variable: maximum wave height
    INTEGER                     :: VarId_flux, VarId_vel             ! variables: flux and mean absolute velocity
    INTEGER                     :: VarId_uplift                      ! variable: topo/bathy uplift
    INTEGER                     :: s       			     ! auxiliari: status counter
    INTEGER, DIMENSION(3)       :: dimarray			     ! auxiliari: array dimension
    INTEGER, DIMENSION(500)     :: stat    			     ! auxiliari: status array
        
    REAL(kind=sp), allocatable  :: aux_loc(:,:)
    REAL(kind=sp), parameter    :: geo_factor = 1./rad

    REAL(KIND=dp)  ::   X0, Y0, D, L, W, TH, DL, HH, RD

    REAL(KIND=dp)  :: att_mw
    integer ::   att_epi_i,  att_epi_j

    ! Set nc output filename

    ncfile = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_inundation.nc'

    ! Reduce the mesh to inundated area, i.e., elements with
    ! nodes on land and mwh >= threshold for at least one node
    ! in an element (we need the neighbouring nodes for the interpolation)

    inun_nod_inv(1:nod2D) = -1
    k=0
    DO el=1,elem2D
       if (all(nodhn_init(elem2D_nodes(1:3,el)) <= 0) .and. &
            any( mwh(elem2D_nodes(1:3,el)) >= mwh_threshold) ) then
          ! collect the element
          k=k+1
          inun_elem(k) = el
          ! mark all nodes as to be collected
          inun_nod_inv(elem2D_nodes(1:3,el)) = 0  
       endif
    ENDDO
    inun_elem2D = k

    if (inun_elem2D == 0) then
       print *,'!!! There is no inundated area in the computational mesh. !!!'
       print *,trim(ncfile), ' will not be written.'
       return
    endif

    ! Now, collect all nodes on land needed for the interpolation
    ! of the inundation
    k = 0
    DO n=1,nod2D
       if ( inun_nod_inv(n) == 0 ) then
          k = k+1
          inun_nod(k) = n
          inun_nod_inv(n) = k
       endif
    END DO
    inun_nod2D = k

    ! Finally, adjust the numbering of the inundated elements' nodes
    DO k = 1,inun_elem2D
       inun_elem_nod(1:3,k) = inun_nod_inv(elem2D_nodes(1:3,inun_elem(k)))
    ENDDO

    ! coordinates and bounding box
    allocate(aux_loc(2,inun_nod2D))
    
    IF (coordinate_type == 1) THEN
       do k = 1, inun_nod2D
          aux_loc(1:2, k) = coord_nod2D(1:2, inun_nod(k)) * geo_factor
       end do
    else
       do k = 1, inun_nod2D
          aux_loc(1:2, k) = coord_nod2D(1:2, inun_nod(k))
       end do
    END IF

    inun_bbox_xmin = aux_loc(1, 1)
    inun_bbox_xmax = aux_loc(1, 1)
    inun_bbox_ymin = aux_loc(2, 1)
    inun_bbox_ymax = aux_loc(2, 1)

    do k = 2, inun_nod2D
       inun_bbox_xmin = min(inun_bbox_xmin, aux_loc(1, k))
       inun_bbox_xmax = max(inun_bbox_xmax, aux_loc(1, k))
       inun_bbox_ymin = min(inun_bbox_ymin, aux_loc(2, k))
       inun_bbox_ymax = max(inun_bbox_ymax, aux_loc(2, k))
    enddo
    
    

!----- open file and write global attributes

    s = 1

    stat(s) = NF_CREATE(TRIM(ncfile), 0, fileid) 
    s = s + 1
    attstr  = 'Tsunami Simulation - Inundation'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'title',       LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1
    attstr  = 'AWI'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'institution', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'Conventions', LEN_TRIM('CF-1.0'),  &
                              'CF-1.0') 


    WRITE(StrMeshID, '(I3.3)') MeshID
    attstr  = 'TsunAWI Revision '//revision//' - Mesh: '//StrMeshID
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'source', LEN_TRIM(attstr),     &
                              TRIM(attstr)) 
    s = s + 1
    attstr  = 'no entry'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'history', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1
    attstr  = 'none'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'references', LEN_TRIM(attstr), &
                              TRIM(attstr)) 
    s = s + 1
    attstr  = 'none'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'comment', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1

    attstr  = 'AWI'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'provider', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1

    attstr  = 'TsunAWI'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'model', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1

    attstr  = revision
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'modelVersion', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1


    attstr  = 'inundation'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'scope', LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1

    
    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'boundingBox_UpperRight_Lon',         NF_DOUBLE, 1, real(inun_bbox_xmax,dp)) 
    s = s + 1

    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'boundingBox_UpperRight_Lat',         NF_DOUBLE, 1,  real(inun_bbox_ymax,dp)) 
    s = s + 1

    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'boundingBox_LowerLeft_Lon',         NF_DOUBLE, 1,  real(inun_bbox_xmin,dp)) 
    s = s + 1

    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'boundingBox_LowerLeft_Lat',         NF_DOUBLE, 1,  real(inun_bbox_ymin,dp)) 

    
    s = s + 1
    
    attval=T_end
    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'integrationTime_in_seconds',         NF_DOUBLE, 1, attval)
    s = s + 1

    attValInt=inun_nod2D
    stat(s) = NF_PUT_ATT_INT(fileid,    NF_GLOBAL,  'degreeOfFreedom',          NF_INT,    1, attValInt)



!----- DEFINE DIMENSIONS ------------------------------

    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'nodes_2D',     inun_nod2D, DimId_n2D)             
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'elements_2D', inun_elem2D, DimId_el2D)        
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'one',   1, dim1)                           
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'two',   2, dim2)                           
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'three', 3, dim3)                         
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'four',  4, dim4)                          

    DO i = 1,  s
       IF (stat(i) /= NF_NOERR) then
          PRINT *, 'nc_inundation, NetCDF error',stat(i),' in dimension definitions, no.', i
          print *, nf_strerror(stat(i))
       end IF
    END DO


!----- DEFINE VARIABLES ---------------------------------

    s = 1

!----- coordinates

    dimarray(1) = dim2
    dimarray(2) = DimId_n2D
    
    stat(s) = NF_DEF_VAR(fileid, 'surface_locations', NF_FLOAT, 2, dimarray(1:2), VarId_loc2D) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'longitude',         NF_FLOAT, 1, DimId_n2D,     VarId_lon) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'latitude',          NF_FLOAT, 1, DimId_n2D,     VarId_lat) 
    s = s + 1

!----- grid

    dimarray(1) = dim3
    dimarray(2) = DimId_el2D

    stat(s) = NF_DEF_VAR(fileid, 'triangles', NF_INT,   2, dimarray(1:2), VarId_trian) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'tri_area',  NF_FLOAT, 1, DimId_el2D,    VarId_tri_area) 
    s = s + 1


!----- topography

    stat(s) = NF_DEF_VAR(fileid, 'initial_topography', NF_FLOAT, 1, DimId_n2D,    VarId_topo_ini) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'topography', NF_FLOAT, 1, DimId_n2D,    VarId_topo) 
    s = s + 1

!----- F I E L D S 


    stat(s) = NF_DEF_VAR(fileid, 'max_wave_height',   NF_FLOAT, 1, DimId_n2D, VarId_mwh) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'max_abs_vel',       NF_FLOAT, 1, DimId_n2D, VarId_vel) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'max_flux',          NF_FLOAT, 1, DimId_n2D, VarId_flux) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'initial_uplift',    NF_FLOAT, 1, DimId_n2D, VarId_uplift) 

    !----- vector variables
    
    DO i = 1,  s
       IF (stat(i) /= NF_NOERR) then
          PRINT *, 'nc_inundation, NetCDF error',stat(i),' in variable definition, no.', i
          print *, nf_strerror(stat(i))
       end IF
    END DO

!----- DEFINE ATTRIBUTES -----------------------------------------

    s = 1

    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_loc2D,    'long_name',   19, 'surface_coordinates') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_loc2D,    'units',       18, 'degrees east/north') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_lon,      'long_name',    9, 'longitude') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_lon,      'units',       12, 'degrees_east') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_lat,      'long_name',    8, 'latitude') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_lat,      'units',       13, 'degrees_north') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tri_area, 'long_name',   13, 'triangle area') 
    s = s + 1
    
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo_ini, 'long_name',   18, 'initial_topography') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo_ini, 'field',       18, 'initial_topography') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo_ini, 'connections', 20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo_ini, 'positions',   17, 'surface_locations') 
    s = s + 1
    
    attstr  = 'meters (before shock bathymetry and topography)'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo_ini, 'units',       LEN_TRIM(attstr), attstr) 
    s = s + 1
    
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'long_name',   10, 'topography') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'field',       10, 'topography') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'connections', 20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'positions',   17, 'surface_locations') 
    s = s + 1
    
    attstr  = 'meters (after shock bathymetry and topography)'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'units',       LEN_TRIM(attstr), attstr) 
    s = s + 1
    
    attstr='initial uplift due to earthquake'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uplift,   'long_name',   LEN_TRIM(attstr), attstr) 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uplift,   'field',       14, 'initial_uplift') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uplift,   'connections', 20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uplift,   'positions',   17, 'surface_locations') 
    s = s + 1
    
    attstr  = 'maximum wave height\n (above after shock topography on land)'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'long_name',   LEN_TRIM(attstr), attstr) 
    s = s + 1
    
    attstr  = 'meters'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'units',       LEN_TRIM(attstr), attstr) 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'field',       15, 'max_wave_height') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'connections', 20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'positions',   17, 'surface_locations') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'long_name',   36, 'maximum nodal velocity approximation') 
    s = s + 1
    
    attstr  = 'm/sec'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'units',       LEN_TRIM(attstr), attstr) 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'field',       11, 'max_abs_vel'); s=s+1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'connections', 20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'positions',   17, 'surface_locations') 
    s = s + 1

    attstr  = 'maximum flux in nodes\n (maximum value of the product of velocity and wave height)'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,      'long_name',   LEN_TRIM(attstr), attstr) 
    s = s + 1
    
    attstr  = 'm^2/sec'
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,     'units',        LEN_TRIM(attstr), attstr) 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,     'field',         8, 'max_flux') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,     'connections',  20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,     'positions',    17, 'surface_locations') 

    DO i = 1, s
       IF (stat(i) /= NF_NOERR) then
          PRINT *, 'nc_inundation, NetCDF error',stat(i),' in attribute assignments, no.', i
          print *, nf_strerror(stat(i))
       endif
    END DO

    s = 1

    stat(s) = NF_ENDDEF(fileid) 
    

!----- write variables

    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_loc2D,   AUX_loc(1:2, 1:inun_nod2D)) 
    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_lon,     AUX_loc(1, 1:inun_nod2D) ) 
    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_lat,     AUX_loc(2, 1:inun_nod2D) ) 
    s = s + 1

!----- triangles

    stat(s) = NF_PUT_VAR_INT(fileid,  VarId_trian,    inun_elem_nod(1:3, 1:inun_elem2D) - 1) 
    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_tri_area, REAL(voltriangle(inun_elem(1:inun_elem2D)), sp)) 
    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_topo,     REAL(nodhn(inun_nod(1:inun_nod2D)), sp)) 
    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_topo,     REAL(nodhn_init(inun_nod(1:inun_nod2D)), sp)) 
    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_uplift,   REAL(ssh_init(inun_nod(1:inun_nod2D)), sp))  
    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_mwh,      mwh(inun_nod(1:inun_nod2D)) ) 

    DO i = 1, s
       IF (stat(i) /= NF_NOERR) then
          PRINT *,'nc_inundation, NetCDF error',stat(i),' in writing variables, no.', i
          print *, nf_strerror(stat(i))
       endif
    END DO

    stat(1) = NF_CLOSE(fileid)

    IF (stat(1) /= NF_NOERR) THEN
        PRINT *, 'nc_inundation, NetCDF error',stat(1),' in closing NetCDF file'
        print *, nf_strerror(stat(1))
    END IF

  END SUBROUTINE nc_inundation
#endif
