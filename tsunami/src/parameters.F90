! Part of TsunAWI
! Copyright (c) 2014
! Alfred Wegener Institut,
! Helmholtz Centre for Polar and Marine Research
! Contact: rakowsky@awi.de


! MODULE: PARAMETERS
! contains all the parameters and variables entering the 
! shallow water equations and handles the namelist.


MODULE PARAMETERS

  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int16   

    IMPLICIT NONE
    SAVE
 
!-------------------------------------------------------------
    !CHARACTER(LEN=5), parameter :: revision='1.9.2'
    CHARACTER(LEN=5), parameter :: revision='r2.0.1'
!------------------------------------------------------------- 

    integer, parameter        :: dp=real64, sp=real32

!  If you want to toggle between single and double precision, change the working
!  precision:

#ifdef WORKING_PRECISION
    integer, parameter        :: wp=WORKING_PRECISION  ! Working precision
#else    
    integer, parameter        :: wp=sp  ! Working precision 
#endif    
    REAL(kind=dp)             :: time, time0

    REAL                      :: mw
    CHARACTER(LEN=27)         :: TSID

    REAL(kind=wp), PARAMETER  :: pi = 3.141592653589793_wp, rad = pi / 180.0_wp
    REAL(kind=wp), PARAMETER  :: omega = 2.0_wp * pi /(24.0_wp * 60.0_wp * 60.0_wp)
    REAL(kind=wp)             :: g = 9.806_wp                          ![m/s^2]
    REAL(kind=wp), PARAMETER  :: r_earth = 6.3675d6                   ![m]
    REAL(kind=wp), PARAMETER  :: density_0 = 1027.7_wp                 ![kg/m^3]
    INTEGER                  :: viscosity_type
    REAL(kind=wp)             :: smag_fact
    REAL(kind=wp)             :: Ah0                           ![m^2/s]
    REAL(kind=wp)             :: alpha, Cd, Dcr
    REAL(kind=WP)             :: min_layer_thickness
    REAL(kind=wp)             :: small = 0.0001_wp

    LOGICAL                   :: cli_eq_params, cli_quakeml
    REAL(kind=wp)             :: eq_mag, eq_epi_x, eq_epi_y
    real(kind=wp)             :: mng_min, mng_max, mng_general

    REAL(kind=dp)             :: dt, T_end, T_chunk, T_out, T_out_tidegauge
    REAL(kind=wp)             :: uvmax, sshmax, sshmin, dmin
    REAL(kind=wp)             :: scaling ! transfer to local cartesian


    LOGICAL                  :: enable_ruptgen_scenario, enable_okada_scenario, enable_aifdr_scenario
    LOGICAL                  :: enable_quakeml
    LOGICAL                  :: enable_read_initial_field, enable_benchmark, enable_idealised
    LOGICAL                  :: enable_dhdt
    INTEGER                  :: benchmark_ident


    LOGICAL                  :: read_restart, write_restart
    INTEGER                  :: snap_int, snap_int_tidegauge, snap_int_raster
  
    INTEGER                  :: momadv_type
    INTEGER                  :: wetting_drying_type
    INTEGER                  :: nmb_iter_smooth_inicond, nmb_iter_smooth_topo

    ! ----------------
    ! verbosity:   0: almost nothing
    !              1: current timestep
    !              2: add sshmin, sshmax, uvmax 
    !              4: 
    !              5: add timer
    !              6: 
    !              7: 
    !              8: 
    !              9: 
    ! ----------------

    INTEGER :: verbosity

        ! ----------------
        ! coordinate types:   1: geographical coordinates(lon/lat) [deg E/deg N]
        !                     2: cartesian coordinates [m]
        ! ----------------
    INTEGER :: coordinate_type 

        ! ----------------  
        ! rotation types:     0: no rotation
        !                     1: constant Coriolis force
        !                     2: latitude dependent Coriolis force (default)
        ! ----------------
    INTEGER :: rotation_type


    CHARACTER(len=100)         :: MeshPath     
    CHARACTER(len=100)         :: TopoFile
    CHARACTER(len=100)         :: OutputPath
    CHARACTER(len=100)         :: output_prefix
    CHARACTER(len=100)         :: IniCondPath
    CHARACTER(len=100)         :: TideGaugeFile
    CHARACTER(len=100)         :: initial_shape 
    CHARACTER(len=100)         :: QuakemlFile
    CHARACTER(len=100)         :: BottomFrictionFile
    
    ! depending on scenario type, initial fields 
    INTEGER                    :: MeshID, Epi_i, Epi_j

    CHARACTER(len=100)         :: okada_parameter_file, okada_fault_list_file
    CHARACTER(len=100)         :: init_filename


    LOGICAL                    :: enable_snapshots, enable_nc_out, enable_bin_out
    LOGICAL                    :: write_diag_in_nc_snapshots
    LOGICAL                    :: enable_ascii_out, enable_raster_out
    REAL                       :: smooth_vel_time
    LOGICAL                    :: write_groundtrack_data  
    LOGICAL                    :: write_ascii_neighbourhood
    LOGICAL                    :: write_velocity, write_ttt, write_cfl
    LOGICAL                    :: write_energy
    LOGICAL                    :: write_tidegauge_data
    LOGICAL                    :: write_final_checksums, write_nc_inundation
    REAL(kind=wp)              :: nc_snapshot_accuracy_ssh

    logical                    :: enable_variable_bottom_friction
    LOGICAL                    :: use_bv_roughness


    ! Declaration of specific variables
    REAL(kind=wp)              :: H_mn, H_crit  !(for GETM)
    
    ! For ascii output of mareograms at tide gauge locations
    INTEGER                         :: n_tidegauges
    CHARACTER(len=100), allocatable :: tidegauge_ascii_filename(:)
    INTEGER, allocatable            :: tidegauge_node_index(:)

    ! parameters for aifdr/GA scenarios
    character(len=100)          :: source_zone
    character(len=2)            :: sz

CONTAINS

  subroutine read_commandline

! The command line is read first. If it does not fit to namelist input, it should win.

    integer            :: num_args, count
    character(len=20)  :: argument, val
    integer            :: stat, pointpos
    LOGICAL            :: have_cli_mag, have_cli_epi_x, have_cli_epi_y
    LOGICAL            :: have_cli_tsid
    
    have_cli_mag     = .false.
    have_cli_epi_x   = .false. 
    have_cli_epi_y   = .false.
    have_cli_tsid    = .false.
    cli_quakeml      = .false.
    enable_quakeml   = .false.

    num_args = COMMAND_ARGUMENT_COUNT()

    if (num_args==0) return  ! No arguments, nothing to be done here
    
    count = 0

    do while (count < num_args)
    
       count = count+1
       call GET_COMMAND_ARGUMENT(count, argument,status=stat)

       if (stat==0) then

          select case(trim(argument))

          case ('--epi_x','--epi_lon','--lon','--x','-epi_x','-epi_lon','-lon','-x')
             have_cli_epi_x = .true.
             count = count+1
             call GET_COMMAND_ARGUMENT(count, val, status=stat)
             if (stat==0) read(val,*,iostat=stat) eq_epi_x

          case ('--epi_y','--epi_lat','--lat','--y','-epi_y','-epi_lat','-lat','-y')
             have_cli_epi_y = .true.
             count = count+1
             call GET_COMMAND_ARGUMENT(count, val, status=stat)
             if (stat==0) read(val,*,iostat=stat) eq_epi_y

          case ('--xy','--lonlat','--epi_xy','--epi_lonlat', &
                 '-xy', '-lonlat', '-epi_xy', '-epi_lonlat')
             
             have_cli_epi_x = .true.
             have_cli_epi_y = .true.

             count = count+1
             call GET_COMMAND_ARGUMENT(count, val, status=stat)

             if (stat==0) then
                read(val,*,iostat=stat) eq_epi_x

                if (stat==0) then
                   count = count+1
                   call GET_COMMAND_ARGUMENT(count, val, status=stat)
                      
                   if (stat==0) read(val,*,iostat=stat) eq_epi_y
                endif
             endif

          case ('--mw','--mag','--magnitude','--m', &
                 '-mw', '-mag', '-magnitude', '-m')
             have_cli_mag = .true.
             count = count+1
             call GET_COMMAND_ARGUMENT(count, val,status=stat)

             if (stat==0) read(val,*,iostat=stat) eq_mag

          case ('--id','--ID','--event','--tsid','--TSID', &
                 '-id', '-ID', '-event', '-tsid', '-TSID')
             have_cli_tsid = .true.
             count = count+1
             call GET_COMMAND_ARGUMENT(count, val,status=stat)

             TSID = trim(val)

             if (len(val) > len(TSID)) print *,"Trunkating scenario ID to",TSID

          case ('--quakeml','-quakeml')

#ifdef NO_FOXY
             print *,"!!! Compiled with -DNO_FOXY - cannot parse quakeML !!!"
             print *,"!!! Specify tsunami in a different way, or rebuild TsunAWI w/o -DNO_FOXY !!!"
             print *,"!!! Stopping here !!!"
             STOP 1
#else
             cli_quakeml = .true.
             count = count+1
             call GET_COMMAND_ARGUMENT(count, val,status=stat)

             QuakemlFile = trim(val)
#endif             
          case default 

             print *,'Error - unknown command line argument:',trim(argument)
             stat = 1

          end select
       endif

       if (stat /=0 )  then
          print *,"!!! Problem parsing the command line options - stopping here !!! "
          STOP 1
       endif

    end do

    !=======================================================
    ! Some checks - far from complete - on reasonable input
    !=======================================================

    if (have_cli_epi_x .and. (eq_epi_x < -180. .or. eq_epi_x > 180.)) then
       print *, "!!! Command line argument EQ epicenter logitude out of range - stopping here !!!"
       STOP 1
    endif

    if (have_cli_epi_y .and. (eq_epi_y < -90. .or. eq_epi_y > 90.))  then
       print *, "!!! Command line argument EQ epicenter latitude out of range - stopping here !!!"
       STOP 1
    endif
                         
    if (have_cli_mag .and. (eq_mag < 0. .or. eq_mag > 10.))  then
       print *,"!!! Command line argument EQ magnitude out of range - stopping here !!!"
       STOP 1
    endif

    

    if (cli_quakeml) then

       ! Final check, then prepare QuakeML input
       
       if (have_cli_mag .or. have_cli_epi_x .or. have_cli_epi_y) then
          print *,"!!! Command line: specify either QuakeML file _or_ EQ parameters !!!"
          print *,"!!! - stopping here                                              !!!"
          STOP 1
       endif       
       
       ! If we don't get a specific ID, we'll take the name of the QuakeML file minus the suffix
       if (.not. have_cli_tsid) then
          pointpos = scan(trim(QuakemlFile),".", BACK=.true.)
          if (pointpos > 1) then
             TSID = QuakemlFile(1:pointpos-1)
          else
             TSID = QuakemlFile
          endif
       endif

       enable_quakeml = .true.

    else

       ! Final check
       
       cli_eq_params = (have_cli_mag .and. have_cli_epi_x .and. have_cli_epi_y) 
 
       if (cli_eq_params .neqv. (have_cli_mag .or. have_cli_epi_x .or. have_cli_epi_y)) then
          print *, "!!! Command line arguments not complete:                        !!!"
          print *, "!!! If EQ parameters are specified here, then                   !!!"
          print *, "!!! all of magnitude and epicenter lon, lat must be given.      !!!" 
          print *, "!!! Stopping here !!!"
          STOP 1
       end if
       
       if (cli_eq_params) then
          if (.not. have_cli_tsid) call derive_TSID
       endif

    endif
    
  end subroutine read_commandline


!================================================================================================

  subroutine read_namelist

    implicit none

    INTEGER                     :: RuptGen_ver, ScenID
    CHARACTER(len=3)            :: StrEpi_i
    CHARACTER(len=2)            :: StrEpi_j
    CHARACTER(len=4)            :: StrMW
    CHARACTER(len=3)            :: StrMeshID
    REAL(kind=wp), DIMENSION(2) :: aux
    integer                     :: n, stat

    NAMELIST /general/ MeshPath, IniCondPath, OutputPath, TopoFile, TideGaugeFile, output_prefix,  &
        coordinate_type, rotation_type, momadv_type, wetting_drying_type, &
        use_bv_roughness,  &
        nmb_iter_smooth_inicond, nmb_iter_smooth_topo, verbosity

    NAMELIST /para/ dt, T_end, T_chunk, T_out, T_out_tidegauge, smooth_vel_time, &
         viscosity_type, smag_fact, Ah0, alpha, Cd, min_layer_thickness, &
         Dcr, H_mn, H_crit

    NAMELIST /model_init/ enable_ruptgen_scenario, enable_okada_scenario, enable_read_initial_field, &
         enable_benchmark, benchmark_ident, enable_idealised, enable_dhdt, initial_shape, enable_aifdr_scenario

    NAMELIST /model_output/ enable_snapshots, enable_nc_out, write_diag_in_nc_snapshots,  &
        enable_bin_out, enable_ascii_out, enable_raster_out,  &
        write_restart, write_groundtrack_data, write_ascii_neighbourhood, write_velocity, write_cfl, write_ttt, write_energy, &
        write_tidegauge_data, write_final_checksums, nc_snapshot_accuracy_ssh, write_nc_inundation

    NAMELIST /model_init_okada_details/ okada_parameter_file, okada_fault_list_file 

    NAMELIST /model_init_ruptgen_details/ RuptGen_ver, MeshID, Epi_i, Epi_j, mw

    NAMELIST /model_init_aifdr_details/ MeshID, source_zone, sz, Epi_i, Epi_j, mw

    NAMELIST /model_init_from_file/ init_filename 

    NAMELIST /bottom_friction/ enable_variable_bottom_friction, BottomFrictionFile, &
                               mng_min, mng_max
    
    integer            :: infile_nml
    integer            :: infile_gfz

! Set namelist parameters to "invalid" 

!    NAMELIST /general/ 
    MeshPath = ""
    IniCondPath = "" 
    OutputPath = "" 
    TopoFile = "" 
    TideGaugeFile= ""
    output_prefix = ""
    coordinate_type = -1 
    rotation_type = -1
    momadv_type = -1
    wetting_drying_type = -1
    use_bv_roughness = .false. 
    nmb_iter_smooth_inicond = -1
    nmb_iter_smooth_topo = -1
    verbosity = -1

!    NAMELIST /para/ 
    dt = -1. 
    T_end = -1. 
    T_chunk = -1. 
    T_out = -1.
    T_out_tidegauge = -1.
    smooth_vel_time = -1.
    viscosity_type = -1 
    smag_fact = -1.
    Ah0 = -1. 
    alpha = -1. 
    Cd = -1.
    min_layer_thickness = -1.
    Dcr = -1.
    H_mn = -1.
    H_crit = -1.
    mng_min = -1. 
    mng_max = -1.

!    NAMELIST /model_init/ 
    enable_ruptgen_scenario =.false. 
    enable_okada_scenario = .false. 
    enable_read_initial_field = .false.
    enable_benchmark = .false.
    benchmark_ident = -1 
    enable_idealised = .false.
    enable_dhdt = .false.
    initial_shape=""
    enable_aifdr_scenario=.false.

!    NAMELIST /model_output/ 
    enable_snapshots = .false.
    enable_nc_out = .false.
    write_diag_in_nc_snapshots = .false. 
    enable_bin_out = .false.
    enable_ascii_out = .false.
    write_restart = .false.
    write_groundtrack_data = .false.
    write_ascii_neighbourhood = .false.
    write_velocity = .false.
    write_cfl = .false.
    write_ttt = .false.
    write_energy = .false.
    write_tidegauge_data = .false.
    write_final_checksums = .false.
    nc_snapshot_accuracy_ssh = -1.
    write_nc_inundation = .false.
    
!    NAMELIST /model_init_okada_details/ 
    okada_parameter_file= ""
    okada_fault_list_file=""


!    NAMELIST /model_init_from_file/ 
    init_filename= ""

!    NAMELIST /bottom_friction/
    enable_variable_bottom_friction = .false.
    BottomFrictionFile = ""
    mng_min = -1.
    mng_max = -1.

    write (*,*) 'Reading  the namelist'
    
    ! Now, read the general parts of the namelist. Open and close inbetween, such that
    ! the order of the blocks does not matter. 
    OPEN(newunit=infile_nml, file = 'namelist.tsunami', iostat=stat)
    READ(infile_nml, NML  = general, iostat=stat)
    CLOSE(infile_nml)
    OPEN(newunit=infile_nml, file = 'namelist.tsunami', iostat=stat)
    READ(infile_nml, NML  = model_output, iostat=stat)
    CLOSE(infile_nml)
    OPEN(newunit=infile_nml, file = 'namelist.tsunami', iostat=stat)
    READ(infile_nml, NML  = para, iostat=stat)
    CLOSE(infile_nml)
   

! Is namelist for variable bottom friction present?
! (Has to be checked first, otherwise, Cd is set to a default)
    OPEN(newunit=infile_nml, file = 'namelist.tsunami', iostat=stat)
    READ(infile_nml, NML = bottom_friction, iostat=stat)
    CLOSE(infile_nml)
    if (stat == 0) then
       if (enable_variable_bottom_friction) then
          if (LEN_TRIM(BottomFrictionFile) == 0) then
             write (*,*) '/bottom_friction/ Setting missing value: BottomFrictionFile="mng2d.out"'
             BottomFrictionFile="mng2d.out"
          endif
          if (mng_min < 0.) then
             write (*,*) '/bottom_friction/ mng_min must be >= 0. Setting mng_min=0.'
             mng_min=0.
          endif
       end if
    end if


    CLOSE(infile_nml)

! and check for unset values. Set those variables to a default. 
    
    if (LEN_TRIM(MeshPath) == 0) then
       write (*,*) '/general/ Setting missing value: MeshPath="./"'
       MeshPath = "./"
    endif

    if ( LEN_TRIM(IniCondPath) == 0) then
       print *, '/general/ Setting missing value: IniCondPath="./"'
       IniCondPath = "./"
    end if
    if (LEN_TRIM(OutputPath) == 0) then
       print *, '/general/ Setting missing value: OutputPath="./"' 
       OutputPath = "./"
    endif

    if (LEN_TRIM(TopoFile) == 0) then
       print *, '/general/ Setting missing value: TopoFile="nodhn.out"'  
       TopoFile="nodhn.out"
    end if
    
    if (LEN_TRIM(output_prefix) == 0) then
       print *, '/general/ Empty output_prefix="" - this is ok'
    endif
    
    if (coordinate_type == -1) then
       print *, '/general/ Setting missing value: coordinate_type=1 (geographical coordinates: degE, degN)' 
       coordinate_type = 1
    endif

    if (rotation_type == -1) then
       if (coordinate_type == 1) then
          print *, '/general/ Setting missing value: rotation_type=2 (latitude dependent Coriolis force)'
          rotation_type = 2
       else
          print *, '/general/ Setting missing value: rotation_type=0 (no Coriolis force)'
          rotation_type = 0
       endif
    end if

    if (momadv_type == -1) then
       print *, '/general/ Setting missing value: momadv_type=1 (non-linear advection, full P1 scheme)'
       momadv_type = 1
    end if

    if (wetting_drying_type == -1) then
       print *, '/general/ Setting missing value: wetting_drying_type=1 (extrapolation)'
       wetting_drying_type = 1
    endif


    if (nmb_iter_smooth_inicond == -1) then
       print *, '/general/ Setting missing value: nmb_iter_smooth_inicond=0  (do not smooth initial elevation)'
       nmb_iter_smooth_inicond=0
    end if

    if (nmb_iter_smooth_topo == -1) then
       print *, '/general/ Setting missing value: nmb_iter_smooth_topo=2  (smooth topography)'
       nmb_iter_smooth_topo = 2
    end if

    if (verbosity == -1) then
       print *, '/general/ Setting missing value: verbosity=1'
       verbosity = 1
    endif

!    NAMELIST /model_output/

#ifdef NO_NETCDF
    IF (enable_nc_out) THEN 
      print *,'/model_output/ enable_nc_out is ignored, because the code was built with -DNO_NETCDF.'
      print *,'               Recompile without -DNO_NETCDF to enable netcdf output.'
   endif
#else
   
    if (nc_snapshot_accuracy_ssh <= 0.) then
       print *,'/model_output/ nc_snapshot_accuracy_ssh unset, zero, or negative: writing SSH snapshots without rounding'
    else
       print *,'Setting SSH snapshot accuracy to:',nc_snapshot_accuracy_ssh
    endif   
#endif    

!    NAMELIST /para/ 
    if (dt < 0. ) then
       print *, '/para/ Setting missing value: dt=1. (timestep 1s - check with grid resolution for CFL!'
       dt = 1.
    end if

    if (T_end < 0.) then
       print *, '/para/ Setting missing value: T_end=3600. (1h modeltime)' 
       T_end = 3600.
    end if

    if ( write_restart .and. T_chunk < 0.) then
       print *, '/para/ Setting missing value: T_chunk=T_end (write restart only at the end)' 
       T_chunk = T_end
    end if

    if (enable_snapshots .and. T_out < 0.) then
       print *, '/para/ Setting missing value: T_out=T_end (write snapshot only at end)'
       T_out = T_end
    end if

    if (write_tidegauge_data) then
       if (T_out_tidegauge < 0.) then
          print *, '/para/ Setting missing value: T_out_tidegauge=T_out (write TG data with snapshot)'
          T_out_tidegauge=T_out
       end if
       if (LEN_TRIM(TideGaugeFile) == 0) then
          print *, '/general/ Setting missing value: TideGaugeFile="./TG_list.txt" (File with TG locations)'
          TideGaugeFile="./TG_list.txt"
       endif
    end if


    if (T_end>12.*3600. .and. smooth_vel_time < 0.) then
       print *, '/para/ Setting missing value: smooth_vel_time=43200. (for long runs: smooth velocity field every 12h)'
       smooth_vel_time = 12.*3600.
    end if

    if (viscosity_type == -1) then
       print *, '/para/ Setting missing value: viscosity_type=1 (linear viscosity)' 
       viscosity_type = 1
    end if

    if (viscosity_type ==2 .and. smag_fact < 0.) then
       print *, '/para/ Setting missing value:  smag_fact=0.3 (factor for Smagorinski viscosity '
       smag_fact = .3
    end if

    if (viscosity_type ==1 .and. Ah0 < 0.) then
       print *, '/para/ Setting missing value: Ah0=.05 (factor for linear viscosity)' 
       Ah0 = .05
    end if

    if (alpha < 0.) then
       print *, '/para/ Setting missing value: alpha=0.05 (Robert-Asselin-Filter for leap-frog timestep)' 
       alpha = 0.05
    end if

    if ((.not. enable_variable_bottom_friction) .and. Cd < 0. ) then
       print *, '/para/ Setting missing value: Cd=0.02 (Manning parameter)'
       Cd = 0.02
    end if

    if (min_layer_thickness < 0.) then
       print *,'/para/ Setting missing value:  min_layer_thickness=', small
    else
       small =  min_layer_thickness
    endif
    if (use_bv_roughness .and. Dcr < 0.) then
       print *, '/para/ Setting missing value: Dcr=0.1 (critical depth for use of Brunt-Vaisala frequency)'
       Dcr = 0.1
    end if

    if (wetting_drying_type==2) then
       if (H_mn < 0.) then
          print *, '/para/ Setting missing value: H_mn=0.02 (min depth for GETM wetting+drying)'
          H_mn=0.02
       end if
       if (H_crit < 0.) then
          print *, '/para/ Setting missing value: H_crit=0.1 (critical depth for GETM wetting+drying) '
          H_crit=0.1
       end if
    endif


!    NAMELIST /model_init/  

    OPEN(newunit=infile_nml, file = 'namelist.tsunami', iostat=stat)
    READ(infile_nml, NML = model_init, iostat=stat)
    CLOSE(infile_nml)

    if (cli_eq_params) then
! /model_init/ is mutually exclusive with specifying parameters via command line or QuakeML
       if (stat==0) then
          print *,"/model_init/ This part will be ignored, as EQ parameters are specified on the command line."
       endif
       enable_idealised = .true.

    elseif (cli_quakeml) then
       if (stat==0) then
          print *,"/model_init/ This part will be ignored, as EQ parameters are specified in QuakeML file"
       endif
    else

       if (.not. (enable_ruptgen_scenario .or. enable_okada_scenario .or. enable_read_initial_field .or.  &
            enable_benchmark .or. enable_idealised .or. enable_aifdr_scenario)) then
          print *,'/model_init/No scenario type given. We take: enable_benchmark=.true.,  benchmark_ident=1'
          enable_benchmark = .true.
          benchmark_ident = 1
          
       elseif (enable_benchmark) then
          if (enable_ruptgen_scenario .or. enable_okada_scenario .or. enable_idealised) then
             print *,'/model_init/ Too many scenario types given. We take: enable_benchmark=.true.'
             enable_ruptgen_scenario = .false.
             enable_okada_scenario = .false.
             enable_idealised = .false.
          endif
          if (benchmark_ident == -1) then
             print *,'/model_init/ Setting missing value benchmark_ident=1 (sloping beach)'
             benchmark_ident = 1
          elseif (benchmark_ident==4 .or. benchmark_ident > 5 ) THEN
             print *, "/model_init/ Specified Benchmark is not implemented!"
             STOP 1
          end if
       END IF



       !========================================================================
       ! read specific parts of the namelist and check for unset values

       if (enable_ruptgen_scenario ) then 
          
          !    NAMELIST /model_init_ruptgen_details/ 
          RuptGen_ver = -1
          MeshID = -1 
          Epi_i = -1 
          Epi_j = -1
          mw = -1.
          
          OPEN(newunit=infile_nml, file = 'namelist.tsunami', iostat=stat)
          read(infile_nml, nml  = model_init_ruptgen_details, iostat=stat)
          CLOSE(infile_nml)
          
          if (MeshID == -1) then
             print *, '/model_init_ruptgen_details/ Setting  missing value: MeshId=0 (for output filename)'
             MeshID = 0
          end if
          
          if (Epi_i == -1 .or. Epi_j == -1 .or. mw < 0.) then
             print *, "/model_init_ruptgen_details/ Sorry, Epi_i, Epi_j, and mw are needed for RuptGen scenarios! Stopping here."
             STOP 1
          end if
          
          write(StrEpi_i,  '(I3.3)') Epi_i
          write(StrEpi_j,  '(I2.2)') Epi_j
          write(StrMW,     '(F4.2)') mw
          write(StrMeshID, '(I3.3)') MeshID
          
          TSID = 'SZ_r'//revision//'_m'//StrMeshID//'_'//StrEpi_i//'_'//StrEpi_j//'_mw'//StrMW//'_'
          ! for example SZ_r018A_m001_000_00_mw8.5
          
          if (enable_read_initial_field) then
             OPEN(newunit=infile_nml, file = 'namelist.tsunami', iostat=stat)
             read(infile_nml, NML = model_init_from_file, iostat=stat)
             CLOSE(infile_nml)
          end if
          
          if ( LEN_TRIM(init_filename) == 0 ) then
             init_filename = 'ssh0_mw'//StrMW//'_'//StrEpi_i//'_'//StrEpi_j//'.out'
          end if
          
          IF (RuptGen_ver==1) THEN
             ScenID = (Epi_i - 1) * 15 + Epi_j
             open(newunit=infile_gfz, file = 'GFZ_epicenters.txt')
          else
             if (RuptGen_ver==-1) then
                print *,'/model_init_ruptgen_details/ Setting missing value: RuptGen_ver=2 (RuptGen 2.0 or 2.1)'
             end if
             ScenID = (Epi_i - 1) * 25 + Epi_j
             OPEN(newunit=infile_gfz, file = 'GFZ_epicenters-2.0.txt')
          END IF
          do n = 1, ScenID
             read(infile_gfz, *) aux(:)
          end DO
          close(infile_gfz)
        
          eq_mag   = mw
          eq_epi_x = aux(1)
          eq_epi_y = aux(2)
          

       else if (enable_aifdr_scenario) then


          !      NAMELIST /model_init_aifdr_scenario/
          MeshID=-1
          source_zone=''
          sz=''
          mw=-1.
          epi_i=-1
          epi_j=-1
          
          OPEN(newunit=infile_nml, file = 'namelist.tsunami', iostat=stat)
          read(infile_nml, nml  = model_init_aifdr_details, iostat=stat)
          CLOSE(infile_nml)
          
          if (Epi_i == -1 .or. Epi_j == -1 .or. mw < 0.) then
             print *, "/model_init_aifdr_details/ Sorry, Epi_i, Epi_j, and mw are needed for aifdr  scenarios! Stopping here."
             STOP 1
          end if
          write(StrEpi_i,  '(I3.3)') Epi_i
          write(StrEpi_j,  '(I2.2)') Epi_j
          write(StrMW,     '(F4.2)') mw
          write(StrMeshID, '(I3.3)') MeshID
          
          TSID = sz//'_'//revision//'_m'//StrMeshID//'_'//StrEpi_i//'_'//StrEpi_j//'_mw'//StrMW//'_'
          ! for example SZ_r018A_m001_000_00_mw8.5
          
          if (enable_read_initial_field) then
             
             OPEN(newunit=infile_nml, file = 'namelist.tsunami', iostat=stat)
             read(infile_nml, NML = model_init_from_file, iostat=stat)
             CLOSE(infile_nml)
             
          end if
          !if ( LEN_TRIM(init_filename) == 0 ) then
          !  init_filename = 'ssh0_mw'//StrMW//'_'//StrEpi_i//'_'//StrEpi_j//'.out'
          !end if
          
          
       else if ( enable_okada_scenario ) then
          TSID='OK_'
          OPEN(newunit=infile_nml, file = 'namelist.tsunami', iostat=stat)
          read(infile_nml, NML  = model_init_okada_details, iostat=stat)
          CLOSE(infile_nml)
          
       else if ( enable_benchmark ) then
          TSID='BM_'
       else
          TSID='TSW_'
       end if
       
       if ((.not. enable_ruptgen_scenario) .and. (.not. enable_aifdr_scenario) .and.  enable_read_initial_field) then
          
          OPEN(newunit=infile_nml, file = 'namelist.tsunami', iostat=stat)
          read(infile_nml, NML = model_init_from_file, iostat=stat)
          CLOSE(infile_nml)
          if ( LEN_TRIM(init_filename) == 0 ) then
             print *, " /model_init_from_file/ enable_read_initial_field=.true., but init_filename is empty!"
             STOP 1
          endif
       end If
    endif

    close(infile_nml)

  end subroutine read_namelist

  subroutine derive_TSID

    ! If no scenario ID is given or can be derived from initial conditions,
    ! we compute one here from the current date and time, as it should be unique 

    integer                 :: datetime(8)
    integer                 :: date_seconds, pos
 
    character(len=62), parameter :: letter="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

    call DATE_AND_TIME(values=datetime)
    
    date_seconds = mod(datetime(1),10) * 366 * 24 * 3600   & ! years in decade     
                          + datetime(2) * 31 * 24 * 3600   & ! month
                          + datetime(3)      * 24 * 3600   & ! day
                          + datetime(5)           * 3600   & ! hour
                          + datetime(6)             * 60   & ! minute
                          + datetime(7)                      ! second       

    pos = mod(date_seconds,62)+1
    TSID = letter(pos:pos)

    do while (date_seconds > 62)
       date_seconds = date_seconds/62
       if (date_seconds < 52) exit
       pos = mod(date_seconds,62)+1
       TSID = letter(pos:pos)//trim(TSID)
    enddo
    pos = mod(date_seconds,52)+1
    TSID = 'TS_'//letter(pos:pos)//trim(TSID)//"_"
    

    print *,'No ID given - generating a tsunami scenario ID based on date+time.'
    print *,'Tsunami scenario ID is set to TSID=',TSID

  end subroutine derive_TSID

END MODULE PARAMETERS
!
!---------------------------------------------------------------------
