!############################################################################
!#                                                                          #
!# Main Program as part of the TsunAWI Code                                 #
!# AWI Tsunami research unit                                                #
!#--------------------------------------------------------------------------#
!#                                                                          #
!#                                                                          #
!# TTTTTTTT   SSSSSS  UU    UU  NN     NN      A      WW       WW  IIIIII   #
!# TTTTTTTT  SSSSSSS  UU    UU  NNN    NN     AAA     WW       WW  IIIIII   #
!#    TT     SS       UU    UU  NNNN   NN    AA AA    WW       WW    II     #
!#    TT     SSSSSS   UU    UU  NN NN  NN   AA   AA   WW       WW    II     #
!#    TT      SSSSSS  UU    UU  NN  NN NN  AAAAAAAAA  WW   W   WW    II     #
!#    TT          SS  UU    UU  NN   NNNN  AAAAAAAAA  WW WW WW WW    II     #
!#    TT     SSSSSSS  UUUUUUUU  NN    NNN  AA     AA  WWW     WWW  IIIIII   #
!#    TT     SSSSSS    UUUUUU   NN     NN  AA     AA  WW       WW  IIIIII   #
!#                                                                          #
!# TsunAWI - A Tsunami model on Unstructured meshes                         #
!#           provided by Alfred Wegener Institute,                          #
!#           Helmholtz Centre for Polar and Marine Research                 #
!#                                                                          #
!# Simulation code for the tsunami propagation and runup.                   #
!#                                                                          #
!# Version 1.9.1                                                            #
!#                                                                          #
!# Copyright (c) 2019                                                       #
!#               Alfred Wegener Institut,                                   #
!#               Helmholtz Centre for Polar and Marine Research             #
!#                                                                          #
!# AUTHORS: N. Rakowsky, A. Androsov, J.Behrens, S. Braune, A. Fuchs,       #
!#          S.Harig, W. Hiller, J. Schroeter, D. Sein, D. Sidorenko,        #
!#          O. Startseva, E.Taguchi                                         #
!#                                                                          #
!# AFFILIATION: Alfred Wegener Institute,                                   #
!#              Helmholtz Centre for Polar and Marine Research              #
!#              Am Handelshafen 12                                          #
!#              27570 Bremerhaven                                           #
!#              Germany                                                     #
!#                                                                          #
!# CONTACT: rakowsky@awi.de                                                 #
!#                                                                          #
!# DATE: 12 August 2019                                                     #
!#                                                                          #
!# Any distribution of this code must contain this complete header and      # 
!# the accompanying file TsunAWILicence.txt containing the complete licence #
!#                                                                          #
!#                                                                          # 
!############################################################################



   !----------------------------------------------------------------------------
   !
   !       ********************   M A I N   *****************************
   !----------------------------------------------------------------------------
PROGRAM MAIN

!$    USE OMP_LIB

    USE PARAMETERS
    USE MESH
    USE ELEMENTS
    USE BENCHMARK
    USE SWE
    USE INITIAL_CONDITIONS
    USE OKADA_FAULT_PARAMETERS
    USE IO_RASTER

    IMPLICIT NONE

    INTEGER               :: i_time, n
    INTEGER               :: cnt, el

    CHARACTER(len=100)    :: fname

    LOGICAL               :: normal_end, snapshot_step = .false.
    INTEGER               :: open_status
    INTEGER               :: smooth_vel_int
    REAL(kind=dp)         :: t0lp, t1lp, t0, t1, t2, t3, t4,t5, t0io
    REAL(kind=dp)         :: t_loop_begin, t_loop_end

    REAL(kind=wp)         :: checksum_ssh
    INTEGER               :: restartfile_nodhn
    
    ! Initialize variables

    normal_end=.TRUE.

    call read_commandline
    if (enable_quakeml) call read_quakeml
    call read_namelist

    !===========================================================
    !================ INITIALIZE MESH ==========================

    ! Set up the mesh and define basis functions
    CALL read_2Dmesh
    CALL standard_element_definition
    CALL mesh_scaling
    CALL basisfunctions

    CALL build_nghbr_arrays

    CALL build_edg_arrays
    CALL alloc_swe_fields

    !&&&&&&&&&&&&&&&

    ! tidal roughness, Brunt-Vaisila frequency
    if (use_bv_roughness) then
       call set_up_bv_roughness
    end if

    ! Manning roughness
    if (enable_variable_bottom_friction) then
       call set_up_mng_roughness
    end if

    call smooth_topography

    !================ END INITIALIZE MESH ======================
    !===========================================================


    !===========================================================
    !================== INITIALIZE FIELDS ======================


    if (enable_raster_out) call io_raster_init()

    ! Determine if there is already a restart file

    restart_file = trim(OutputPath)//'RESTART_'//trim(TSID)//trim(output_prefix)
    
    fname=TRIM(restart_file)//'_nodhn'
    OPEN(newunit=restartfile_nodhn,file=fname,STATUS="OLD",form='unformatted',IOSTAT=open_status)
    CLOSE(restartfile_nodhn)
    
    if (open_status==0) then
       read_restart=.TRUE.
       PRINT *,'LOADING RESTART FILES'

       IF (enable_okada_scenario) call okada_read_faultlist

       call bin_restart_in

    else
       ! Start "from scratch", no restart

       read_restart=.FALSE.
       PRINT *,'STARTING NEW MODEL RUN'

       IF (enable_okada_scenario) THEN
          call okada_read_faultlist
          call initial_conditions_param
       END IF
       
       if ( enable_idealised) then
          call start_point_init
       elseif ( enable_ruptgen_scenario .or. enable_read_initial_field ) then
          call start_point_init_file 
       end if
   
       ! Initial ssh for benchmarks, if not read from file.
       if (enable_benchmark) call benchmark_init
          

       sshmin = 1000._wp
       sshmax = -1000._wp
!$OMP PARALLEL DO  REDUCTION(MIN:sshmin) REDUCTION(MAX:sshmax)
       do n=1,nod2D
          ssh0(n) = max(ssh0(n),min(0._wp,-nodhn(n)))
          ssh1(n) = ssh0(n)
          ssh2(n) = ssh0(n)          

          sshmin = min(sshmin,ssh0(n))
          sshmax = max(sshmax,ssh0(n))
       end do
!$OMP END PARALLEL DO

       PRINT *,'INIT MESH, initial SSH: minimum=', sshmin,', maximum=',sshmax

       time0 = 0.
       time=time0

       ! compute estimated travel time?
       if (write_ttt) then
          call compute_ttt
          if (enable_raster_out) call io_raster_ttt 
       endif
        
#ifndef NO_NETCDF
        IF (enable_nc_out) THEN
          call nc_init
          call nc_out(0)
       END IF
#endif

       if (enable_raster_out) call io_raster_out ! write initial ssh

       IF (enable_benchmark) call benchmark_output

       if (write_ttt) deallocate(ttt)
       
    END IF


    !================== END INITIALIZE FIELDS ==================
    !===========================================================

    CALL build_index_bnd

    
    if (enable_snapshots) then
       snap_int = max(1,INT(T_out / dt))
       WRITE(*,*) 'snap_int=',snap_int
    endif
    if (write_tidegauge_data) then
       snap_int_tidegauge = max(1,INT(T_out_tidegauge / dt))

       call ascii_out_tidegauge_init
    endif

    
    smooth_vel_int = INT(smooth_vel_time/dt)

    i_time = 0

    save_vel_at_nodes=.FALSE.
    smooth_vel=.FALSE.

    PRINT *,'INITIAL TIME :',time,time0

!$OMP PARALLEL DO 
    do n=1,nod2D
       depth(n) = max(ssh2(n)+nodhn(n), 0._wp)
    end do
!$OMP END PARALLEL DO


    !======================================================
    !===================== MAIN LOOP ======================
    !======================================================
!$  t_loop_begin=omp_get_wtime()

    DO WHILE ( time < min(time0+T_chunk,T_end) )

!$      if (verbosity >= 1) t0lp=omp_get_wtime()

        i_time = i_time + 1
        if (enable_snapshots) snapshot_step = (MOD(i_time,snap_int) == 0) 

!$      if (verbosity >= 2) t0=omp_get_wtime()

!$OMP PARALLEL PRIVATE(save_vel_at_nodes)
         if (wetting_drying_type==1) then
           CALL compute_gradient_EXTRAPOLATE
        else
           CALL compute_gradient_GETM
        end if

        IF (enable_okada_scenario) THEN
!$OMP MASTER
           IF ( time <= ttl_rupt_time .and. &             
                     time > 0.0 .AND. ANY(ABS(rupt_delay-time) < 1.e-4) ) THEN
              CALL update_initial_ssh
           END IF
!$OMP END MASTER
!$OMP BARRIER
        ENDIF

!$OMP MASTER
!$      if (verbosity >= 2) t1=omp_get_wtime()
!$OMP END MASTER

        save_vel_at_nodes = (mod(i_time,smooth_vel_int)==smooth_vel_int-2) 

        CALL compute_velocity_at_nodes

        if (write_energy .and. (snapshot_step .or. verbosity > 0)) &
             call compute_energy
!$OMP END PARALLEL

!$      t2=omp_get_wtime()

        smooth_vel = (mod(i_time,smooth_vel_int)==0)

        if (wetting_drying_type==1) then
           CALL compute_velocity_EXTRAPOLATE
        else
           CALL compute_velocity_GETM
        end if

!$      if (verbosity >= 2) t3=omp_get_wtime()

        if (wetting_drying_type==1) then
           CALL compute_ssh_EXTRAPOLATE
        else
           CALL compute_ssh_GETM
        end if        

!$      if (verbosity >= 2) t4=omp_get_wtime()

        if (enable_benchmark) call benchmark_adjust_ssh

        IF (verbosity >= 1 .and. write_energy) &
           print *,'potential and kinetic Energy, Epot, Ekin=', Epot, Ekin
        
           
        IF (verbosity==1) THEN
           WRITE(*,'(I5,A3,f6.3)') i_time," t=",time/3600.
        ELSE IF (verbosity>=2) THEN
           checksum_ssh=0.

!$OMP PARALLEL DO REDUCTION(+:checksum_ssh)
           DO n = 1, nod2D
              checksum_ssh = checksum_ssh + ssh2(n)
           ENDDO
!$OMP END PARALLEL DO
           WRITE(*,'(I5,A3,f6.3,A8,E10.3,A8,E10.3,A7,E10.3,A6,E10.3)') &
                i_time," t=",time/3600.," smin=",sshmin," smax=",sshmax," uvmax=",uvmax," dmin=",dmin
           print *, " check=",checksum_ssh
        END IF ! verbosity
              
        IF ( uvmax > 900. ) THEN 
            WRITE(*,*) 'Explosion, uvmax=', uvmax
            normal_end=.FALSE.
#ifndef NO_NETCDF
            IF (enable_nc_out) CALL nc_out(i_time)
#endif
            exit 
        END IF
        time = time0 + dt * real(i_time,wp)

        ! binary, netCDF, and/or restart output
        IF ( snapshot_step) THEN
!$         t0io = omp_get_wtime()

#ifndef NO_NETCDF
           IF (enable_nc_out) CALL nc_out(i_time)
#endif
           if (enable_ascii_out) call ascii_out

           if (enable_benchmark) call benchmark_output

!$         if (verbosity>=5) then
!$            print *,'MAIN LOOP - output           : ',omp_get_wtime() - t0io
!$         endif
        END IF

        IF (write_tidegauge_data) then 
           if (mod(i_time,snap_int_tidegauge)==0) call ascii_out_tidegauge
        END IF
        
        IF (enable_snapshots .and. enable_raster_out) then
           if (mod(i_time,snap_int_raster)==0) call io_raster_out
        END IF
        
!$      if (verbosity == 1) then
!$         print *,'TOTAL TIME LOOP :',omp_get_wtime() - t0lp
!$      elseif (verbosity>=2) then

!$         print *,'MAIN LOOP - gradient         : ',t1-t0
!$         print *,'MAIN LOOP - velocity_in_nodes: ',t2-t1
!$         print *,'MAIN LOOP - velocity         : ',t3-t2
!$         print *,'MAIN LOOP - ssh, radiate     : ',t4-t3

!$         print *,'TOTAL TIME LOOP :',omp_get_wtime() - t0lp
!$
!$      end if

    END DO

!$  t_loop_end=omp_get_wtime()

!$  print *,'Time stepping took ',real(t_loop_end-t_loop_begin,sp),'seconds'
    
    !======================================================
    !================== END OF MAIN LOOP ===================
    !======================================================


    ! Prepare restart output
    if (time<T_end .AND. normal_end .AND. write_restart) then
       call bin_restart_out
    end if

#ifndef NO_NETCDF
    IF ( enable_nc_out ) then
       if (write_nc_inundation) call nc_inundation
       call nc_finalize
    endif
#endif

    if (enable_benchmark) call benchmark_finalize

    if (enable_ascii_out) CALL ascii_out_finalize
    
    if (write_final_checksums) call write_checksums

    if (enable_raster_out) call io_raster_finalize

    if (normal_end) then
       stop 0
    else
       stop 1
    end if

END PROGRAM MAIN


