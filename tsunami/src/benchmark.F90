! Part of TsunAWI
! Copyright (c) 2014
! Alfred Wegener Institut,
! Helmholtz Centre for Polar and Marine Research
! Contact: rakowsky@awi.de


! MODULE: BENCHMARK -- define benchmark specific code here
! 
!  benchmark_ident is set in the namelist (parameters.F90)
!
!  benchmark_ident = 1  ! Run up on a sloping beach, grid bm1
!                  = 2  ! Okushiri Channel experiment, grid bm2

MODULE BENCHMARK

  use parameters
  use mesh
  use swe
  implicit none
  save
  private
  public :: benchmark_init, benchmark_adjust_ssh, benchmark_finalize, &
            benchmark_output

  integer, parameter :: bm_channel_slopingbeach=1
  integer, parameter :: bm_channel_okushiri=2

  ! For benchmark_ident==2, bm_channel_okushiri
  REAL(kind=wp), allocatable, dimension(:) :: toku, boku
  INTEGER                                 :: max_okushiri_timesteps
  INTEGER                                 :: n_okushiri_timesteps
  INTEGER                                 :: idx_ch5, idx_ch7, idx_ch9
  INTEGER                                 :: outfile_oku

  ! For benchmark_ident==5
    REAL(kind=wp), parameter :: T_frc=2.02
    REAL(kind=wp), parameter :: amp_frc=0.01

  
CONTAINS

subroutine benchmark_init
  
  integer :: n
  integer :: infile_boundary

  ! For benchmark_ident==2, bm_channel_okushiri
  real(kind=wp), dimension(2), parameter :: loc_ch5=(/ 4.521, 1.196 /)
  real(kind=wp), dimension(2), parameter :: loc_ch7=(/ 4.521, 1.696 /)
  real(kind=wp), dimension(2), parameter :: loc_ch9=(/ 4.521, 2.196 /)
  real(kind=wp)                          :: d5, d7, d9, d
  integer                                :: ierr

  select case(benchmark_ident)

  case(bm_channel_okushiri)

     OPEN(newunit=infile_boundary, file=TRIM(MeshPath)//'init_boundary.out',  status='old')

     max_okushiri_timesteps = 9000
     allocate(toku(max_okushiri_timesteps), boku(max_okushiri_timesteps))
     do n=1,max_okushiri_timesteps
        read(infile_boundary,*,iostat=ierr) toku(n), boku(n)
        if (ierr /= 0) then
           n_okushiri_timesteps = n-1
        endif
     end do
     
     
     close(infile_boundary)

     do n=1,nbnd
        ssh0(nod_in_bnd(n))=boku(1)
        ssh1(nod_in_bnd(n))=boku(1)
        ssh2(nod_in_bnd(n))=boku(1)
     end do
 
     ! Initialize output of channel mareograms
     idx_ch5=1
     idx_ch7=1 
     idx_ch9=1
     d5 = (coord_nod2D(1,1) - loc_ch5(1))*(coord_nod2D(1,1) - loc_ch5(1)) &
        + (coord_nod2D(2,1) - loc_ch5(2))*(coord_nod2D(2,1) - loc_ch5(2))
     d7 = (coord_nod2D(1,1) - loc_ch7(1))*(coord_nod2D(1,1) - loc_ch7(1)) &
        + (coord_nod2D(2,1) - loc_ch7(2))*(coord_nod2D(2,1) - loc_ch7(2))
     d9 = (coord_nod2D(1,1) - loc_ch9(1))*(coord_nod2D(1,1) - loc_ch9(1)) &
        + (coord_nod2D(2,1) - loc_ch9(2))*(coord_nod2D(2,1) - loc_ch9(2)) 
     do n=2,nod2D
        d = (coord_nod2D(1,n) - loc_ch5(1))*(coord_nod2D(1,n) - loc_ch5(1)) &
          + (coord_nod2D(2,n) - loc_ch5(2))*(coord_nod2D(2,n) - loc_ch5(2))
        if (d < d5) then
           d5 = d
           idx_ch5 = n
        endif

        d = (coord_nod2D(1,n) - loc_ch7(1))*(coord_nod2D(1,n) - loc_ch7(1)) &
          + (coord_nod2D(2,n) - loc_ch7(2))*(coord_nod2D(2,n) - loc_ch7(2))
        if (d < d7) then
           d7 = d
           idx_ch7 = n
        endif

        d = (coord_nod2D(1,n) - loc_ch9(1))*(coord_nod2D(1,n) - loc_ch9(1)) &
          + (coord_nod2D(2,n) - loc_ch9(2))*(coord_nod2D(2,n) - loc_ch9(2))
        if (d < d9) then
           d9 = d
           idx_ch9 = n
        endif
     end do

  open(newunit=outfile_oku, file=trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_ch579.out',status='new')

  end select

end subroutine benchmark_init

!-----------------------------------------------------------------

subroutine benchmark_adjust_ssh
  
  integer       :: i, n
  integer, save :: i0=2
  real(kind=wp)  :: boku0, boku1, boku2, boku_m1
  real(kind=wp)  :: t0, t1, t2, tm1

  select case(benchmark_ident)

  case(bm_channel_okushiri)
     
     if (time<22.5 ) then

        tm1 = max(time-3*dt,0.)
        t0 = max(time-2*dt,0.)
        t1 = max(time-dt,0.)
        t2 = time
        
        do i = i0,n_okushiri_timesteps
           if (tm1 >= toku(i-1) .and. tm1 < toku(i)) then
              boku_m1 = boku(i-1) + (tm1-toku(i-1))*(boku(i)-boku(i-1))/(toku(i)-toku(i-1))
              i0=i  
           end if
           if (t0 >= toku(i-1) .and. t0 < toku(i)) then
              boku0 = boku(i-1) + (t0-toku(i-1))*(boku(i)-boku(i-1))/(toku(i)-toku(i-1))  
           end if
           if (t1 >= toku(i-1) .and. t1 < toku(i)) then
              boku1 = boku(i-1) + (t1-toku(i-1))*(boku(i)-boku(i-1))/(toku(i)-toku(i-1))  
           end if
           if (t2 >= toku(i-1) .and. t2 < toku(i)) then
              boku2 = boku(i-1) + (t2-toku(i-1))*(boku(i)-boku(i-1))/(toku(i)-toku(i-1))  
              exit
           end if
        end do
     
        if (verbosity>0) print *,time,i, boku_m1, boku0, boku1, boku2
        do n=1,nbnd
           ssh0(nod_in_bnd(n))=boku0 + alpha*(boku1 - 2.*boku0 + boku_m1)
           ssh1(nod_in_bnd(n))=boku1
           ssh2(nod_in_bnd(n))=boku2
        enddo
     end if

     case(5)

        do n=1,nod2D
           if (index_nod2D(n)==5)  then
              ssh2(n)=amp_frc * sin(2.*pi*time/T_frc)
           end if
        end do

  end select

end subroutine benchmark_adjust_ssh

!-----------------------------------------------------------------
subroutine benchmark_output
 
  select case(benchmark_ident)

  case(bm_channel_okushiri)
     write(outfile_oku,*)   real(time,sp), real(ssh2(idx_ch5),sp), &
          real(ssh2(idx_ch7),sp), real(ssh2(idx_ch9),sp)
  end select

end subroutine benchmark_output

!-----------------------------------------------------------------
subroutine benchmark_finalize

  select case(benchmark_ident)

  case(bm_channel_okushiri)
     deallocate(toku,boku)
     close(outfile_oku)
  end select

end subroutine benchmark_finalize


END MODULE BENCHMARK
