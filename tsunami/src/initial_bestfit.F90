subroutine initial_bestfit

  use PARAMETERS
  implicit none

  real(kind=sp), parameter :: no_value=-999.
  real(kind=sp)            :: epi_xy(2,150,25), lon, lat
  integer                  :: rupt_min_i=150, rupt_max_i=0
  integer                  :: i, j, i0, j0
  integer                  :: infile_epi, i_stat
  real(kind=sp)            :: dist, mindist
  CHARACTER(len=3)         :: StrEpi_i
  CHARACTER(len=2)         :: StrEpi_j
  CHARACTER(len=3)         :: StrMW

! ============================
! First, check for a RuptGen scenario
! We assume them to be located in MeshPath
! ============================


! In ruptgen_ijxy.dat, all epicenters are stored
  OPEN(newunit=infile_epi,   file= TRIM(MeshPath) // 'RuptGen2.1/ruptgen_ijxy.out',  status= 'old', iostat=i_stat)

  if (i_stat==0) then

     if (verbosity >0) print *,'Bestfit: Searching for closest RuptGen scenario'

     epi_xy(:,:,:) = no_value

     do while (i_stat==0)
        read(infile_epi,*, iostat=i_stat) i,j,lon,lat
        epi_xy(1,i,j) = lon
        epi_xy(2,i,j) = lat
        rupt_min_i = min(rupt_min_i, i)
        rupt_max_i = max(rupt_max_i, i)
     end do
     
     CLOSE(infile_epi)

! Now: brute force search for closest epicenter
! This can be done much more sophisticated and faster...
     mindist = 180.*180.
     do j=1,25
        do i=rupt_min_i, rupt_max_i
           dist = (epi_xy(1,i,j)-eq_epi_x)*(epi_xy(1,i,j)-eq_epi_x) &
                + (epi_xy(2,i,j)-eq_epi_y)*(epi_xy(2,i,j)-eq_epi_y) 
           if (dist < mindist) then
              mindist = dist
              i0 = i
              j0 = j
           end if
        end do
     end do

     if (mindist <= 0.2) then

        enable_read_initial_field = .true.

        write(StrEpi_i,  '(I3.3)') i0
        write(StrEpi_j,  '(I2.2)') j0
        write(StrMW,     '(F3.1)') eq_mag
        
        IniCondPath   = MeshPath
        init_filename = 'RuptGen2.1/ssh0_mw'//StrMW//'_'//StrEpi_i//'_'//StrEpi_j//'.out'

        if (verbosity > 0) then
           print *,'Bestfit: RuptGen scenario found with epicenter i,j=',i0,j0
           print *,'         RuptGen epicenter lon, lat=',epi_xy(1,i0,j0),epi_xy(2,i0,j0)
           print *,'         Initial conditions will be taken from'
           print *,'        ',trim(IniCondPath)//trim(init_filename)
        endif

        return
     endif
  end if

! ===========================
! Last option: a cosine bell
! =========================== 
  if (verbosity > 0) print *,'Bestfit: No sources found, falling back to cosine bell' 
  enable_idealised = .true.
  initial_shape    = 'cosine'

end subroutine initial_bestfit
