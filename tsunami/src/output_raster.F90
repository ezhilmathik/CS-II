MODULE io_raster

! This module interpolates values from the unstructured mesh to a raster
! and writes directly to the simple raster format Golden Surfer Grid.
! The setup, well, could be better. The possibilities to toggle the output
! by the namelist are rather limited: resolution, bounding box, accuracy, timestep.
! Anything else has to be hard coded, e.g. on land, in the ocean, both; linear
! interpolation or maximum value; fields to be written. 

  use parameters
  use mesh
  use swe
  
  implicit none

  private
  public :: io_raster_init, io_raster_out, io_raster_finalize, io_raster_ttt

  INTERFACE interpolate
     module procedure interpolate_r8
     module procedure interpolate_r4
  END INTERFACE interpolate
  
  real(kind=sp)      :: raster_ssh_accuracy  ! accuracy of MWH/SSH output 
  real(kind=sp)      :: raster_eta_accuracy  ! accuracy of ETA output 
  real(kind=sp)      :: raster_dx, raster_dy ! pixel size for the raster
                                             ! degree for geo coord., meter for cartesian
  real(kind=sp)      :: raster_ssh_T_out     ! time step [sec] for SSH slices
  real(kind=sp)      :: raster_boundingbox(4) ! bounding box: [xmin, xmax, ymin, ymax] 
  logical            :: raster_ascii         ! Ascii or binary output?
  
  real(kind=sp)      :: no_value = -999999.
 

  
  type rastertype
     real(kind=sp)               :: dx, dy        ! pixel size
     integer                     :: M, N          ! raster grid 
     real(kind=sp), allocatable  :: x(:), y(:)    ! raster coordinates, middle of pixel
     real(kind=sp), allocatable  :: A(:,:)        ! Raster data
     real(kind=sp)               :: Amin          ! data range
     real(kind=sp)               :: Amax          !   ''
     integer                     :: mstart        ! All values /=0 and /= undefined
     integer                     :: mend          ! are inside
     integer                     :: nstart        ! A(mstart:mend, nstart:nend) 
     integer                     :: nend          !
     character(len=3)            :: interpolation_type ! 'lin', 'max'
     
  ! For raster data and the interpolation from the triangular mesh
     real(kind=sp), allocatable  :: w(:,:,:)      ! interpolation weights
     integer, allocatable        :: nodi(:,:,:)   ! (element) nodes from which to interpolate

     integer, allocatable        :: nodesInPixel_idx(:)
     integer, allocatable        :: nodesInPixel(:)
  end type rastertype
  
  ! For the start, we need only one raster data type with lin. interpolation,
  ! for one resolution
  type(rastertype) :: rg    ! r=raster, g=general
  
  ! And for MWH on land, the worst case raster
  type(rastertype) :: rmwh  ! r=raster, for mwh
  
! For the Surfer 6 Grid format: parameter for "no value"
  real(kind=sp), parameter      :: r_notdefined = 1.70141e+038
  
contains

SUBROUTINE io_raster_init

  call read_raster_namelist
! The "general" interpolation for ssh, ttt: in water, linear
  call setup_raster_and_linearInterpolation(raster_dx, raster_dy,rg)

! An the worst case setup for MWH   
!NR Skipped, as it is buggy, and done only once for MWH. No need to
!NR set up auilliary fields and such.
!  call setup_raster_and_maxVal(raster_dx, raster_dy,rmwh)

END SUBROUTINE io_raster_init

SUBROUTINE io_raster_ttt

  character(len=200) :: grdfile_ttt

  call interpolate(ttt, raster_eta_accuracy, rg)

  grdfile_ttt = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_ttt.grd'
  call write_raster(grdfile_ttt,rg)


END SUBROUTINE io_raster_ttt


!=====================================================================

SUBROUTINE setup_raster_and_linearInterpolation(dx,dy,r)

  use mesh
  
  implicit none
  
  real(kind=sp), intent(in)     :: dx, dy
  type(rastertype), intent(out) :: r

  real(kind=sp), parameter :: cyclic_length    = 360.
  
  real(kind=sp)  :: xi_min, xi_max, yi_min, yi_max 
  real(kind=sp)  :: x, x1, x2, x3, y, y1, y2, y3, w1, w2, w3
  real(kind=sp)  :: xmax, xmin, ymin, ymax, detinv
  real(kind=sp)  :: invdx, invdy
  integer        :: m, n, el, mstart, mend, nstart, nend
  real(kind=sp)  :: geo_factor

  r%interpolation_type='lin'
  
  r%dx = dx
  r%dy = dy

  invdx = 1./dx
  invdy = 1./dy

  geo_factor=1.  
  if (coordinate_type==1) geo_factor = 1./rad

  ! Determine raster
  xi_min= real(  floor(raster_boundingbox(1)*invdx))*dx
  xi_max= real(ceiling(raster_boundingbox(2)*invdx))*dx
  yi_min= real(  floor(raster_boundingbox(3)*invdy))*dy
  yi_max= real(ceiling(raster_boundingbox(4)*invdy))*dy
  
  r%M = ceiling((xi_max-xi_min)*invdx) +1
  r%N = ceiling((yi_max-yi_min)*invdy) +1

  allocate(r%x(r%M), r%y(r%N))
  allocate(r%A(r%M, r%N))
  allocate(r%w(3,r%M,r%N), r%nodi(3,r%M,r%N))
  
  ! Coordinates of the raster
!$OMP target teams & !!!!$OMP PARALLEL default(shared) &
!$OMP           private(m, x1, x2, x3, xmin, xmax, mstart, mend, &
!$OMP                   n, y1, y2, y3, ymin, ymax, nstart, nend, &
!$OMP                      w1, w2, w3, detinv, el)
!$OMP distribute parallel DO simd
  do m=1,r%M
     r%x(m) = xi_min + real(m-1)*dx
  enddo
!$OMP END distribute parallel DO simd !!!!!$OMP END DO NOWAIT
!$OMP distribute parallel DO simd
  do n=1,r%N
     r%y(n) = yi_min + real(n-1)*dy
  enddo
!$OMP END distribute parallel DO simd !!!!!!$OMP END DO NOWAIT
!$OMP distribute parallel DO simd  
  ! Initialize interpolation weights
  do n=1,r%N
     do m=1,r%M
        r%nodi(1:3,m,n) = 0
        r%w(   1:3,m,n) = 0.
     enddo
  enddo
!$OMP END distribute parallel DO simd
!$OMP END target teams 
  
  !NR To be on the save side, no OpenMP for this loop.
  !NR It is possible that raster nodes are handled multiple times
  !NR by two elements (located on the edge) or even by a full patch
  !NR if raster node == mesh node.
  do el =1,elem2D

!NR We also need the values on land for LEXIS
     if (.true.) then
!NR    if (minval(nodhn_init(elem2D_nodes(1:3,el))) >= 0.) then
!NR        ! consider only wet nodes here
!NR         ! and consider to implement a namelist switch
!NR         ! to toggle the output at land
        
        
        y1 = coord_nod2D(2, elem2D_nodes(1,el)) * geo_factor
        y2 = coord_nod2D(2, elem2D_nodes(2,el)) * geo_factor
        y3 = coord_nod2D(2, elem2D_nodes(3,el)) * geo_factor
        ymin= min(min(y1,y2),y3)
        ymax= max(max(y1,y2),y3)
        nstart = max(0, ceiling((ymin-r%y(1))*invdy)) +1
        nend   = min(r%N-1, floor((ymax-r%y(1))*invdy)) +1    
        
        if (nstart > nend) cycle


        x1 = coord_nod2D(1, elem2D_nodes(1,el)) * geo_factor
        x2 = coord_nod2D(1, elem2D_nodes(2,el)) * geo_factor
        x3 = coord_nod2D(1, elem2D_nodes(3,el)) * geo_factor
        
        xmin= min(min(x1,x2),x3)
        xmax= max(max(x1,x2),x3)
        
        
        if (coordinate_type==1) then
           ! Does this element cross the 180E/W-meridian?
           if (xmax - xmin > 180.) then
              ! We have to regard this element twice, once for -180 degr at the
              ! left boarder of the the raster, once for +180degr at the right
              ! boarder.

              ! 1. case, shift all x-values to around -180degr
              if (x1 > 0.) x1 = x1 - cyclic_length
              if (x2 > 0.) x2 = x2 - cyclic_length
              if (x3 > 0.) x3 = x3 - cyclic_length

              xmin= min(min(x1,x2),x3)
              xmax= max(max(x1,x2),x3)
 
              mstart = max(0, ceiling((xmin-r%x(1))*invdx)) +1
              mend   = min(r%M-1, floor((xmax-r%x(1))*invdx)) +1    
        
              if (mstart <= mend)  call interpolation_weight_for_el
              
              ! 2. case, from around -180degr shift to +180degr
              x1 = x1 + cyclic_length
              x2 = x2 + cyclic_length
              x3 = x3 + cyclic_length
              xmin= xmin + cyclic_length
              xmax= xmax + cyclic_length
           endif
        endif

        mstart = max(0, ceiling((xmin-r%x(1))*invdx)) +1
        mend   = min(r%M-1, floor((xmax-r%x(1))*invdx)) +1    
        
        if (mstart <= mend) call interpolation_weight_for_el
        
     end if ! element wet?
  end do ! el

contains

  subroutine interpolation_weight_for_el
    
    detinv=1._8/((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
        
    do n = nstart, nend
       y = r%y(n)
       do m = mstart, mend
          x = r%x(m)
          w3=((x1-x)*(y2-y) - (x2-x)*(y1-y))
          
          if (w3 >= 0.) then
             w2=((x3-x)*(y1-y) - (x1-x)*(y3-y))
             if (w2>=0.) then
                w1=((x2-x)*(y3-y) - (x3-x)*(y2-y))
                if (w1 >= 0.) then
                   r%w(1,m,n) = w1*detinv
                   r%w(2,m,n) = w2*detinv
                   r%w(3,m,n) = w3*detinv
                   r%nodi(1:3,m,n) = elem2D_nodes(1:3,el)
                end if ! w1
             end if ! w2
          end if ! w3
       end do ! m
    end do ! n
        
  end subroutine interpolation_weight_for_el
end SUBROUTINE setup_raster_and_linearInterpolation

!=====================================================================

SUBROUTINE setup_raster_and_maxVal(dx,dy,r)

!NR  Not yet ready!!! This routine seems to be buggy, and even if fixed 
!NR  this approach would only work correctly for pixels >> triangles.

  use mesh
  
  implicit none
  
  real(kind=sp), intent(in)     :: dx, dy
  type(rastertype), intent(out) :: r

  real(kind=sp), parameter :: cyclic_length    = 360.
  
  real(kind=sp)  :: xi_min, xi_max, yi_min, yi_max 
  real(kind=sp)  :: xmax, xmin, ymin, ymax
  real(kind=sp)  :: invdx, invdy
  integer        :: m, n, el, mstart, mend, nstart, nend
  real(kind=sp)  :: geo_factor  
  integer        :: max_nodesInPixel, i, j, k, iA, iB
  
  r%interpolation_type='max'
  
  r%dx = dx
  r%dy = dy

  invdx = 1./dx
  invdy = 1./dy

  geo_factor=1.  
  if (coordinate_type==1) geo_factor = 1./rad

  ! Determine raster
  xi_min= real(  floor(raster_boundingbox(1)*invdx))*dx
  xi_max= real(ceiling(raster_boundingbox(2)*invdx))*dx
  yi_min= real(  floor(raster_boundingbox(3)*invdy))*dy
  yi_max= real(ceiling(raster_boundingbox(4)*invdy))*dy
  
  r%M = ceiling((xi_max-xi_min)*invdx) +1
  r%N = ceiling((yi_max-yi_min)*invdy) +1

  allocate(r%x(r%M), r%y(r%N))
  allocate(r%A(r%M, r%N))
  allocate(r%nodesInPixel_idx(0:r%M*r%N))
  
  ! Coordinates of the raster
!$OMP target teams & !!!!$OMP PARALLEL default(shared) &
!$OMP           private(m, n, k)
!$OMP distribute parallel DO simd
  do m=1,r%M
     r%x(m) = xi_min + real(m-1)*dx
  enddo
!$OMP END distribute parallel DO simd !!!!$OMP END DO NOWAIT
!$OMP distribute parallel DO simd 
  do n=1,r%N
     r%y(n) = yi_min + real(n-1)*dy
  enddo
!$OMP END distribute parallel DO simd !!!!$OMP END DO NOWAIT
!$OMP distribute parallel DO simd
  do k=0,r%M*r%N
     r%nodesInPixel_idx(k) = 0
  end do
!$OMP END distribute parallel DO simd
!$OMP END target teams 

 
  max_nodesInPixel = 0
  do n=1,nod2D
     ! Only nodes on land
     if (nodhn_init(n) <= 0._wp) then

        i = floor((coord_nod2D(1,n)*geo_factor - xi_min +.5*dx))*invdx +1
        j = floor((coord_nod2D(2,n)*geo_factor - yi_min +.5*dy))*invdy +1
        
        if (i>0 .and. i<=r%M .and. j>0 .and. j<=r%N) &
             r%nodesInPixel_idx(i+ (j-1)*r%M) = r%nodesInPixel_idx(i+(j-1)*r%M) + 1
        
     end if ! node on land?
  enddo


  do k = 2, r%M*r%N
     r%nodesInPixel_idx(k) = r%nodesInPixel_idx(k) + r%nodesInPixel_idx(k-1) 
  end do

  allocate(r%nodesInPixel(r%nodesInPixel_idx(r%M*r%N)))  ! should be nod2D at most, but nodes on a raster edge might be counted twice

  r%nodesInPixel(:) = -1

  do n=1,nod2D
     ! Only nodes on land
     if (nodhn_init(n) <= 0._wp) then

        i = floor((coord_nod2D(1,n)*geo_factor - xi_min +.5*dx))*invdx +1
        j = floor((coord_nod2D(2,n)*geo_factor - yi_min +.5*dy))*invdy +1
        
        if (i>0 .and. i<=r%M .and. j>0 .and. j<=r%N) then

           iA = r%nodesInPixel_idx(i+(j-1)*r%M-1)+1 
           iB = r%nodesInPixel_idx(i+(j-1)*r%M)
           do k=iA, iB
              if (r%nodesInPixel(k) == -1) then
                 r%nodesInPixel(k) = n
                 exit ! k-loop
              endif
           enddo
        endif ! node in bounding box?
     end if ! node on land?
  enddo
  
  
end SUBROUTINE setup_raster_and_maxVal
!======================================================================

SUBROUTINE read_raster_namelist

  implicit none
  
  NAMELIST /output_raster/ raster_dx, raster_dy, raster_ssh_T_out, &
       raster_ssh_accuracy, raster_eta_accuracy, raster_boundingbox, raster_ascii

  character(len=7) :: coordunit
  
  integer  :: infile_nml

  if (coordinate_type==1) then
     coordunit = 'degrees'
  else
     coordunit = 'meters '
  endif
  
  ! Initialize values as "unset" 

  raster_dx = -1.  ! pixel size (degree) for the raster
  raster_dy = -1.  ! assume "square" pixels 

  raster_ssh_T_out = -1.  ! write timesteps every 10 minutes

  raster_ssh_accuracy = -1.  ! round ssh to given value [m], or do not round for <= 0.
  raster_eta_accuracy = -1.  ! round eta to given value [min], or do not round for <= 0.

  raster_boundingbox=(/1.,-1.,1.,-1./)

  raster_ascii = .false.

  OPEN(newunit=infile_nml, file = 'namelist.tsunami')
  READ(infile_nml, NML = output_raster)
  CLOSE(infile_nml)

  if (verbosity > 0) then
     print *,'Raster data (surfer grid) will be written for ssh and mwh'
  endif

  ! Check input and set default values
  
  if (raster_dx <= 0. ) then
     if (raster_dy <= 0.) then
        raster_dx = .02
        raster_dy = raster_dx
        write(*,*) '/output_raster/ Setting missing values raster_dx= 0.2',coordunit
        write(*,*) '/output_raster/                    and raster_dy= 0.2',coordunit
     else
        raster_dx = raster_dy
        write(*,*) '/output_raster/ Setting missing value raster_dx to raster_dy (',raster_dx,coordunit,')'
     endif
  else
     if (raster_dy <= 0.) then
        raster_dy = raster_dx
        write(*,*) '/output_raster/ Setting missing value raster_dy to raster_dx (',raster_dy,coordunit,')'
     end if
  endif

  if (raster_ssh_T_out <= 0.) then
     raster_ssh_T_out = -1
     write(*,*) '/output_raster/ raster_ssh_T_out is not set or less than zero'
     write(*,*) '                There will be no output of SSH raster data.'

     snap_int_raster = 2*T_end / dt + 100  ! that's never reached during the tsunami run
                                           ! 2* and +100 to be on the very safe side
  else
     snap_int_raster = max(1, int(raster_ssh_T_out/dt))
  endif
  
  if ( raster_boundingbox(1) >= raster_boundingbox(2) .or. &
       raster_boundingbox(3) >= raster_boundingbox(4) .or. &
       raster_boundingbox(1) >= BoundingBox_xmax      .or. &
       raster_boundingbox(2) <= BoundingBox_xmin      .or. &
       raster_boundingbox(3) >= BoundingBox_ymax      .or. &
       raster_boundingbox(4) <= BoundingBox_ymin        ) then


     write(*,*) '/output_raster/ No or invalid bounding box specified,'
     write(*,*) '                the whole computional area will be used.'

     raster_boundingbox(1) = BoundingBox_xmin
     raster_boundingbox(2) = BoundingBox_xmax
     raster_boundingbox(3) = BoundingBox_ymin
     raster_boundingbox(4) = BoundingBox_ymax

  endif

END SUBROUTINE read_raster_namelist

!====================================
SUBROUTINE io_raster_out

  character(len=6)   :: timestring  ! that's enough for 999999s = 11 days
  character(len=200) :: grdfile_ssh

  write(timestring,'(i6.6)') int(time)
  
  call interpolate(ssh2,raster_ssh_accuracy, rg)

  grdfile_ssh = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_ssh_'//timestring//'.grd'
  call write_raster(grdfile_ssh,rg)
  
END SUBROUTINE io_raster_out

SUBROUTINE write_raster(grdfile, r)

  character(len=200), intent(in) :: grdfile
  type(rastertype), intent(in)   :: r
  integer                        :: outfile_grd
  integer                        :: m,n
  
  if (verbosity >0) then
     print *,'Writing raster output to', trim(grdfile)
  endif
  
  if (.not. raster_ascii) then
     ! Write binary file

     open(newunit=outfile_grd,file=trim(grdfile), form='unformatted',access='stream')
     write(outfile_grd) 'DSBB'
     write(outfile_grd) int(r%mend-r%mstart+1,2), int(r%nend-r%nstart+1,2)
     write(outfile_grd) real(r%x(r%mstart),dp),   real(r%x(r%mend),dp)
     write(outfile_grd) real(r%y(r%nstart),dp),   real(r%y(r%nend),dp)
     write(outfile_grd) real(r%Amin,dp),          real(r%Amax,dp)
     write(outfile_grd) r%A(r%mstart:r%mend, r%nstart:r%nend) 
     close(outfile_grd)

  else
     ! Write ascii file

     open(newunit=outfile_grd,file=trim(grdfile))
     write(outfile_grd,*) 'DSBA'
     write(outfile_grd,*) int(r%mend-r%mstart+1,2), int(r%nend-r%nstart+1,2)
     write(outfile_grd,*) real(r%x(r%mstart),dp),   real(r%x(r%mend),dp)
     write(outfile_grd,*) real(r%y(r%nstart),dp),   real(r%y(r%nend),dp)
     write(outfile_grd,*) real(r%Amin,dp),          real(r%Amax,dp)
     do n=r%nstart,r%nend
        do m=r%mstart,r%mend
           if (r%A(m,n) == r_notdefined) then
              write(outfile_grd,fmt="(a)",advance="no") '1.70141e+038'
           else 
              write(outfile_grd,fmt="(f0.3)",advance="no") r%A(m,n)
           endif
           write(outfile_grd,fmt="(a)",advance="no") ' '
        end do
     end do
     close(outfile_grd)
     
  end if
  
END SUBROUTINE write_raster

!====================================
SUBROUTINE io_raster_finalize

  character(len=200) :: grdfile_eta, grdfile_mwh, grdfile_maxMWH_land
  real :: dx, dy
  
  ! Here, write MWH [and ETA - but w/o quality control, always prefer TTT!)

!  grdfile_eta = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_eta.grd'
  grdfile_mwh = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_mwh.grd'
  
  call interpolate(mwh, raster_ssh_accuracy, rg)
  call write_raster(grdfile_mwh,rg)
  
!  call interpolate(arrival_time, raster_eta_accuracy, rg)
!  call write_raster(grdfile_eta,rg)

  grdfile_mwh = trim(OutputPath)//trim(TSID)//trim(output_prefix)//'_mwh_maxLand.grd'

  call interpolate_maxOnLand(mwh, raster_ssh_accuracy, rg)
  call write_raster(grdfile_mwh,rg)


END SUBROUTINE io_raster_finalize

!====================================

SUBROUTINE interpolate_r8(val, eps, r)

  real(kind=dp), intent(in)       :: val(nod2D)
  real(kind=sp), intent(in)       :: eps
  type(rastertype), intent(inout) :: r
  real(kind=sp)                   :: epsinv
  integer                         :: m, n
  integer                         :: mstart, mend
  integer                         :: nstart, nend
  real(kind=sp)                   :: Amin, Amax

    mstart = r%M
    mend = 1
    nstart = r%N
    nend = 1
    
    Amin=0.
    Amax=0.

    if (eps > 0.) then
       epsinv=1./eps     ! with rounding

!$OMP target teams distribute PARALLEL DO & !!!!$OMP PARALLEL DO default(shared) &
!$OMP             private(n,m) &
!$OMP             reduction(min:mstart, nstart, Amin) &       
!$OMP             reduction(max:mend,   nend,   Amax)       
       do n=1,r%N
          do m=1,r%M
             r%A(m,n) = r_notdefined
             
             if (r%nodi(1,m,n)>0) then             
                
                r%A(m,n) = eps*int(epsinv*sum(r%w(1:3,m,n)*real(val(r%nodi(1:3,m,n)),4)))
                
                if (abs(r%A(m,n)) > 0.) then
                   mstart = min(m, mstart)
                   mend   = max(m, mend)
                   nstart = min(n, nstart)
                   nend   = max(n, nend)
                endif
                Amin = min(Amin,r%A(m,n))
                Amax = max(Amax,r%A(m,n))
             endif
          end do
       end do
!$OMP END target teams distribute PARALLEL DO       
    else
       ! No rounding
       
!$OMP target teams distribute PARALLEL DO default(shared) &
!$OMP             private(n,m) &
!$OMP             reduction(min:mstart, nstart, Amin) &       
!$OMP             reduction(max:mend,   nend,   Amax)       
       do n=1,r%N
          do m=1,r%M
             r%A(m,n) = r_notdefined
             
             if (r%nodi(1,m,n)>0) then             
                
                r%A(m,n) = sum(r%w(1:3,m,n)*real(val(r%nodi(1:3,m,n)),4))
                
                if (abs(r%A(m,n)) > 0.) then
                   mstart = min(m, mstart)
                   mend   = max(m, mend)
                   nstart = min(n, nstart)
                   nend   = max(n, nend)
                endif
                Amin = min(Amin,r%A(m,n))
                Amax = max(Amax,r%A(m,n))
             endif
          end do
       end do
!$OMP END target teams distribute PARALLEL DO       
    end if

    r%mstart = mstart
    r%mend   = mend 
    r%nstart = nstart
    r%nend   = nend
    
    r%Amin   = Amin
    r%Amax   = Amax


  END SUBROUTINE interpolate_r8

!====================================
SUBROUTINE interpolate_r4(val, eps, r)

  real(kind=sp), intent(in)       :: val(nod2D)
  real(kind=sp), intent(in)       :: eps
  type(rastertype), intent(inout) :: r
  real(kind=sp)                   :: epsinv
  integer                         :: m, n, iA, iB
  

  r%mstart = r%M
  r%mend = 1
  r%nstart = r%N
  r%nend = 1
    
  r%Amin=0.
  r%Amax=0.

  select case(r%interpolation_type)

  case('lin')

     if (eps > 0.) then
        epsinv=1./eps     ! with rounding

        do n=1,r%N
           do m=1,r%M
              r%A(m,n) = r_notdefined

              if (r%nodi(1,m,n)>0) then             

                 r%A(m,n) = eps*int(epsinv*sum(r%w(1:3,m,n)*val(r%nodi(1:3,m,n))))

                 if (abs(r%A(m,n)) > 0.) then
                    r%mstart = min(m, r%mstart)
                    r%mend   = max(m, r%mend)
                    r%nstart = min(n, r%nstart)
                    r%nend   = max(n, r%nend)
                 endif
                 r%Amin = min(r%Amin,r%A(m,n))
                 r%Amax = max(r%Amax,r%A(m,n))
              endif
           end do
        end do

     else
        ! No rounding
        do n=1,r%N
           do m=1,r%M
              r%A(m,n) = r_notdefined

              if (r%nodi(1,m,n)>0) then             

                 r%A(m,n) = sum(r%w(1:3,m,n)* val(r%nodi(1:3,m,n)) )

                 if (abs(r%A(m,n)) > 0.) then
                    r%mstart = min(m, r%mstart)
                    r%mend   = max(m, r%mend)
                    r%nstart = min(n, r%nstart)
                    r%nend   = max(n, r%nend)
                 endif
                 r%Amin = min(r%Amin,r%A(m,n))
                 r%Amax = max(r%Amax,r%A(m,n))
              endif
           end do
        end do
     end if


  case('max')

!NR  Not yet ready!!! This maximum value extraction seems to be buggy, and 
!NR  even if fixed this approach would only work correctly for pixels >> triangles.

     if (eps > 0.) then
        epsinv=1./eps     ! with rounding

        do n=1,r%N
           do m=1,r%M
              r%A(m,n) = r_notdefined
              
              iA = r%nodesInPixel_idx(m + (n-1)*r%M -1) +1
              iB = r%nodesInPixel_idx(m + (n-1)*r%M)
              
              if (iB > iA) then
              
                 r%A(m,n) = eps*int(epsinv*maxval(val(r%nodesInPixel(iA:iB))))
              
                 if (abs(r%A(m,n)) > 0.) then
                    r%mstart = min(m, r%mstart)
                    r%mend   = max(m, r%mend)
                    r%nstart = min(n, r%nstart)
                    r%nend   = max(n, r%nend)
                 endif
                 r%Amin = min(r%Amin,r%A(m,n))
                 r%Amax = max(r%Amax,r%A(m,n))
              endif

              r%A(m,n) = real(iB-iA)
           end do
        end do

     else ! no rounding
        do n=1,r%N
           do m=1,r%M
              r%A(m,n) = r_notdefined
              
              iA = r%nodesInPixel_idx(m + (n-1)*r%M -1) +1
              iB = r%nodesInPixel_idx(m + (n-1)*r%M)
              
              if (iB > iA) then
              
                 r%A(m,n) = maxval(val(r%nodesInPixel(iA:iB)))
              
                 if (abs(r%A(m,n)) > 0.) then
                    r%mstart = min(m, r%mstart)
                    r%mend   = max(m, r%mend)
                    r%nstart = min(n, r%nstart)
                    r%nend   = max(n, r%nend)
                 endif
                 r%Amin = min(r%Amin,r%A(m,n))
                 r%Amax = max(r%Amax,r%A(m,n))
              endif
           end do
        end do
     endif ! rounding

   end select


  END SUBROUTINE interpolate_r4

  SUBROUTINE finalize_raster(r)
    type(rastertype), intent(inout) :: r

    deallocate(r%A)
    if (allocated(r%x)) deallocate(r%x, r%y)
    if (allocated(r%w)) deallocate(r%w, r%nodi)
    if (allocated(r%nodesInPixel)) deallocate(r%nodesInPixel, r%nodesInPixel_idx)
    
  END SUBROUTINE finalize_raster

!=====================================================================

SUBROUTINE  interpolate_maxOnLand(val, eps, r)

  use mesh
  
  implicit none
  
  real(kind=sp), intent(in)     :: val(nod2D)
  real(kind=sp), intent(in)     :: eps
  type(rastertype), intent(inout) :: r

  
  real(kind=sp)  :: invdx, invdy, Amin, Amax, dx, dy
  integer        :: i, j, m, n, mstart, mend, nstart, nend, nland
  real(kind=sp)  :: geo_factor, value
  real(kind=sp), parameter :: mwh_threshold=0.001

  dx    = r%dx
  dy    = r%dy
  invdx = 1./dx
  invdy = 1./dy

  geo_factor=1.  
  if (coordinate_type==1) geo_factor = 1./rad

  mstart = r%M
  mend = 1
  nstart = r%N
  nend = 1
  
  Amin=0.
  Amax=0.

  ! We need a linear interpolation at those pixels which do not contain
  ! mesh nodes.
  ! However, only nodes on land must be used for this interpolation.
  ! Otherwise, the inundation height can be overestimated, in particular
  ! in regions with a very steep coast.

  do n=1,r%N
     do m=1,r%M
        r%A(m,n) = r_notdefined

        if (r%nodi(1,m,n)>0) then             
           
           nland = count(nodhn_init(r%nodi(1:3,m,n)) <= 0.)
           value = 0.
           
           if (nland == 3) then
              value = sum(r%w(1:3,m,n)* val(r%nodi(1:3,m,n)) )
           elseif (nland >= 1) then
              ! interpolate bathymetry to determine if it is a pixel
              ! on land
              
              if (sum(r%w(1:3,m,n)* nodhn_init(r%nodi(1:3,m,n)) ) <= 0.) then
                 ! Now, only take values at land nodes into consideration
                 if (nodhn_init(r%nodi(1,m,n)) <=0.) value = val(r%nodi(1,m,n))
                 if (nodhn_init(r%nodi(2,m,n)) <=0.) value = max(value, val(r%nodi(2,m,n)))
                 if (nodhn_init(r%nodi(3,m,n)) <=0.) value = max(value, val(r%nodi(3,m,n)))
              endif
           endif

             
           if (value > mwh_threshold) then
              mstart = min(m, mstart)
              mend   = max(m, mend)
              nstart = min(n, nstart)
              nend   = max(n, nend)

              
              Amax = max(Amax,r%A(m,n))
              r%A(m,n) = value
           endif
        endif
     end do
  end do
  
  !NR consider to switch the two blocks above and below: First, maximum value of all nodes,
  !NR fill gaps second.
  
  ! Now, we collect the maximum value at all nodes inside a pixel

  do n=1,nod2D
     ! Only wet nodes on land 
     if (nodhn_init(n)<=0._wp .and. val(n)>=mwh_threshold ) then

        i = floor((coord_nod2D(1,n)*geo_factor - r%x(1) + .5*dx)*invdx) +1
        j = floor((coord_nod2D(2,n)*geo_factor - r%y(1) + .5*dy)*invdy) +1

        if (i>0 .and. i<=r%M .and. j>0 .and. j<r%N) then

           if (r%A(i,j) == r_notdefined) then
              r%A(i,j) = val(n)
                
              mstart = min(i, mstart)
              mend   = max(i, mend)
              nstart = min(j, nstart)
              nend   = max(j, nend)
              
              Amax = max(Amax,r%A(i,j))
           else
              r%A(i,j) = max(r%A(i,j), val(n))
              Amax = max(Amax,r%A(i,j))
           endif 
        endif ! node in bounding box?
     end if ! node on land and mwh>0?
  enddo

  r%interpolation_type='max'

  r%mstart = mstart
  r%mend   = mend
  r%nstart = nstart
  r%nend   = nend
  
!  print *,'rastersize:',mstart,mend,nstart,nend  
!  print *,'rastersize:',mend-mstart+1, nend-nstart+1

  r%Amin = 0.
  r%Amax = Amax

end SUBROUTINE interpolate_maxOnLand
END MODULE io_raster
