! Part of TsunAWI
! Copyright (c) 2014
! Alfred Wegener Institut,
! Helmholtz Centre for Polar and Marine Research
! Contact: rakowsky@awi.de

! MODULE: MESH --  setup 2D mesh used for the tsunami simulations and the nodal properties of the mesh

MODULE MESH

  USE DATA_TYPES
  USE PARAMETERS
  
  IMPLICIT NONE

  SAVE  


  INTEGER                                           :: elem2D
  INTEGER,   ALLOCATABLE, DIMENSION(:,:)            :: elem2D_nodes

  INTEGER,   ALLOCATABLE, DIMENSION(:,:)            :: edge_nod2D
  INTEGER                                           :: edg2D
  INTEGER,   ALLOCATABLE, DIMENSION(:,:)            :: edge_in_elem2D
  INTEGER,   ALLOCATABLE, DIMENSION(:,:)            :: elem2D_edges
  REAL(kind=wp),      ALLOCATABLE, DIMENSION(:,:)   :: edg_nxy
  REAL(kind=wp),      ALLOCATABLE, DIMENSION(:)     :: edg_length

  INTEGER                                           :: nod2D, nbnd, nedgbnd
 
  REAL(kind=wp),      ALLOCATABLE, DIMENSION(:,:)   :: coord_nod2D
  real(kind=wp)                                     :: BoundingBox_xmin 
  real(kind=wp)                                     :: BoundingBox_xmax 
  real(kind=wp)                                     :: BoundingBox_ymin 
  real(kind=wp)                                     :: BoundingBox_ymax 
  INTEGER,           ALLOCATABLE, DIMENSION(:)      :: index_nod2D
  REAL(kind=wp),      ALLOCATABLE, DIMENSION(:)     :: cos_nod
  INTEGER,           ALLOCATABLE, DIMENSION(:)      :: nod_in_bnd, edg_in_bnd
  INTEGER,           ALLOCATABLE, DIMENSION(:,:)    :: nodngh_in_bnd
  REAL(kind=wp),      ALLOCATABLE, DIMENSION(:,:)   :: nodedgfac_in_bnd

  TYPE(addresstype), ALLOCATABLE, DIMENSION(:)      :: nghbr_nod2D

  REAL(kind=wp),      ALLOCATABLE, DIMENSION(:)     :: nodhn, nodhn_aux
  REAL(kind=sp),      ALLOCATABLE, DIMENSION(:)     :: nodhn_init
  REAL(kind=wp),      ALLOCATABLE, DIMENSION(:)     :: rgh_edg

! To streamline compute_velocity_at_nodes
! Compact storage for edges adjacent to nodes
  integer,           ALLOCATABLE, DIMENSION(:)      :: nodedge_ptr, nodedge_addr
  real(kind=wp),      ALLOCATABLE, DIMENSION(:)     :: nodedge_wght

! Compact storage for elements adjacent to nodes
  integer,           ALLOCATABLE, DIMENSION(:)      :: nodelem_ptr, nodelem_addr
  real(kind=wp),      ALLOCATABLE, DIMENSION(:,:)   :: nodelem_wght ! To streamline compute_ssh

  real(kind=wp)                                     :: shortestEdge, longestEdge

  CHARACTER(len=100)                                :: netcdf_file
  CHARACTER(len=100)                                :: restart_file

  
!NR Remainer of hard-coded extraction of satellite tracks for 2004
!NR Indian Ocean Tsunami
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:) :: J1_data, TP_data
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:) :: TrBF_J1, TrBF_TP
  INTEGER, ALLOCATABLE, DIMENSION(:) :: TrElems_J1, TrElems_TP

  ! bottom friction, variable bottom roughness
  real(KIND=wp), allocatable, dimension(:) :: mng_edg2D


CONTAINS
  

SUBROUTINE mesh_scaling  

  ! Transforms degrees in rad in coord_nod2D(2,nod2D)     
  ! Constructs the arrays cos_nod(nod2D)

  !use ELEMENTS

  IMPLICIT NONE

  INTEGER      :: n
  

! Set scaling for earth or plane,
! define element wise field of mean cosine(latitude) in triangles.

    if ( coordinate_type == 1 ) then
! on the sphere 
       scaling = r_earth 
       ALLOCATE(cos_nod(nod2D))
!$OMP PARALLEL DO default(shared) private(n) 
       DO n=1, nod2D         
          ! transform coordinates from deg to rad
          coord_nod2D(1:2,n) = coord_nod2D(1:2,n) * rad

          ! scaling according to latitude
          cos_nod(n) = COS(coord_nod2D(2,n))
       END DO
!$OMP END PARALLEL DO      

    elseif ( coordinate_type == 2 ) then
! in the plane 
      scaling=1._wp
    end if

END SUBROUTINE mesh_scaling


SUBROUTINE read_2Dmesh()
  !
  !--------------------------------------------------------------------------
  ! READS 2D MESH AND VARIABLES CONCERNING THE 2D MESH    
  !
  ! OUTPUT:
  !   Number of 2D-Nodes: nod2D
  !   Number of 2D-Elements: elem2D
  !   Number of 2D-element-neighbours is always 3
  !
  !   Mesh-Arrays: 
  !     coord_nod2D(2,nod2D)
  !     elem2D_nodes(3,elem2D)
  !     index_nod2D(nod2D)
  !     nodhn(nod2D)
  !--------------------------------------------------------------------------
  !
  !use ELEMENTS
  !
!$  use OMP_LIB

  IMPLICIT NONE
  !
  INTEGER            :: infile_nod, infile_elem, infile_depth
  INTEGER            :: nq, ind
  REAL(kind=wp)      :: x1, x2
  INTEGER            :: n, i_err
  character(len=200) :: a_dummy
  
  OPEN(newunit=infile_nod,   file= TRIM(MeshPath) // 'nod2d.out',  status= 'old')
  OPEN(newunit=infile_elem,  file= TRIM(MeshPath) // 'elem2d.out', status= 'old')
  OPEN(newunit=infile_depth, file= TRIM(MeshPath)//TRIM(TopoFile),  status= 'old')
  WRITE(*,*) '2D mesh is opened'

  READ(infile_nod, *) nod2D       

  ALLOCATE(coord_nod2D(2,nod2D), index_nod2D(nod2D))
!$OMP PARALLEL DO 
  DO n=1, nod2D
     coord_nod2D(:,n) = 0.
     index_nod2D(n)   = 0.
  END DO
!$OMP END PARALLEL DO

  ! Test, if nod2d.out contains 3 or 4 columns.
  ! The first column, the number of the node, may be skipped.
  READ(infile_nod,'(a)') a_dummy
  READ(a_dummy,*,IOSTAT=i_err) nq, x1, x2, ind

  if (i_err==0) then
     coord_nod2D(1,nq) = x1
     coord_nod2D(2,nq) = x2
     index_nod2D(nq)   = ind
     DO n=2, nod2D
        READ(infile_nod, *) nq, x1, x2, ind
        coord_nod2D(1,nq) = x1
        coord_nod2D(2,nq) = x2
        index_nod2D(nq)   = ind
     END DO
  else
     READ(a_dummy,*) coord_nod2D(1:2,1),  index_nod2D(1)
     DO n=2, nod2D
        READ(infile_nod, *) coord_nod2D(1:2,n),  index_nod2D(n)
     end DO
  end if
  CLOSE(infile_nod)

  ! Determine Bounding Box for later use

  BoundingBox_xmin = coord_nod2D(1,1)
  BoundingBox_xmax = coord_nod2D(1,1)
  BoundingBox_ymin = coord_nod2D(2,1)
  BoundingBox_ymax = coord_nod2D(2,1)
!$OMP PARALLEL DO &
!$OMP       reduction(min: BoundingBox_xmin, BoundingBox_ymin) &  
!$OMP       reduction(max: BoundingBox_xmax, BoundingBox_ymax)  
  do n=2,nod2D
     BoundingBox_xmin = min(BoundingBox_xmin,coord_nod2D(1,n))
     BoundingBox_xmax = max(BoundingBox_xmax,coord_nod2D(1,n))
     BoundingBox_ymin = min(BoundingBox_ymin,coord_nod2D(2,n))
     BoundingBox_ymax = max(BoundingBox_ymax,coord_nod2D(2,n))
  end do
!$OMP END PARALLEL DO  
  !
  !  READS THE NODE-NUMBERS OF EACH ELEMENT 
  !  
  READ(infile_elem, *)  elem2D      
  ALLOCATE(elem2D_nodes(3,elem2D))

  !NR For OpenMP, we need data locality
!$OMP PARALLEL DO  
!$  DO n=1, elem2D
!$     elem2D_nodes(:,n) = 0
!$  END DO
!$OMP END PARALLEL DO

  read(infile_elem, *) elem2D_nodes(1:3,1:elem2D)
  CLOSE(infile_elem)
  !
  ALLOCATE(nodhn(nod2D), nodhn_init(nod2D))
 
!$OMP PARALLEL DO
!$  DO n=1,nod2D
!$     nodhn(n) = 0._wp
!$  END DO
!$OMP END PARALLEL DO
  
  READ(infile_depth, *) nodhn(1:nod2D)

  ! The initial topography is kept for reference, e.g. to determine land nodes
  ! before the rupture lifts or lowers nodes. The main purpose is to retain the
  ! original land-sea-mask for output.
  ! Furthermore, nodhn will be smoothed.
!$OMP PARALLEL DO
!$  DO n=1,nod2D
!$     nodhn_init(n) = real(nodhn(n), sp)
!$  END DO
!$OMP END PARALLEL DO
 

  CLOSE(infile_depth)

END SUBROUTINE read_2Dmesh

!--------- --------- --------- --------- --------- --------- ---------

SUBROUTINE set_up_bv_roughness

  IMPLICIT NONE

  integer                          :: edg, nd1, nd2, n
  REAL(kind=wp), DIMENSION(nod2D)  :: rgh_nod2D
  integer                          :: infile_rgh

  OPEN(newunit=infile_rgh,file= TRIM(MeshPath) // 'rgh.out',  status= 'old')
  DO n=1, nod2D
     read(infile_rgh,*) rgh_nod2D(n)
  END DO
  CLOSE(infile_rgh)     
  

  ALLOCATE(rgh_edg(edg2D))
  
!$OMP PARALLEL DO 
  do edg=1,edg2D

     nd1=edge_nod2D(1,edg)
     nd2=edge_nod2D(2,edg)

     if ( nodhn(nd1)<150._wp .OR. nodhn(nd2)<150._wp ) then
        rgh_edg(edg)=0. 
     else
        rgh_edg(edg)=0.5_wp*(rgh_nod2D(nd1)+rgh_nod2D(nd2))
     endif

  end do
!$OMP END PARALLEL DO

END SUBROUTINE set_up_bv_roughness



subroutine set_up_mng_roughness

  implicit none

  integer                         :: i, n, nmb
  integer                         :: edg
  real(kind=wp), dimension(nod2D) :: mng_aux, mng_nod2D
  integer                         :: infile_mng

  print *,'Setting up Manning roughness'

!SH Manning roughness field

    open(newunit=infile_mng, file= trim(MeshPath)//trim(BottomFrictionFile), status= 'old')
    read(infile_mng,*) mng_nod2D(1:nod2D)
    
    close(infile_mng)

if (mng_max > 0.) then
!$OMP PARALLEL DO
  do n=1,nod2D
     mng_nod2D(n) = min(max(mng_nod2D(n),mng_min),  mng_max) 
  end do
!$OMP END PARALLEL DO
else  
!$OMP PARALLEL DO
  do n=1,nod2D
     mng_nod2D(n) = max(mng_nod2D(n),mng_min) 
  end do
!$OMP END PARALLEL DO
endif

  print *,'Before smoothing:'
  print *,'Minimal Manning roughness : ',minval(mng_nod2D)
  print *,'Maximal Manning roughness : ',maxval(mng_nod2D)

  do i = 1,2  
!$OMP PARALLEL 
!$OMP DO  private(nmb)
    do n = 1, nod2D
      nmb = nghbr_nod2D(n)%nmb
      mng_aux(n) = sum(mng_nod2D(nghbr_nod2d(n)%addr(:))) / real(nmb,wp)
    end do
!$OMP END DO 
!$OMP DO 
    do n = 1, nod2D
       mng_nod2D(n)=mng_aux(n)
    end do
!$OMP END DO 
!$OMP END PARALLEL
  end do

  print *,'After smoothing:'
  print *,'Minimal Manning roughness : ',minval(mng_nod2D)
  print *,'Maximal Manning roughness : ',maxval(mng_nod2D)
  

  allocate(mng_edg2D(edg2D))

!$OMP PARALLEL DO 
  do edg=1,edg2D
    mng_edg2D(edg)=0.5_wp*sum(mng_nod2D(edge_nod2D(1:2,edg)))
  end do  
!$OMP END PARALLEL DO

  print *,'Values in edges :'
  print *,'Minimal Manning roughness : ',minval(mng_edg2D)
  print *,'Maximal Manning roughness : ',maxval(mng_edg2D)



end subroutine set_up_mng_roughness


SUBROUTINE build_nghbr_arrays()
  !
  !use ELEMENTS
  !
  IMPLICIT NONE
  
  INTEGER                 :: j,k,m,a,b,c, count, el,ml(1)
  INTEGER, DIMENSION(100) :: aux=0
  INTEGER                 :: ind(nod2D)
  integer                 :: n, iA, iZ
  integer                 :: outfile_nghb
  !
  !--------------- 2D mesh:
  ! Builds sparse storage for patches: nodelem_ptr, nodelem_addr
  !   

  ind(:) = 0
  DO j = 1, elem2D
     ind(elem2D_nodes(1:3,j)) = ind(elem2D_nodes(1:3,j)) + 1
  END DO

  allocate(nodelem_ptr(nod2D+1))
! first touch initialization for data locality in OpenMP
!$OMP PARALLEL DO 
!$  DO n = 1, nod2D+1
!$     nodelem_ptr(n) = 0
!$  END DO
!$OMP END PARALLEL DO

  nodelem_ptr(1) = 1

  DO n = 1, nod2D
     nodelem_ptr(n+1) = nodelem_ptr(n) + ind(n)
  END DO

  allocate(nodelem_addr(nodelem_ptr(nod2D+1)-1))

! first touch initialization for data locality in OpenMP 
!$OMP PARALLEL DO  PRIVATE(n,j) 
  DO n = 1, nod2D  
     do j = nodelem_ptr(n), nodelem_ptr(n+1)-1
        nodelem_addr(j) = 0
     enddo
  END DO
!$OMP END PARALLEL DO      

  DO el = 1, elem2D 
     DO k = 1, 3  
        a = elem2D_nodes(k,el)
        jLoop: do j = nodelem_ptr(a), nodelem_ptr(a+1)-1
           if (nodelem_addr(j) == 0 ) then
              nodelem_addr(j) = el
              exit jLoop
           endif
        enddo jLoop
     end DO
  END DO


  ! The list of elements is ordered, and no sorting is needed
  !
 
  ALLOCATE(nghbr_nod2D(nod2D))
  
  DO n = 1, nod2D
     count = 1
     aux(1) = n
     DO m = nodelem_ptr(n), nodelem_ptr(n+1)-1 
        el = nodelem_addr(m)
        DO k = 1, 3
           a = elem2D_nodes(k,el)       
           IF(a/=n .and. all(aux(1:count) /= a)) THEN  
              count = count + 1         
              aux(count) = a
           END IF
        END DO
     END DO
     nghbr_nod2D(n)%nmb = count
     ALLOCATE(nghbr_nod2D(n)%addr(count))
     
     ! we need to sort array aux(1:count)
     DO m = count, 1, -1
        ml = MAXLOC(aux(1:count))
        nghbr_nod2D(n)%addr(m) = aux(ml(1))
        aux(ml(1)) = -999
     END DO
  END DO

    ! perform auxiliary output
    if (write_ascii_neighbourhood) then
       open(newunit=outfile_nghb, file="neighbourhood.out")
       do n=1,nod2D
          WRITE(outfile_nghb,*) nghbr_nod2D(n)%nmb
          do j=1,nghbr_nod2D(n)%nmb
             WRITE(outfile_nghb,*) nghbr_nod2D(n)%addr(j)
          end do
       end do
       close(outfile_nghb)
    end if


END SUBROUTINE build_nghbr_arrays


SUBROUTINE build_edg_arrays
  ! build_edg_arrays
  !OUTPUT:
  !edg2D
  !edge_nod2D     (2, edg2D)                       
  !edge_in_elem2D (2, edg2D)
  !edg_nxy        (2, edg2D)
  !edg_length     (   edg2D)
  !elem2D_edges   (3, elem2D)
  !--------------------------------------------------------------------------
  !edg2D
  !the total number of edges
  !
  !edge_nod2D(:, edg) 
  !keeps two nodes which form the edge edg. edge_nod2D(edg,2)>edge_nod2D(edg, 1). 
  !
  !elem2D_edges(:, elem)
  !keeps the  edges which belong to the elem. (/1, 2/), (/2, 3/), (/3, 1/) are 
  !the nodes (elem2D_nodes(:, elem)) of a triangle elem which form these edges.
  !
  !edge_in_elem2D(:, edg)
  !keeps two triangle which have the common edge edg. If it is only one triangle 
  !the second value will be zero. edge_in_elem2D(1, edg) is on the left side of 
  !a direction (edge_nod2D(1, edg) edge_nod2D(2, edg)). 
  !
  !edg_nxy(:, edg)
  !is the normal (nx, ny) to the edge always points from the left triangle to the 
  !right one. If the edge belongs to the open boundary than it will be the outer 
  !normal to the edge.
  !
  !edg_length(edg)
  ! length of the edge [m]
  !
  IMPLICIT NONE
  !
  INTEGER                                           :: n, el, nmb, cnt
  INTEGER                                           :: elnmb, nghbrnmb
  INTEGER, ALLOCATABLE, DIMENSION(:)                :: nghbrnodes
  INTEGER                                           :: i, j, k, elnodes(3), iA, iZ, edg
  REAL(kind=wp)                                     :: vec1(3), vec2(3), nvec(3), r
  
  edg2D = 0
  cnt = 1
  vec1 = 0.0
  vec2 = 0.0
  nvec = 0.0
  DO n = 1, nod2D
     nghbrnmb = nghbr_nod2D(n)%nmb
     edg2D = edg2D + COUNT(nghbr_nod2D(n)%addr > n)
  END DO
  
  if (verbosity >= 5) then 
     print *,'Number of nodes    in 2D mesh: ', nod2D
     print *,'          elements in 2D mesh: ', elem2D
     print *,'          edges    in 2D mesh: ', edg2D
  endif

  ALLOCATE(edge_nod2D(2,edg2D))
  ALLOCATE(edge_in_elem2D(2,edg2D))
  ALLOCATE(elem2D_edges(3,elem2D))
  ALLOCATE(edg_nxy(2,edg2D))
  ALLOCATE(edg_length(edg2D))
!$OMP PARALLEL 
!$OMP DO 
  do n=1,  edg2D
     edge_nod2D(:,n)     = 0
     edge_in_elem2D(:,n) = 0
  enddo
!$OMP END DO NOWAIT
!$OMP DO 
  do n=1, elem2D
     elem2D_edges(:,n)   = 0
  enddo
!$OMP END DO
!$OMP END PARALLEL 

  DO n = 1, nod2D
     nghbrnmb = nghbr_nod2D(n)%nmb
     ALLOCATE(nghbrnodes(nghbrnmb))
     nghbrnodes = nghbr_nod2D(n)%addr
     DO i = 1, nghbrnmb
        IF(nghbrnodes(i) <= n) CYCLE
        edge_nod2D(:,cnt) = (/ n, nghbrnodes(i) /)
        vec1(1:2) = coord_nod2D(:,nghbrnodes(i)) - coord_nod2D(:,n)
        
       if (coordinate_type == 1) then
          !correct for global grids for jumps across -180/180
          if (abs(vec1(1)) > .5_wp*pi) vec1(1) = -sign(1._wp,vec1(1))*(2._wp*pi-abs(vec1(1)))
          ! and scale according to latitude
          vec1(1) = vec1(1) * .5_WP*(cos_nod(n)+cos_nod(nghbrnodes(i)))
       endif


        DO j = nodelem_ptr(n),nodelem_ptr(n+1)-1 
           el = nodelem_addr(j) 
           elnodes = elem2D_nodes(:,el)
           !
           !Finds whether edge cnt belongs to the element el or not
           !
           IF(COUNT((elnodes == n).OR.(elnodes == nghbrnodes(i))) == 2) THEN
              DO k = 1, 3
                 IF((elnodes(k) /= n).AND.(elnodes(k) /= nghbrnodes(i))) THEN
                    vec2(1:2) = coord_nod2D(:,elnodes(k)) - coord_nod2D(:,n)
                    EXIT
                 END IF
              END DO
              nvec(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
              !
              !Finds whether element el is on the left or on the right of vec1
              !
              IF(nvec(3) > 0) THEN
                 edge_in_elem2D(1,cnt) = el
              ELSE
                 edge_in_elem2D(2,cnt) = el
              END IF
              !
              !Finds which nodes of the triangle form the edge cnt
              !
              IF(COUNT((elnodes(1:2) == n).OR.(elnodes(1:2) == nghbrnodes(i))) == 2) THEN
                 elem2D_edges(1, el)=cnt
              ELSEIF(COUNT((elnodes(2:3) == n).OR.(elnodes(2:3) == nghbrnodes(i))) == 2) THEN
                 elem2D_edges(2, el) = cnt
              ELSE
                 elem2D_edges(3, el) = cnt
              END IF
           END IF
        END DO

        
        if (coordinate_type == 1) then
           !correct for global grids for jumps across -180/180
           if (abs(vec2(1)) > .5_wp*pi) vec2(1) = -sign(1._wp,vec2(1))*(2._wp*pi-abs(vec2(1)))
           ! and scale according to latitude
           vec2(1) = vec2(1) * .5_WP*(cos_nod(n)+cos_nod(nghbrnodes(i)))
        endif
        !
        !Now construct the normal to the edge
        !
        edg_nxy(1,cnt) =  vec1(2)
        edg_nxy(2,cnt) = -vec1(1)      
        edg_length(cnt) = SQRT(SUM(vec1*vec1)) * scaling
        r = SQRT(SUM(edg_nxy(1:2,cnt)*edg_nxy(1:2,cnt)))
        edg_nxy(1:2,cnt) = edg_nxy(1:2,cnt) / r 


        IF((MINVAL(edge_in_elem2D(:,cnt)) == 0).AND.(nvec(3) < 0)) THEN
           edg_nxy(1,cnt) = -edg_nxy(1,cnt)
           edg_nxy(2,cnt) = -edg_nxy(2,cnt)
           edge_in_elem2D(1,cnt) = edge_in_elem2D(2,cnt)
           edge_in_elem2D(2,cnt) = 0
	   edge_nod2D(:,cnt) = (/ edge_nod2D(2,cnt), edge_nod2D(1,cnt) /)
        END IF
        !
        !Come to the next edge
        !
        cnt = cnt + 1
     END DO
     DEALLOCATE(nghbrnodes)
  END DO


  shortestEdge=minval(edg_length)
  longestEdge=maxval(edg_length)

  if (allocated(cos_nod)) deallocate(cos_nod)

! Build compact storage: node with surrounding edges 
  allocate(nodedge_ptr(nod2D+1))
! first touch initialization for data locality in OpenMP
!$OMP PARALLEL DO 
!$  DO n = 1, nod2D+1
!$     nodedge_ptr(n) = 0
!$  END DO
!$OMP END PARALLEL DO

  nodedge_ptr(1) = 1
  do n=1,nod2D
    nodedge_ptr(n+1) = nodedge_ptr(n) +  nghbr_nod2D(n)%nmb -1
 enddo

 allocate(nodedge_addr(nodedge_ptr(nod2D+1)-1))

! first touch initialization for data locality in OpenMP 
!$OMP PARALLEL DO  PRIVATE(n,j) 
  DO n = 1, nod2D  
     do j = nodedge_ptr(n), nodedge_ptr(n+1)-1
        nodedge_addr(j) = 0
     enddo
  END DO
!$OMP END PARALLEL DO      


 do edg=1,edg2D

! insert edg as adjacent to it's 1st and 2nd node
    do j=1,2
       n = edge_nod2D(j,edg)
       iA = nodedge_ptr(n)
       iZ = nodedge_ptr(n+1)-1
       do i=iA,iZ
          if (nodedge_addr(i)==0) then
             nodedge_addr(i) = edg
             exit
          endif
       enddo
    enddo ! j=1,2

 end do ! edg


END SUBROUTINE build_edg_arrays

SUBROUTINE smooth_topography
  
  INTEGER                 :: i, n
  REAL(kind=wp)           :: depval, topo_min
  REAL(kind=wp)           :: nodhn_min, nodhn_max  


  !============== PREPARING TOPOGRAPHY =========


  !setting topography in coastal nodes to zero
  i = 0
  nodhn_min = 99999.
  nodhn_max =-99999.
!$OMP PARALLEL DO  REDUCTION(+:i) &
!$OMP             REDUCTION(MIN:nodhn_min) &
!$OMP             REDUCTION(MAX:nodhn_max) 
  DO n = 1, nod2D
     IF ( index_nod2D(n) == 3 .AND. nodhn(n) /= 0.0 ) THEN
        nodhn(n) = 0.0
        i = i + 1
     END IF
     nodhn_min = min(nodhn_min,nodhn(n))
     nodhn_max = max(nodhn_max,nodhn(n))
  END DO
!$OMP END PARALLEL DO

  PRINT *,'Set topograhy in',i,' coastal nodes to Zero'

  ! smoothing	the topography 

  topo_min=-.1_wp   ! coastal land nodes should remain _land_ 
  if (enable_benchmark) topo_min=0._wp

  if (nmb_iter_smooth_topo>0) then
     allocate(nodhn_aux(nod2D))
     print *,'Topography before smoothing: ',nodhn_min, nodhn_max

     do i = 1, nmb_iter_smooth_topo
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,depval)
!$OMP DO  
        do n = 1, nod2D

           !Check for (seldom) case of degenerated grid: node without neighbours
           if (nghbr_nod2D(n)%nmb > 0 ) then
              depval = sum(nodhn(nghbr_nod2d(n)%addr(:))) / real(nghbr_nod2D(n)%nmb,wp)

              if ( all(nodhn(nghbr_nod2d(n)%addr(:)) >= 0.) .OR.   &
                   all(nodhn(nghbr_nod2d(n)%addr(:)) <= 0.) ) then

                 nodhn_aux(n) = depval

              else
                 if (nodhn(n) >= 0.) then
                    nodhn_aux(n) = max(depval,0._wp)
                 else
                    nodhn_aux(n) = min(depval,topo_min)
                 end if
              end if

           else
              nodhn_aux(n) = nodhn(n)
           endif

        end do
!$OMP END DO
!$OMP DO 
        do n = 1, nod2D
           nodhn(n) = nodhn_aux(n)
        end do
!$OMP END DO          
!$OMP END PARALLEL
     end do
     deallocate(nodhn_aux)


     !setting topography in coastal nodes *AGAIN* to zero
     i = 0
     nodhn_min = 99999.
     nodhn_max =-99999.
!$OMP PARALLEL DO REDUCTION(MIN:nodhn_min) &
!$OMP             REDUCTION(MAX:nodhn_max) &
!$OMP             REDUCTION(+:i) 
     do n = 1, nod2D
        if ( index_nod2D(n) == 3 .and. nodhn(n) /= 0.0 ) then
           nodhn(n) = 0.0
           i = i + 1
        end if
        nodhn_min = min(nodhn_min,nodhn(n))
        nodhn_max = max(nodhn_max,nodhn(n))
     end do
!$OMP END PARALLEL DO
     print *,'topography after smoothing : ',nodhn_min, nodhn_max
     PRINT *,'Set topograhy in',i,' coastal nodes to Zero again'

  end if

END SUBROUTINE smooth_topography


END MODULE MESH

