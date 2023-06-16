! Part of TsunAWI
! Copyright (c) 2014
! Alfred Wegener Institut,
! Helmholtz Centre for Polar and Marine Research
! Contact: rakowsky@awi.de


MODULE SWE

! Description:
!  This module computes the velocity field and the sea surface elevation. 
!  This module contains shallow water equations related subroutines
!

  use PARAMETERS
  implicit none

  REAL(kind=wp), ALLOCATABLE, DIMENSION(:)   :: ssh0, ssh1, ssh2
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:) :: grad
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:)   :: ssh_init
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:)   :: depth 
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:)   :: dhdt
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:) :: uv0, uv1, uv2
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:) :: uv_node, uv_node_old
  REAL(kind=sp), ALLOCATABLE, DIMENSION(:)   :: MAX_ABS_vel, MAX_flux
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:)   :: coriolis_param
  INTEGER,      ALLOCATABLE, DIMENSION(:)   :: i_wet

! arrival time, max wave height
  REAL(KIND=sp), ALLOCATABLE, DIMENSION(:)   :: arrival_time
  REAL(KIND=sp), ALLOCATABLE, DIMENSION(:)   :: ttt
  REAL(KIND=sp), ALLOCATABLE, DIMENSION(:)   :: mwh

! Energy
  REAL(kind=wp)                              :: Ekin, Epot  

  LOGICAL :: smooth_vel, save_vel_at_nodes	

! For GETM inundation scheme
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:)   :: mom_fact
  REAL(kind=wp) :: H_mult

! Some auxiliary fields that can be precomputed
  REAL(kind=wp), allocatable, dimension(:,:)   :: bafuvel_edgelem


CONTAINS

  !
  !----------------------------------------------------------------------------
  !               V E L O C I T Y
  !----------------------------------------------------------------------------
  !
! Computes the velocity field



subroutine compute_velocity_EXTRAPOLATE

  use MESH
  use ELEMENTS
  !
  implicit none
  integer          :: i, el, edg,  el1, el2
  real(kind=wp)     :: vol, r, uvabs
  real(kind=wp)     :: Fx1, Fx2, Fy1, Fy2
  real(kind=wp)     :: hu, tmp, tmpu, tmpv, volsum
  real(kind=wp)     :: dxy_uv(8,elem2D)
  real(kind=wp)     :: u_edg, v_edg, abs_vel

  real(kind=wp)     :: hu_4_3rdrt, u_rhs, v_rhs
  
  real(kind=wp)     :: A_vel, A_vel_max, mng, beta


!  momadv_type = 1 -------> Full P1 scheme for ADVECTION
!
!  momadv_type = 2 -------> P1--P_NC scheme for ADVECTION

  A_vel_max=-99999.  
  uvmax = 0.
  beta = 1._wp-2._wp*alpha
  
!$OMP   PARALLEL DEFAULT(SHARED) &
!$OMP&           PRIVATE(edg,vol,volsum,Fx1,Fy1,Fx2,Fy2,el1,el2,el,i, &
!$OMP&                  u_edg,v_edg,tmpu,tmpv,A_vel,u_rhs,v_rhs, &
!$OMP&                  hu,hu_4_3rdrt,r,abs_vel, mng) & 
!$OMP&         REDUCTION(MAX:A_vel_max,uvmax)

  !Apply filter for computational mode in Leap frog scheme

!$OMP DO 
  do edg=1,edg2D
     uv0(1,edg)    = beta*uv1(1,edg) + alpha*(uv2(1,edg) + uv0(1,edg))
     uv0(2,edg)    = beta*uv1(2,edg) + alpha*(uv2(2,edg) + uv0(2,edg))
     uv1(1,edg)    = uv2(1,edg)
     uv1(2,edg)    = uv2(2,edg)
     uv2(1,edg)    = 0._wp
     uv2(2,edg)    = 0._wp
  end do
!$OMP END DO

  if (smooth_vel) then
!$OMP MASTER
     print *,'SMOOTHING VELOCITY FIELD'
!$OMP END MASTER
!$OMP DO  
     do edg=1,edg2D
        uv2(1,edg) = 0.5*sum(uv_node_old(1,edge_nod2D(:,edg)))
        uv2(2,edg) = 0.5*sum(uv_node_old(2,edge_nod2D(:,edg)))
     end do
!$OMP END DO

  else

     if (viscosity_type==2) then
!$OMP DO 
        do el =1, elem2D
           if ( i_wet(el)>0 ) then 
              
              ! derivatives in element 'el' derived from nodal projections
              ! bafuzeta

              dxy_uv(5,el) = uv_node(1,elem2D_nodes(1,el))*bafu_xy(1,el)  &
                        + uv_node(1,elem2D_nodes(2,el))*bafu_xy(2,el)  & 
                        + uv_node(1,elem2D_nodes(3,el))*bafu_xy(3,el)   
           
              dxy_uv(6,el) = uv_node(1,elem2D_nodes(1,el))*bafu_xy(4,el)  &
                        + uv_node(1,elem2D_nodes(2,el))*bafu_xy(5,el)  & 
                        + uv_node(1,elem2D_nodes(3,el))*bafu_xy(6,el) 
                        
              dxy_uv(7,el) = uv_node(2,elem2D_nodes(1,el))*bafu_xy(1,el)  &
                        + uv_node(2,elem2D_nodes(2,el))*bafu_xy(2,el)  & 
                        + uv_node(2,elem2D_nodes(3,el))*bafu_xy(3,el) 
           
              dxy_uv(8,el) = uv_node(2,elem2D_nodes(1,el))*bafu_xy(4,el)  &
                        + uv_node(2,elem2D_nodes(2,el))*bafu_xy(5,el)  & 
                        + uv_node(2,elem2D_nodes(3,el))*bafu_xy(6,el) 

              A_vel = 3._wp*smag_fact * &
                   sqrt(dxy_uv(5,el)**2 + 0.5_wp*(dxy_uv(6,el) + dxy_uv(7,el))**2 + dxy_uv(8,el)**2)

              if (verbosity >= 5) A_vel_max = max(A_vel_max, A_vel)

              ! derivatives in element 'el' 
              ! bafu_vel = -2* permutation (rotation) of bafu_xy by construction
              dxy_uv(1,el) = (uv0(1,elem2D_edges(1,el))*bafu_xy(3,el)  &
                         + uv0(1,elem2D_edges(2,el))*bafu_xy(1,el)  & 
                         + uv0(1,elem2D_edges(3,el))*bafu_xy(2,el)) * A_vel *(-2._wp) 
           
              dxy_uv(2,el) = (uv0(1,elem2D_edges(1,el))*bafu_xy(6,el)  &
                         + uv0(1,elem2D_edges(2,el))*bafu_xy(4,el)  & 
                         + uv0(1,elem2D_edges(3,el))*bafu_xy(5,el)) * A_vel *(-2._wp)       
                        
              dxy_uv(3,el) = (uv0(2,elem2D_edges(1,el))*bafu_xy(3,el)  &
                         + uv0(2,elem2D_edges(2,el))*bafu_xy(1,el)  & 
                         + uv0(2,elem2D_edges(3,el))*bafu_xy(2,el)) * A_vel *(-2._wp)      
           
              dxy_uv(4,el) = (uv0(2,elem2D_edges(1,el))*bafu_xy(6,el)  &
                         + uv0(2,elem2D_edges(2,el))*bafu_xy(4,el)  & 
                         + uv0(2,elem2D_edges(3,el))*bafu_xy(5,el)) * A_vel *(-2._wp)       
           end if
        end do
!$OMP END DO 
     elseif (momadv_type > 0) then

        A_vel = 3._wp*Ah0  /(2._wp*dt)
!$OMP DO 
        do el =1, elem2D
           if ( i_wet(el) > 0 ) then 
        

              ! derivatives in element 'el' 
              ! bafu_vel = -2* permutation (rotation) of bafu_xy by construction
              dxy_uv(1,el) = (uv0(1,elem2D_edges(1,el))*bafu_xy(3,el)  &
                         + uv0(1,elem2D_edges(2,el))*bafu_xy(1,el)  & 
                         + uv0(1,elem2D_edges(3,el))*bafu_xy(2,el)) * A_vel *(-2._wp)   
           
              dxy_uv(2,el) = (uv0(1,elem2D_edges(1,el))*bafu_xy(6,el)  &
                         + uv0(1,elem2D_edges(2,el))*bafu_xy(4,el)  & 
                         + uv0(1,elem2D_edges(3,el))*bafu_xy(5,el)) * A_vel *(-2._wp)    
                        
              dxy_uv(3,el) = (uv0(2,elem2D_edges(1,el))*bafu_xy(3,el)  &
                         + uv0(2,elem2D_edges(2,el))*bafu_xy(1,el)  & 
                         + uv0(2,elem2D_edges(3,el))*bafu_xy(2,el)) * A_vel  *(-2._wp)  
           
              dxy_uv(4,el) = (uv0(2,elem2D_edges(1,el))*bafu_xy(6,el)  &
                         + uv0(2,elem2D_edges(2,el))*bafu_xy(4,el)  & 
                         + uv0(2,elem2D_edges(3,el))*bafu_xy(5,el)) * A_vel *(-2._wp)   

              ! derivatives in element 'el' derived from nodal projections

              dxy_uv(5,el) = uv_node(1,elem2D_nodes(1,el))*bafu_xy(1,el)  &
                        + uv_node(1,elem2D_nodes(2,el))*bafu_xy(2,el)  & 
                        + uv_node(1,elem2D_nodes(3,el))*bafu_xy(3,el)   
           
              dxy_uv(6,el) = uv_node(1,elem2D_nodes(1,el))*bafu_xy(4,el)  &
                        + uv_node(1,elem2D_nodes(2,el))*bafu_xy(5,el)  & 
                        + uv_node(1,elem2D_nodes(3,el))*bafu_xy(6,el) 
                        
              dxy_uv(7,el) = uv_node(2,elem2D_nodes(1,el))*bafu_xy(1,el)  &
                        + uv_node(2,elem2D_nodes(2,el))*bafu_xy(2,el)  & 
                        + uv_node(2,elem2D_nodes(3,el))*bafu_xy(3,el) 
           
              dxy_uv(8,el) = uv_node(2,elem2D_nodes(1,el))*bafu_xy(4,el)  &
                        + uv_node(2,elem2D_nodes(2,el))*bafu_xy(5,el)  & 
                        + uv_node(2,elem2D_nodes(3,el))*bafu_xy(6,el) 

           end if
        end do
!$OMP END DO 
     else
        
        A_vel = 3._wp*Ah0  /(2._wp*dt)
!$OMP DO 
        do el =1, elem2D
           if ( i_wet(el)>0 ) then 

              ! derivatives in element 'el' 
              ! bafu_vel = -2* permutation (rotation) of bafu_xy by construction
              dxy_uv(1,el) = (uv0(1,elem2D_edges(1,el))*bafu_xy(3,el)  &
                         + uv0(1,elem2D_edges(2,el))*bafu_xy(1,el)  & 
                         + uv0(1,elem2D_edges(3,el))*bafu_xy(2,el)) * A_vel  *(-2._wp)   
           
              dxy_uv(2,el) = (uv0(1,elem2D_edges(1,el))*bafu_xy(6,el)  &
                         + uv0(1,elem2D_edges(2,el))*bafu_xy(4,el)  & 
                         + uv0(1,elem2D_edges(3,el))*bafu_xy(5,el)) * A_vel *(-2._wp)     
                        
              dxy_uv(3,el) = (uv0(2,elem2D_edges(1,el))*bafu_xy(3,el)  &
                         + uv0(2,elem2D_edges(2,el))*bafu_xy(1,el)  & 
                         + uv0(2,elem2D_edges(3,el))*bafu_xy(2,el)) * A_vel *(-2._wp)    
           
              dxy_uv(4,el) = (uv0(2,elem2D_edges(1,el))*bafu_xy(6,el)  &
                         + uv0(2,elem2D_edges(2,el))*bafu_xy(4,el)  & 
                         + uv0(2,elem2D_edges(3,el))*bafu_xy(5,el)) * A_vel *(-2._wp)     
           endif
        end do
!$OMP END DO
     end if

! most important case: handle ifs outside the loop
  if ((momadv_type == 1) .and. (.not. enable_variable_bottom_friction) &
                         .and. (.not. use_bv_roughness)) then
!$OMP DO  
     do edg=1,edg2D

        if ( i_wet(edge_in_elem2D(1,edg)) > 0  .or.  &
             i_wet(edge_in_elem2D(2,edg)) > 0 )then

           ! Full P1 scheme for ADVECTION
           u_edg = 0.5_wp*(uv_node(1,edge_nod2D(1,edg)) + uv_node(1,edge_nod2D(2,edg)))
           v_edg = 0.5_wp*(uv_node(2,edge_nod2D(1,edg)) + uv_node(2,edge_nod2D(2,edg)))
            
           abs_vel = sqrt(uv1(1,edg)**2 + uv1(2,edg)**2)  

           el1=edge_in_elem2D(1,edg)
           el2=edge_in_elem2D(2,edg)

           uvmax = max(uvmax,abs_vel)

           !Bottom friction (Manning)
                 
           r = 1._wp/(1._wp + 2._wp*dt * g * Cd*Cd * abs_vel * &
                max(0.5_wp*(depth(edge_nod2D(1,edg))+depth(edge_nod2D(2,edg))),Dcr)**(-4._wp/3._wp))
           
          
           if (i_wet(el1) > 0 .and. i_wet(el2) > 0) then

              volsum = 1._wp/ (voltriangle(el1) + voltriangle(el2))

              u_rhs = voltriangle(el1) * ( g*grad(1,el1)       & ! Gradient
                  + dxy_uv(1,el1) * bafuvel_edgelem(1,edg)     & ! Viscosity
                  + dxy_uv(2,el1) * bafuvel_edgelem(2,edg)     & ! Viscosity
                  + dxy_uv(5,el1)*u_edg + dxy_uv(6,el1)*v_edg) & ! Advection 
                    +                                          &
                      voltriangle(el2) * ( g*grad(1,el2)       &
                  + dxy_uv(1,el2) * bafuvel_edgelem(3,edg)     &
                  + dxy_uv(2,el2) * bafuvel_edgelem(4,edg)     &
                  + dxy_uv(5,el2)*u_edg + dxy_uv(6,el2)*v_edg)   

              v_rhs = voltriangle(el1) * ( g*grad(2,el1)       &
                  + dxy_uv(3,el1) * bafuvel_edgelem(1,edg)     &
                  + dxy_uv(4,el1) * bafuvel_edgelem(2,edg)     &
                  + dxy_uv(7,el1)*u_edg + dxy_uv(8,el1)*v_edg) &
                    +                                          & 
                      voltriangle(el2) * ( g*grad(2,el2)       &
                  + dxy_uv(3,el2) * bafuvel_edgelem(3,edg)     &
                  + dxy_uv(4,el2) * bafuvel_edgelem(4,edg)     &
                  + dxy_uv(7,el2)*u_edg + dxy_uv(8,el2)*v_edg )
                                  

              ! ("-" due to integration by parts)
              uv2(1,edg) = (uv0(1,edg) + 2._wp*dt*(-u_rhs*volsum + coriolis_param(edg) *uv1(2,edg))) * r
              uv2(2,edg) = (uv0(2,edg) + 2._wp*dt*(-v_rhs*volsum - coriolis_param(edg) *uv1(1,edg))) * r      

           elseif (i_wet(el2) > 0) then

              u_rhs = g*grad(1,el2)                          &  ! Gradient
                    + dxy_uv(1,el2) * bafuvel_edgelem(3,edg) &  ! Viscosity
                    + dxy_uv(2,el2) * bafuvel_edgelem(4,edg) &
                    + dxy_uv(5,el2)*u_edg + dxy_uv(6,el2)*v_edg ! Advection

              v_rhs = g*grad(2,el2) &
                    + dxy_uv(3,el2) * bafuvel_edgelem(3,edg) &
                    + dxy_uv(4,el2) * bafuvel_edgelem(4,edg) &
                    + dxy_uv(7,el2)*u_edg + dxy_uv(8,el2)*v_edg
               
              ! ("-" due to integration by parts)
              uv2(1,edg) = (uv0(1,edg) + 2._wp*dt*(-u_rhs + coriolis_param(edg) *uv1(2,edg))) * r
              uv2(2,edg) = (uv0(2,edg) + 2._wp*dt*(-v_rhs - coriolis_param(edg) *uv1(1,edg))) * r      

           else ! only i_wet(el1) > 0

              u_rhs =  g*grad(1,el1)                         &  ! Gradient
                    + dxy_uv(1,el1) * bafuvel_edgelem(1,edg) &  ! Viscosity
                    + dxy_uv(2,el1) * bafuvel_edgelem(2,edg) &
                    + dxy_uv(5,el1)*u_edg + dxy_uv(6,el1)*v_edg ! Advection

              v_rhs =  g*grad(2,el1) &
                    + dxy_uv(3,el1) * bafuvel_edgelem(1,edg) &
                    + dxy_uv(4,el1) * bafuvel_edgelem(2,edg) &
                    + dxy_uv(7,el1)*u_edg + dxy_uv(8,el1)*v_edg 
                               
           ! ("-" due to integration by parts)
              uv2(1,edg) = (uv0(1,edg) + 2._wp*dt*(- u_rhs + coriolis_param(edg) *uv1(2,edg))) * r
              uv2(2,edg) = (uv0(2,edg) + 2._wp*dt*(- v_rhs - coriolis_param(edg) *uv1(1,edg))) * r   
           
              !  Boundary Condition? Possible only here.
           
              if ( ( el2 == 0 .and. &
                   (index_nod2D(edge_nod2D(1,edg))==1 .or. index_nod2D(edge_nod2D(2,edg))==1 )) .or. &
                   (index_nod2D(edge_nod2D(1,edg))==4 .or. index_nod2D(edge_nod2D(2,edg))==4 ) ) then

                 tmpu = edg_nxy(2,edg) * (edg_nxy(2,edg) * uv2(1,edg) - edg_nxy(1,edg) * uv2(2,edg))
                 tmpv = edg_nxy(1,edg) * (edg_nxy(1,edg) * uv2(2,edg) - edg_nxy(2,edg) * uv2(1,edg))
                 uv2(1,edg) = tmpu
                 uv2(2,edg) = tmpv
                
              end if

           end if ! wet element? No boundary edge: el==0?!

        end if

     end do       ! (edg)
     
!$OMP END DO
else ! all other cases, ifs remain in loop body

!$OMP DO  
     do edg=1,edg2D

        if ( i_wet(edge_in_elem2D(1,edg)) > 0  .or.  &
             i_wet(edge_in_elem2D(2,edg)) > 0 )then

           if (momadv_type == 1) then      ! Full P1 scheme for ADVECTION
              u_edg = 0.5*(uv_node(1,edge_nod2D(1,edg)) + uv_node(1,edge_nod2D(2,edg)))
              v_edg = 0.5*(uv_node(2,edge_nod2D(1,edg)) + uv_node(2,edge_nod2D(2,edg)))
           elseif (momadv_type == 2) then  ! P1--P_NC scheme for ADVECTION
              u_edg = uv1(1,edg)
              v_edg = uv1(2,edg)
           endif
 
           el1=edge_in_elem2D(1,edg)
           el2=edge_in_elem2D(2,edg)

           !Bottom friction (Manning)
           mng = Cd
           if (enable_variable_bottom_friction)  mng = mng_edg2D(edg)
           
           abs_vel = sqrt(uv1(1,edg)**2 + uv1(2,edg)**2)  
           uvmax = max(uvmax,abs_vel)

           if (.not. use_bv_roughness) then  
              r = 1._wp/(1._wp + 2._wp*dt * g * mng*mng * abs_vel * &
                max(0.5_wp*(depth(edge_nod2D(1,edg))+depth(edge_nod2D(2,edg))),Dcr)**(-4._wp/3._wp))
           else
              hu = max(0.5_wp*(depth(edge_nod2D(1,edg))+depth(edge_nod2D(2,edg))),Dcr) 

              r = 1._wp/(1._wp + 2._wp*dt * g * mng*mng * abs_vel * (hu**(-4._wp/3._wp)) &
                   + 20._wp*dt *rgh_edg(edg) / hu)                       
           end if

           if (i_wet(el1) > 0 ) then

              volsum = voltriangle(el1)                 

              Fx1 = g*grad(1,el1) &  ! Gradient 
                  + dxy_uv(1,el1) * bafuvel_edgelem(1,edg) & ! Viscosity
                  + dxy_uv(2,el1) * bafuvel_edgelem(2,edg)

              Fy1 = g*grad(2,el1) &
                  + dxy_uv(3,el1) * bafuvel_edgelem(1,edg) &
                  + dxy_uv(4,el1) * bafuvel_edgelem(2,edg)

              
              if (momadv_type>0) then     !Advection
                 Fx1 = Fx1 + dxy_uv(5,el1)*u_edg + dxy_uv(6,el1)*v_edg 
                 Fy1 = Fy1 + dxy_uv(7,el1)*u_edg + dxy_uv(8,el1)*v_edg 
              endif

              Fx1 = Fx1*voltriangle(el1)
              Fy1 = Fy1*voltriangle(el1)
           else
              volsum = 0._wp
              Fx1    = 0._wp
              Fy1    = 0._wp
           endif
           if (i_wet(el2) > 0 ) then

              volsum = volsum + voltriangle(el2)

              Fx2 = g*grad(1,el2) &
                  + dxy_uv(1,el2) * bafuvel_edgelem(3,edg) &
                  + dxy_uv(2,el2) * bafuvel_edgelem(4,edg)

              Fy2 = g*grad(2,el2) &
                  + dxy_uv(3,el2) * bafuvel_edgelem(3,edg) &
                  + dxy_uv(4,el2) * bafuvel_edgelem(4,edg)
                 
              !Advection
              if (momadv_type>0) then
                 Fx2 = Fx2 + dxy_uv(5,el2)*u_edg + dxy_uv(6,el2)*v_edg 
                 Fy2 = Fy2 + dxy_uv(7,el2)*u_edg + dxy_uv(8,el2)*v_edg 
              end if
              
              Fx2 = Fx2*voltriangle(el2)
              Fy2 = Fy2*voltriangle(el2)
           else
              Fx2    = 0._wp
              Fy2    = 0._wp                 
           endif
                 
           u_rhs = (Fx1 + Fx2)/volsum
           v_rhs = (Fy1 + Fy2)/volsum

           ! ("-" due to integration by parts)
           uv2(1,edg) = (uv0(1,edg) + 2._wp*dt*(-u_rhs + coriolis_param(edg) *uv1(2,edg))) * r
           uv2(2,edg) = (uv0(2,edg) + 2._wp*dt*(-v_rhs - coriolis_param(edg) *uv1(1,edg))) * r  

           !  Boundary Condition?
           
           if (el2 == 0 .and. &
              (index_nod2D(edge_nod2D(1,edg))==1 .or. index_nod2D(edge_nod2D(2,edg))==1)) then

              tmpu = edg_nxy(2,edg) * (edg_nxy(2,edg) * uv2(1,edg) - edg_nxy(1,edg) * uv2(2,edg))
              tmpv = edg_nxy(1,edg) * (edg_nxy(1,edg) * uv2(2,edg) - edg_nxy(2,edg) * uv2(1,edg))
              uv2(1,edg) = tmpu
              uv2(2,edg) = tmpv
           endif
                    
        end if
        
     end do       ! (edg)
     
!$OMP END DO

  end if
end if
 
!$OMP END PARALLEL
  
  if (verbosity>=5) then
     print *,'Compute_velocity - max. viscosity  : ',A_vel_max
  end if


end subroutine compute_velocity_EXTRAPOLATE



! !ROUTINE: SUBROUTINE compute_velocity_GETM
!
! !INTERFACE:
!
  !
  !----------------------------------------------------------------------------
  !               V E L O C I T Y
  !----------------------------------------------------------------------------
  !
! Computes the velocity field

subroutine compute_velocity_GETM

  use MESH
  use ELEMENTS
  !
  implicit none
  integer          :: el, edg,  elin
  real(kind=wp)    :: vol, r, vol2, uvabs
  real(kind=wp)    :: Fx, Fy , ax, ay
  real(kind=wp)    :: hu, nx, ny, tmp, tmpu, tmpv, volsum
  integer          :: nod1, nod2, elnod(3), edges(3)
  real(kind=wp)    :: dx_u_nod, dy_u_nod, dx_v_nod, dy_v_nod
  real(kind=wp)    :: u_edg, v_edg, abs_vel

  real(kind=wp)    :: hu_4_3rdrt
  
  real(kind=wp)    :: A_vel, A_vel_max, mng
  real(kind=wp)    :: mom_fact_edg



!$  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
!$  real(kind=wp) :: omp_get_wtime, t1, t2, t3, t4


! NEW
!
!  momadv_type = 1 -------> Full P1 scheme for ADVECTION
!
!  momadv_type = 2 -------> P1--P_NC scheme for ADVECTION
!

! END NEW

  A_vel_max=-99999.  
  


!$OMP   PARALLEL DEFAULT(SHARED) &
!$OMP&           PRIVATE(edg,vol,volsum,Fx,Fy,elin,el, ax, ay, &
!$OMP&                  dx_u_nod,dy_u_nod,dx_v_nod,dy_v_nod,u_edg,v_edg,nx,ny,tmpu,tmpv,A_vel, &
!$OMP&                  hu,hu_4_3rdrt,r,abs_vel, nod1, nod2, edges, elnod, mng, mom_fact_edg) & 
!$OMP&         REDUCTION(MAX:A_vel_max,uvmax)

  !Apply filter for computational mode in Leap frog scheme

!$OMP DO  
  do edg=1,edg2D
     uv0(1,edg)    = uv1(1,edg) + alpha*(uv2(1,edg) -2.*uv1(1,edg) +uv0(1,edg))
     uv0(2,edg)    = uv1(2,edg) + alpha*(uv2(2,edg) -2.*uv1(2,edg) +uv0(2,edg))
     uv1(1,edg)    = uv2(1,edg)
     uv1(2,edg)    = uv2(2,edg)
  end do
!$OMP END DO

  uvmax = 0.

  if (smooth_vel) then
!$OMP MASTER
     print *,'SMOOTHING VELOCITY FIELD'
!$OMP END MASTER
!$OMP DO  
     do edg=1,edg2D

        uv2(1,edg) = 0.5*sum(uv_node_old(1,edge_nod2D(:,edg)))
        uv2(2,edg) = 0.5*sum(uv_node_old(2,edge_nod2D(:,edg)))

     end do
!$OMP END DO

  else


!$OMP DO  

     do edg=1,edg2D

        uv2(1,edg)    = 0.
        uv2(2,edg)    = 0.

        nod1 = edge_nod2D(1,edg)
        nod2 = edge_nod2D(2,edg)


        if ( depth(nod1) > small  .and. &
             depth(nod2) > small) then

           volsum=0.

           !SH
           mom_fact_edg = 0.5 * (mom_fact(nod1) + mom_fact(nod2))

           if (momadv_type == 1) then      ! Full P1 scheme for ADVECTION
              u_edg    = 0.5*(uv_node(1,nod1) + uv_node(1,nod2))
              v_edg    = 0.5*(uv_node(2,nod1) + uv_node(2,nod2))
           elseif (momadv_type == 2) then  ! P1--P_NC scheme for ADVECTION
              u_edg = uv1(1,edg)
              v_edg = uv1(2,edg)
           endif
 
           do elin=1, 2

              el=edge_in_elem2D(elin, edg)

              if (el>0) then

                 elnod = elem2D_nodes(:, el)

                 !if ( all(ssh2(elnod)+nodhn(elnod)>small) &
                 !    .and. max(abs(grad(1,el)),abs(grad(2,el)))>0. ) then
                 !if ( max(abs(grad(1,el)),abs(grad(2,el)))>0. ) then

                    vol=voltriangle(el)
                    volsum = volsum + vol

                    edges = elem2D_edges(:,el)

                    if (momadv_type>0 .or. viscosity_type==2) then
                       
                       ! derivatives in element 'el' derived from nodal projections
                       dx_u_nod = sum(uv_node(1,elnod(:))*bafu_xy(1:3,el))
                       dy_u_nod = sum(uv_node(1,elnod(:))*bafu_xy(4:6,el))
                       dx_v_nod = sum(uv_node(2,elnod(:))*bafu_xy(1:3,el))
                       dy_v_nod = sum(uv_node(2,elnod(:))*bafu_xy(4:6,el))
                       
                    end if

                    if (viscosity_type==2) then
                       ! Smagorinsky Diffusion
                 
                       A_vel = smag_fact *vol *sqrt(dx_u_nod**2. + 0.5*(dx_v_nod + dy_u_nod)**2. + dy_v_nod**2.)
                    else 
                       A_vel = Ah0 *vol/(2.*dt)
                    end if
              
                    if (verbosity >= 5) A_vel_max = max(A_vel_max, A_vel)
              
                    !Viscosity terms 
                    if (edges(1) == edg) then

                       ax = 3.*A_vel* bafu_xy(3,el)  
                       ay = 3.*A_vel* bafu_xy(6,el)
                                           
                    elseif (edges(2) == edg) then
                       ax = 3.*A_vel* bafu_xy(1,el)  
                       ay = 3.*A_vel* bafu_xy(4,el)

                    else  ! edges(3) = edg
                       ax = 3.*A_vel* bafu_xy(2,el)  
                       ay = 3.*A_vel* bafu_xy(5,el)
                    endif

                       Fx = ax * (uv0(1,edges(1))*bafu_xy(3,el) &
                                + uv0(1,edges(2))*bafu_xy(1,el) &
                                + uv0(1,edges(3))*bafu_xy(2,el))&
                          + ay * (uv0(1,edges(1))*bafu_xy(6,el) &
                                + uv0(1,edges(2))*bafu_xy(4,el) &
                                + uv0(1,edges(3))*bafu_xy(5,el))
                    
                       
                       Fy = ax * (uv0(2,edges(1))*bafu_xy(3,el) &
                                + uv0(2,edges(2))*bafu_xy(1,el) &
                                + uv0(2,edges(3))*bafu_xy(2,el))&
                          + ay * (uv0(2,edges(1))*bafu_xy(6,el) &
                                + uv0(2,edges(2))*bafu_xy(4,el) &
                                + uv0(2,edges(3))*bafu_xy(5,el))
                    !Advection
                    
                    if (momadv_type>0) then
                       Fx = Fx + (dx_u_nod*u_edg + dy_u_nod*v_edg )
                       Fy = Fy + (dx_v_nod*u_edg + dy_v_nod*v_edg )
                    end if
                 
                    !SH                 
                    !pressure gradient terms
                    !Fx = Fx +g*grad(1,el)
                    !Fy = Fy +g*grad(2,el)

                    !pressure gradient terms
                    Fx = mom_fact_edg*(Fx +g*grad(1,el))
                    Fy = mom_fact_edg*(Fy +g*grad(2,el))

                 
                    ! ("-" due to integration by parts)
                    uv2(1,edg) = uv2(1,edg) - Fx*vol
                    uv2(2,edg) = uv2(2,edg) - Fy*vol

                 !endif ! wet element 
              end if ! boundary edge: el==0?!
           end do   ! elin
     
           if (volsum > 0.) then


!GETM EXPERIMENT
              if (rotation_type==2) then
              uv2(1,edg) = uv0(1,edg) + 2.*dt*( uv2(1,edg)/volsum + mom_fact_edg*coriolis_param(edg) *uv1(2,edg))
              uv2(2,edg) = uv0(2,edg) + 2.*dt*( uv2(2,edg)/volsum - mom_fact_edg*coriolis_param(edg) *uv1(1,edg))       
!                 uv2(1,edg) = uv0(1,edg) + 2.*dt*mom_fact_edg*(uv2(1,edg)/volsum + coriolis_param(edg) *uv1(2,edg))
!                 uv2(2,edg) = uv0(2,edg) + 2.*dt*mom_fact_edg*(uv2(2,edg)/volsum - coriolis_param(edg) *uv1(1,edg))       
              else
                 uv2(1,edg) = uv0(1,edg) + 2.*dt*uv2(1,edg)/volsum 
                 uv2(2,edg) = uv0(2,edg) + 2.*dt*uv2(2,edg)/volsum 
!                 uv2(1,edg) = uv0(1,edg) + 2.*dt* mom_fact_edg*uv2(1,edg)/volsum 
!                 uv2(2,edg) = uv0(2,edg) + 2.*dt* mom_fact_edg*uv2(2,edg)/volsum 
              endif
                                
              if (enable_variable_bottom_friction) then  
                 mng = mng_edg2D(edg)
              else
                 mng = Cd
              end if


              !Bottom friction (Manning)
              hu=0.5* ( depth(nod1) + depth(nod2) )
              hu_4_3rdrt = exp(-(4._wp/3._wp)*log(max(hu,Dcr)))
              abs_vel = sqrt(uv1(1,edg)**2 + uv1(2,edg)**2)
              uvmax            = max(uvmax,abs_vel)

              if (.not. use_bv_roughness) then 

                 r = 1._wp/(1._wp + 2._wp *dt * mom_fact_edg * &
                      ( g * mng*mng * sqrt(uv1(1,edg)**2 + uv1(2,edg)**2) *hu_4_3rdrt))
              else
                 r = 1._wp/(1._wp + 2._wp *dt * mom_fact_edg * &
                      ( g * mng*mng * sqrt(uv1(1,edg)**2 + uv1(2,edg)**2) *hu_4_3rdrt + &
                      10.*rgh_edg(edg)/max(hu,Dcr) ))

!GETM EXPERIMENT
!              r = 1. + 2.*dt * g * mng*mng * abs_vel *hu_4_3rdrt
!                 r = r + 20.*dt*rgh_edg(edg)/max(hu,Dcr) 
              end if

              uv2(1,edg) = uv2(1,edg)*r
              uv2(2,edg) = uv2(2,edg)*r
             
              !  Boundary Condition
           
              if (edge_in_elem2D(2,edg) == 0) then 
        
                 nx = edg_nxy(1,edg)
                 ny = edg_nxy(2,edg)
            
                 tmpu = ny * (ny * uv2(1,edg) - nx * uv2(2,edg))
                 tmpv = nx * (nx * uv2(2,edg) - ny * uv2(1,edg))
                 uv2(1,edg) = tmpu
                 uv2(2,edg) = tmpv
                 
              end if

           end if
        end if
        
     end do       ! (edg)
     
!$OMP END DO

  end if
  
!$OMP END PARALLEL
  


  if (verbosity>=5) then
     print *,'Compute_velocity - max. viscosity  : ',A_vel_max
  end if


end subroutine compute_velocity_GETM



SUBROUTINE compute_ssh_EXTRAPOLATE

! Description:
!   This subroutine computes the sea surface elevation.


  USE MESH
  USE ELEMENTS
  !
  IMPLICIT NONE
  INTEGER       :: n, i, j, el, iA, iZ
  REAL(kind=wp)  :: vol, rhs_ssh, ssh_SIGN, ssh_extr
  REAL(kind=wp)  :: el_uv(2,elem2D)
  LOGICAL        :: l_wetnod(nod2D)

!$ real(kind=wp) :: omp_get_wtime, t1, t2

  REAL(kind=wp)  :: cnt, sshlimit, beta
  INTEGER :: nd

  ! Declaration of specific variables

  sshmax = -1000. 
  sshmin =  1000. 
  dmin   =  1000. 
  beta = 1._wp - 2._wp*alpha

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&         PRIVATE(vol,i,j,n,el,rhs_ssh, ssh_extr, &
!$OMP&                 ssh_sign,nd,cnt,sshlimit, iA, iZ) &
!$OMP&         REDUCTION(MIN:sshmin,dmin) &
!$OMP&         REDUCTION(MAX:sshmax)

!$OMP DO  
  do n = 1,nod2D
     ssh0(n) = beta*ssh1(n) + alpha*(ssh2(n) + ssh0(n))
     ssh1(n) = ssh2(n)
     l_wetnod(n) = (ssh0(n) + nodhn(n)) > small
     ssh2(n) = 0._wp
  end do
!$OMP END DO
!$OMP DO 
  do el = 1,elem2D

     if (l_wetnod(elem2D_nodes(1,el)) .or. l_wetnod(elem2D_nodes(2,el))  .or. &
                                           l_wetnod(elem2D_nodes(3,el)) ) then

        el_uv(1,el) = uv1(1,elem2D_edges(1,el)) * (depth(elem2D_nodes(1,el))+depth(elem2D_nodes(2,el))) &
                    + uv1(1,elem2D_edges(2,el)) * (depth(elem2D_nodes(2,el))+depth(elem2D_nodes(3,el))) &
                    + uv1(1,elem2D_edges(3,el)) * (depth(elem2D_nodes(1,el))+depth(elem2D_nodes(3,el)))

        el_uv(2,el) = uv1(2,elem2D_edges(1,el)) * (depth(elem2D_nodes(1,el))+depth(elem2D_nodes(2,el))) &
                    + uv1(2,elem2D_edges(2,el)) * (depth(elem2D_nodes(2,el))+depth(elem2D_nodes(3,el))) &
                    + uv1(2,elem2D_edges(3,el)) * (depth(elem2D_nodes(1,el))+depth(elem2D_nodes(3,el)))

        if (abs(el_uv(1,el)) > 0. .or. abs(el_uv(2,el)) > 0.) then
           ! Consider element in all computations
           i_wet(el)=2
        else
           ! Only relevant for extrapolation at dry nodes
           i_wet(el) = 1
        endif

     else
        ! Dry element, can be skipped in all computations
        i_wet(el) = 0
     endif
  enddo
!$OMP END DO

!$OMP DO  

  do n = 1,nod2D

     iA = nodelem_ptr(n)
     iZ = nodelem_ptr(n+1)-1

     if (l_wetnod(n)) then

        ! Result for ssh 

           ssh2(n) = ssh0(n) + dt * sum(nodelem_wght(1,iA:iZ) * el_uv(1,nodelem_addr(iA:iZ)) &
                                      + nodelem_wght(2,iA:iZ) * el_uv(2,nodelem_addr(iA:iZ)), &
                                       mask = (i_wet(nodelem_addr(iA:iZ))==2)) 

        ! update maximum wave height
        if (index_nod2D(n) /= 2) mwh(n) = max(mwh(n),real(ssh2(n) + min(nodhn(n),0.),sp))
        
     else

        cnt      = 0._wp
        sshlimit = 0._wp
        ssh_extr = 0._wp

        do i = iA, iZ
 
           if (i_wet(nodelem_addr(i)) > 0) then
              do j=1,3
                 nd = elem2D_nodes(j, nodelem_addr(i))

                 if (l_wetnod(nd)) then   
                                                           
                    ssh_extr = ssh_extr &
                         + grad(1,nodelem_addr(i))*(coord_nod2D(1,n)-coord_nod2D(1,nd)) &
                         + grad(2,nodelem_addr(i))*(coord_nod2D(2,n)-coord_nod2D(2,nd))

                    sshlimit = sshlimit+ssh0(nd)
                    cnt = cnt+1._wp
                 end if
              end do
           end if
        end do

        if (cnt > 0._wp) then
           ssh2(n) = (sshlimit + scaling*min(ssh_extr,0._wp)) /cnt

           ! update maximum wave height
           if (index_nod2D(n) /= 2)  mwh(n) = max(mwh(n),real(ssh2(n) + min(nodhn(n),0.),sp))

           if (verbosity>1) then
              if (ssh2(n)>50.) print *,'SSH2 TOO LARGE : ',n,el,&
                   ( (coord_nod2D(1,n)-coord_nod2D(1,nd))*grad(1,el)+ &
                   (coord_nod2D(2,n)-coord_nod2D(2,nd))*grad(2,el)  ) * scaling, &
                   ssh0(nd),ssh0(n),nodhn(n),ssh2(n)
              
              if (ssh2(n)<-30.) print *,'SSH2 TOO SMALL : ',n,el,&
                   ssh0(n),nodhn(n),ssh2(n)
           endif
        endif

     end if
     
     depth(n) = max(ssh2(n) + nodhn(n), 0._wp)

     IF (verbosity>1) THEN

        if (nodhn(n) > 0._wp) then
           sshmax = max(sshmax, ssh2(n))
        else
           sshmax = max(sshmax, depth(n))  
        end if

        if (ssh2(n) >= -nodhn(n)) then
           sshmin = min(sshmin, ssh2(n))
        else
           dmin = min(dmin, ssh2(n))
        end if

     END IF

     ! arrival time
     if (depth(n)> 0.0001_wp .and. ssh2(n) > 0.01_wp &
                            .and. arrival_time(n) < -.5 ) then
           arrival_time(n) = time
        
     end if
  end do   ! in

!$OMP END DO

! Moving bottom
  if (enable_dhdt) CALL readdhdt

! Open boundary 
  call radiate

!$OMP END PARALLEL


END SUBROUTINE compute_ssh_EXTRAPOLATE

  !----------------------------------------------------------------------------
  !               S S H
  !----------------------------------------------------------------------------
SUBROUTINE compute_ssh_GETM

! Description:
!   This subroutine computes the sea surface elevation.
!

  USE MESH
  USE ELEMENTS
  !
  IMPLICIT NONE
  INTEGER                              :: n, i, j, el, iA, iZ
  INTEGER                              :: edges(3), elnodes(3)
  REAL(kind=wp)                       :: vol, gsshx, gsshy, factor, ssh_SIGN
  REAL(kind=wp)                       :: a
  REAL                                 :: vol2, ssh1mean,nodhnmean

  REAL(kind=wp)                       :: hu(3), sshelmax
  INTEGER                            :: imax
  
!$ real(kind=wp) :: omp_get_wtime, t1, t2

  sshmax = -1000. 
  sshmin =  1000. 
  dmin   =  1000. 


!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&         PRIVATE(vol,i,j,n,el,edges,elnodes,gsshx,gsshy,factor,a, &
!$OMP&                 ssh_sign,vol2, ssh1mean,nodhnmean, iA, iZ) &
!$OMP&         REDUCTION(MIN:sshmin,dmin) &
!$OMP&         REDUCTION(MAX:sshmax)

!$OMP DO  
  do n = 1,nod2D
     ssh0(n) = ssh1(n) + alpha*(ssh2(n) - 2. * ssh1(n) + ssh0(n))
     ssh1(n) = ssh2(n)
  end do
!$OMP END DO

!$OMP DO  

  do n = 1,nod2D

     iA = nodelem_ptr(n)
     iZ = nodelem_ptr(n+1)-1
     vol=0.
     ssh2(n) = 0.

!!$     if (ssh1(n)+nodhn(n)>small) then

     if (iZ >= iA ) then

        do i = iA, iZ
           el       = nodelem_addr(i)
           vol      = vol + voltriangle(el)
           edges    = elem2D_edges(:,el)
           elnodes  = elem2D_nodes(:,el)
           
           if (elnodes(1) == n) then
              gsshx = bafu_xy(1, el)
              gsshy = bafu_xy(4, el)
           elseif (elnodes(2) == n) then
              gsshx = bafu_xy(2, el)
              gsshy = bafu_xy(5, el)
           else
              gsshx = bafu_xy(3, el)
              gsshy = bafu_xy(6, el)
           end if
           
        
           a      = 0.
           !factor = ((uv1(1,edges(3))+uv1(1,edges(1))) *gsshx + (uv1(2,edges(3))+uv1(2,edges(1))) *gsshy)   &
           !     * max(nodhn(elnodes(1)) + ssh1(elnodes(1)), a)                          &
           !     + ((uv1(1,edges(1))+uv1(1,edges(2))) *gsshx + (uv1(2,edges(1))+uv1(2,edges(2))) *gsshy)   &
           !     *max(nodhn(elnodes(2)) + ssh1(elnodes(2)), a)                          &
           !     + ((uv1(1,edges(2))+uv1(1,edges(3))) *gsshx + (uv1(2,edges(2))+uv1(2,edges(3))) *gsshy)   &
           !     *max(nodhn(elnodes(3)) + ssh1(elnodes(3)), a)         

           factor=0.
           if (nodhn(elnodes(1)) + ssh1(elnodes(1)) > small) factor = factor+ &
                ((uv1(1,edges(3))+uv1(1,edges(1))) *gsshx + (uv1(2,edges(3))+uv1(2,edges(1))) *gsshy) * &
                (nodhn(elnodes(1)) + ssh1(elnodes(1)))

           if (nodhn(elnodes(2)) + ssh1(elnodes(2)) > small) factor = factor+ &
                ((uv1(1,edges(1))+uv1(1,edges(2))) *gsshx + (uv1(2,edges(1))+uv1(2,edges(2))) *gsshy) * &
                (nodhn(elnodes(2)) + ssh1(elnodes(2)))

           if (nodhn(elnodes(3)) + ssh1(elnodes(3)) > small) factor = factor+ &
                ((uv1(1,edges(2))+uv1(1,edges(3))) *gsshx + (uv1(2,edges(2))+uv1(2,edges(3))) *gsshy) * &
                (nodhn(elnodes(3)) + ssh1(elnodes(3)))
           
           ssh2(n) = ssh2(n)+factor*voltriangle(el)
           
        end do

        ! Result for ssh   
        
           ssh2(n) = ssh0(n) + dt * ssh2(n) / vol

     endif
           
        !if (index_nod2D(n)==5 ) !SH

!!$     else
!!$
!!$        vol2=0.
!!$        ssh1mean=0.
!!$        nodhnmean=0.
!!$
!!$
!!$        do i = 1, nod_in_elem2D(n)%nmb
!!$           el       = nod_in_elem2D(n)%addr(i)
!!$           vol      = vol + voltriangle(el)
!!$           edges    = elem2D_edges(:,el)
!!$           elnodes  = elem2D_nodes(:,el)
!!$           
!!$           if (elnodes(1) == n) then
!!$              gsshx = bafu_xy(1, el)
!!$              gsshy = bafu_xy(4, el)
!!$           elseif (elnodes(2) == n) then
!!$              gsshx = bafu_xy(2, el)
!!$              gsshy = bafu_xy(5, el)
!!$           else
!!$              gsshx = bafu_xy(3, el)
!!$              gsshy = bafu_xy(6, el)
!!$           end if
!!$           
!!$        
!!$           factor=0.
!!$           if (nodhn(elnodes(1)) + ssh1(elnodes(1)) > small) factor=factor+ &
!!$                ((uv1(1,edges(3))+uv1(1,edges(1))) *gsshx + (uv1(2,edges(3))+uv1(2,edges(1))) *gsshy) * &
!!$                (nodhn(elnodes(1)) + ssh1(elnodes(1)))
!!$           if (nodhn(elnodes(2)) + ssh1(elnodes(2)) > small) factor=factor+ &
!!$                ((uv1(1,edges(1))+uv1(1,edges(2))) *gsshx + (uv1(2,edges(1))+uv1(2,edges(2))) *gsshy) * &
!!$                (nodhn(elnodes(2)) + ssh1(elnodes(2)))
!!$           if (nodhn(elnodes(3)) + ssh1(elnodes(3)) > small)factor=factor+ &
!!$                ((uv1(1,edges(2))+uv1(1,edges(3))) *gsshx + (uv1(2,edges(2))+uv1(2,edges(3))) *gsshy) * &
!!$                (nodhn(elnodes(3)) + ssh1(elnodes(3)))
!!$           
!!$           a      = 0.
!!$           !factor = ((uv1(1,edges(3))+uv1(1,edges(1))) *gsshx + (uv1(2,edges(3))+uv1(2,edges(1))) *gsshy)   &
!!$           !     * max(nodhn(elnodes(1)) + ssh1(elnodes(1)), a)                          &
!!$           !     + ((uv1(1,edges(1))+uv1(1,edges(2))) *gsshx + (uv1(2,edges(1))+uv1(2,edges(2))) *gsshy)   &
!!$           !     * max(nodhn(elnodes(2)) + ssh1(elnodes(2)), a)                          &
!!$           !     + ((uv1(1,edges(2))+uv1(1,edges(3))) *gsshx + (uv1(2,edges(2))+uv1(2,edges(3))) *gsshy)   &
!!$           !     * max(nodhn(elnodes(3)) + ssh1(elnodes(3)), a)         
!!$           
!!$           ssh2(n) = ssh2(n)+1./6.*factor*voltriangle(el)
!!$           
!!$           if (count(ssh1(elnodes)+nodhn(elnodes)>small)==2) then
!!$              ssh1mean=ssh1mean+(sum(ssh1(elnodes))-ssh1(n))*vol
!!$              vol2=vol2+vol
!!$              nodhnmean=nodhnmean+(sum(nodhn(elnodes))-nodhn(n))*vol
!!$           end if
!!$
!!$        end do
!!$
!!$        ! Result for ssh       
!!$           ssh2(n) = ssh0(n) + dt * ssh2(n) / vol
!!$
!!$        if (vol2>0.) then
!!$
!!$           ssh1mean=0.5*(ssh1mean/vol2)
!!$           nodhnmean=0.5*(nodhnmean/vol2)
!!$
!!$           if (index_nod2D(n)==6 .and. ssh1mean+nodhnmean+ssh2(n)+nodhn(n)>small) then
!!$              ssh2(n)=ssh2(n)-nodhn(n)+nodhnmean
!!$              index_nod2D(n)=5
!!$           else if (index_nod2D(n)==5 .and. ssh1mean+nodhnmean+ssh2(n)-nodhn(n)<=small) then
!!$              ssh2(n)=ssh2(n)+nodhn(n)-nodhnmean
!!$              index_nod2D(n)=6
!!$           end if
!!$
!!$        end if
!!$
!!$     end if

     ! correct if below ground
!     if (ssh2(n)<-nodhn(n)) ssh2(n) = -nodhn(n)
     ssh2(n) = max(min(-nodhn(n),0._wp) , ssh2(n))

     ! update maximum wave height
     mwh(n) = max(mwh(n), ssh2(n) + max(nodhn(n),0.))
     

     ! arrival time

     IF (verbosity>1) THEN

        depth(n) = ssh2(n) + nodhn(n)

        if (nodhn(n)>0.) then
           sshmax = max(sshmax, ssh2(n))
        else
           if ( depth(n) > small) sshmax = max(sshmax, depth(n))  
        end if

        if (ssh2(n)>=-nodhn(n)) then
           sshmin=min(sshmin, ssh2(n))
        else
           dmin=min(dmin, ssh2(n))
        end if

     END IF

     if (arrival_time(n) < -.5 .and. nodhn(n) >= 0. ) then
        if (ssh2(n) > 0.01) arrival_time(n) = time
     end if
  end do   ! in


!$OMP END DO




!!$!$OMP DO 
!!$
!!$  ! find elements with less than three wet nodes and adjust for wetting
!!$
!!$  do el = 1,elem2D
!!$
!!$     elnodes = elem2D_nodes(:, el)
!!$     hu      = ssh2(elnodes) + nodhn(elnodes)
!!$
!!$     if (any(hu(:) >= small) .and. any(hu(:) < small)) then
!!$
!!$        do i=1,3
!!$
!!$           if (nodhn(elnodes(i))<0. &
!!$                .and. -nodhn(elnodes(i))<sshelmax) then
!!$
!!$              !if (ssh2(imax)>=ssh1(imax)) then
!!$
!!$              ssh2(elnodes(i))=max(ssh2(elnodes(i)),-nodhn(elnodes(i)))
!!$              index_nod2D(elnodes(i))=5
!!$
!!$              !else
!!$              !   ssh2(elnodes(i))=min(ssh2(elnodes(i)),-nodhn(elnodes(i)))
!!$              !end if
!!$
!!$           end if
!!$        end do
!!$
!!$     end if
!!$
!!$  end do
!!$!$OMP END DO 

! Moving bottom
  if (enable_dhdt) CALL readdhdt

! Open boundary 
  call radiate

!$OMP END PARALLEL



END SUBROUTINE compute_ssh_GETM

!
!Androsov 13.07.2020--------------------------------------------------
! subroutine energy
!
! Compute Kinetik and potential energy
!
subroutine compute_energy

  use MESH
  use ELEMENTS
  use PARAMETERS
!
  implicit none
  integer                  :: el, elnodes(3), depH, a, b
  real(kind=wp), parameter :: one_third=1._wp/3._wp, plop0=1026.0_wp
  
!$OMP MASTER  
   Epot = 0.0_wp
   Ekin = 0.0_wp   
!$OMP END MASTER
!$OMP BARRIER

!NR orphaned loop, local variables are private by default
   
!$OMP DO REDUCTION(+:Epot,Ekin) 
   do el=1,elem2d

      elnodes = elem2D_nodes(:,el)

      depH = sum(ssh2(elnodes) + nodhn(elnodes))*one_third
      a    = sum(ssh2(elnodes))*one_third
      b    = sum(sqrt(uv_node(1,elnodes)*uv_node(1,elnodes) &
                    + uv_node(2,elnodes)*uv_node(2,elnodes)))*one_third

      Epot = Epot + voltriangle(el)      *a*a
      Ekin = Ekin + voltriangle(el)*depH *b*b

   enddo
!$OMP END DO
   
!$OMP MASTER
      Epot = 0.5_wp*g*plop0* Epot
      Ekin = 0.5_wp*  plop0* Ekin
!$OMP END MASTER   
   
 end subroutine compute_energy


SUBROUTINE readdhdt

! Description:
!   This subroutine computes the change of the bathmetry dh/dt.


  USE MESH
  USE ELEMENTS
 
  IMPLICIT NONE
  INTEGER (kind=sp)              :: i,n, idt
  character*12                   :: cfile 
  REAL    (kind=wp)              :: raux
  integer                        :: infile_dhdt

  raux = 15._wp

  IF ( (MOD(time,raux) < smALL) .AND. (INT(time/15.) < 42.1) ) THEN
!$OMP MASTER
     cfile = 'dhdt_xxx.out'
     WRITE(cfile(6:8),'(I3.3)') INT(time/15.)
     WRITE(*,*) 'READING dh/dt: ', cfile
     OPEN(newunit=infile_dhdt, file = TRIM(MESHPATH)//cfile,status='old')

     DO i = 1,nod2D
        READ(infile_dhdt,*) dhdt(i)
     END DO

     CLOSE(infile_dhdt)
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO  
      do n=1,nod2D
         dhdt(n) = -30. * dhdt(n)      !sh bathymetry is measured DOwnwards
         nodhn(n) = nodhn(n) + dhdt(n) * dt
         ssh2(n) = ssh2(n) - 2._wp*dt*dhdt(n)
      enddo
!$OMP END  DO
  ENDIF

  RETURN

END SUBROUTINE readdhdt

SUBROUTINE build_index_bnd

! Description:
!   This subroutine builds the index of the boundary nodes.

  USE MESH
  USE ELEMENTS
  !
  IMPLICIT NONE
  INTEGER             :: i, j, k, n, edg, iA, iZ
  REAL(kind=wp)      :: volfac 
 
  nbnd = 0
  DO  n = 1,nod2D 
     IF ( index_nod2D(n) == 2 ) THEN
        nbnd = nbnd + 1
     END IF
  END DO
  
  ALLOCATE (nod_in_bnd(nbnd))
  ALLOCATE (nodngh_in_bnd(2,nbnd))
  ALLOCATE (nodedgfac_in_bnd(2,nbnd))
  !
  WRITE(*,*) "nbnd=", nbnd
 
  nbnd = 0
  DO  n = 1,nod2D
     IF ( index_nod2D(n) == 2 ) THEN


        iA = nodedge_ptr(n)
        iZ = nodedge_ptr(n+1)-1
        if (iZ <= iA) cycle  ! Skip degenerated, isolated nodes
                             ! Should not happen in a good mesh...
        
        nbnd             = nbnd + 1
        nod_in_bnd(nbnd) = n

        volfac = 1._wp/sum(voltriangle(nodelem_addr(nodelem_ptr(n): &
                                                   nodelem_ptr(n+1)-1)))

        j=0
        ! Check the edges adjacent to nod n if they are boundary edges.
        ! When two boundary edges are found, it is fine.
        do i = iA, iZ
           edg = nodedge_addr(i)
           if (edge_in_elem2D(2,edg) == 0 ) then

              j = j+1

              nodedgfac_in_bnd(j,nbnd) = edg_length(edg)*volfac

              if (n /= edge_nod2D(1,edg)) then
                 nodngh_in_bnd(j,nbnd) = edge_nod2D(1,edg)
              else
                 nodngh_in_bnd(j,nbnd) = edge_nod2D(2,edg)
              endif

              if (j==2) exit ! Assume that each boundary node has two 
                             ! adjacent boundary edges. It could be 4 if
                             ! two pieces of sea touch in this node
                             ! Should not happen in a good mesh...
           endif
        enddo
        
        if (j<2) then 
           print *,'Check grid: Open boundary node ',n,' has no boundary edges?!'
           nbnd = nbnd-1
        endif
     END IF
  END DO

  if (.not. write_cfl) then
     deallocate(edg_length)
  else
     edg_length(:)=1._wp/edg_length(:)
  endif
  
  RETURN
END SUBROUTINE build_index_bnd

SUBROUTINE radiate

! Description:
!   This subroutine computes the radiation boundary conditions that are used
!   for the correction in calculation of the sea elevation.

  USE MESH
  USE ELEMENTS
  !
  IMPLICIT NONE
  INTEGER                              :: i, n1, n2, n
  REAL(kind=wp)                       :: radiation, hu1, hu2, hu
  
!$OMP DO schedule(auto) 
!DIR$ IVDEP  ! OMP SIMD has more implications and seems to slow more loops down than speeding them up.
             ! But it might help the compiler to know that all n values are different and there is
             ! no dependcy among loop iterations.
  do i=1,nbnd
      
     n  = nod_in_bnd(i)

     if (abs(ssh2(n)) > 0._wp) then
        hu1 = ssh1(nodngh_in_bnd(1,i)) + nodhn(nodngh_in_bnd(1,i))
        hu2 = ssh1(nodngh_in_bnd(2,i)) + nodhn(nodngh_in_bnd(2,i))
        hu  = ssh1(n)                  + nodhn(n)
     
        radiation = nodedgfac_in_bnd(1,i)*sqrt(.125_wp*g*(max(hu1+hu,0._wp))) &
                  + nodedgfac_in_bnd(2,i)*sqrt(.125_wp*g*(max(hu2+hu,0._wp))) 

        ssh2(n) = ssh2(n) / (1._wp + 6._wp*dt*radiation)

        ! update maximum wave height
        mwh(n) = max(mwh(n),ssh2(n) + min(nodhn(n),0.))
        
     endif
  enddo
!$OMP END DO 

  RETURN	 
END SUBROUTINE radiate  



SUBROUTINE compute_gradient_EXTRAPOLATE

! Description:
!   This subroutine computes the surface gradient.

  USE MESH
  USE ELEMENTS

  IMPLICIT NONE

  INTEGER             :: i,j, k, n, el, el1, elaround(50), n1, n2 ,n3
  REAL(kind=wp)        :: volsum, hu(3), grad_x, grad_y

! Orphaned OpenMP region, local variables are private by default

!$OMP DO  

  ! find elements with three wet nodes and determine ssh gradient

  do el = 1,elem2D    
    
     i_wet(el) = count(depth(elem2D_nodes(1:3,el)) >= small) -3

     ! If all nodes are wet, i_wet(el) == 0 and we will calculate the gradient.         
     ! If all nodes are dry, i_wet(el) == -3 and we forget about el for the time being         
     ! Otherwise (i_wet(el) == -1, -2), the gradient will be extrapolated from                 
     ! neighbouring wet elements in a second loop.                        

     if ( i_wet(el) == 0 .and. any(abs(ssh2(elem2D_nodes(1:3,el)))>0._wp)) then

        grad(1,el) = bafu_xy(1,el)*ssh2(elem2D_nodes(1,el)) &
                   + bafu_xy(2,el)*ssh2(elem2D_nodes(2,el)) &
                   + bafu_xy(3,el)*ssh2(elem2D_nodes(3,el))

        grad(2,el) = bafu_xy(4,el)*ssh2(elem2D_nodes(1,el)) &
                   + bafu_xy(5,el)*ssh2(elem2D_nodes(2,el)) &
                   + bafu_xy(6,el)*ssh2(elem2D_nodes(3,el))

        ! We want to know it more precisely - only if the gradient is not zero, 
        ! we need it in some further computations. 
        i_wet(el) = count(abs(grad(1:2,el))>0._wp) 

     else

        grad(1,el)=0._wp
        grad(2,el)=0._wp
        
     end if

  end do  ! el          
!$OMP END DO

!$OMP DO 
  do el=1,elem2D
     if (i_wet(el) == -1 .or. i_wet(el) == -2) then

        k  = 0

! cycle through patch around the 1st node (if it is wet)
        if (depth(elem2D_nodes(1,el)) >= small) then
           do j = nodelem_ptr(elem2D_nodes(1,el)), nodelem_ptr(elem2D_nodes(1,el)+1)-1
              if ( (i_wet(nodelem_addr(j)) >= 0) ) then
                 k           = k + 1
                 elaround(k) = nodelem_addr(j)  
              end if
           end do
        end if
! cycle through patch around the 2nd node (if it is wet)
        if (depth(elem2D_nodes(2,el)) >= small) then
           do j = nodelem_ptr(elem2D_nodes(2,el)), nodelem_ptr(elem2D_nodes(2,el)+1)-1
              if ( (i_wet(nodelem_addr(j)) >= 0) .and. all(elaround(1:k) /= nodelem_addr(j))) then
                 k           = k + 1
                 elaround(k) = nodelem_addr(j)  
              end if
           end do
        end if
! cycle through patch around the 3rd node (if it is wet)
        if (depth(elem2D_nodes(3,el)) >= small) then
           do j = nodelem_ptr(elem2D_nodes(3,el)), nodelem_ptr(elem2D_nodes(3,el)+1)-1
              if ( (i_wet(nodelem_addr(j)) >= 0) .and. all(elaround(1:k) /= nodelem_addr(j))) then
                 k           = k + 1
                 elaround(k) = nodelem_addr(j)  
              end if
           end do
        end if
        
        if ( k == 1) then           
           grad(1,el) = grad(1,elaround(1))
           grad(2,el) = grad(2,elaround(1))
        elseif ( k > 0 ) then
           volsum = 1._wp/sum(voltriangle(elaround(1:k)))
           grad(1,el) = sum(grad(1,elaround(1:k))*voltriangle(elaround(1:k))) * volsum
           grad(2,el) = sum(grad(2,elaround(1:k))*voltriangle(elaround(1:k))) * volsum
        endif


     endif ! i_wet(el) = -1, partly wet, extrapolate
  end do  ! el
!$OMP END DO NOWAIT


END SUBROUTINE compute_gradient_EXTRAPOLATE


SUBROUTINE compute_gradient_GETM

! Description:
!   This subroutine computes the surface gradient.

  USE MESH
  USE ELEMENTS

  IMPLICIT NONE

  INTEGER                            :: i, j, k, n, el, el1, elnodes(3), elin, elaround(50)
  REAL(kind=wp)                       :: vol, volsum, x, y
  REAL(kind=wp)                       :: hu(3), sshelmax


!$OMP DO  

  ! find elements with three wet nodes and determine ssh gradient

  do el = 1,elem2D


     if ( depth(elem2D_nodes(1,el)) >= small .and. &
          depth(elem2D_nodes(2,el)) >= small .and. &
          depth(elem2D_nodes(3,el)) >= small ) then

        grad(1,el) = bafu_xy(1,el)*ssh2(elem2D_nodes(1, el)) &
                   + bafu_xy(2,el)*ssh2(elem2D_nodes(2, el)) &
                   + bafu_xy(3,el)*ssh2(elem2D_nodes(3, el))

        grad(2,el) = bafu_xy(4,el)*ssh2(elem2D_nodes(1, el)) &
                   + bafu_xy(5,el)*ssh2(elem2D_nodes(2, el)) &
                   + bafu_xy(6,el)*ssh2(elem2D_nodes(3, el))
     else
        grad(1,el)=0.
        grad(2,el)=0.
     end if

  end do  ! el          
!$OMP END DO NOWAIT

!$OMP DO 
  do n=1,nod2D

     mom_fact(n) = max( min( (depth(n)-H_mn)*H_mult, 1._wp), 0._wp)

  end do
!$OMP END DO NOWAIT



END SUBROUTINE compute_gradient_GETM


!
subroutine compute_velocity_at_nodes
!
  use MESH
  use ELEMENTS
!  
  implicit none
  integer                      :: el,  k, nod, ie
  integer                      :: n, iA, iZ
  real(kind=wp)                 :: volsum
  real(kind=sp)                 :: abs_vel
  
! Orphaned OpenMP region, local variables are private by default 

 if (.not. save_vel_at_nodes) then
!$OMP DO  
  do n = 1,nod2D

    if (depth(n) > small) then

        iA = nodedge_ptr(n)
        iZ = nodedge_ptr(n+1)-1

        uv_node(1,n) = sum(uv1(1,nodedge_addr(iA:iZ))*nodedge_wght(iA:iZ))
        uv_node(2,n) = sum(uv1(2,nodedge_addr(iA:iZ))*nodedge_wght(iA:iZ))

        if (abs(uv_node(1,n)) > 0._wp .or. abs(uv_node(2,n)) > 0._wp) then
           max_flux(n)    = max(max_flux(n), max_abs_vel(n) * real(depth(n),sp) )
           ! update max_abs_vel as a second step - this is not exact, but good enough,
           ! and faster (allows to overlap computations)
           max_abs_vel(n) = max(max_abs_vel(n),  &
                sqrt(real(uv_node(1,n),sp)**2 + real(uv_node(2,n),sp)**2))
        endif
     end if
  enddo
!$OMP END DO NOWAIT
else
!$OMP DO  
   do n=1,nod2D
      if (depth(n) > small) then

         iA = nodedge_ptr(n)
         iZ = nodedge_ptr(n+1)-1

         uv_node(1,n) = sum(uv1(1,nodedge_addr(iA:iZ))*nodedge_wght(iA:iZ))
         uv_node(2,n) = sum(uv1(2,nodedge_addr(iA:iZ))*nodedge_wght(iA:iZ))

         if (abs(uv_node(1,n)) > 0._wp .or. abs(uv_node(2,n)) > 0._wp) then
            max_flux(n)    = max(max_flux(n), max_abs_vel(n) * real(depth(n),sp) )
            max_abs_vel(n) = max(max_abs_vel(n),  &
                 sqrt(real(uv_node(1,n),sp)**2 + real(uv_node(2,n),sp)**2))
         endif
         uv_node_old(1,n) = uv_node(1,n)
         uv_node_old(2,n) = uv_node(2,n)
         
      else
        uv_node_old(1,n) = 0._wp
        uv_node_old(2,n) = 0._wp
     end if
  end do
!$OMP END DO NOWAIT
end if
  

end subroutine compute_velocity_at_nodes
 


subroutine compute_ttt

  USE MESH
  USE ELEMENTS
  
  implicit none

  integer      :: status(nod2D)
  integer      :: nod_todo(nod2D),nod_todonext(nod2D), nod_count(nod2D)
  integer      :: n, i, j, ngh, l_ini, nd, ngh_next, ii, k, nnew
  real(kind=wp) :: topo_fact, X, Y, Xn, Yn, dist, lc, tnow

  allocate(ttt(nod2D))
  
  ! Init fields. Arrival time =0 for all nodes
  ! with inital ssh <= 80% of maximum initial ssh 
  l_ini = 0
  
  sshmax = 0.
!$OMP DO  REDUCTION(MAX:sshmax)
  do n = 1,nod2D
     if (nodhn(n) >= 0.) then
        sshmax = max(sshmax, ssh0(n))
        status(n) = 0
     else
        status(n) = 3
     endif
     ttt(n) = -0.000001
     nod_todo(n) = 0.
     nod_todonext(n) = 0
     nod_count(n) = 0
  end do
!$OMP END DO
  print *,'ttt with ssh-max= ',sshmax

  do n = 1,nod2D
     
     if (nodhn(n) >= 0. .and. ssh0(n) > .8*sshmax) then
        ttt(n) = 0.
        status(n) = 3
        l_ini = l_ini+1
        nod_todonext(l_ini) = n
     endif

  end do
  print *,'l_ini ',l_ini

  
  ! Starting from these nodes calculate travel time
  ! from the model depth values
  
  ! Determine the nodes on the rim of the initial field
  ngh = 0

  do i=1,l_ini
     nd = nod_todonext(i)

     do j = 1,nghbr_nod2D(nd)%nmb         ! Loop over al neighbouring nodes 
        n = nghbr_nod2D(nd)%addr(j)  
        if (status(n) == 0) then
           ngh = ngh+1
           nod_todo(ngh) = nd
           exit
        endif
     end do
  end do
  
  print *,'ngh = ', ngh
  lc = 1.  
  tnow = 60.

  Iterate: do ii=1,100000 
     ngh_next = 0
     nnew =0

     do k = 1, ngh
        nd = nod_todo(k)

        if (ttt(nd) <= tnow) then

           topo_fact = 1./sqrt(9.81*max(5.d0, nodhn(nd)))

           X = coord_nod2D(1, nd)
           Y = coord_nod2D(2, nd)
           if (coordinate_type == 1) lc = cos(Y)

           do j = 1,nghbr_nod2D(nd)%nmb 
              n = nghbr_nod2D(nd)%addr(j)  

              if (status(n) <= 2) then
                 Xn = coord_nod2D(1, n)
                 Yn = coord_nod2D(2, n)
                 dist = scaling * sqrt((lc*(X-Xn))**2 + (Y-Yn)**2)
                 
                 ttt(n) = ttt(n) + ttt(nd) + dist*topo_fact
                 nod_count(n) = nod_count(n)+1

                 if (status(n) == 0) then
                    ngh_next = ngh_next + 1;
                    nod_todonext(ngh_next) = n;
                    status(n) = 2;
                    nnew = nnew+1
                 endif
              endif
           end do
           
        else  ! ttt(nd) > tnow 
           ngh_next = ngh_next + 1;
           nod_todonext(ngh_next) = nd;  
        endif

     end do

     if (nnew==0) then
        tnow = tnow+60.
        if (tnow > 60.*60.*24.) exit Iterate
        if (ngh_next==0)  exit Iterate
     else

        do k=1,ngh_next
           n = nod_todonext(k)
           if (status(n)==2) ttt(n) = ttt(n) / real(nod_count(n),dp)
           nod_todo(k)  = n
           nod_todonext(k)  = 0
           status(n)    = 3
        end do

        ngh = ngh_next
     endif
  end do Iterate

  do n=1,nod2D
     if (ttt(n) <  0) ttt(n) = -1.
  enddo


end subroutine compute_ttt

!---------------------------------------------------------------------------


subroutine alloc_swe_fields

  USE MESH
  USE ELEMENTS
  implicit none
  integer      :: i, n, elin,edg,el, edgA, edgB, iA, iZ, ie, iA_el, iZ_el
  real(kind=wp) :: volfac

    allocate(ssh1(nod2D), ssh2(nod2D), ssh0(nod2D), ssh_init(nod2D),depth(nod2D))
    allocate(uv0(2,edg2D), uv1(2,edg2D), uv2(2,edg2D))
    allocate(max_abs_vel(nod2D), max_flux(nod2D))
    allocate(grad(2,elem2d), i_wet(0:elem2D))
    allocate(uv_node_old(2,nod2D))
    allocate(uv_node(2,nod2D))

    i_wet(0) = 0

!$OMP PARALLEL
!$OMP DO 
    do n=1,nod2D
       ssh1(n)=0.
       ssh2(n)=0.
       ssh0(n)=0.
       depth(n) = 0.
       max_abs_vel(n)=0.
       max_flux(n)=0.       
    enddo
!$OMP END DO NOWAIT
!$OMP DO 
    do i=1,edg2D
       uv0(1,i)=0.
       uv1(1,i)=0.
       uv2(1,i)=0.

       uv0(2,i)=0.
       uv1(2,i)=0.
       uv2(2,i)=0.
    enddo
!$OMP END DO NOWAIT
!$OMP DO 
    do i=1,elem2D
       i_wet(i) = 0 !.false.
    enddo
!$OMP END DO
!$OMP END PARALLEL
 
    if (wetting_drying_type==2) then
       allocate(mom_fact(nod2D))
 !$OMP PARALLEL DO 
       do n=1,nod2D
          mom_fact(n)=0.
       end do
!$OMP END PARALLEL DO
       H_mult=1./(H_crit-H_mn) ! auxiliary ratio to be multiplied
       print *,'mom_fact allocated, H_mn= ',H_mn,' H_crit= ',H_crit 
    endif

    ALLOCATE(coriolis_param(edg2D))
    IF ( rotation_type == 2 ) THEN
!$OMP PARALLEL DO 
       DO i = 1, edg2D
          coriolis_param(i) = 2._wp * omega * SIN(0.5_wp * SUM(coord_nod2D(2, edge_nod2D(:,i))))
       END DO
!$OMP END PARALLEL DO
       else
!$OMP PARALLEL DO 
       DO i = 1, edg2D
          coriolis_param(i) = 0._wp
       END DO
!$OMP END PARALLEL DO          
    END IF

    if (enable_dhdt) then
       ALLOCATE (dhdt(nod2D))
!$OMP PARALLEL DO  
       do n=1,nod2D
          dhdt(n)=0.
       end do
!$OMP END PARALLEL DO
    end if

!---- Diagnostic fields
    ! set up fields for arrival time and max. wave height

    ALLOCATE(arrival_time(nod2D))
    ALLOCATE(mwh(nod2D))
    
!$OMP PARALLEL DO 
    do i=1, nod2D
       arrival_time(i) = -1.
       mwh(i)          =  0.
    enddo
!$OMP END PARALLEL DO

    
! Fields with precomputed values to speed up computation later on
    allocate(bafuvel_edgelem(4,edg2D)) ! To be used in compute_velocity, viscosity-computation

    allocate(nodedge_wght(nodedge_ptr(nod2D+1)-1))   ! To streamline velocity_at_nodes
    
    allocate(nodelem_wght(2,nodelem_ptr(nod2D+1)-1)) ! To streamline compute_ssh

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elin,edg,el,iA,iZ,iA_el, iZ_el,volfac,edgA,edgB,ie)
!$OMP DO 
    do edg=1,edg2D
       do elin=1, 2

          el=edge_in_elem2D(elin, edg)
          if (el == 0) cycle

          if (elem2D_edges(1,el) == edg) then
             
             bafuvel_edgelem(2*elin-1, edg) = -2._wp*bafu_xy(3,el)*voltriangle(el)
             bafuvel_edgelem(2*elin,   edg) = -2._wp*bafu_xy(6,el)*voltriangle(el)
             
          elseif (elem2D_edges(2,el) == edg) then
             
             bafuvel_edgelem(2*elin-1, edg) = -2._wp*bafu_xy(1,el)*voltriangle(el)
             bafuvel_edgelem(2*elin,   edg) = -2._wp*bafu_xy(4,el)*voltriangle(el)                    
                 
          else !if (edges(3) == edg) then
                    
             bafuvel_edgelem(2*elin-1, edg) = -2._wp*bafu_xy(2,el)*voltriangle(el)
             bafuvel_edgelem(2*elin,   edg) = -2._wp*bafu_xy(5,el)*voltriangle(el)                    
             
          endif
       enddo
    enddo
!$OMP END DO NOWAIT
!$OMP DO 
  do n = 1,nod2D

        iA = nodedge_ptr(n)
        iZ = nodedge_ptr(n+1)-1
        if (iZ <= iA) cycle
        
        nodedge_wght(iA:iZ) = 0.

        iA_el = nodelem_ptr(n)
        iZ_el = nodelem_ptr(n+1)-1
        volfac = 0.5_wp/sum(voltriangle(nodelem_addr(iA_el:iZ_el)))

        do ie = iA_el, iZ_el

           el = nodelem_addr(ie)
           
           if (n == elem2D_nodes(1,el)) then
              edgA = elem2D_edges(1,el)
              edgB = elem2D_edges(3,el)

           elseif (n == elem2D_nodes(2,el)) then
              edgA = elem2D_edges(1,el)
              edgB = elem2D_edges(2,el)

           else !if (n == elem2D_nodes(3,el)) then
              edgA = elem2D_edges(2,el)
              edgB = elem2D_edges(3,el)
           end if
           
           do i=iA, iZ
              if (nodedge_addr(i) == edgA .or. nodedge_addr(i) == edgB) then
                 nodedge_wght(i) = nodedge_wght(i) + voltriangle(el)*volfac
              endif
           enddo
        enddo
  enddo    
!$OMP END DO NOWAIT
!$OMP DO 
  do n = 1,nod2D

     iA = nodelem_ptr(n)
     iZ = nodelem_ptr(n+1)-1

     if (iZ >= iA) then
        volfac = 1._wp/sum(voltriangle(nodelem_addr(iA:iZ)))

        do i = iA, iZ
           el    = nodelem_addr(i) 
        
           if (elem2D_nodes(1,el) == n) then

              nodelem_wght(1,i) = voltriangle(el)* bafu_xy(1,el)*volfac
              nodelem_wght(2,i) = voltriangle(el)* bafu_xy(4,el)*volfac
        
           elseif (elem2D_nodes(2,el) == n) then 
              nodelem_wght(1,i) = voltriangle(el)* bafu_xy(2,el)*volfac
              nodelem_wght(2,i) = voltriangle(el)* bafu_xy(5,el)*volfac
        
           else
              nodelem_wght(1,i) = voltriangle(el)* bafu_xy(3,el)*volfac
              nodelem_wght(2,i) = voltriangle(el)* bafu_xy(6,el)*volfac
        
           endif
                      
        end do
     endif
  enddo
!$OMP END DO
!$OMP END PARALLEL

  end subroutine alloc_swe_fields

! !Description: This module contains shallow water equations related subroutines.

END MODULE SWE
