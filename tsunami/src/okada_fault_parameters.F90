! Part of TsunAWI
! Copyright (c) 2014
! Alfred Wegener Institut,
! Helmholtz Centre for Polar and Marine Research
! Contact: rakowsky@awi.de


! !MODULE: OKADA_FAULT_PARAMETERS - contains fault parameters andvariables nee-
!	ded for the determination of initial conditions from these parameters
! !DESCRIPTION: This module contains parameters related to initial conditions 
!               based on Okada parameters. These are read from a file with all
!               sets of parameters. Starting from these parameters the 
!               initial conditions are determined in a regular auxiliar grid
!               by a tool developed in the IUGG/IOC TIME Project.
!               In a second step these initial conditions are interpolated
!               onto the model grid
!
! Definitions are:\\
! \begin{tabular}{l|lr}\hline
! Name & Description                    & Unit \\ \hline
! pi   & $\pi=3.141592653589793$        & -- \\
! rad  & radians: $rad=\frac{\pi}{180}$ & -- \\
! X0   & Fault origin longitude         & deg west \\
! Y0   & Fault Origin latitude          & deg north\\
! D    & Slip (Dislocation) & m\\
! L    & fault length & m \\
! W    & fault width & m\\
! TH   & Strike Angle & deg \\
! DL   & Dip Angle & deg\\
! HH   & fault depth & m\\
! RD   & Direction of dislocation & deg\\
! XCD  & origin of computational domain lon & deg west \\
! YCD  & origin of computational domain lat & deg north\\
! Z    & Array with initial conditions &\\
! dr   & resolution of initial cond in regular auxiliar mesh & deg\\
! nfault & number of fault planes in simulation& -- \\
! ttl\_nfault & total number of fault planes & -- \\
! faults  & array with list of faults in simulation (size: nfault) & -- \\
! AllFaults & array with all fault parameters (size: ttl\_nfault x 11) & -- \\ \hline
! \end{tabular}


MODULE OKADA_FAULT_PARAMETERS

    use PARAMETERS
    IMPLICIT NONE
    SAVE


    INTEGER                                   :: nfault, ttl_nfault

    INTEGER,parameter                         :: ig=501, jg=501

    REAL(KIND=dp), DIMENSION(ig,jg)            :: Z
    
    !  [min] resolution in spherical coordinates
    REAL(KIND=dp), PARAMETER                   :: dr = 2.0     
    REAL(KIND=dp), PARAMETER                   :: Z_thresh = 0.000001  !0.30 0.05
  
    !  [m]   resolution in cartesian coordinates (not used)
    REAL(KIND=dp), PARAMETER                   :: dx = 15.0e3  

    !  List of faults from file
    INTEGER,      ALLOCATABLE, DIMENSION(:)   :: faults 
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: AllFaults
    REAL(KIND=dp)                              :: ttl_rupt_time
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)   :: rupt_delay
    INTEGER                                    :: infile_okada, infile_list

CONTAINS

 SUBROUTINE okada_read_faultlist

  implicit none

  integer      :: iostatus, open_status, n
  real(kind=dp) :: dummy  

  open(newunit=infile_okada, file = okada_parameter_file, form = 'formatted')
  iostatus=0
  ttl_nfault=0
  do while(iostatus==0)
     read (infile_okada,*,iostat=iostatus) dummy
     ttl_nfault=ttl_nfault+1
  end do
  close(infile_okada)
  ttl_nfault=ttl_nfault-1
  print *,'Total number of faults in file: ',ttl_nfault
  
  open(newunit=infile_list, file = okada_fault_list_file, form = 'formatted')
  read(infile_list, *) nfault
  allocate(faults(nfault), AllFaults(11,ttl_nfault), rupt_delay(nfault))
  read(infile_list, *) faults
  read(infile_list, *) rupt_delay
  ttl_rupt_time = rupt_delay(nfault)
  close(infile_list)

END SUBROUTINE okada_read_faultlist

SUBROUTINE initial_conditions_param

    IMPLICIT NONE

    INTEGER :: i, n
    INTEGER :: infile_okada

    REAL(KIND=dp)                              :: X0, Y0, D, L, W, TH, DL, HH, RD
    !  Origin computational domain                    
    REAL(KIND=dp)                              :: XCD, YCD             

    Z = 0.0
  
    OPEN(newunit=infile_okada, file=okada_parameter_file, form='formatted')
    DO i = 1 , ttl_nfault
        READ(infile_okada,*) AllFaults(:,i)
    END DO
    close(infile_okada)

    DO i = 1 , nfault
        IF (rupt_delay(i) == 0.0) THEN

            PRINT *,'Processing Fault Nmb. ',i,faults(i)

            X0  = AllFaults(1,faults(i))
            Y0  = AllFaults(2,faults(i))
            D   = AllFaults(3,faults(i))  ! Slip
            L   = AllFaults(4,faults(i))
            W   = AllFaults(5,faults(i))
            TH  = AllFaults(6,faults(i))  ! Strike angle
            DL  = AllFaults(7,faults(i))  ! Dip angle
            HH  = AllFaults(8,faults(i))  ! Depth
            RD  = AllFaults(9,faults(i))  ! Rake angle
            XCD = AllFaults(10,faults(i))
            YCD = AllFaults(11,faults(i))

            CALL deform(X0, Y0, D, L, W, TH, DL, HH, RD, XCD, YCD)             
            PRINT *,i,minval(Z),maxval(Z)
        END IF
    END DO

    CALL interpolate_Z(XCD,YCD)

    

END SUBROUTINE initial_conditions_param


SUBROUTINE update_initial_ssh
!  Update ssh if the rupture takes several time steps

    IMPLICIT NONE

    INTEGER :: i

    REAL(KIND=dp)                              :: X0, Y0, D, L, W, TH, DL, HH, RD
    !  Origin computational domain                    
    REAL(KIND=dp)                              :: XCD, YCD        

     
    Z = 0.0_dp
  
    DO i = 1 , nfault

        IF (ABS(rupt_delay(i)-time) < 1.e-5) THEN

           PRINT *,'Processing Fault Nmb. ',i,faults(i)

           X0  = AllFaults(1 ,faults(i))
           Y0  = AllFaults(2 ,faults(i))
           D   = AllFaults(3 ,faults(i))
           L   = AllFaults(4 ,faults(i))
           W   = AllFaults(5 ,faults(i))
           TH  = AllFaults(6 ,faults(i))
           DL  = AllFaults(7 ,faults(i))
           HH  = AllFaults(8 ,faults(i))
           RD  = AllFaults(9 ,faults(i))
           XCD = AllFaults(10,faults(i))
           YCD = AllFaults(11,faults(i))

           CALL deform(X0, Y0, D, L, W, TH, DL, HH, RD, XCD, YCD)   
           PRINT *,i,minval(Z),maxval(Z)

        END IF

     END DO

     CALL interpolate_Z_again(XCD,YCD)

END SUBROUTINE update_initial_ssh


REAL(KIND=dp) FUNCTION atn(ax, ay)

    IMPLICIT NONE

    REAL(KIND=dp)            :: ax, ay

    REAL(KIND=dp), PARAMETER :: gx = 1.0e-6
    REAL(KIND=dp)            :: aax, aay, p, sr

    aax = ABS(ax) ; aay = ABS(ay)
    p = ax*ay
  
    IF (aax <= gx .AND. aay <= gx) THEN
        WRITE(6,'("ATAN -- ax=",E15.7,2X,"ay=",E15.7)') ax,ay
        atn = 0.2
    ELSE
        sr  = atan2(aax,aay)
        atn = sign(sr,p)
    END IF

END FUNCTION atn

SUBROUTINE deform(X0, Y0, D, L, W, TH, DL, HH, RD, XCD, YCD)
  
    IMPLICIT NONE

    REAL(KIND=dp)                              :: X0, Y0, D, L, W, TH, DL, HH, RD
    !  Origin computational domain                    
    REAL(KIND=dp)                              :: XCD, YCD  

    INTEGER      :: i, j
    REAL(KIND=dp) :: XL, YL, H1, H2, DS, DD
    REAL(KIND=dp) :: XX, YY                         ! nod(x,y)
    REAL(KIND=dp) :: X1, X2, X3
    REAL(KIND=dp) :: F1, F2, F3, F4, G1, G2, G3, G4
    REAL(KIND=dp) :: US, UD
    INTEGER      :: iif, jjf

    XL = pi * r_earth * (X0-XCD) * cos(rad*YCD)/180.0
    YL = pi * r_earth * (Y0-YCD) / 180.0
    H1 = HH / sin(rad*DL)
    H2 = HH / sin(rad*DL) + W
    DS =  D * cos(rad*RD)
    DD =  D * sin(rad*RD)
    !write(6,*) XL, YL, H1, H2, DS,  DD
    
    DO j = 1 , jg
       YY = pi *r_earth*DR*(j-1)/(60.0*180.0)
       DO i = 1 , ig
        
            XX = pi *r_earth*DR*(i-1) * cos(rad*(YCD+DR*(j-1)/60.))/(60.*180.)
            X1 = (XX-XL) * sin(rad*TH) + (YY-YL)*cos(rad*TH) - L/2.0
            X2 = (XX-XL) * cos(rad*TH) - (YY-YL)*sin(rad*TH) + HH/tan(rad*DL)
            X3 = 0.0_dp
        
            CALL uscal(X1, X2, X3,  L/2., H2, rad*DL, F1)
            CALL uscal(X1, X2, X3,  L/2., H1, rad*DL, F2)
            CALL uscal(X1, X2, X3, -L/2., H2, rad*DL, F3)
            CALL uscal(X1, X2, X3, -L/2., H1, rad*DL, F4)
            CALL udcal(X1, X2, X3,  L/2., H2, rad*DL, G1)
            CALL udcal(X1, X2, X3,  L/2., H1, rad*DL, G2)
            CALL udcal(X1, X2, X3, -L/2., H2, rad*DL, G3)
            CALL udcal(X1, X2, X3, -L/2., H1, rad*DL, G4)

            US  = (F1-F2-F3+F4) * DS/(12.0*pi)
            UD  = (G1-G2-G3+G4) * DD/(12.0*pi)
          
            Z(i,j) = Z(i,j)+US+UD
        
            IF(ABS(Z(i,j)) <= Z_thresh) THEN
                Z(i,j) = 0.0_dp
            END IF
        
        END DO
    
    END DO

!.........................................................................
END SUBROUTINE deform
!.........................................................................

SUBROUTINE uscal(X1, X2, X3, C, CC, DPL, F)
! computation of the vertical displacement due to strike slip

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(IN)  :: X1, X2, X3, C, CC, DPL
    REAL(KIND=dp), INTENT(OUT) :: F

    REAL(KIND=dp) :: sn, cs, c1, c2, c3
    REAL(KIND=dp) :: r, r2, r3, q, q2, q3 
    REAL(KIND=dp) :: h, k
    REAL(KIND=dp) :: a1, a2, a3
    REAL(KIND=dp) :: b1, b2, b3, b4, b5, b6, b7, b8, b9
    REAL(KIND=dp) :: b10, b11, b12, b13, b14

    sn  = sin(DPL)
    cs  = cos(DPL)
    c1  = C; c2 = CC*cs; c3 = CC*sn
    r   = sqrt((X1-c1)**2+(X2-c2)**2+(X3-c3)**2)
    q   = sqrt((X1-c1)**2+(X2-c2)**2+(X3+c3)**2)
    r2  = X2*sn - X3*cs
    r3  = X2*cs - X3*sn
    q2  = X2*sn + X3*cs
    q3  =-X2*cs + X3*sn
    h   = sqrt(q2**2+(q3+CC)**2)
    k   = sqrt((X1-c1)**2+q2**2)
    a1  = log(r+r3-CC)
    a2  = log(q+q3+CC)
    a3  = log(q+X3+c3)
    b1  = 1. + 3.0*(tan(DPL))**2
    b2  = 3.0 * tan(DPL)/cs
    b3  = 2.0 * r2 * sn
    b4  = q2 + x2*sn
    b5  = 2.0 * r2**2*cs
    b6  = r*(r+r3-CC)
    b7  = 4.0 * q2 * X3 * sn**2
    b8  = 2.0 * (q2+X2*sn) * (x3+q3*sn)
    b9  = q * (q+q3+CC)
    b10 = 4.0 * q2 * x3 * sn
    b11 = (x3+c3) - q3*sn
    b12 = 4.0 * q2**2 * q3*x3*cs*sn
    b13 = 2.0*q+q3+cc
    b14 = q**3 * (q+q3+CC)

    F   = cs * (a1 + b1*a2 -b2*a3)                  &
          + b3/r + 2.0*sn*b4/q - b5/b6              &
          + (b7-b8)/b9 + b10*b11/q**3 - b12*b13/b14

END SUBROUTINE uscal


SUBROUTINE udcal(X1, X2, X3, C, CC, DPL, F)
! computation of the vertical displacement due to dip slip

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(IN)  :: X1, X2, X3, C, CC, DPL
    REAL(KIND=dp), INTENT(OUT) :: F

    REAL(KIND=dp) :: sn, cs, c1, c2, c3
    REAL(KIND=dp) :: r, r2, r3, q, q2, q3 
    REAL(KIND=dp) :: h, k
    REAL(KIND=dp) :: a1, a2
    REAL(KIND=dp) :: b1, b2, b3, d1, d2, d3, d4, d5, d6
    REAL(KIND=dp) :: t1, t2, t3
!    REAL(KIND=dp), EXTERNAL :: atn

    sn = sin(DPL)
    cs = cos(DPL)
    c1 = C; c2 = CC*cs; c3 = CC*sn
    r  = sqrt((X1-c1)**2+(X2-c2)**2+(X3-c3)**2)
    q  = sqrt((X1-c1)**2+(X2-c2)**2+(X3+c3)**2)
    r2 = X2*sn - X3*cs
    r3 = X2*cs - X3*sn
    q2 = X2*sn + X3*cs
    q3 =-X2*cs + X3*sn
    h  = sqrt(q2**2+(q3+CC)**2)
    k  = sqrt((X1-c1)**2+q2**2)
    a1 = log(r+X1-c1)
    a2 = log(q+X1-c1)
    b1 = q*(q+x1-c1)
    b2 = r*(r+x1-c1)
    b3 = q*(q+q3+CC)
    d1 = X1-c1 ; d2 = X2-c2 ; d3 = X3-c3
    d4 = X3+c3 ; d5 = r3-CC ; d6 = q3+CC
    t1 = atn(d1*d2, (h+d4)*(q+k))
    t2 = atn(d1*d5, r2*r)
    t3 = atn(d1*d6, q2*q)
  
    F  = sn*(d2*(2*d3/b2+4*d3/b1-4*c3*X3*d4*(2*q+d1)/(b1**2*q))      &
         - 6*t1 + 3*t2 -6*t3)                                        &
         + cs * (a1-a2-2*(d3**2)/b2                                  &
         - 4*(d4**2-c3*X3)/b1 - 4*c3*X3*d4**2*(2*q+X1-c1)/(b1**2*q)) &
         + 6*X3 * (cs*sn*(2*d6/b1+d1/b3)-q2*(sn**2-cs**2)/b1)

END SUBROUTINE udcal



SUBROUTINE interpolate_Z(XCD,YCD)

    USE ELEMENTS
    USE MESH
    USE SWE
    
    
    IMPLICIT NONE

    REAL(KIND=dp)   :: lonZ(jg), latZ(ig),XCD,YCD
    INTEGER        :: i, j
    REAL(KIND=dp)   :: xi, yi, xx0, xx1, yy0, yy1
    INTEGER        :: xind, yind
    REAL(KIND=dp)   :: xres, yres 

    REAL(KIND=dp)   :: Z_t(jg,ig)

    INTEGER        :: n, nmb


    print *,'INIT DATA',MINVAL(Z), MAXVAL(Z)

    DO i = 1,ig
        latZ(i) = real(i-1,dp)*dr/60.
    END DO

    DO i = 1,jg
        lonZ(i) = real(i-1,dp)*dr/60.
    END DO

    lonZ = lonZ + XCD
    latZ = latZ + YCD

    lonZ = lonz * rad
    latZ = latZ * rad

    xres = lonZ(2) - lonZ(1)
    yres = latZ(2) - latZ(1)


    DO i = 1,ig
        DO j = 1,jg
            Z_t(j,i) = Z(i,j)
        END DO
    END DO

    DO i=1,nod2D

        xi = coord_nod2D(1,i)
        yi = coord_nod2D(2,i)

        IF (xi>lonZ(1) .AND. xi<lonz(jg)            &
            .AND. yi>latZ(1) .AND. yi<latZ(ig)) THEN

            xind = INT((xi-lonZ(1))/xres)+1
            yind = INT((yi-latZ(1))/yres)+1

            xx0 = lonZ(xind)
            yy0 = latZ(yind)
            xx1 = lonZ(xind+1)
            yy1 = latZ(yind+1)

            !print *,xind,yind,xx0,xi,xx1,yy0,yi,yy1

            IF ((yi-yy0)/yres <= 1-(xi-xx0)/xres) THEN
                ssh0(i) = Z_t(yind,xind) &
                          + (Z_t(yind,xind+1)-Z_t(yind,xind))*(xi-xx0)/xres                     &
                          + (Z_t(yind+1,xind)-Z_t(yind,xind))*(yi-yy0)/yres 
            ELSE
                ssh0(i) = Z_t(yind+1,xind+1) &
                          + (Z_t(yind+1,xind)-Z_t(yind+1,xind+1))*(xx1-xi)/xres                  &  
                          + (Z_t(yind,xind+1)-Z_t(yind+1,xind+1))*(yy1-yi)/yres
            END IF
        
        END IF

    END DO

  ! smoothing initial conditions
    DO i = 1,nmb_iter_smooth_inicond
        DO n = 1,nod2D
            nmb = nghbr_nod2D(n)%nmb
            ssh_init(n) = SUM(ssh0(nghbr_nod2D(n)%addr(:)))/nmb
        END DO
        ssh0 = ssh_init
    END DO

  ! Adjust topography
    DO i = 1,nod2D
        nodhn(i) = nodhn(i) - ssh0(i)
    END DO


    ssh_init = ssh0
    ssh0 = max(ssh0,-nodhn) !sh

    DO i = 1,nod2D
        IF (nodhn(i) < 0.0 .AND. ssh0(i) > ABS(nodhn(i))) then
            ssh0(i) = nodhn(i)
        END IF
    END DO

    ssh1 = ssh0
    ssh2 = ssh0


END SUBROUTINE interpolate_Z


SUBROUTINE interpolate_Z_again(XCD,YCD)

    USE MESH
    USE ELEMENTS
    USE SWE

    IMPLICIT NONE

    REAL(KIND=dp)   :: lonZ(jg), latZ(ig),XCD,YCD
    INTEGER        :: i, j
    REAL(KIND=dp)   :: xi, yi, xx0, xx1, yy0, yy1
    INTEGER        :: xind, yind
    REAL(KIND=dp)   :: xres, yres 

    REAL(KIND=dp)   :: Z_t(jg,ig)
    REAL(KIND=dp)   :: dssh

    
    
    DO i = 1,ig
        latZ(i) = (i-1)*dr/60.
    END DO

    DO i = 1,jg
        lonZ(i) = (i-1)*dr/60.
    END DO

    lonZ = lonZ + XCD
    latZ = latZ + YCD

    lonZ = lonz * rad
    latZ = latZ * rad

    xres = lonZ(2)-lonZ(1)
    yres = latZ(2)-latZ(1)

    
    DO i = 1,ig
        DO j = 1,jg
            Z_t(j,i) = Z(i,j)
        END DO
    END DO


    DO i = 1,nod2D

        xi = coord_nod2D(1,i)
        yi = coord_nod2D(2,i)

        IF (xi > lonZ(1) .AND. xi < lonz(jg)            &
            .AND. yi > latZ(1) .AND. yi < latZ(ig)) THEN

            xind = INT((xi-lonZ(1))/xres) + 1
            yind = INT((yi-latZ(1))/yres) + 1

            xx0 = lonZ(xind)
            yy0 = latZ(yind)
            xx1 = lonZ(xind+1)
            yy1 = latZ(yind+1)

            !print *,xind,yind,xx0,xi,xx1,yy0,yi,yy1

            IF ((yi-yy0)/yres <= 1-(xi-xx0)/xres) THEN
                dssh = Z_t(yind,xind)  &
                       + (Z_t(yind,xind+1)-Z_t(yind,xind))*(xi-xx0)/xres     &
                       + (Z_t(yind+1,xind)-Z_t(yind,xind))*(yi-yy0)/yres
            ELSE
                dssh = Z_t(yind+1,xind+1) &
                       + (Z_t(yind+1,xind)-Z_t(yind+1,xind+1))*(xx1-xi)/xres  &
                       + (Z_t(yind,xind+1)-Z_t(yind+1,xind+1))*(yy1-yi)/yres
            END IF

            ssh0(i)  = ssh0(i) + dssh
            ssh1(i)  = ssh1(i) + dssh
            ssh2(i)  = ssh2(i) + dssh

            nodhn(i) = nodhn(i) - dssh

        END IF

    END DO


END SUBROUTINE interpolate_Z_again
  


END MODULE OKADA_FAULT_PARAMETERS
