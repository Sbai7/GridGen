!******************************************************************************
!*
!*    GridGen: Command-line Grid Generation Program and Much More.
!*    Copyright (C) 2016 BRGM (French Geological Survey) 
!*
!*    This file is part of GwMove hydrogeological software authored by
!*    Dr. M. A. Sbai
!*
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    any later version.
!*
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!*
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*
!*    Redistribution and use in source and binary  forms, with or without
!*    modification, are permitted provided that the following conditions are met:
!*
!*  - Redistributions of  source code must retain the above copyright notice,
!*    this list of conditions and the disclaimer below.
!*
!*  - Redistributions in binary form must reproduce the above copyright notice,
!*    this list of  conditions and the disclaimer (as noted below) in the
!*    documentation and/or other materials provided with the distribution.
!*
!******************************************************************************

SUBROUTINE InterpolateSoils()
! Interpolate soil ID's with a natural neighbor algorithm into an existing 3D 
! grid from a list of hydrogeological cross-sections of the study area. 

   USE GridGlobals

   IMPLICIT NONE
   
   DOUBLE PRECISION, PARAMETER :: big = 1.d+10
   
   CHARACTER(len=80) :: titlep
   INTEGER :: i, j, k, nr, ii, jj ! some indices
   INTEGER :: nn            ! number of nodes 
   INTEGER :: ne            ! number of elements/cells
   INTEGER :: lcount        ! line counter in input file 
   INTEGER :: zeroCount     ! empty cells (no soil) counter 
   INTEGER :: Orient        ! cross-section orientation flag
   INTEGER :: n1, n2        ! cross-section grid size
   DOUBLE PRECISION :: s    ! distance from the plane
   DOUBLE PRECISION :: smax ! max influence distance 
   DOUBLE PRECISION :: DistanceToPLane
                            ! Distance between a point & a plane in 3D
   LOGICAL :: ret           ! return flag
   LOGICAL :: GetPlaneCoeffs ! function
   DOUBLE PRECISION, DIMENSION(3) :: x0, x1, x2     ! points in a plane
   DOUBLE PRECISION, DIMENSION(3) :: x_proj         ! projected point 
   DOUBLE PRECISION, DIMENSION(2) :: Xmin, Xmax     ! cross-section extents
   DOUBLE PRECISION, DIMENSION(4) :: plane_eq       ! Coeffs of plane equation
   DOUBLE PRECISION, DIMENSION(2) :: x_buff         ! 2D point buffer
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xn_ ! packed nodal coordinates
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xe  !!!! make it a global var !!!!
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: smin  ! history of min distances to cells
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: kp       ! cross-section soil types
   
   
   nn = nx*ny*nz
   ne = (nx-1)*(ny-1)*(nz-1)
   
   ! Pack coordinates into xn_ array 
   IF (.not. ALLOCATED(xn_)) ALLOCATE (xn_(nn,3))
   nr = 0
   DO k = 1, nz
      DO j = 1, ny
         DO i = 1, nx
            nr       = nr + 1
            xn_(nr,1) = x(i,j)
            xn_(nr,2) = y(i,j)
            xn_(nr,3) = z(i,j,k)
         END DO
      END DO
   END DO 

   ! Calculate coordinates of cells centers 
   IF (.not. ALLOCATED(xe)) ALLOCATE (xe(ne,3))
   xe = 0.d0
   DO i = 1, ne
      DO j = 1, 3
         DO k = 1, 8
            xe(i,j) = xe(i,j) + xn_(node(i,k),j)
         END DO
         xe(i,j) = xe(i,j)/8.d0
      END DO
   END DO 
   
   ! Initializations
   IF (.not. ALLOCATED(smin))   ALLOCATE(smin(ne))
   smin = 0.99d0 * big

   ! Open expected cross sections file with '.sec' extension
   lcount = 1
   OPEN (sec_unit,FILE=CHAR(project_pathname)//'.sec',STATUS='old',ERR=20)
   
   WRITE (con_unit,'(20X,a)') "Soil types interpolation ..."
   
   ! Process information from gridded material cross-sections
   i = 0
   cs: DO
      
      !*** Read data for this plane
      i = i + 1 
      
      ! What's we're doing here?
      
      READ(sec_unit,'(a)',ERR=30,END=10) titlep ! Soil profile title 
      lcount = lcount + 1
      WRITE(con_unit,*) "Processing data plane: ", TRIM(titlep)
      
      READ(sec_unit,*,ERR=30) Orient          ! Profile orientation (1=X_orth,2=Y,3=Z)
      lcount = lcount + 1
      READ(sec_unit,*,ERR=30) (x0(j),j=1,3)   ! First  point x0 in that plane 
      lcount = lcount + 1
      READ(sec_unit,*,ERR=30) (x1(j),j=1,3)   ! Second point x1 ... 
      lcount = lcount + 1
      READ(sec_unit,*,ERR=30) (x2(j),j=1,3)   ! Third  point x2 ...
      lcount = lcount + 1
      READ(sec_unit,*,ERR=30) Xmin(1),Xmax(1), & 
                              Xmin(2),Xmax(2) ! Profile extents 
      lcount = lcount + 1
      READ(sec_unit,*,ERR=30) smax            ! Maximal influence distance
                                              ! for NN interpolation 
      lcount = lcount + 1
      READ(sec_unit,*,ERR=30) n1,n2           ! Grid size of the soil material profile
      lcount = lcount + 1
      
      IF (ALLOCATED(kp)) DEALLOCATE(kp)
      ALLOCATE (kp(n1,n2))                    ! Read a unifrom grid 
      READ(sec_unit,*,ERR=30) ((kp(j,k),j=1,n1),k=1,n2)
      lcount = lcount + 1
      
      !*** End of reading records section
        
      ! Calculate plane equation 
      ret = GetPlaneCoeffs(x0,x1,x2,plane_eq)
      ! ... if an error occurs 
      IF (RET == .FALSE.) THEN 
         WRITE(con_unit,*) "WARNING: Points of cross-section plane",i," are colinear!"
         CYCLE cs
      END IF
      
      ! Detect the set of nearby candidate elements and 
      ! fill them with the closest material index
      elem: DO k = 1, ne
         
         ! Calculate distance between this cell center and the plane 
         s = DistanceToPLane(xe(k,:),plane_eq)
      
         ! Calculate the coordinates of the projected point on the plane 
         CALL PlaneProjCoord(xe(k,:),plane_eq,x_proj)
         
         ! Depending on the plane orientation select corresponding 2D point 
         IF (Orient == 1 )     THEN    ! X-orthogonal soil profile 
            x_buff(1) = x_proj(2) ! y 
            x_buff(2) = x_proj(3) ! z 
         ELSE IF (Orient == 2) THEN    ! Y-orthogonal soil profile  
            x_buff(1) = x_proj(1) ! x 
            x_buff(2) = x_proj(3) ! z 
         ELSE IF (Orient == 3) THEN    ! Z-orthogonal soil profile 
            x_buff(1) = x_proj(1) ! x 
            x_buff(2) = x_proj(2) ! y 
         END IF 
      
         ! Get local indices on the soil grid to identify 
         ! the candidate material ID
         CALL LocalIndices(x_buff,Xmin,Xmax,n1,n2,ii,jj)
                      
         ! If within influence area select that soil type
         IF (s < smin(k) .and. s < smax) THEN 
            smin(k) = s
            soilId(k) = kp(ii,jj)
         END IF
      
      END DO elem  ! elements loop

   END DO cs       ! cross-section loop 

10 CONTINUE
   
   ! it's done.
   CLOSE (sec_unit)

   ! Count the number of empty cells (zero soil ID) 
   zeroCount = 0
   DO i = 1, ne 
      IF (smin(i) >= big) THEN
         soilId(i) = 0
         zeroCount = zeroCount + 1
      END IF
   END DO
   
   IF (zeroCount > 0) &
   WRITE(con_unit,*) "WARNING: Number of empty cells = ", zeroCount 
   
20 CONTINUE

   ! Free local vectors 
   IF (ALLOCATED(xn_))      DEALLOCATE (xn_)
   IF (ALLOCATED(xe))       DEALLOCATE (xe)
   IF (ALLOCATED(smin))     DEALLOCATE (smin)
   IF (ALLOCATED(kp))       DEALLOCATE (kp)

   RETURN 

30 CONTINUE
   WRITE(con_unit,*) "ERROR: Unexpected input in file: ", CHAR(project_pathname)//'.sec', & 
                     " at line number: ", lcount
   STOP 
   
END SUBROUTINE InterpolateSoils


LOGICAL FUNCTION GetPlaneCoeffs(pt1,pt2,pt3,a)
! Calculate the four coefficients of plane's equation from given three points.

   IMPLICIT NONE
   ! Coordinates of points P1, P2, P3 
   DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: pt1,pt2,pt3
   DOUBLE PRECISION, DIMENSION(4), INTENT(OUT) :: a

   DOUBLE PRECISION, PARAMETER :: eps = 1.D-12
   INTEGER :: i 
   
   ! a(1:3) are computed from the cross-product of vectors P1P2 & P1P3
   a(1) = (pt2(2)-pt1(2))*(pt3(3)-pt1(3))-(pt3(2)-pt1(2))*(pt2(3)-pt1(3)) 
   a(2) = (pt2(3)-pt1(3))*(pt3(1)-pt1(1))-(pt3(3)-pt1(3))*(pt2(1)-pt1(1)) 
   a(3) = (pt2(1)-pt1(1))*(pt3(2)-pt1(2))-(pt3(1)-pt1(1))*(pt2(2)-pt1(2)) 
   
   ! a(4) given by simple application of the plane equation to point P1
   a(4) = -(a(1)*pt1(1)+a(2)*pt1(2)+a(3)*pt1(3)) 
   
   ! Check if the three points are colinear?
   IF (abs(a(1)) < eps .AND. abs(a(2)) < eps .AND. abs(a(3)) < eps) THEN 
      GetPlaneCoeffs = .FALSE. 
   ELSE 
      GetPlaneCoeffs = .TRUE. 
   END IF  
      
END FUNCTION GetPlaneCoeffs


DOUBLE PRECISION FUNCTION DistanceToPLane(x,a) 
! Calculate the shortest distance between point x and the plane whose equation 
! coefficients are given in array a. 

   IMPLICIT NONE 
   DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: x 
   DOUBLE PRECISION, DIMENSION(4), INTENT(IN) :: a 

   ! This is the formula used to calculate the normal distance, d, between 
   ! point (x',y',z') and the plane a*x + b*y + c*z + d = 0 :
   !              | a*x' + b*y' + c*z' + d |
   !       d = --------------------------------
   !               (a^2 + b^2 + c^2)^(1/2)
   DistanceToPLane = ABS(a(1)*x(1)+a(2)*x(2)+a(3)*x(3)+a(4))
   DistanceToPLane = DistanceToPLane/SQRT(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

END FUNCTION DistanceToPLane


SUBROUTINE PlaneProjCoord(xc,a,x)
! Calculate spatial coordinates of the projection of point xc, which is given 
! array x, on the plane whose equation coefficients are given in array a.

   IMPLICIT NONE 
   DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: xc 
   DOUBLE PRECISION, DIMENSION(4), INTENT(IN)  :: a
   DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: x
   
   DOUBLE PRECISION :: t ! parametric coordinate of point xc
   
   ! The following formula results from insertion of the three parametric 
   ! equations (forming the line joining the point xc and its projection x)
   ! into the equation of the plane. This simple manipulation will determine 
   ! the parameter t, such that this point will be, at the same time, on the 
   ! line and the plane.  
   t = -(a(1)*xc(1)+a(2)*xc(2)+a(3)*xc(3))
   t = t/(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
   
   ! Finally, use the parametric equations to get the coordinates 
   ! of projected point in the plane.
   x(1:3) = xc(1:3)+a(1:3)*t 

   RETURN 
END SUBROUTINE PlaneProjCoord


SUBROUTINE LocalIndices(x,Xmin,Xmax,ni,nj,i,j) 
! Returns local indices of a given point in a gridded cross-sectional 
! soil types plane. 

   IMPLICIT NONE 
   DOUBLE PRECISION, DIMENSION(2), INTENT(IN)  :: x
   DOUBLE PRECISION, DIMENSION(2), INTENT(IN)  :: Xmin
   DOUBLE PRECISION, DIMENSION(2), INTENT(IN)  :: Xmax
   INTEGER, INTENT(IN)  :: ni,nj 
   INTEGER, INTENT(OUT) :: i,j
   
   DOUBLE PRECISION :: dx, dy 
   
   dx = (Xmax(1)-Xmin(1))/ni    ! uniform spacing is assumed
   dy = (Xmax(2)-Xmin(2))/nj    ! idem for y-spacing 
   
   i = int(x(1)/dx) + 1         ! local i index 
   j = int(x(2)/dy) + 1         ! local j index 
   
   RETURN 
END SUBROUTINE LocalIndices