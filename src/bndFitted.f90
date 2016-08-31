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

SUBROUTINE BndFittedMesh()
! GridGen: A 2D/3D GEOHYDROLOGICAL GRIDDING PACKAGE

   USE GridGlobals

   IMPLICIT NONE
   CHARACTER(LEN=5)  :: test
   INTEGER, ALLOCATABLE, DIMENSION(:) :: nlayers
   INTEGER :: i, j, k, nr, nxy, ne
   DOUBLE PRECISION, PARAMETER :: tol = 0.01d0


   WRITE (con_unit,'(20X,a)') "Generation of the finite element grid ..."

!**** XY-GRID **********************************************************

   ! Read general input
   WRITE (con_unit,*) "Reading input file..."
   OPEN (grd_unit,file=CHAR(project_pathname)//'.grd',status='old')
   DO i = 1, 3
      READ (grd_unit,'(a)') title(i) ! replace this by variable length strings
   END DO
   READ (grd_unit,*) nx, ny, nz
   IF (nx < 2 .or. ny < 2 .or. nz < 2) THEN
      WRITE(con_unit,*) "ERROR in GRIDGEN inputs: bad dimensions!"
      STOP
   END IF

   ! Initialization
   nxy = nx*ny
   ne  = (nx-1)*(ny-1)*(nz-1)

   ! Dynamic Allocation
   IF (.not.ALLOCATED(in))   ALLOCATE(in(nxy))
   IF (.not.ALLOCATED(jn))   ALLOCATE(jn(nxy))
   IF (.not.ALLOCATED(xn))   ALLOCATE(xn(nxy))
   IF (.not.ALLOCATED(yn))   ALLOCATE(yn(nxy))
   IF (.not.ALLOCATED(join)) ALLOCATE(join(nxy))
   
   IF (.not.ALLOCATED(node))    ALLOCATE(node(ne,8))
   IF (.not. ALLOCATED(soilId)) ALLOCATE(soilId(ne))

   ! Input of fixed xy-coordinates
   WRITE(con_unit,*) "Gridding in the XY plane"
   DO i = 1, nxy + 1
      READ (grd_unit,'(a5)',END=10) test
      BACKSPACE grd_unit
      IF (test == "layer") GOTO 10
      READ (grd_unit,*) in(i),jn(i),xn(i),yn(i),join(i)
      IF (in(i) > nx) THEN
         WRITE(con_unit,*) "ERROR in GRIDGEN input: bad i-index " ! must also give the line number where the error occurs
      END IF
      IF (jn(i) > ny) THEN
         WRITE(con_unit,*) "ERROR in GRIDGEN input: bad j-index"  ! idem ...
      END IF
   END DO
   WRITE(con_unit,*) "Warning: too much xy-nodes specified"
10 nr = i-1

   ! Calculation of fixed xy-coordinates
   IF (.not.ALLOCATED(x))    ALLOCATE(x(nx,ny))
   IF (.not.ALLOCATED(y))    ALLOCATE(y(nx,ny))
   IF (.not.ALLOCATED(kode)) ALLOCATE(kode(nx,ny))

   ! Initializations
   x(1:nx,1:ny) = 0.d0
   y(1:nx,1:ny) = 0.d0
   kode(1:nx,1:ny) = 0
   
   soilId(1:ne) = 0

   CALL xyFixed(nr)

   ! Calculation of boundary xy-coordinates
   CALL xyBord(tol)

   ! Calculation of free coordinates in the xy-plane
   CALL xyFree(tol)

!**** Projection of Z-LAYERS *******************************************

   WRITE(con_unit,*) "Projection to Z-layers"
   IF (.not.ALLOCATED(nlayers)) ALLOCATE(nlayers(nz))
   nlayers(1:nz) = 0

   if (.not.ALLOCATED(xp)) ALLOCATE(xp(nxy,nz))
   if (.not.ALLOCATED(yp)) ALLOCATE(yp(nxy,nz))
   if (.not.ALLOCATED(zp)) ALLOCATE(zp(nxy,nz))
   if (.not.ALLOCATED(np)) ALLOCATE(np(nz))

   ! Input of fixed z-coordinates for each plane
   DO i = 1, nz
      READ (grd_unit,'(5x,i10)',END=30,ERR=30) nlayers(i)
      DO j = 1, nxy
         READ(grd_unit,'(a5)',end=20) test
         BACKSPACE grd_unit
         IF (test == "layer") GOTO 20
         READ (grd_unit,*) xp(j,i),yp(j,i),zp(j,i) 
      END DO
      WRITE(con_unit,*) "Warning: too much z-data specified"
20    np(i) = j-1
   END DO
30 nr = i-1

   IF (nlayers(1) /= 1) THEN
      WRITE(con_unit,*) "ERROR: first layer not specified"
      STOP
   END IF
   IF (nlayers(nr) /= nz) THEN
      WRITE(con_unit,*) "ERROR: last layers not specified"
      STOP
   END IF
   CLOSE (grd_unit)

   ! Calculation of z-coordinates
   IF (.not.ALLOCATED(z)) ALLOCATE(z(nx,ny,nz))
   CALL zFixed(1,nlayers(1))
   DO i = 2, nr
      CALL zFixed(i,nlayers(i))
      CALL zFree(nlayers(i-1),nlayers(i))
   END DO
   
   ! Builds FE connectivity table
   CALL BuildConnectivity()
   
   ! Interpolate soils from cross section planes 
   ! when the associated files exists 
   CALL InterpolateSoils()

!**** OUTPUT SECTION ***************************************************

   WRITE (con_unit,'(20X,a)') "Generation of output files ..."
   
   WRITE(con_unit,*) "Writing native mesh file"
   CALL WriteMeshFile()

   WRITE(con_unit,*) "Writing mesh to Tecplot file ..."
   CALL WriteTecplotFile()
   
   ! Free allocated memory -- global variables --
   IF (ALLOCATED(x))    DEALLOCATE(x)
   IF (ALLOCATED(y))    DEALLOCATE(y)
   IF (ALLOCATED(z))    DEALLOCATE(z)
   IF (ALLOCATED(kode)) DEALLOCATE(kode)
   IF (ALLOCATED(in))   DEALLOCATE(in)
   IF (ALLOCATED(jn))   DEALLOCATE(jn)
   IF (ALLOCATED(xn))   DEALLOCATE(xn)
   IF (ALLOCATED(yn))   DEALLOCATE(yn)
   IF (ALLOCATED(join)) DEALLOCATE(join)
   IF (ALLOCATED(xp))   DEALLOCATE(xp)
   IF (ALLOCATED(yp))   DEALLOCATE(yp)
   IF (ALLOCATED(zp))   DEALLOCATE(zp)
   IF (ALLOCATED(np))   DEALLOCATE(np)
   
   IF (ALLOCATED(node))     DEALLOCATE(node)
   IF (ALLOCATED(soilId))   DEALLOCATE (soilId)

   ! -- local variables --
   IF (ALLOCATED(nlayers)) DEALLOCATE(nlayers)

   RETURN 
END SUBROUTINE BndFittedMesh