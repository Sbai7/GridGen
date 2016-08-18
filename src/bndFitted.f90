!******************************************************************************
!*
!* This file is part of GwMove hydrogeological software developed by
!* Dr. M. A. Sbai
!*
!* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
!* modification, are permitted provided that the following conditions are met:
!*
!*  - Redistributions of  source code must  retain the above  copyright notice,
!*    this list of conditions and the disclaimer below.
!*  - Redistributions in binary form must reproduce the above copyright notice,
!*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
!*    documentation and/or other materials provided with the distribution.
!*
!******************************************************************************

SUBROUTINE BndFittedMesh()
! GridGen: A 2D/3D GEOHYDROLOGICAL GRIDDING PACKAGE

   USE GridGlobals

   IMPLICIT NONE
   CHARACTER(LEN=5)  :: test
   INTEGER, ALLOCATABLE, DIMENSION(:) :: nlayers
   INTEGER :: i, j, k, nr, nxy
   DOUBLE PRECISION, PARAMETER :: tol = 0.01d0


   WRITE (6,'(20X,a)') "Generation of the finite element grid ..."

!**** XY-GRID **********************************************************

   ! Read general input
   WRITE (6,*) "Reading input file..."
   OPEN (2,file=CHAR(project_pathname)//'.grd',status='old')
   DO i = 1, 3
      READ (2,'(a)') title(i) ! replace this by variable length strings
   END DO
   READ (2,*) nx, ny, nz
   IF (nx < 2 .or. ny < 2 .or. nz < 2) THEN
      WRITE(6,*) "ERROR in GRIDGEN inputs: bad dimensions!"
      STOP
   END IF

   ! Initialization
   nxy = nx*ny

   ! Dynamic Allocation
   IF (.not.ALLOCATED(in))   ALLOCATE(in(nxy))
   IF (.not.ALLOCATED(jn))   ALLOCATE(jn(nxy))
   IF (.not.ALLOCATED(xn))   ALLOCATE(xn(nxy))
   IF (.not.ALLOCATED(yn))   ALLOCATE(yn(nxy))
   IF (.not.ALLOCATED(join)) ALLOCATE(join(nxy))

   ! Input of fixed xy-coordinates
   WRITE(6,*) "Gridding in the XY plane"
   DO i = 1, nxy + 1
      READ (2,'(a5)',END=10) test
      BACKSPACE 2
      IF (test == "layer") GOTO 10
      READ (2,*) in(i),jn(i),xn(i),yn(i),join(i)
      IF (in(i) > nx) THEN
         WRITE(6,*) "ERROR in GRIDGEN input: bad i-index " ! must also give the line number where the error occurs
      END IF
      IF (jn(i) > ny) THEN
         WRITE(6,*) "ERROR in GRIDGEN input: bad j-index"  ! idem ...
      END IF
   END DO
   WRITE(6,*) "Warning: too much xy-nodes specified"
10 nr = i-1

   ! Calculation of fixed xy-coordinates
   IF (.not.ALLOCATED(x))    ALLOCATE(x(nx,ny))
   IF (.not.ALLOCATED(y))    ALLOCATE(y(nx,ny))
   IF (.not.ALLOCATED(kode)) ALLOCATE(kode(nx,ny))

   ! initialization
   x(1:nx,1:ny) = 0.d0
   y(1:nx,1:ny) = 0.d0
   kode(1:nx,1:ny) = 0

   CALL xyFixed(nr)

   ! Calculation of boundary xy-coordinates
   CALL xyBord(tol)

   ! Calculation of free coordinates in the xy-plane
   CALL xyFree(tol)

!**** Projection of Z-LAYERS *******************************************

   WRITE(6,*) "Projection to Z-layers"
   IF (.not.ALLOCATED(nlayers)) ALLOCATE(nlayers(nz))
   nlayers(1:nz) = 0

   if (.not.ALLOCATED(xp)) ALLOCATE(xp(nxy,nz))
   if (.not.ALLOCATED(yp)) ALLOCATE(yp(nxy,nz))
   if (.not.ALLOCATED(zp)) ALLOCATE(zp(nxy,nz))
   if (.not.ALLOCATED(np)) ALLOCATE(np(nz))

   ! Input of fixed z-coordinates for each plane
   DO i = 1, nz
      READ (2,'(5x,i10)',end=30,err=30) nlayers(i)
      DO j = 1, nxy
         READ(2,'(a5)',end=20) test
         BACKSPACE 2
         IF (test == "layer") GOTO 20
         READ (2,*) xp(j,i),yp(j,i),zp(j,i) 
      END DO
      WRITE(6,*) "Warning: too much z-data specified"
20    np(i) = j-1
   END DO
30 nr = i-1

   IF (nlayers(1) /= 1) THEN
      WRITE(6,*) "ERROR: first layer not specified"
      STOP
   END IF
   IF (nlayers(nr) /= nz) THEN
      WRITE(6,*) "ERROR: last layers not specified"
      STOP
   END IF
   CLOSE (2)

   ! Calculation of z-coordinates
   IF (.not.ALLOCATED(z)) ALLOCATE(z(nx,ny,nz))
   CALL zFixed(1,nlayers(1))
   DO i = 2, nr
      CALL zFixed(i,nlayers(i))
      CALL zFree(nlayers(i-1),nlayers(i))
   END DO

!**** OUTPUT SECTION ***************************************************

   WRITE(6,*) "Writing native mesh file"
   CALL WriteMeshFile()

   WRITE(6,*) "Writing mesh to Tecplot file ..."
   ! Write the tecplot file header      
   OPEN (9,file=CHAR(project_pathname)//'.plt',status='unknown')
   WRITE(9,*) 'TITLE = "',title(1),'"'
   ! Write the ordered zone records
   WRITE(9,*) 'VARIABLES = "x", "y", "z"'
   WRITE(9,*) 'ZONE T="ordered zone", I=', nx,', J=',ny,', K=',nz,', F=POINT'
   DO k = 1, nz
      DO j = 1, ny
         DO i = 1, nx
            WRITE(9,'(3EN16.3)') x(i,j), y(i,j), z(i,j,k)
         END DO 
      END DO 
   END DO 
   CLOSE (9)


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

   ! -- local variables --
   IF (ALLOCATED(nlayers)) DEALLOCATE(nlayers)

   RETURN 
END SUBROUTINE BndFittedMesh
