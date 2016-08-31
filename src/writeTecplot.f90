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

Subroutine WriteTecplotFile()

   USE GridGlobals

   IMPLICIT NONE 
   INTEGER :: i,j,k,nr,ne
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: XX
   
   OPEN (tec_unit,file=CHAR(project_pathname)//'.plt',status='unknown')
   
   ! Write the tecplot file header      
   WRITE(tec_unit,*) 'TITLE = "',title(1),'"'
   WRITE(tec_unit,*) 'VARIABLES = "X", "Y", "Z", "soilId"'
   WRITE(tec_unit,*) 'ZONE T="',title(2),'"'
   WRITE(tec_unit,*) 'I=', nx,', J=',ny,', K=',nz
   WRITE(tec_unit,*) 'ZONETYPE = ORDERED, DATAPACKING = BLOCK'
   WRITE(tec_unit,*) 'VARLOCATION = ([4] = CELLCENTERED)'
   
   ! Write x,y,z sequences 
   CALL WRITE_XY_SEQ('X')
   CALL WRITE_XY_SEQ('Y')
   CALL WRITE_XY_SEQ('Z')
   
   ! Write soil ID's in cells
   ne = (nx-1)*(ny-1)*(nz-1)
   nr = 1 
   l1: DO 
      IF (ne-nr >= 2) THEN 
         WRITE(tec_unit,'(3I6)') soilId(nr), soilId(nr+1), soilId(nr+2) 
         nr = nr + 3 
      ELSE IF  (ne-nr == 1) THEN 
         WRITE(tec_unit,'(3I6)') soilId(nr), soilId(nr+1)
         EXIT l1
      ELSE
         WRITE(tec_unit,'(3I6)')  soilId(nr)
         EXIT l1 
      END IF
   END DO l1 
      
   ! close the file
   CLOSE (tec_unit)
   
RETURN 
END SUBROUTINE WriteTecplotFile


SUBROUTINE WRITE_XY_SEQ(var)

   USE GridGlobals
   
   IMPLICIT NONE 
   CHARACTER(LEN=1), INTENT(IN) :: var
   
   INTEGER :: i,j,k,nr
   DOUBLE PRECISION :: v 
   

   nr = 0
   DO k = 1, nz
   DO j = 1, ny 
   DO i = 1, nx 
      
      nr = nr + 1
      
      IF (var .EQ. 'X') v = x(i,j) 
      IF (var .EQ. 'Y') v = y(i,j) 
      IF (var .EQ. 'Z') v = z(i,j,k) 
      
      IF (nr == 3) THEN 
         WRITE(tec_unit,'(EN16.3)') v 
         nr = 0 
      ELSE
         WRITE(tec_unit,'(EN16.3)',ADVANCE='NO') v
      END IF 
   
   END DO
   END DO
   END DO

RETURN 
END SUBROUTINE WRITE_XY_SEQ  