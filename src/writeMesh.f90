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

SUBROUTINE WriteMeshFile()
! OUTPUT OF THE MESH TO A DATA FILE
  
   USE ISO_VARYING_STRING
   USE GridGlobals

   IMPLICIT NONE

   INTEGER :: i,j,k,kg,nn,ne,nul
   DOUBLE PRECISION :: zero

!**** MAIN DATA ********************************************************

   OPEN (msh_unit,file=CHAR(project_pathname)//'.dat')

   DO i = 1, 3
      WRITE(dat_unit,'(a,i2,t30,a,a)') "Title",i,": ",title(i)
   END DO
   nn = nx*ny*nz
   WRITE (dat_unit,'(a,t30,a,4i10)') "Number of nodes",":",nn,nx,ny,nz
   ne = (nx-1)*(ny-1)*(nz-1)
   WRITE (dat_unit,'(a,t30,a,i10)') "Number of elements",":",ne      
   nul = 0
   zero = 0.d0
   WRITE (dat_unit,'(a,t30,a,i10)') "Number of soil types",":",nul ! UPDATE
   !AS WRITE (dat_unit,'(a,t30,a,3g10.3)') "Time parameters",":",zero,zero,zero

   CLOSE(dat_unit)

!**** NODES AND ELEMENT DATA *******************************************

   open (msh_unit,FILE=CHAR(project_pathname)//'.fem',STATUS="UNKNOWN")

   WRITE(msh_unit,'(i8,a40)') nx, "# Number of nodes in x-direction"
   WRITE(msh_unit,'(i8,a40)') ny, "# Number of nodes in y-direction"
   WRITE(msh_unit,'(i8,a40)') nz, "# Number of nodes in z-direction"
   WRITE(msh_unit,'(i8,a40)') nn, "# Number of nodes"
   WRITE(msh_unit,'(i8,a40)') ne, "# Number of cells"
   
   ! Nodal coordinates 
   DO k = 1, nz
      DO j = 1, ny
         DO i = 1, nx
            WRITE(msh_unit,'(3EN16.3)') x(i,j), y(i,j), z(i,j,k)
         END DO 
      END DO 
   END DO 
   
   ! Mesh FE connectivity table + element properties ...  
   DO i = 1, ne
      WRITE(msh_unit,'(8i8,i8)') node(i,1:8), soilId(i) 
   END DO 
   
   CLOSE(msh_unit)

   RETURN
END SUBROUTINE WriteMeshFile
