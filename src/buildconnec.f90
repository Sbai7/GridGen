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

SUBROUTINE BuildConnectivity()
! Builds connectivity table of the FE mesh 

   USE GridGlobals

   IMPLICIT NONE
   INTEGER :: i,j,k, nr

   nr = 0
   DO k = 1, nz-1
      DO j = 1, ny-1
         DO i = 1, nx-1
            nr = nr + 1
            node(nr,1) = i + (j-1)*nx + (k-1)*nx*ny
            node(nr,2) = node(nr,1) + nx
            node(nr,3) = node(nr,1) + nx + 1
            node(nr,4) = node(nr,1) + 1
            node(nr,5) = node(nr,1) + nx*ny
            node(nr,6) = node(nr,5) + nx
            node(nr,7) = node(nr,5) + nx + 1
            node(nr,8) = node(nr,5) + 1
         END DO
      END DO
   END DO

RETURN 
END SUBROUTINE BuildConnectivity