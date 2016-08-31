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

SUBROUTINE zFixed(n,nl) 
! SUBROUTINE FOR CALCULATING FIXED Z-LAYERS 

   USE GridGlobals

   IMPLICIT NONE
   integer, intent(in) :: n
   integer, intent(in) :: nl

   integer :: i,j,k
   double precision :: up,down,s


   do i = 1, nx
      jloop: do j = 1, ny
         up = 0.d0
         down = 0.d0
         do k = 1, np(n)
            s = sqrt((x(i,j)-xp(k,n))**2+(y(i,j)-yp(k,n))**2)
            if (s == 0.d0) then
               z(i,j,nl) = zp(k,n)
               cycle jloop
            else
               up = up + zp(k,n)/s
               down = down+1/s
            end if
         end do
         z(i,j,nl) = up/down
      end do jloop
   end do

return
END SUBROUTINE zFixed
