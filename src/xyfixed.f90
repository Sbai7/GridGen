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

SUBROUTINE xyFixed(nr)
! SUBROUTINE FOR CALCULATING FIXED XY-COORDINATES

   USE GridGlobals, only : x,y,kode,in,jn,xn,yn,join,nx,ny

   IMPLICIT NONE
   integer, intent(in) :: nr

   integer :: l,id,jd,ida,jda,k,ij

   
   l = 0
   loop: do
      l = l + 1
      if (l > nr) exit loop

      if (join(l) == 1) then
         ! Join the nodes
         id = in(l+1)-in(l)
         jd = jn(l+1)-jn(l)
         ida = abs(id)
         jda = abs(jd)
         do k = 1, max(ida,jda)
            if (mod(ida,k) == 0 .and. mod(jda,k) == 0) ij = k
         end do
         do k = 0, ij
            x(in(l)+id*k/ij,jn(l)+jd*k/ij) = xn(l) + k*(xn(l+1)-xn(l))/ij
            y(in(l)+id*k/ij,jn(l)+jd*k/ij) = yn(l) + k*(yn(l+1)-yn(l))/ij
            kode(in(l)+id*k/ij,jn(l)+jd*k/ij) = 1
         end do 
      else
         ! Isolated node
         x(in(l),jn(l)) = xn(l)
         y(in(l),jn(l)) = yn(l)
         kode(in(l),jn(l)) = 1
      end if
   end do loop

return
END SUBROUTINE xyFixed
