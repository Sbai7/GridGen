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

SUBROUTINE xyBord(tolerance)
! SUBROUTINE FOR CALCULATING BOUNDARY XY-COORDINATES
! We may also send the number of iterations to solve the Laplacian problem
! which might be useful to identify possible problems (bad coordinates, etc)

   USE GridGlobals, only : x,y,kode,nx,ny

   IMPLICIT NONE
   double precision, intent(in) :: tolerance

   double precision, allocatable, dimension(:) :: xseq, yseq 
   integer, allocatable, dimension(:) :: kodeseq 

   integer :: i
   double precision :: emax, er, eabs

   
   ! Allocate local dynamic arrays
   if (.not.allocated(xseq))    allocate(xseq(2*(nx+ny)))
   if (.not.allocated(yseq))    allocate(yseq(2*(nx+ny)))
   if (.not.allocated(kodeseq)) allocate(kodeseq(2*(nx+ny)))

   xseq(1:nx) = x(1:nx,1)
   yseq(1:nx) = y(1:nx,1)
   kodeseq(1:nx) = kode(1:nx,1)

   xseq(nx+1:nx+ny-1) = x(nx,2:ny)
   yseq(nx+1:nx+ny-1) = y(nx,2:ny)
   kodeseq(nx+1:nx+ny-1) = kode(nx,2:ny)

   do i = nx-1, 1, -1
      xseq(nx-1+ny+(nx-i)) = x(i,ny)
      yseq(nx-1+ny+(nx-i)) = y(i,ny)
      kodeseq(nx-1+ny+(nx-i)) = kode(i,ny)
   end do

   do i = ny-1, 1, -1
      xseq(nx-1+ny-1+nx+(ny-i)) = x(1,i)
      yseq(nx-1+ny-1+nx+(ny-i)) = y(1,i)
      kodeseq(nx-1+ny-1+nx+(ny-i)) = kode(1,i)
   end do

   xseq(2*(nx+ny-1)) = x(2,1)
   yseq(2*(nx+ny-1)) = y(2,1)
   kodeseq(2*(nx+ny-1)) = kode(2,1)

   do
      emax = 0.d0
      do i = 2, 2*(nx+ny)-3
         if (kodeseq(i) /= 1) then
            er = (xseq(i-1)+xseq(i+1))/2. - xseq(i)
            eabs = abs(er)
            if (emax < eabs) emax = eabs
            xseq(i) = xseq(i) + 1.85*er  
            er = (yseq(i-1)+yseq(i+1))/2. - yseq(i)
            eabs = abs(er)
            if (emax < eabs) emax = eabs
            yseq(i) = yseq(i) + 1.85*er  
         end if
      end do

      xseq(1) = xseq(2*nx+2*ny-3)   
      yseq(1) = yseq(2*nx+2*ny-3)   
      kodeseq(1) = kodeseq(2*nx+2*ny-3)   
      xseq(2*(nx+ny-1)) = xseq(2)
      yseq(2*(nx+ny-1)) = yseq(2)
      kodeseq(2*(nx+ny-1)) = kodeseq(2)

      if (emax < tolerance) exit
   end do

   x(1:nx,1) = xseq(1:nx)
   y(1:nx,1) = yseq(1:nx)
   kode(1:nx,1) = kodeseq(1:nx)

   x(nx,2:ny) = xseq(nx+1:nx+ny-1)
   y(nx,2:ny) = yseq(nx+1:nx+ny-1)
   kode(nx,2:ny) = kodeseq(nx+1:nx+ny-1)

   do i = nx-1, 1, -1
      x(i,ny) = xseq(nx-1+ny+(nx-i))
      y(i,ny) = yseq(nx-1+ny+(nx-i))
      kode(i,ny) = kodeseq(nx-1+ny+(nx-i))
   end do

   do i = ny-1, 1, -1
      x(1,i) = xseq(nx-1+ny-1+nx+(ny-i))
      y(1,i) = yseq(nx-1+ny-1+nx+(ny-i))
      kode(1,i) = kodeseq(nx-1+ny-1+nx+(ny-i))
   end do

   x(2,1) = xseq(2*(nx+ny-1))
   y(2,1) = yseq(2*(nx+ny-1))
   kode(2,1) = kodeseq(2*(nx+ny-1))

   ! Deallocation
   if (allocated(xseq))     deallocate(xseq)
   if (allocated(yseq))     deallocate(yseq)
   if (allocated(kodeseq))  deallocate(kodeseq)

return
END SUBROUTINE xyBord
