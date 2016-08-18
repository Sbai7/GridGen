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
