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

SUBROUTINE xyFree(tolerance)
! SUBROUTINE FOR CALCULATING FREE XY-COORDINATES
! We may return the number of iterations required for convergence.
! Another important parameter for the algorithm is the relaxation 
! factor (now fixed to 1.85), we can experiment different values
! to find an optimum (for each grid size).
      
   USE GridGlobals, only : x,y,kode,nx,ny

   IMPLICIT NONE
   double precision, intent(in) :: tolerance

   integer :: i,j
   double precision :: emax,er,eabs


   do
      emax = 0.d0

      do i = 2, nx-1
         do j = 2, ny-1
            if (kode(i,j) /= 1) then
               er = (x(i-1,j)+x(i,j-1)+x(i+1,j)+x(i,j+1))/4.d0 - x(i,j)
               eabs = abs(er)
               if (emax < eabs) emax = eabs
               x(i,j) = x(i,j) + 1.85d0*er  
               er = (y(i-1,j)+y(i,j-1)+y(i+1,j)+y(i,j+1))/4.d0 - y(i,j)
               eabs = abs(er)
               if (emax < eabs) emax = eabs
               y(i,j) = y(i,j) + 1.85d0*er
            end if 
         end do
      end do

      if (emax < tolerance) exit 
   end do

return      
END SUBROUTINE xyFree
