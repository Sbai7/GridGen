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

SUBROUTINE zFree(nl1,nl2)
! SUBROUTINE FOR CALCULATING FREE Z-LAYERS 

   USE GridGlobals, only : nx,ny,z

   IMPLICIT NONE
   integer, intent(in) :: nl1
   integer, intent(in) :: nl2

   integer :: nl,i,j
   double precision :: a


   do nl = nl1, nl2
      do i = 1, nx
         do j = 1, ny
            a = real(nl-nl1)/(nl2-nl1)
            z(i,j,nl) = z(i,j,nl1)*(1-a)+z(i,j,nl2)*a
         end do
      end do
   end do

return
END SUBROUTINE zFree
