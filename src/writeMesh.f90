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

SUBROUTINE WriteMeshFile()
! OUTPUT OF THE MESH TO A DATA FILE
  
   USE ISO_VARYING_STRING
   USE GridGlobals, only : x,y,z,nx,ny,nz,project_pathname,title

   IMPLICIT NONE

   integer, allocatable, dimension(:,:) :: node

   integer :: i,j,k,kg,nr,nn,ne,nul
   double precision :: zero

!**** MAIN DATA ********************************************************

   open (3,file=CHAR(project_pathname)//'.dat')

   do i = 1, 3
      write(3,'(a,i2,t30,a,a)') "Title",i,": ",title(i)
   end do
   nn = nx*ny*nz
   write (3,'(a,t30,a,4i10)') "Number of nodes",":",nn,nx,ny,nz
   ne = (nx-1)*(ny-1)*(nz-1)
   write (3,'(a,t30,a,i10)') "Number of elements",":",ne      
   nul = 0
   zero = 0.d0
   write (3,'(a,t30,a,i10)') "Number of soil types",":",nul
   write (3,'(a,t30,a,3g10.3)') "Time parameters",":",zero,zero,zero

   close(3)

!**** NODES AND ELEMENT DATA *******************************************

   open (3,FILE=CHAR(project_pathname)//'.fem',STATUS="UNKNOWN")

   WRITE(3,'(i8,a40)') nx, "# Number of nodes in x-direction"
   WRITE(3,'(i8,a40)') ny, "# Number of nodes in y-direction"
   WRITE(3,'(i8,a40)') nz, "# Number of nodes in z-direction"
   WRITE(3,'(i8,a40)') nn, "# Number of nodes"
   WRITE(3,'(i8,a40)') ne, "# Number of cells"
   
   DO k = 1, nz
      DO j = 1, ny
         DO i = 1, nx
            WRITE(3,'(3EN16.3)') x(i,j), y(i,j), z(i,j,k)
         END DO 
      END DO 
   END DO 
   
   IF (.not.ALLOCATED(node)) ALLOCATE(node((nx-1)*(ny-1)*(nz-1),8))
   kg = 0
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
   
   DO i = 1, ne
      WRITE(3,'(8i8,i8)') node(i,1:8), kg 
   END DO 
   CLOSE(3)
   IF (ALLOCATED(node)) DEALLOCATE(node)

return
END SUBROUTINE WriteMeshFile
