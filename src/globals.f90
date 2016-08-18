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

MODULE GridGlobals

   USE ISO_VARYING_STRING

   TYPE(VARYING_STRING) :: project_pathname
   CHARACTER(LEN=48) :: title(3)

   ! Needed Grid dimensions along each axis
   integer :: nx
   integer :: ny
   integer :: nz

   ! Allocatable Arrays
   double precision, allocatable, dimension(:,:) :: x
   double precision, allocatable, dimension(:,:) :: y
   integer, allocatable, dimension(:,:) :: kode

   integer, allocatable, dimension(:) :: in
   integer, allocatable, dimension(:) :: jn
   integer, allocatable, dimension(:) :: join
   double precision, allocatable, dimension(:) :: xn
   double precision, allocatable, dimension(:) :: yn

   double precision, allocatable, dimension(:,:,:) :: z

   double precision, allocatable, dimension(:,:) :: xp
   double precision, allocatable, dimension(:,:) :: yp
   double precision, allocatable, dimension(:,:) :: zp
   integer, allocatable, dimension(:) :: np

END MODULE GridGlobals
