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

MODULE GridGlobals

   USE ISO_VARYING_STRING
   
   CHARACTER(LEN=5), PARAMETER :: GridGen_version = "2.0.0"

   TYPE(VARYING_STRING) :: project_pathname
   CHARACTER(LEN=48) :: title(3)

   ! Needed Grid dimensions along each axis
   INTEGER :: nx
   INTEGER :: ny
   INTEGER :: nz

   ! Allocatable Arrays
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: y
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: kode

   INTEGER, ALLOCATABLE, DIMENSION(:) :: in
   INTEGER, ALLOCATABLE, DIMENSION(:) :: jn
   INTEGER, ALLOCATABLE, DIMENSION(:) :: join
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xn
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: yn

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: z

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xp
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: yp
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: zp
   INTEGER, ALLOCATABLE, DIMENSION(:) :: np
   
   ! Code for elmeents soil types 
   INTEGER, ALLOCATABLE, DIMENSION(:) :: soilId
   
   ! ***
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: node
    
    ! Files units numbers 
    INTEGER, PARAMETER :: con_unit = 6
    INTEGER, PARAMETER :: grd_unit = 20
    INTEGER, PARAMETER :: msh_unit = 12
    INTEGER, PARAMETER :: sec_unit = 13
    INTEGER, PARAMETER :: dat_unit = 14
    INTEGER, PARAMETER :: tec_unit = 15

END MODULE GridGlobals