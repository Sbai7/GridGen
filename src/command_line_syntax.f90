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

SUBROUTINE COMMAND_LINE_SYNTAX()
! Prints a list of program supported command line switches
   
   USE GridGlobals
   IMPLICIT NONE 
   
   WRITE (con_unit,*) &
   "GRIDGEN: Command-line Grid Generation Program and Much More - Version ", &
   GridGen_version
   
   WRITE (con_unit,*) "GRIDGEN [-ih] <input-file>"
   WRITE (con_unit,*)
   WRITE (con_unit,*) " <input-file> a string identifier of the project datafile"
   WRITE (con_unit,*) "   -i  To specify the input project file"
   WRITE (con_unit,*) "   -h  To get this help screen"
   
   STOP
END SUBROUTINE COMMAND_LINE_SYNTAX
