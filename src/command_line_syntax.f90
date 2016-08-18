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

SUBROUTINE COMMAND_LINE_SYNTAX()
! syntax: Prints a list of command line switches
   write (6,*) "GRIDGEN: Boundary-fitted FEM grid generator for groundwater flow models - Version 1.0"
   write (6,*) "GRIDGEN [-ih] <input-file>"
   write (6,*)
   write (6,*) " <input-file> a string identifier of the project datafile"
   write (6,*) "   -i  To specify the input project file"
   write (6,*) "   -h  To get this help screen"
   STOP
END SUBROUTINE COMMAND_LINE_SYNTAX
