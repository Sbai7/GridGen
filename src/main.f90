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

PROGRAM GRID_GENERATOR
! Main program for the grid generation pre-processor program 

  USE GridGlobals
  
  IMPLICIT NONE
  
  CALL PARSE_COMMAND_LINE(project_pathname)
  WRITE(6,*) "Project = ", CHAR(project_pathname)
  
  CALL BndFittedMesh()
  
END PROGRAM GRID_GENERATOR