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

SUBROUTINE PARSE_COMMAND_LINE(project_name)
! PARSE_COMMAND_LINE: READS THE COMMAND LINE, IDENTIFY SWITCHES,
! AND SETUP COMMAND LINE OPTIONS.
  USE ISO_VARYING_STRING
  IMPLICIT NONE
  TYPE(VARYING_STRING), INTENT(OUT) :: project_name 
  CHARACTER(len=1000) :: name
  INTEGER :: argc
  INTEGER*2 :: i, icount
  LOGICAL :: found
  REAL :: NARGS
  INTRINSIC GETARG
  INTRINSIC ADJUSTL
  INTRINSIC TRIM
  argc = COMMAND_ARGUMENT_COUNT ()
! gets command-line arguments
! for instance there is only one, which is
! the project filename, including the path ... 
  IF (argc .LT. 1) CALL COMMAND_LINE_SYNTAX()
  icount = 0
  found = .false.
arguments:DO i=1,argc
    CALL GET_COMMAND_ARGUMENT(i,name)
    IF (name(1:1) .EQ. '-') THEN
! command line switch
      icount = icount + 1
      IF (name(2:2) .EQ. 'i') THEN
        found = .true.
      ELSE IF (name(2:2) .EQ. 'h' .OR. name(2:2) .EQ. '') THEN
        CALL COMMAND_LINE_SYNTAX()
      END IF
    ELSE IF (found) THEN
      project_name = TRIM(ADJUSTL(name))
    END IF
  END DO arguments
  IF (LEN(TRIM(project_name)) .EQ. 0 .OR. (found .EQV. .false.)) CALL &
&   COMMAND_LINE_SYNTAX()
  RETURN
END SUBROUTINE PARSE_COMMAND_LINE
