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
!* Author:
!*
!*   Robert Renka,
!*   University of North Texas,
!*   (817) 565-2767.
!*
!* Reference:
!*
!*   Robert Renka,
!*   Algorithm 661: QSHEP3D, Quadratic Shepard method for trivariate
!*   interpolation of scattered data,
!*   ACM Transactions on Mathematical Software,
!*   Volume 14, 1988, pages 151-152.
!*
!******************************************************************************

MODULE TOMS661

   IMPLICIT NONE

   ! The number of data points to be used in
   ! the least squares fit for coefficients defining the nodal functions Q(K).  
   ! A recommended value is NQ = 17.  9 <= NQ <= MIN ( 40, N-1 ).
   INTEGER, PARAMETER, PRIVATE :: NQ = 17
   
   ! the number of nodes within (and defining) 
   ! the radii of influence R(K) which enter into the weights W(K).  For N 
   ! sufficiently large, a recommended value is NW = 32.  
   ! 1 <= NW <= min(40,N-1).
   INTEGER, PARAMETER, PRIVATE :: NW = 20
   
   ! The number of rows, columns, and planes in the cell grid (as defined in 
   ! subroutine store3).  A box containing the nodes is partitioned into cells
   ! in order to increase search efficiency. The cell density (average number 
   ! of nodes per cell) is D = N/(NR**3).  A recommended value, based on 
   ! empirical evidence, is D = 3, so NR = (N/3)**(1/3).  (1 <= NR).
   INTEGER, PARAMETER, PRIVATE :: NR = 3
   
   ! A cell array such that LCELL(I,J,K) contains the index (for X, Y, and Z)
   ! of the first node (node with smallest index) in cell (I,J,K), or 
   ! LCELL(I,J,K) = 0 if no nodes are contained in the cell.  The orner of cell
   ! (I,J,K) which is farthest from the box corner defined by XYZMIN has 
   ! coordinates (XMIN+I*DX,YMIN+J*DY,ZMIN+K*DZ), where (XMIN,YMIN,ZMIN) are 
   ! the elements of XYZMIN.
   INTEGER, DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: lcell
   
   ! Next-node indices such that lnext(l) contains the index of the next node
   ! in the cell which contains node l, or lnext(l) = l if l is the last node 
   ! in the cell for l = 1,...,n.  (the nodes contained in a cell are ordered 
   ! by their indices.) if, for example, cell (i,j,k) contains nodes 2, 3, and 
   ! 5 (and no others), then lcell(i,j,k) = 2, lnext(2) = 3, lnext(3) = 5, and 
   ! lnext(5) = 5.  lnext is not defined if ier /= 0.
   INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: lnext
   
   ! The minimum nodal coordinates xmin, ymin, and zmin (in that order) unless 
   ! ier = 1.  the opposite corner of the box defined by the nodes is
   ! (xmin+nr*dx,ymin+nr*dy,zmin+nr*dz).
   DOUBLE PRECISION, DIMENSION(3), PRIVATE :: xyzmin
   
   ! The dimensions of the cells unless ier = 1. xyzdel(1) = (xmax-xmin)/nr, 
   ! xyzdel(2) = (ymax-ymin)/nr, and xyzdel(3) = (zmax-zmin)/nr, where xmin,
   ! xmax, ymin, ymax, zmin, and zmax are the extrema of x, y, and z.
   DOUBLE PRECISION, DIMENSION(3), PRIVATE :: xyzdel
   
   ! Square root of the largest element in RSQ, maximum radius R(K).
   DOUBLE PRECISION, PRIVATE :: rmax
   
   ! The squares of the radii r(k) which enter into the weights W(K).
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, PRIVATE :: rsq
   
   ! The coefficients for quadratic nodal function Q(K) in column K.
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: a

   ! Public access members
   PUBLIC  :: init3, del3, interpolate3, gradient3

   ! Private access memebers
   PRIVATE :: QShepard3d, getnp3, givens, rotate, setup3, store3

   !*** Generic interfaces

   interface init3
      module procedure Init_QShepard3d
   end interface 

   interface del3
      module procedure Del_QShepard3d
   end interface
   
   interface interpolate3
      module procedure QShepard3d_Eval
   end interface
   
   interface gradient3
      module procedure QShepard3d_Grad
   end interface


CONTAINS


SUBROUTINE Init_QShepard3d( n, x, y, z, f, ierr )

!*****************************************************************************80
!
! Init_QShepard3d initializes a smooth trivariate interpolant of scattered 3D 
! data.
!

   ! The number of nodes of data values (10 <= N)
   INTEGER ( kind = 4 ), INTENT(IN)  :: n
   
   ! Coordinates of the scatter points
   REAL ( kind = 8 ), INTENT(IN)     :: x(n)
   REAL ( kind = 8 ), INTENT(IN)     :: y(n)
   REAL ( kind = 8 ), INTENT(IN)     :: z(n)
   
   ! Data values at scatter points
   REAL ( kind = 8 ), INTENT(IN)     :: f(n)
   
  ! Error indicator:
  !  0, if no errors were encountered
  !  1, if N, NQ, NW, or NR is out of range
  !  2, if duplicate nodes were encountered
  !  3, if all scatter points are coplanar
   INTEGER ( kind = 4 ), INTENT(OUT) :: ierr
   
   
   ! Allocate memory for module's private allocatable arrays.
   IF (.not.ALLOCATED(lcell)) ALLOCATE (lcell(nr,nr,nr))
   IF (.not.ALLOCATED(lnext)) ALLOCATE (lnext(n))
   IF (.not.ALLOCATED(rsq))   ALLOCATE (rsq(n)) 
   IF (.not.ALLOCATED(a))     ALLOCATE (a(9,n))
   
   ! Prepare the trivariate interpolation scheme defined by Q.
   ! (See QShepard3d Subroutine).
   CALL QShepard3d (n, x, y, z, f, ierr)
   
   RETURN
END SUBROUTINE Init_QShepard3d


SUBROUTINE Del_QShepard3d()

!*****************************************************************************80
!
! Del_QShepard3d frees the memory allocated to trivariate interpolation of 
! scattered 3D data.
!
   
   IF (ALLOCATED(lcell)) DEALLOCATE (lcell) 
   IF (ALLOCATED(lnext)) DEALLOCATE (lnext)
   IF (ALLOCATED(rsq))   DEALLOCATE (rsq) 
   IF (ALLOCATED(a))     DEALLOCATE (a)
   
   RETURN
END SUBROUTINE Del_QShepard3d


SUBROUTINE QShepard3d (n, x, y, z, f, ier)

!*****************************************************************************80
!
! QShepard3d defines a smooth trivariate interpolant of scattered 3D data.
!
!  Discussion:
!
!    This subroutine computes a set of parameters A and RSQ
!    defining a smooth (once continuously differentiable) trivariate 
!    function Q(X,Y,Z) which interpolates data values
!    F at scattered nodes (X,Y,Z).  The interpolant Q may be
!    evaluated at an arbitrary point by function QShepard3d_Eval, 
!    and its first derivatives are computed by subroutine QShepard3Grad.
!
!    The interpolation scheme is a modified quadratic Shepard method:
!
!      Q = (W(1)*Q(1)+W(2)*Q(2)+..+W(N)*Q(N))/(W(1)+W(2)+..+W(N))
!
!    for trivariate functions W(K) and Q(K).  The nodal functions are 
!    given by
!
!      Q(K)(X,Y,Z) = 
!          A(1,K) * DX**2 
!        + A(2,K) * DX * DY 
!        + A(3,K) * DY**2
!        + A(4,K) * DX * DZ
!        + A(5,K) * DY * DZ 
!        + A(6,K) * DZ**2
!        + A(7,K) * DX 
!        + A(8,K) * DY 
!        + A(9,K) * DZ 
!        + F(K)
!
!    where DX = (X-X(K)), DY = (Y-Y(K)), and DZ = (Z-Z(K)).
!
!    Thus, Q(K) is a quadratic function which interpolates the
!    data value at node K.  Its coefficients A(*,K) are obtained
!    by a weighted least squares fit to the closest NQ data
!    points with weights similar to W(K).  Note that the radius
!    of influence for the least squares fit is fixed for each
!    K, but varies with K.
!
!    The weights are taken to be
!
!      W(K)(X,Y,Z) = ( (R(K)-D(K))+ / R(K)*D(K) )**2
!
!    where (R(K)-D(K))+ = 0 if R(K) <= D(K), and D(K)(X,Y,Z)
!    is the euclidean distance between (X,Y,Z) and node K.  The
!    radius of influence R(K) varies with K and is chosen so
!    that NW nodes are within the radius.  Note that W(K) is
!    not defined at node (X(K),Y(K),Z(K)), but Q(X,Y,Z) has
!    limit F(K) as (X,Y,Z) approaches (X(K),Y(K),Z(K)).
!

  ! The number of nodes of data values (10 <= N)
  INTEGER ( kind = 4 ), INTENT(IN)  :: n
  
  ! Coordinates of the scatter points
  REAL ( kind = 8 ), INTENT(IN)     :: x(n)
  REAL ( kind = 8 ), INTENT(IN)     :: y(n)
  REAL ( kind = 8 ), INTENT(IN)     :: z(n)
  
  ! Data values at scatter points
  REAL ( kind = 8 ), INTENT(IN)     :: f(n)
  
  ! Error indicator:
  !  0, if no errors were encountered
  !  1, if N, NQ, NW, or NR is out of range
  !  2, if duplicate nodes were encountered
  !  3, if all nodes are coplanar
  INTEGER ( kind = 4 ), INTENT(OUT) :: ier

  ! root-mean-square distance between k and the
  ! nodes in the least squares system (unless
  ! additional nodes are introduced for stabil-
  ! ity).  the first 6 columns of the matrix
  ! are scaled by 1/avsq, the last 3 by 1/av
  real ( kind = 8 ) av
  
  ! av*av
  real ( kind = 8 ) avsq
  
  ! transpose of the augmented regression matrix
  real ( kind = 8 ) b(10,10)
  
  ! first component of the plane rotation used to
  ! zero the lower triangle of b**t -- computed
  ! by subroutine givens
  real ( kind = 8 ) c
  
  ! minimum of the magnitudes of the diagonal
  ! elements of the regression matrix after
  ! zeros are introduced below the diagonal
  real ( kind = 8 ) dmin
  
  ! tolerance for detecting an ill-conditioned
  ! system.  the system is accepted when DTOL <= DMIN.
  real ( kind = 8 ), parameter :: dtol = 0.01D+00
  
  ! data value at node k -- f(k)
  real ( kind = 8 ) fk
  
  ! index for a, b, npts, xyzmin, xyzmn, xyzdel, and xyzdl
  integer ( kind = 4 ) i
  
  ! do-loop index for back solve
  integer ( kind = 4 ) ib
  
  ! error flag for the call to store3
  integer ( kind = 4 ) ierr
  
  ! irow-1
  integer ( kind = 4 ) irm1
  
  ! row index for b
  integer ( kind = 4 ) irow
  
  ! index for a and b
  integer ( kind = 4 ) j
  
  ! nodal function index and column index for a
  integer ( kind = 4 ) k
  
  ! maximum number of npts elements (must be con-
  ! sistent with the dimension statement above)
  integer ( kind = 4 ) lmax
  
  ! current length of npts
  integer ( kind = 4 ) lnp
  
  ! number of equations in the least squares fit
  integer ( kind = 4 ) neq
  
  ! local copy of n
  integer ( kind = 4 ) nn
  
  ! local copy of nq
  integer ( kind = 4 ) nnq
  
  ! local copy of nr
  integer ( kind = 4 ) nnr
  
  ! local copy of nw
  integer ( kind = 4 ) nnw
  
  ! npts element
  integer ( kind = 4 ) np
  
  ! array containing the indices of a sequence of
  ! nodes to be used in the least squares fit
  ! or to compute rsq.  the nodes are ordered
  ! by distance from k and the last element
  ! (usually indexed by lnp) is used only to
  ! determine rq, or rsq(k) if NQ < NW.
  integer ( kind = 4 ) npts(40)
  
  ! max(nq,nw)
  integer ( kind = 4 ) nqwmax
  
  ! radius of influence which enters into the
  ! weights for q(k) (see subroutine setup3)
  real ( kind = 8 ) rq
  
  ! squared distance between k and npts(lnp) --
  ! used to compute rq and rsq(k)
  real ( kind = 8 ) rs
  
  ! maximum rsq element encountered
  real ( kind = 8 ) rsmx
  
  ! squared distance between k and npts(lnp-1) --
  ! used to compute a relative change in rs
  ! between succeeding npts elements
  real ( kind = 8 ) rsold
  
  ! tolerance for detecting a sufficiently large
  ! relative change in rs.  if the change is
  ! not greater than rtol, the nodes are
  ! treated as being the same distance from k
  real ( kind = 8 ), parameter :: rtol = 1.0D-05
  
  ! current value of rsq(k)
  real ( kind = 8 ) rws
  
  ! second component of the plane givens rotation
  real ( kind = 8 ) s
  
  ! marquardt stabilization factor used to damp
  ! out the first 6 solution components (second
  ! partials of the quadratic) when the system
  ! is ill-conditioned.  as sf increases, the
  ! fitting function approaches a linear
  real ( kind = 8 ), parameter :: sf = 1.0D+00
  
  ! sum of squared euclidean distances between
  ! node k and the nodes used in the least
  ! squares fit (unless additional nodes are
  ! added for stability)
  real ( kind = 8 ) sum2
  
  ! temporary variable for accumulating a scalar
  ! product in the back solve
  real ( kind = 8 ) t
  
  ! coordinates of node k -- x(k), y(k), z(k)
  real ( kind = 8 ) xk
  real ( kind = 8 ) yk
  real ( kind = 8 ) zk
  
  ! local variables for xyzdel
  real ( kind = 8 ) xyzdl(3)
  
  ! local variables for xyzmin
  real ( kind = 8 ) xyzmn(3)


  nn = n
  nnq = nq
  nnw = nw
  nnr = nr
  nqwmax = max ( nnq, nnw )
  lmax = min ( 40, nn - 1 )

  if ( nnq < 9 ) then
    ier = 1
    return
  end if

  if ( nnw < 1 ) then
    ier = 1
    return
  end if

  if ( lmax < nqwmax ) then
    ier = 1
    return
  end if

  if ( nnr < 1 ) then
    ier = 1
    return
  end if
!
!  Create the cell data structure, and initialize RSMX.
!
  call store3 ( nn, x, y, z, ierr )

  if ( ierr /= 0 ) then
    xyzmin(1:3) = xyzmn(1:3)
    xyzdel(1:3) = xyzdl(1:3)
    ier = 3
    return
  end if

  rsmx = 0.0D+00
!
!  Outer loop on node K.
!
  do k = 1, nn

    xk = x(k)
    yk = y(k)
    zk = z(k)
    fk = f(k)
!
!  Mark node K to exclude it from the search for nearest neighbors.
!
    lnext(k) = -lnext(k)
!
!  Initialize for loop on NPTS.
!
    rs = 0.0D+00
    sum2 = 0.0D+00
    rws = 0.0D+00
    rq = 0.0D+00
    lnp = 0
!
!  Compute NPTS, LNP, rws, neq, rq, and avsq.
!
1   continue

    sum2 = sum2 + rs

    if ( lnp == lmax ) then
      go to 3
    end if

    lnp = lnp + 1
    rsold = rs
    call getnp3 ( xk, yk, zk, x, y, z, np, rs )

    if ( rs == 0.0D+00 ) then
      go to 21
    end if

    npts(lnp) = np

    if ( ( rs - rsold ) / rs < rtol ) then
      go to 1
    end if

    if ( rws == 0.0D+00 .and. nnw < lnp ) then
      rws = rs
    end if

    if ( rq /= 0.0D+00 .or. lnp <= nnq ) then
      go to 2
    end if
!
!  RQ = 0 (not yet computed) and NQ < LNP.  rq =
!  sqrt(rs) is sufficiently large to (strictly) include
!  nq nodes.  The least squares fit will include neq =
!  lnp-1 equations for 9 <= nq <= neq < lmax <= n-1.
!
    neq = lnp - 1
    rq = sqrt ( rs )
    avsq = sum2 / real ( neq, kind = 8 )
!
!  Bottom of loop -- test for termination.
!
2   continue

    if ( nqwmax < lnp ) then
      go to 4
    else
      go to 1
    end if
!
!  All LMAX nodes are included in npts.  rws and/or rq**2 is
!  (arbitrarily) taken to be 10 percent larger than the
!  distance RS to the last node included.
!
3   continue

    if ( rws == 0.0D+00 ) then
      rws = 1.1D+00 * rs
    end if

    if ( rq == 0.0D+00 ) then
      neq = lmax
      rq = sqrt ( 1.1D+00 * rs )
      avsq = sum2 / real ( neq,  kind = 8 )
    end if
!
!  Store RSQ(K), update RSMX if necessary, and compute av.
!
4   continue

    rsq(k) = rws
    rsmx = max ( rsmx, rs )
    av = sqrt ( avsq )
!
!  Set up the augmented regression matrix (transposed) as the
!  columns of B, and zero out the lower triangle (upper
!  triangle of b) with givens rotations -- qr decomposition
!  with orthogonal matrix q not stored.
!
    i = 0

5   continue

    i = i + 1
    np = npts(i)
    irow = min ( i, 10 )
    call setup3 ( xk, yk, zk, fk, x(np), y(np), z(np), f(np), av, avsq, &
      rq, b(1,irow) )

    if ( i == 1 ) then
      go to 5
    end if

    irm1 = irow-1

    do j = 1, irow-1
      call givens ( b(j,j), b(j,irow), c, s )
      call rotate ( 10-j, c, s, b(j+1,j), b(j+1,irow) )
    end do

    if ( i < neq ) then
      go to 5
    end if
!
!  Test the system for ill-conditioning.
!
    dmin = min ( abs(b(1,1)), abs(b(2,2)), abs(b(3,3)), &
                 abs(b(4,4)), abs(b(5,5)), abs(b(6,6)), &
                 abs(b(7,7)), abs(b(8,8)), abs(b(9,9)) )

    if ( dtol <= dmin * rq ) then
      go to 13
    end if

    if ( neq == lmax ) then
      go to 10
    end if
!
!  Increase RQ and add another equation to the system to
!  improve the conditioning.  The number of NPTS elements
!  is also increased if necessary.
!
7   continue

    rsold = rs
    neq = neq + 1
    if ( neq == lmax ) go to 9
    if ( neq == lnp )  go to 8
!
!  NEQ < LNP
!
    np = npts(neq+1)
    rs = (x(np)-xk)**2 + (y(np)-yk)**2 + (z(np)-zk)**2
    if ( ( rs - rsold ) / rs < rtol ) go to 7
    rq = sqrt ( rs )
    go to 5
!
!  Add an element to NPTS.
!
8   continue

    lnp = lnp + 1
    call getnp3 ( xk, yk, zk, x, y, z, np, rs )

    if ( np == 0 ) then
      go to 21
    end if

    npts(lnp) = np

    if ( ( rs - rsold ) / rs < rtol ) then
      go to 7
    end if

    rq = sqrt ( rs )
    go to 5

9   continue

    rq = sqrt ( 1.1D+00 * rs )
    go to 5
!
!  Stabilize the system by damping second partials.  Add
!  multiples of the first six unit vectors to the first
!  six equations.
!
10  continue

    do i = 1, 6

      b(i,10) = sf
      b(i+1:10,10) = 0.0D+00

      do j = i, 9
        call givens ( b(j,j), b(j,10), c, s )
        call rotate ( 10-j, c, s, b(j+1,j), b(j+1,10) )
      end do

    end do
!
!  Test the stabilized system for ill-conditioning.
!
    dmin = min ( abs(b(1,1)),abs(b(2,2)),abs(b(3,3)), &
                 abs(b(4,4)),abs(b(5,5)),abs(b(6,6)), &
                 abs(b(7,7)),abs(b(8,8)),abs(b(9,9)) )
!
!  No unique solution due to collinear nodes.
!
    if ( dmin * rq < dtol ) then
      xyzmin(1:3) = xyzmn(1:3)
      xyzdel(1:3) = xyzdl(1:3)
      ier = 3
      return
    end if
!
!  Solve the 9 by 9 triangular system for the coefficients
!
13  continue

    do ib = 1, 9
      i = 10-ib
      t = 0.0D+00
      do j = i+1, 9
        t = t + b(j,i)*a(j,k)
      end do
      a(i,k) = (b(10,i)-t)/b(i,i)
    end do
!
!  Scale the coefficients to adjust for the column scaling.
!
    a(1:6,k) = a(1:6,k) / avsq
    a(7,k) = a(7,k) / av
    a(8,k) = a(8,k) / av
    a(9,k) = a(9,k) / av
!
!  Unmark K and the elements of NPTS.
!
    lnext(k) = -lnext(k)
    do i = 1, lnp
      np = npts(i)
      lnext(np) = -lnext(np)
    end do

  end do
!
!  No errors encountered.
!
!as  xyzmin(1:3) = xyzmn(1:3)
!as  xyzdel(1:3) = xyzdl(1:3)
  rmax = sqrt ( rsmx )
  ier = 0
  return
!
!  N, NQ, NW, or NR is out of range.
!
20 continue

  ier = 1
  return
!
!  Duplicate nodes were encountered by GETNP3.
!
21 ier = 2
  RETURN
END SUBROUTINE QShepard3d


FUNCTION QShepard3d_Eval ( px, py, pz, n, x, y, z, f )

!*****************************************************************************80
!
! QShepard3d_Eval evaluates the interpolant function Q(X,Y,Z) created by 
! QShepard3d.
!
!  Discussion:
!
!    This function returns the value Q(PX,PY,PZ) where Q is
!    the weighted sum of quadratic nodal functions defined in
!    subroutine QShepard3d.  QShepard3Grad may be called to compute
!    a gradient of Q along with the value, or to test for errors.
!
!    This function should not be called if a nonzero error flag was
!    returned by QShepard3d.
!

  ! The point P at which Q is to be evaluated.
  REAL ( kind = 8 ), INTENT(IN)    :: px
  REAL ( kind = 8 ), INTENT(IN)    :: py
  REAL ( kind = 8 ), INTENT(IN)    :: pz
  
  ! The number of nodes and data values defining Q. (10 <= N).
  INTEGER ( kind = 4 ), INTENT(IN) :: n
  
  ! The node coordinates and data values interpolated by Q.
  REAL ( kind = 8 ), INTENT(IN)    :: x(n)
  REAL ( kind = 8 ), INTENT(IN)    :: y(n)
  REAL ( kind = 8 ), INTENT(IN)    :: z(n)
  REAL ( kind = 8 ), INTENT(IN)    :: f(n)

  ! The function value Q(PX,PY,PZ) unless N, NR, XYZDEL, or RMAX is invalid,
  ! in which case the value 0 is returned.
  real ( kind = 8 ) QShepard3d_Eval
  
  real ( kind = 8 ) delx
  real ( kind = 8 ) dely
  real ( kind = 8 ) delz
  real ( kind = 8 ) dxsq
  real ( kind = 8 ) dysq
  real ( kind = 8 ) dzsq
  real ( kind = 8 ) ds
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) dz
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmin
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) kmin
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lp
  real ( kind = 8 ) rd
  real ( kind = 8 ) rds
  real ( kind = 8 ) rs
  real ( kind = 8 ) sw
  real ( kind = 8 ) swq
  real ( kind = 8 ) w
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xp
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yp
  real ( kind = 8 ) zmax
  real ( kind = 8 ) zmin
  real ( kind = 8 ) zp


  xp = px
  yp = py
  zp = pz
  xmin = xyzmin(1)
  ymin = xyzmin(2)
  zmin = xyzmin(3)
  dx = xyzdel(1)
  dy = xyzdel(2)
  dz = xyzdel(3)

  if ( n < 10 ) then
    QShepard3d_Eval = 0.0D+00
    return
  end if

  if ( nr < 1  .or.  dx <= 0.0 &
         .or.  dy <= 0.0  .or.  dz <= 0.0  .or. &
         rmax < 0.0 ) then
    QShepard3d_Eval = 0.0D+00
    return
  end if
!
!  Set IMIN, imax, jmin, jmax, kmin, and kmax to cell indices
!  defining the range of the search for nodes whose radii
!  include P.  The cells which must be searched are those
!  intersected by (or contained in) a sphere of radius rmax
!  centered at P.
!
  imin = int((xp-xmin-rmax)/dx) + 1
  imin = max ( imin, 1 )

  imax = int((xp-xmin+rmax)/dx) + 1
  imax = min ( imax, nr )

  jmin = int((yp-ymin-rmax)/dy) + 1
  jmin = max ( jmin, 1 )

  jmax = int((yp-ymin+rmax)/dy) + 1
  jmax = min ( jmax, nr )

  kmin = int((zp-zmin-rmax)/dz) + 1
  kmin = max ( kmin, 1 )

  kmax = int((zp-zmin+rmax)/dz) + 1
  kmax = min ( kmax, nr )
!
!  Test for no cells within the sphere of radius RMAX.
!
  if ( imax < imin .or. &
       jmax < jmin .or. &
       kmax < kmin ) then
    QShepard3d_Eval = 0.0D+00
    return
  end if
!
!  Accumulate weight values in SW and weighted nodal function
!  values in SWQ.  The weights are w(l) = ((r-d)+/(r*d))**2
!  for r**2 = rsq(l) and d = distance between P and node L.
!
  sw = 0.0D+00
  swq = 0.0D+00
!
!  Outer loop on cells (i,j,k).
!
  do k = kmin, kmax
    do j = jmin, jmax
      do i = imin, imax

        l = lcell(i,j,k)

        if ( l == 0 ) then
          cycle
        end if
!
!  Inner loop on nodes L.
!
        do

          delx = xp - x(l)
          dely = yp - y(l)
          delz = zp - z(l)

          dxsq = delx * delx
          dysq = dely * dely
          dzsq = delz * delz

          ds = dxsq + dysq + dzsq
          rs = rsq(l)

          if ( ds < rs ) then

            if ( ds == 0.0D+00 ) then
              QShepard3d_Eval = f(l)
              return
            end if

            rds = rs * ds
            rd = sqrt ( rds )
            w = ( rs + ds - rd - rd ) / rds
            sw = sw + w

            swq = swq + w *( a(1,l) * dxsq + a(2,l)*delx*dely + &
              a(3,l) * dysq + a(4,l)*delx*delz + &
              a(5,l) * dely*delz + a(6,l)*dzsq + &
              a(7,l) * delx + a(8,l)*dely + &
              a(9,l) * delz + f(l) )

          end if

          lp = l
          l = lnext(lp)

          if ( l == lp ) then
            exit
          end if

        end do

      end do
    end do
  end do
!
!  SW = 0 iff P is not within the radius R(L) for any node L.
!
  if ( sw == 0.0D+00 ) then
    QShepard3d_Eval = 0.0D+00
  else
    QShepard3d_Eval = swq / sw
  end if

  RETURN
END FUNCTION QShepard3d_Eval


SUBROUTINE QShepard3d_Grad ( px, py, pz, n, x, y, z, f, q, qx, qy, qz, ier )

!*****************************************************************************80
!
!  QShepard3d_Grad computes the value and gradient of the interpolant function.
!
!  Discussion:
!
!    This subroutine computes the value and gradient at (PX,PY,PZ) of 
!    the interpolatory function Q defined in subroutine QShepard3d.  
!
!    Q(X,Y,Z) is a weighted sum of quadratic nodal functions.
!

  ! The point P at which Q and its partials are to be evaluated.
  REAL ( kind = 8 ), INTENT(IN)     :: px
  REAL ( kind = 8 ), INTENT(IN)     :: py
  REAL ( kind = 8 ), INTENT(IN)     :: pz
  
  ! The number of nodes and data values defining Q. (10 <= N)
  INTEGER ( kind = 4 ), INTENT(IN)  :: n
  
  ! The node coordinates and data values interpolated by Q.
  REAL ( kind = 8 ), INTENT(IN)     :: x(n)
  REAL ( kind = 8 ), INTENT(IN)     :: y(n)
  REAL ( kind = 8 ), INTENT(IN)     :: z(n)
  REAL ( kind = 8 ), INTENT(IN)     :: f(n)
  
  ! The value of Q at (PX,PY,PZ) unless IER == 1, in which case no values 
  ! are returned.
  REAL ( kind = 8 ), INTENT(OUT)    :: q
  
  ! The first partial derivatives of Q at (PX,PY,PZ) unless IER == 1.
  REAL ( kind = 8 ), INTENT(OUT)    :: qx
  REAL ( kind = 8 ), INTENT(OUT)    :: qy
  REAL ( kind = 8 ), INTENT(OUT)    :: qz
  
  ! IER, error indicator
  !  0, if no errors were encountered.
  !  1, if N, NR, XYZDEL, or RMAX is invalid.
  !  2, if no errors were encountered but (PX.PY.PZ) is not within the 
  !     radius R(K) for any node K (and thus Q = QX = QY = QZ = 0).
  INTEGER ( kind = 4 ), INTENT(OUT) :: ier

  real ( kind = 8 ) delx
  real ( kind = 8 ) dely
  real ( kind = 8 ) delz
  real ( kind = 8 ) ds
  real ( kind = 8 ) dx
  real ( kind = 8 ) dxsq
  real ( kind = 8 ) dy
  real ( kind = 8 ) dysq
  real ( kind = 8 ) dz
  real ( kind = 8 ) dzsq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmin
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) kmin
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lp
  real ( kind = 8 ) ql
  real ( kind = 8 ) qlx
  real ( kind = 8 ) qly
  real ( kind = 8 ) qlz
  real ( kind = 8 ) rd
  real ( kind = 8 ) rds
  real ( kind = 8 ) rs
  real ( kind = 8 ) sw
  real ( kind = 8 ) swq
  real ( kind = 8 ) swqx
  real ( kind = 8 ) swqy
  real ( kind = 8 ) swqz
  real ( kind = 8 ) sws
  real ( kind = 8 ) swx
  real ( kind = 8 ) swy
  real ( kind = 8 ) swz
  real ( kind = 8 ) t
  real ( kind = 8 ) w
  real ( kind = 8 ) wx
  real ( kind = 8 ) wy
  real ( kind = 8 ) wz
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xp
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yp
  real ( kind = 8 ) zmax
  real ( kind = 8 ) zmin
  real ( kind = 8 ) zp


  xp = px
  yp = py
  zp = pz
  xmin = xyzmin(1)
  ymin = xyzmin(2)
  zmin = xyzmin(3)
  dx = xyzdel(1)
  dy = xyzdel(2)
  dz = xyzdel(3)

  if ( n < 10  .or.  nr < 1 &
       .or.  dx <= 0.       &
       .or.  dy <= 0.0D+00  &
       .or.  dz <= 0.0D+00  &
       .or.  rmax < 0.0D+00 ) then
    ier = 1
    return
  end if
!
!  Set IMIN, IMAX, jmin, jmax, kmin, and kmax to cell indices
!  defining the range of the search for nodes whose radii
!  include P.  The cells which must be searched are those
!  intersected by or contained in a sphere of radius RMAX
!  centered at P.
!
  imin = int((xp-xmin-rmax)/dx) + 1
  imin = max ( imin, 1 )

  imax = int((xp-xmin+rmax)/dx) + 1
  imax = min ( imax, nr )

  jmin = int((yp-ymin-rmax)/dy) + 1
  jmin = max ( jmin, 1 )

  jmax = int((yp-ymin+rmax)/dy) + 1
  jmax = min ( jmax, nr )

  kmin = int((zp-zmin-rmax)/dz) + 1
  kmin = max ( kmin, 1 )

  kmax = int((zp-zmin+rmax)/dz) + 1
  kmax = min ( kmax, nr )
!
!  Test for no cells within the sphere of radius RMAX.
!
  if ( imax < imin .or. &
       jmax < jmin .or. &
       kmax < kmin ) then
    q = 0.0D+00
    qx = 0.0D+00
    qy = 0.0D+00
    qz = 0.0D+00
    ier = 2
    return
  end if
!
!  Q = swq/sw = sum(w(l)*q(l))/sum(w(l)) where the sum is
!  from l = 1 to N, q(l) is the quadratic nodal function,
!  and w(l) = ((r-d)+/(r*d))**2 for radius r(l) and dist-
!  ance d(l).  thus
!
!    qx = ( swqx * sw - swq * swx ) / sw**2
!    qy = ( swqy * sw - swq * swy ) / sw**2
!    qz = ( swqz * sw - swq * swz ) / sw**2
!
!  where swqx and swx are partial derivatives with respect
!  to X of swq and sw, respectively.  swqy, swy, swqz, and
!  swz are defined similarly.
!
  sw   = 0.0D+00
  swx  = 0.0D+00
  swy  = 0.0D+00
  swz  = 0.0D+00
  swq  = 0.0D+00
  swqx = 0.0D+00
  swqy = 0.0D+00
  swqz = 0.0D+00
!
!  Outer loop on cells (i,j,k).
!
  do k = kmin, kmax
    do j = jmin, jmax
      do i = imin, imax

        l = lcell(i,j,k)
        if ( l == 0 ) then
          cycle
        end if
!
!  Inner loop on nodes L.
!
        do

          delx = xp - x(l)
          dely = yp - y(l)
          delz = zp - z(l)
          dxsq = delx * delx
          dysq = dely * dely
          dzsq = delz * delz
          ds = dxsq + dysq + dzsq
          rs = rsq(l)

          if ( ds < rs ) then

            if ( ds == 0.0D+00 ) then
              q = f(l)
              qx = a(7,l)
              qy = a(8,l)
              qz = a(9,l)
              ier = 0
              return
            end if

            rds = rs * ds
            rd = sqrt ( rds )
            w = ( rs + ds - rd - rd ) / rds
            t = 2.0D+00 * ( rd - rs ) / ( ds * rds )
            wx = delx * t
            wy = dely * t
            wz = delz * t

            qlx = 2.0D+00 *a(1,l)*delx + a(2,l)*dely + a(4,l)*delz
            qly = a(2,l)*delx + 2.0D+00 * a(3,l)*dely + a(5,l)*delz
            qlz = a(4,l)*delx + a(5,l)*dely + 2.0D+00 * a(6,l)*delz
            ql = (qlx*delx + qly*dely + qlz*delz) / 2.0D+00 + &
              a(7,l)*delx + a(8,l)*dely + a(9,l)*delz + f(l)
            qlx = qlx + a(7,l)
            qly = qly + a(8,l)
            qlz = qlz + a(9,l)

            sw = sw + w
            swx = swx + wx
            swy = swy + wy
            swz = swz + wz
            swq = swq + w*ql
            swqx = swqx + wx*ql + w*qlx
            swqy = swqy + wy*ql + w*qly
            swqz = swqz + wz*ql + w*qlz

          end if

          lp = l
          l = lnext(lp)

          if ( l == lp ) then
            exit
          end if

        end do

      end do
    end do
  end do
!
!  SW = 0 iff P is not within the radius R(L) for any node L.
!
  if ( sw /= 0.0D+00 ) then

    q   = swq / sw
    sws = sw * sw
    qx  = ( swqx * sw - swq * swx ) / sws
    qy  = ( swqy * sw - swq * swy ) / sws
    qz  = ( swqz * sw - swq * swz ) / sws
    ier = 0
!
!  No cells contain a point within RMAX of P, or
!  SW = 0 and thus RSQ(L) <= DS for all L.
!
  else

    q   = 0.0D+00
    qx  = 0.0D+00
    qy  = 0.0D+00
    qz  = 0.0D+00
    ier = 2

  end if

  RETURN
END SUBROUTINE QShepard3d_Grad


SUBROUTINE getnp3 ( px, py, pz, x, y, z, np, dsq )

!*****************************************************************************80
!
!  GETNP3 finds the closest node to a given point.
!
!  Discussion:
!
!    Given a set of N nodes and the data structure defined in
!    subroutine STORE3, this subroutine uses the cell method to
!    find the closest unmarked node NP to a specified point P.
!
!    NP is then marked by setting LNEXT(NP) to -LNEXT(NP).  (A
!    node is marked if and only if the corresponding lnext element 
!    is negative.  The absolute values of LNEXT elements,
!    however, must be preserved.)  Thus, the closest M nodes to
!    P may be determined by a sequence of M calls to this routine.  
!    Note that if the nearest neighbor to node K is to
!    be determined (PX = X(K), PY = Y(K), and PZ = Z(K)), then
!    K should be marked before the call to this routine.
!
!    The search is begun in the cell containing (or closest
!    to) P and proceeds outward in box-shaped layers until all
!    cells which contain points within distance R of P have
!    been searched, where R is the distance from P to the first
!    unmarked node encountered (infinite if no unmarked nodes
!    are present).
!

  ! Coordinates of the point P whose nearest unmarked neighbor is to be found.
  REAL ( kind = 8 ), INTENT(IN)     :: px
  REAL ( kind = 8 ), INTENT(IN)     :: py
  REAL ( kind = 8 ), INTENT(IN)     :: pz
  
  ! Coordinates of the nodes.
  REAL ( kind = 8 ), INTENT(IN)     :: x(*)
  REAL ( kind = 8 ), INTENT(IN)     :: y(*)
  REAL ( kind = 8 ), INTENT(IN)     :: z(*)
  
  ! Index of the nearest unmarked node to P, or 0 if all nodes are marked or 
  ! NR < 1 or an element of XYZDEL is not positive.  LNEXT(NP) < 0 if NP /= 0.
  INTEGER ( kind = 4 ), INTENT(OUT) :: np
  
  ! Squared euclidean distance between P and node NP, or 0 if NP = 0.
  REAL ( kind = 8 ), INTENT(OUT)    :: dsq

  real ( kind = 8 ) delx
  real ( kind = 8 ) dely
  real ( kind = 8 ) delz
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) dz
  logical first
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmin
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k0
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) kmin
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lmin
  integer ( kind = 4 ) ln
  real ( kind = 8 ) r
  real ( kind = 8 ) rsmin
  real ( kind = 8 ) rsq
  real ( kind = 8 ) xp
  real ( kind = 8 ) yp
  real ( kind = 8 ) zp


  xp = px
  yp = py
  zp = pz
  dx = xyzdel(1)
  dy = xyzdel(2)
  dz = xyzdel(3)
!
!  Test for invalid input parameters.
!
  if ( nr < 1 .or. dx <= 0.0D+00 .or. dy <= 0.0D+00 .or. dz <= 0.0D+00 ) then
    np = 0
    dsq = 0.0D+00
    return
  end if
!
!  Initialize parameters --
!
!  first = true iff the first unmarked node has yet to be encountered,
!  imin,...,kmax = cell indices defining the range of the search,
!  delx,dely,delz = px-xyzmin(1), py-xyzmin(2), and pz-xyzmin(3),
!  i0,j0,k0 = cell containing or closest to p,
!  i1,...,k2 = cell indices of the layer whose intersection
!  with the range defined by imin,...,kmax is
!  currently being searched.
!
  first = .true.
  imin  = 1
  imax  = nr
  jmin  = 1
  jmax  = nr
  kmin  = 1
  kmax  = nr
  delx  = xp - xyzmin(1)
  dely  = yp - xyzmin(2)
  delz  = zp - xyzmin(3)
  i0    = int(delx/dx) + 1

  if ( i0 < 1 ) then
    i0 = 1
  end if

  if ( nr < i0 ) then
    i0 = nr
  end if

  j0 = int ( dely / dy ) + 1
  if ( j0 < 1 ) j0 = 1
  if ( nr < j0 ) then
    j0 = nr
  end if

  k0 = int(delz/dz) + 1
  if ( k0 < 1 ) k0 = 1
  if ( nr < k0 ) then
    k0 = nr
  end if

  i1 = i0
  i2 = i0
  j1 = j0
  j2 = j0
  k1 = k0
  k2 = k0
!
!  Outer loop on layers, inner loop on layer cells, excluding
!  those outside the range (imin,imax) x (jmin,jmax) x (kmin,kmax).
!
1 continue

  do k = k1, k2

    if ( kmax < k ) then
      exit
    end if

    if ( k < kmin ) then
      cycle
    end if

    do j = j1, j2

      if ( jmax < j ) then
        exit
      end if

      if ( j < jmin ) then
        cycle
      end if

      do i = i1, i2

        if ( imax < i ) then
          exit
        end if

        if ( i < imin ) then
          cycle
        end if

        if ( k /= k1  .and.  k /= k2  .and. j /= j1  .and.  &
          j /= j2 .and. i /= i1 .and. i /= i2 ) then
          cycle
        end if
!
!  Search cell (i,j,k) for unmarked nodes L.
!
        l = lcell(i,j,k)

        if ( l == 0 ) then
          cycle
        end if
!
!  Loop on nodes in cell (i,j,k).
!
2       continue

        ln = lnext(l)

        if ( ln < 0 ) then
          go to 4
        end if
!
!  Node L is not marked.
!
        rsq = (x(l)-xp)**2 + (y(l)-yp)**2 + (z(l)-zp)**2

        if ( .not. first ) then
          go to 3
        end if
!
!  Node L is the first unmarked neighbor of P encountered.
!  initialize LMIN to the current candidate for NP, and
!  rsmin to the squared distance from P to lmin.  imin,
!  imax, jmin, jmax, kmin, and kmax are updated to define
!  the smallest rectangle containing a sphere of radius
!  r = sqrt(rsmin) centered at P, and contained in (1,nr)
!  x (1,nr) x (1,nr) (except that, if P is outside the
!  box defined by the nodes, it is possible that NR < imin
!  or imax < 1, etc.).  FIRST is reset to
!  false.
!
        lmin = l
        rsmin = rsq
        r = sqrt(rsmin)
        imin = int((delx-r)/dx) + 1
        if ( imin < 1 ) imin = 1
        imax = int((delx+r)/dx) + 1
        if ( nr < imax ) imax = nr
        jmin = int((dely-r)/dy) + 1
        if ( jmin < 1 ) jmin = 1
        jmax = int((dely+r)/dy) + 1
        if ( nr < jmax ) jmax = nr
        kmin = int((delz-r)/dz) + 1
        if ( kmin < 1 ) kmin = 1
        kmax = int((delz+r)/dz) + 1
        if ( nr < kmax ) kmax = nr
        first = .false.
        go to 4
!
!  Test for node L closer than LMIN to P.
!
3       continue
!
!  Update LMIN and RSMIN.
!
        if ( rsq < rsmin ) then
          lmin = l
          rsmin = rsq
        end if
!
!  Test for termination of loop on nodes in cell (i,j,k).
!
4       continue
 
        if ( abs(ln) == l ) then
          cycle
        end if

        l = abs ( ln )
        go to 2

      end do

    end do

  end do
!
! Test for termination of loop on cell layers.
!
  if ( i1 <= imin .and. imax <= i2 .and. &
       j1 <= jmin .and. jmax <= j2 .and. &
       k1 <= kmin .and. kmax <= k2 ) go to 9
  i1 = i1 - 1
  i2 = i2 + 1
  j1 = j1 - 1
  j2 = j2 + 1
  k1 = k1 - 1
  k2 = k2 + 1
  go to 1
!
!  Unless no unmarked nodes were encountered, LMIN is the
!  closest unmarked node to P.
!
9 continue

  if ( .not. first ) then
    np = lmin
    dsq = rsmin
    lnext(lmin) = -lnext(lmin)
  else
    np = 0
    dsq = 0.0D+00
  end if

  RETURN
END SUBROUTINE getnp3


SUBROUTINE givens ( a, b, c, s )

!*****************************************************************************80
!
!  GIVENS constructs a Givens plane rotation.
!
!  Discussion:
!
!    This routine constructs the Givens plane rotation
!
!          ( c  s)
!      g = (     ) 
!          (-s  c)
!
!    where c*c + s*s = 1, which zeros the second entry of the 2-vector 
!    (a b)-transpose.  A call to GIVENS is normally followed by a call 
!    to ROTATE which applies the transformation to a 2 by N matrix.  
!    This routine was taken from LINPACK.
!

  ! On input, the first component of the vector.  
  ! On output, overwritten by r = +/-sqrt(a*a + b*b)
  real ( kind = 8 ), INTENT(INOUT) :: a
  
  ! On input, the second component of the vector.  
  ! On output, overwritten by a value z which allows c and s to be 
  ! recovered as follows --
  !    c = sqrt(1-z*z), s=z     if abs(z) <= 1.
  !    c = 1/z, s = sqrt(1-c*c) if 1 < abs(z).
  real ( kind = 8 ), INTENT(INOUT) :: b
  
  ! c = +/-(a/r)
  real ( kind = 8 ), INTENT(OUT)   :: c
  
  ! s = +/-(b/r)
  real ( kind = 8 ), INTENT(OUT)   :: s
  
  ! local copies of a and b
  real ( kind = 8 ) aa
  real ( kind = 8 ) bb
  
  ! r = c*a + s*b = +/-sqrt(a*a+b*b)
  real ( kind = 8 ) r
  
  ! variables used to scale a and b for computing r
  real ( kind = 8 ) u
  real ( kind = 8 ) v


  aa = a
  bb = b
!
!  Abs(b) < abs(a)
!  Note that R has the sign of A, 0 < C, and S has sign(a)*sign(b).
!
  if ( abs ( bb ) < abs ( aa ) ) then
    u = aa + aa
    v = bb / u
    r = sqrt ( 0.25D+00 + v * v ) * u
    c = aa / r
    s = v * ( c + c )
    b = s
    a = r
    return
  end if
!
!  abs(a) <= abs(b)
!
  if ( bb == 0.0D+00 ) then
    c = 1.0D+00
    s = 0.0D+00
    return
  end if

  u = bb + bb
  v = aa / u
!
!  Store R in A.
!
  a = sqrt ( 0.25D+00 + v * v ) * u
  s = bb / a
  c = v * ( s + s )
!
!  Note that R has the sign of b, 0 < S, and c has sign(a)*sign(b).
!
  b = 1.0D+00
  if ( c /= 0.0D+00 ) then
    b = 1.0D+00 / c
  end if

  return
END SUBROUTINE givens


SUBROUTINE rotate ( n, c, s, x, y )

!*****************************************************************************80
!
!  ROTATE applies a Givens rotation to two vectors.
!
!  Discussion:
!
!    This routine applies the Givens rotation
!
!      ( c  s)
!      (-s  c)
! 
!    to the 2 by n matrix
!
!      (x(1) ... x(n))
!      (y(1) ... y(n))
!

  ! The number of columns to be rotated.
  INTEGER ( kind = 4 ), INTENT(IN) :: n
  
  ! The elements of the Givens rotation.
  ! These may be determined by subroutine GIVENS.
  REAL ( kind = 8 ), INTENT(IN)    :: c
  REAL ( kind = 8 ), INTENT(IN)    :: s
  
  ! The vectors to be rotated.
  REAL ( kind = 8 ), INTENT(INOUT) :: x(n)
  REAL ( kind = 8 ), INTENT(INOUT) :: y(n)

  integer ( kind = 4 ) i
  real ( kind = 8 ) xi
  real ( kind = 8 ) yi


  if ( c == 1.0D+00 .and. s == 0.0D+00 ) then
    return
  end if

  do i = 1, n
    xi = x(i)
    yi = y(i)
    x(i) =  c * xi + s * yi
    y(i) = -s * xi + c * yi
  end do

  return
END SUBROUTINE rotate


SUBROUTINE setup3 ( xk, yk, zk, fk, xi, yi, zi, fi, s1, s2, r, row )

!*****************************************************************************80
!
!  SETUP3 sets up the weighted least-squares fit of the data.
!
!  Discussion:
!
!    This routine sets up the I-th row of an augmented regression matrix 
!    for a weighted least-squares fit of a quadratic function Q(X,Y,Z) 
!    to a set of data values F, where Q(XK,YK,ZK) = FK.  
!
!    The first 6 columns (quadratic terms) are scaled by 1/S2, and columns 
!    7, 8, and 9 (linear terms) are scaled by 1/S1.  The weight is 
!    (R-D)/(R*D) if D < R, and 0 if R <= D, where D is the distance 
!    between nodes I and K.
!

  ! Coordinates and data value at node K (interpolated by q).
  REAL ( kind = 8 ), INTENT(IN)  :: xk
  REAL ( kind = 8 ), INTENT(IN)  :: yk
  REAL ( kind = 8 ), INTENT(IN)  :: zk
  REAL ( kind = 8 ), INTENT(IN)  :: fk
  
  ! coordinates and data value at node I.
  REAL ( kind = 8 ), INTENT(IN)  :: xi
  REAL ( kind = 8 ), INTENT(IN)  :: yi
  REAL ( kind = 8 ), INTENT(IN)  :: zi
  REAL ( kind = 8 ), INTENT(IN)  :: fi
  
  ! reciprocals of the scale factors.
  REAL ( kind = 8 ), INTENT(IN)  :: s1
  REAL ( kind = 8 ), INTENT(IN)  :: s2
  
  ! Radius of influence about node K defining the weight.
  REAL ( kind = 8 ), INTENT(IN)  :: r
  
  ! A row of the augmented regression matrix.
  REAL ( kind = 8 ), INTENT(OUT) :: row(10)

  real ( kind = 8 ) d
  real ( kind = 8 ) dx
  real ( kind = 8 ) dxsq
  real ( kind = 8 ) dy
  real ( kind = 8 ) dysq
  real ( kind = 8 ) dz
  real ( kind = 8 ) dzsq
  integer ( kind = 4 ) i
  real ( kind = 8 ) w
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2


  dx = xi - xk
  dy = yi - yk
  dz = zi - zk

  d = sqrt ( dx**2 + dy**2 + dz**2 )

  if ( d <= 0.0D+00 .or. r <= d ) then
    row(1:10) = 0.0D+00
    return
  end if

  w = ( r - d ) / r / d
  w1 = w / s1
  w2 = w / s2

  row(1) = dx * dx * w2
  row(2) = dx * dy * w2
  row(3) = dy * dy * w2
  row(4) = dx * dz * w2
  row(5) = dy * dz * w2
  row(6) = dz * dz * w2

  row(7) = dx * w1
  row(8) = dy * w1
  row(9) = dz * w1

  row(10) = ( fi - fk ) * w

  RETURN
END SUBROUTINE setup3


SUBROUTINE store3 ( n, x, y, z, ier )

!*****************************************************************************80
!
!  STORE3 sets up a data structure for N scattered nodes in 3D.
!
!  Discussion:
!
!    Given a set of N arbitrarily distributed nodes in three-space, 
!    this subroutine creates a data structure for a cell-based method of 
!    solving closest-point problems.  The smallest box containing the nodes 
!    is partitioned into an NR by NR by NR uniform grid of cells, and 
!    nodes are associated with cells.  In particular, the data structure
!    stores the indices of the nodes contained in each cell.  For a 
!    uniform random distribution of nodes, the nearest node to an 
!    arbitrary point can be determined in constant expected time.
!

  ! The number of nodes. (2 <= N)
  INTEGER ( kind = 4 ), INTENT(IN)  :: n
  
  ! Coordinates of the nodes.
  REAL ( kind = 8 ), INTENT(IN)     :: x(n)
  REAL ( kind = 8 ), INTENT(IN)     :: y(n)
  REAL ( kind = 8 ), INTENT(IN)     :: z(n)
  
  ! Error indicator:
  !   0, if no errors were encountered.
  !   1, if n < 2 or nr < 1.
  !   2, if a component of xyzdel is not positive.
  INTEGER ( kind = 4 ), INTENT(OUT) :: ier

  real ( kind = 8 ) delx
  real ( kind = 8 ) dely
  real ( kind = 8 ) delz
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nnr
  integer ( kind = 4 ) np1
  real ( kind = 8 ) xmn
  real ( kind = 8 ) xmx
  real ( kind = 8 ) ymx
  real ( kind = 8 ) ymn
  real ( kind = 8 ) zmx
  real ( kind = 8 ) zmn


  ier = 0
  nn = n
  nnr = nr

  if ( nn < 2 .or. nnr < 1 ) then
    ier = 1
    return
  end if
!
!  Compute the dimensions of the box containing the nodes.
!
  xmn = minval ( x(1:nn) )
  xmx = maxval ( x(1:nn) )
  ymn = minval ( y(1:nn) )
  ymx = maxval ( y(1:nn) )
  zmn = minval ( z(1:nn) )
  zmx = maxval ( z(1:nn) )

  xyzmin(1) = xmn
  xyzmin(2) = ymn
  xyzmin(3) = zmn
!
!  Compute cell dimensions and test for zero area.
!
  delx = ( xmx - xmn ) / real ( nnr, kind = 8 )
  dely = ( ymx - ymn ) / real ( nnr, kind = 8 )
  delz = ( zmx - zmn ) / real ( nnr, kind = 8 )

  xyzdel(1) = delx
  xyzdel(2) = dely
  xyzdel(3) = delz

  if ( delx == 0.0D+00 .or. dely == 0.0D+00 .or. delz == 0.0D+00 ) then
    ier = 2
    return
  end if
!
!  Initialize LCELL.
!
  lcell(1:nnr,1:nnr,1:nnr) = 0
!
!  Loop on nodes, storing indices in LCELL and LNEXT.
!
  np1 = nn + 1

  do ll = 1, nn

    lb = np1 - ll

    i = int ( ( x(lb) - xmn ) / delx ) + 1
    if ( nnr < i ) then
      i = nnr
    end if

    j = int ( ( y(lb) - ymn ) / dely ) + 1
    if ( nnr < j ) then
      j = nnr
    end if

    k = int ( ( z(lb) - zmn ) / delz ) + 1
    if ( nnr < k ) then
      k = nnr
    end if

    l = lcell(i,j,k)

    if ( l == 0 ) then
      lnext(lb) = lb
    else
      lnext(lb) = l
    end if

    lcell(i,j,k) = lb

  end do

  RETURN
END SUBROUTINE store3


END MODULE TOMS661