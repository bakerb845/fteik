      PURE SUBROUTINE fteik_dummy(i, j, k, nz, nx, ny, nzx,   &
                                  ttimes, tt1, tt2, tt3, tt4, &
                                  tt5, tt6, tt7, tt8)         &
                 BIND(C, NAME='fteik_dummy')
      !$OMP DECLARE SIMD(fteik_dummy) UNIFORM(nz, nx, ny, nzx)
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE 
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
      INTEGER(C_INT), PARAMETER :: sgntz =  1 
      INTEGER(C_INT), PARAMETER :: sgntx =  1 
      INTEGER(C_INT), PARAMETER :: sgnty =  1 
      tt7 = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      tt5 = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      tt6 = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      tt3 = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      tt4 = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      tt2 = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      tt1 = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      tt8 = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      RETURN
      END SUBROUTINE

!>    @brief  Extracts the travel-times for the 1'th sweep.
!>            Note, that this is machine generated code.
!>
!>    @param[in] i        iz'th grid point.  This is Fortran indexed.
!>    @param[in] j        ix'th grid point.  This is Fortran indexed.
!>    @param[in] k        iy'th grid point.  This is Fortran indexed.
!>    @param[in] ttimes   The travel times at all points in the grid.
!>                        this is a vector of dimension [nz x nx x ny].
!>
!>    @param[out] tt1     Travel time at node (i-sgntz,j,k) i.e, tv.
!>    @param[out] tt2     Travel time at node (i,j-sgntx,k) i.e. te.
!>    @param[out] tt3     Travel time at node (i,j,k-sgnty) i.e. tn.
!>    @param[out] tt4     Travel time at node (i-sgntz,j-sgntx,k) i.e. tev.
!>    @param[out] tt5     Travel time at node (i,j-sgntx,k-sgnty) i.e. ten.
!>    @param[out] tt6     Travel time at node (i-sgntz,j,k-sgnty) i.e. tnv.
!>    @param[out] tt7     Travel time at node (i-sgntz,j-sgntx,k-sgnty) i.e. tnve.
!>    @param[out] tt8     travel time at (i,j,k).
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep1TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                     ttimes, ttLoc)     &
                 BIND(C, NAME='fteik_prefetchSweep1TravelTimes64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: ttLoc(8)
      INTEGER(C_INT), PARAMETER :: sgntz =  1
      INTEGER(C_INT), PARAMETER :: sgntx =  1
      INTEGER(C_INT), PARAMETER :: sgnty =  1
      ttLoc(7) = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      ttLoc(5) = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      ttLoc(6) = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      ttLoc(3) = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      ttLoc(4) = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      ttLoc(2) = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      ttLoc(1) = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      ttLoc(8) = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      RETURN
      END SUBROUTINE
!>    @brief  Extracts the travel-times for the 2'th sweep.
!>            Note, that this is machine generated code.
!>
!>    @param[in] i        iz'th grid point.  This is Fortran indexed.
!>    @param[in] j        ix'th grid point.  This is Fortran indexed.
!>    @param[in] k        iy'th grid point.  This is Fortran indexed.
!>    @param[in] ttimes   The travel times at all points in the grid.
!>                        this is a vector of dimension [nz x nx x ny].
!>
!>    @param[out] tt1     Travel time at node (i-sgntz,j,k) i.e, tv.
!>    @param[out] tt2     Travel time at node (i,j-sgntx,k) i.e. te.
!>    @param[out] tt3     Travel time at node (i,j,k-sgnty) i.e. tn.
!>    @param[out] tt4     Travel time at node (i-sgntz,j-sgntx,k) i.e. tev.
!>    @param[out] tt5     Travel time at node (i,j-sgntx,k-sgnty) i.e. ten.
!>    @param[out] tt6     Travel time at node (i-sgntz,j,k-sgnty) i.e. tnv.
!>    @param[out] tt7     Travel time at node (i-sgntz,j-sgntx,k-sgnty) i.e. tnve.
!>    @param[out] tt8     travel time at (i,j,k).
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep2TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                     ttimes, ttLoc)     &
                 BIND(C, NAME='fteik_prefetchSweep2TravelTimes64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: ttLoc(8)
      INTEGER(C_INT), PARAMETER :: sgntz =  1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty =  1
      ttLoc(6) = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      ttLoc(3) = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      ttLoc(7) = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      ttLoc(5) = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      ttLoc(1) = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      ttLoc(8) = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      ttLoc(4) = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      ttLoc(2) = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      RETURN
      END SUBROUTINE
!>    @brief  Extracts the travel-times for the 3'th sweep.
!>            Note, that this is machine generated code.
!>
!>    @param[in] i        iz'th grid point.  This is Fortran indexed.
!>    @param[in] j        ix'th grid point.  This is Fortran indexed.
!>    @param[in] k        iy'th grid point.  This is Fortran indexed.
!>    @param[in] ttimes   The travel times at all points in the grid.
!>                        this is a vector of dimension [nz x nx x ny].
!>
!>    @param[out] tt1     Travel time at node (i-sgntz,j,k) i.e, tv.
!>    @param[out] tt2     Travel time at node (i,j-sgntx,k) i.e. te.
!>    @param[out] tt3     Travel time at node (i,j,k-sgnty) i.e. tn.
!>    @param[out] tt4     Travel time at node (i-sgntz,j-sgntx,k) i.e. tev.
!>    @param[out] tt5     Travel time at node (i,j-sgntx,k-sgnty) i.e. ten.
!>    @param[out] tt6     Travel time at node (i-sgntz,j,k-sgnty) i.e. tnv.
!>    @param[out] tt7     Travel time at node (i-sgntz,j-sgntx,k-sgnty) i.e. tnve.
!>    @param[out] tt8     travel time at (i,j,k).
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep3TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                     ttimes, ttLoc)     &
                 BIND(C, NAME='fteik_prefetchSweep3TravelTimes64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: ttLoc(8)
      INTEGER(C_INT), PARAMETER :: sgntz =  1
      INTEGER(C_INT), PARAMETER :: sgntx =  1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      ttLoc(4) = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      ttLoc(2) = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      ttLoc(1) = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      ttLoc(8) = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      ttLoc(7) = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      ttLoc(5) = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      ttLoc(6) = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      ttLoc(3) = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      RETURN
      END SUBROUTINE
!>    @brief  Extracts the travel-times for the 4'th sweep.
!>            Note, that this is machine generated code.
!>
!>    @param[in] i        iz'th grid point.  This is Fortran indexed.
!>    @param[in] j        ix'th grid point.  This is Fortran indexed.
!>    @param[in] k        iy'th grid point.  This is Fortran indexed.
!>    @param[in] ttimes   The travel times at all points in the grid.
!>                        this is a vector of dimension [nz x nx x ny].
!>
!>    @param[out] tt1     Travel time at node (i-sgntz,j,k) i.e, tv.
!>    @param[out] tt2     Travel time at node (i,j-sgntx,k) i.e. te.
!>    @param[out] tt3     Travel time at node (i,j,k-sgnty) i.e. tn.
!>    @param[out] tt4     Travel time at node (i-sgntz,j-sgntx,k) i.e. tev.
!>    @param[out] tt5     Travel time at node (i,j-sgntx,k-sgnty) i.e. ten.
!>    @param[out] tt6     Travel time at node (i-sgntz,j,k-sgnty) i.e. tnv.
!>    @param[out] tt7     Travel time at node (i-sgntz,j-sgntx,k-sgnty) i.e. tnve.
!>    @param[out] tt8     travel time at (i,j,k).
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep4TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                     ttimes, ttLoc)     &
                 BIND(C, NAME='fteik_prefetchSweep4TravelTimes64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: ttLoc(8)
      INTEGER(C_INT), PARAMETER :: sgntz =  1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      ttLoc(1) = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      ttLoc(8) = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      ttLoc(4) = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      ttLoc(2) = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      ttLoc(6) = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      ttLoc(3) = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      ttLoc(7) = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      ttLoc(5) = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      RETURN
      END SUBROUTINE
!>    @brief  Extracts the travel-times for the 5'th sweep.
!>            Note, that this is machine generated code.
!>
!>    @param[in] i        iz'th grid point.  This is Fortran indexed.
!>    @param[in] j        ix'th grid point.  This is Fortran indexed.
!>    @param[in] k        iy'th grid point.  This is Fortran indexed.
!>    @param[in] ttimes   The travel times at all points in the grid.
!>                        this is a vector of dimension [nz x nx x ny].
!>
!>    @param[out] tt1     Travel time at node (i-sgntz,j,k) i.e, tv.
!>    @param[out] tt2     Travel time at node (i,j-sgntx,k) i.e. te.
!>    @param[out] tt3     Travel time at node (i,j,k-sgnty) i.e. tn.
!>    @param[out] tt4     Travel time at node (i-sgntz,j-sgntx,k) i.e. tev.
!>    @param[out] tt5     Travel time at node (i,j-sgntx,k-sgnty) i.e. ten.
!>    @param[out] tt6     Travel time at node (i-sgntz,j,k-sgnty) i.e. tnv.
!>    @param[out] tt7     Travel time at node (i-sgntz,j-sgntx,k-sgnty) i.e. tnve.
!>    @param[out] tt8     travel time at (i,j,k).
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep5TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                     ttimes, ttLoc)     &
                 BIND(C, NAME='fteik_prefetchSweep5TravelTimes64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: ttLoc(8)
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx =  1
      INTEGER(C_INT), PARAMETER :: sgnty =  1
      ttLoc(5) = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      ttLoc(7) = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      ttLoc(3) = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      ttLoc(6) = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      ttLoc(2) = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      ttLoc(4) = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      ttLoc(8) = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      ttLoc(1) = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      RETURN
      END SUBROUTINE
!>    @brief  Extracts the travel-times for the 6'th sweep.
!>            Note, that this is machine generated code.
!>
!>    @param[in] i        iz'th grid point.  This is Fortran indexed.
!>    @param[in] j        ix'th grid point.  This is Fortran indexed.
!>    @param[in] k        iy'th grid point.  This is Fortran indexed.
!>    @param[in] ttimes   The travel times at all points in the grid.
!>                        this is a vector of dimension [nz x nx x ny].
!>
!>    @param[out] tt1     Travel time at node (i-sgntz,j,k) i.e, tv.
!>    @param[out] tt2     Travel time at node (i,j-sgntx,k) i.e. te.
!>    @param[out] tt3     Travel time at node (i,j,k-sgnty) i.e. tn.
!>    @param[out] tt4     Travel time at node (i-sgntz,j-sgntx,k) i.e. tev.
!>    @param[out] tt5     Travel time at node (i,j-sgntx,k-sgnty) i.e. ten.
!>    @param[out] tt6     Travel time at node (i-sgntz,j,k-sgnty) i.e. tnv.
!>    @param[out] tt7     Travel time at node (i-sgntz,j-sgntx,k-sgnty) i.e. tnve.
!>    @param[out] tt8     travel time at (i,j,k).
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep6TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                     ttimes, ttLoc)     &
                 BIND(C, NAME='fteik_prefetchSweep6TravelTimes64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: ttLoc(8)
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty =  1
      ttLoc(3) = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      ttLoc(6) = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      ttLoc(5) = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      ttLoc(7) = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      ttLoc(8) = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      ttLoc(1) = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      ttLoc(2) = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      ttLoc(4) = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      RETURN
      END SUBROUTINE
!>    @brief  Extracts the travel-times for the 7'th sweep.
!>            Note, that this is machine generated code.
!>
!>    @param[in] i        iz'th grid point.  This is Fortran indexed.
!>    @param[in] j        ix'th grid point.  This is Fortran indexed.
!>    @param[in] k        iy'th grid point.  This is Fortran indexed.
!>    @param[in] ttimes   The travel times at all points in the grid.
!>                        this is a vector of dimension [nz x nx x ny].
!>
!>    @param[out] tt1     Travel time at node (i-sgntz,j,k) i.e, tv.
!>    @param[out] tt2     Travel time at node (i,j-sgntx,k) i.e. te.
!>    @param[out] tt3     Travel time at node (i,j,k-sgnty) i.e. tn.
!>    @param[out] tt4     Travel time at node (i-sgntz,j-sgntx,k) i.e. tev.
!>    @param[out] tt5     Travel time at node (i,j-sgntx,k-sgnty) i.e. ten.
!>    @param[out] tt6     Travel time at node (i-sgntz,j,k-sgnty) i.e. tnv.
!>    @param[out] tt7     Travel time at node (i-sgntz,j-sgntx,k-sgnty) i.e. tnve.
!>    @param[out] tt8     travel time at (i,j,k).
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep7TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                     ttimes, ttLoc)     &
                 BIND(C, NAME='fteik_prefetchSweep7TravelTimes64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: ttLoc(8)
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx =  1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      ttLoc(2) = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      ttLoc(4) = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      ttLoc(8) = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      ttLoc(1) = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      ttLoc(5) = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      ttLoc(7) = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      ttLoc(3) = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      ttLoc(6) = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      RETURN
      END SUBROUTINE
!>    @brief  Extracts the travel-times for the 8'th sweep.
!>            Note, that this is machine generated code.
!>
!>    @param[in] i        iz'th grid point.  This is Fortran indexed.
!>    @param[in] j        ix'th grid point.  This is Fortran indexed.
!>    @param[in] k        iy'th grid point.  This is Fortran indexed.
!>    @param[in] ttimes   The travel times at all points in the grid.
!>                        this is a vector of dimension [nz x nx x ny].
!>
!>    @param[out] tt1     Travel time at node (i-sgntz,j,k) i.e, tv.
!>    @param[out] tt2     Travel time at node (i,j-sgntx,k) i.e. te.
!>    @param[out] tt3     Travel time at node (i,j,k-sgnty) i.e. tn.
!>    @param[out] tt4     Travel time at node (i-sgntz,j-sgntx,k) i.e. tev.
!>    @param[out] tt5     Travel time at node (i,j-sgntx,k-sgnty) i.e. ten.
!>    @param[out] tt6     Travel time at node (i-sgntz,j,k-sgnty) i.e. tnv.
!>    @param[out] tt7     Travel time at node (i-sgntz,j-sgntx,k-sgnty) i.e. tnve.
!>    @param[out] tt8     travel time at (i,j,k).
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep8TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                     ttimes, ttLoc)     &
                 BIND(C, NAME='fteik_prefetchSweep8TravelTimes64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: ttLoc(8)
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      ttLoc(8) = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      ttLoc(1) = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      ttLoc(2) = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      ttLoc(4) = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      ttLoc(3) = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      ttLoc(6) = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      ttLoc(5) = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      ttLoc(7) = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      RETURN
      END SUBROUTINE
!>    @brief Fetches the slowness corresponding to the stencil at grid
!>           point (i,j,k) for the 1'th sweep.  Note, that this is machine
!>           generated code.
!>
!>    @param[in] i           iz'th grid point.  This is Fortran indexed.
!>    @param[in] j           ix'th grid point.  This is Fortran indexed.
!>    @param[in] k           iy'th grid point.  This is Fortran indexed.
!>    @param[in] nzm1        nz - 1.
!>    @param[in] nzm1_nxm1   (nz - 1)*(nx - 1).
!>    @param[in] slow        This is the slowness (s/m) in each cell.  This is a
!>                           a vector of dimension [(nz-1) x (nx-1) x (ny-1)].
!>
!>    @param[out] slowLoc    This is the slowness stencil for the (iz,ix,iy)'th
!>                           grid point for use in the local solver.  It is a vector
!>                           dimension [7].
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep1Slowness64fF(i, j, k, nz, nx, ny, &
                                                       nzm1, nzm1_nxm1, slow, slowLoc) &
                 BIND(C, NAME='fteik_prefetchSweep1Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
      REAL(C_DOUBLE) sloc(19)
      !DIR ATTRIBUTES ALIGN: 64:: sloc
      INTEGER(C_INT) i1, j1, k1,                            &
                     max_im1_1, max_jm1_1, max_km1_1,       &
                     min_im1_nzm1, min_jm1_nxm1, min_km1_nym1
      INTEGER(C_INT), PARAMETER :: sgnvz =  1
      INTEGER(C_INT), PARAMETER :: sgnvx =  1
      INTEGER(C_INT), PARAMETER :: sgnvy =  1
      i1 = i - sgnvz
      j1 = j - sgnvx
      k1 = k - sgnvy
      ! Mins and maxes; i think this is superfluous
      max_im1_1    = MAX(i - 1, 1)
      max_jm1_1    = MAX(j - 1, 1)
      max_km1_1    = MAX(k - 1, 1)
      min_im1_nzm1 = MIN(i, nz - 1)
      min_jm1_nxm1 = MIN(j, nx - 1)
      min_km1_nym1 = MIN(k, ny - 1)
      sloc( 1) = slow(velGrid2indexF(i1,           max_jm1_1,    max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 5) = slow(velGrid2indexF(max_im1_1,    j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 9) = slow(velGrid2indexF(max_im1_1,    max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc(13) = slow(velGrid2indexF(i1, j1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc(15) = slow(velGrid2indexF(i1, max_jm1_1,    k1, nzm1, nzm1_nxm1))
      sloc(17) = slow(velGrid2indexF(max_im1_1,    j1, k1, nzm1, nzm1_nxm1))
      sloc(19) = slow(velGrid2indexF(i1, j1, k1, nzm1, nzm1_nxm1))
      sloc( 6) = slow(velGrid2indexF(min_im1_nzm1, j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc(11) = slow(velGrid2indexF(min_im1_nzm1, max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc(18) = slow(velGrid2indexF(min_im1_nzm1, j1, k1, nzm1, nzm1_nxm1))
      sloc( 3) = slow(velGrid2indexF(i1,           min_jm1_nxm1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc(10) = slow(velGrid2indexF(max_im1_1,    min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc(16) = slow(velGrid2indexF(i1, min_jm1_nxm1, k1, nzm1, nzm1_nxm1))
      sloc(12) = slow(velGrid2indexF(min_im1_nzm1, min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc( 2) = slow(velGrid2indexF(i1,           max_jm1_1,    min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 7) = slow(velGrid2indexF(max_im1_1,    j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(14) = slow(velGrid2indexF(i1, j1, min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 8) = slow(velGrid2indexF(min_im1_nzm1, j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 4) = slow(velGrid2indexF(i1,           min_jm1_nxm1, min_km1_nym1, nzm1, nzm1_nxm1))
      slowLoc(1) = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc(2) = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc(3) = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc(4) = MIN(sloc(13), sloc(14))
      slowLoc(5) = MIN(sloc(15), sloc(16))
      slowLoc(6) = MIN(sloc(17), sloc(18))
      slowLoc(7) = sloc(19)
      RETURN
      END SUBROUTINE
!>    @brief Fetches the slowness corresponding to the stencil at grid
!>           point (i,j,k) for the 2'th sweep.  Note, that this is machine
!>           generated code.
!>
!>    @param[in] i           iz'th grid point.  This is Fortran indexed.
!>    @param[in] j           ix'th grid point.  This is Fortran indexed.
!>    @param[in] k           iy'th grid point.  This is Fortran indexed.
!>    @param[in] nzm1        nz - 1.
!>    @param[in] nzm1_nxm1   (nz - 1)*(nx - 1).
!>    @param[in] slow        This is the slowness (s/m) in each cell.  This is a
!>                           a vector of dimension [(nz-1) x (nx-1) x (ny-1)].
!>
!>    @param[out] slowLoc    This is the slowness stencil for the (iz,ix,iy)'th
!>                           grid point for use in the local solver.  It is a vector
!>                           dimension [7].
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep2Slowness64fF(i, j, k, nz, nx, ny, &
                                                       nzm1, nzm1_nxm1, slow, slowLoc) &
                 BIND(C, NAME='fteik_prefetchSweep2Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
      REAL(C_DOUBLE) sloc(19)
      !DIR ATTRIBUTES ALIGN: 64:: sloc
      INTEGER(C_INT) i1, j1, k1,                            &
                     max_im1_1, max_jm1_1, max_km1_1,       &
                     min_im1_nzm1, min_jm1_nxm1, min_km1_nym1
      INTEGER(C_INT), PARAMETER :: sgnvz =  1
      INTEGER(C_INT), PARAMETER :: sgnvx =  0
      INTEGER(C_INT), PARAMETER :: sgnvy =  1
      i1 = i - sgnvz
      j1 = j - sgnvx
      k1 = k - sgnvy
      ! Mins and maxes; i think this is superfluous
      max_im1_1    = MAX(i - 1, 1)
      max_jm1_1    = MAX(j - 1, 1)
      max_km1_1    = MAX(k - 1, 1)
      min_im1_nzm1 = MIN(i, nz - 1)
      min_jm1_nxm1 = MIN(j, nx - 1)
      min_km1_nym1 = MIN(k, ny - 1)
      sloc( 1) = slow(velGrid2indexF(i1,           max_jm1_1,    max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 9) = slow(velGrid2indexF(max_im1_1,    max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc(15) = slow(velGrid2indexF(i1, max_jm1_1,    k1, nzm1, nzm1_nxm1))
      sloc(11) = slow(velGrid2indexF(min_im1_nzm1, max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc( 3) = slow(velGrid2indexF(i1,           min_jm1_nxm1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 5) = slow(velGrid2indexF(max_im1_1,    j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc(10) = slow(velGrid2indexF(max_im1_1,    min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc(13) = slow(velGrid2indexF(i1, j1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc(16) = slow(velGrid2indexF(i1, min_jm1_nxm1, k1, nzm1, nzm1_nxm1))
      sloc(17) = slow(velGrid2indexF(max_im1_1,    j1, k1, nzm1, nzm1_nxm1))
      sloc(19) = slow(velGrid2indexF(i1, j1, k1, nzm1, nzm1_nxm1))
      sloc( 6) = slow(velGrid2indexF(min_im1_nzm1, j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc(12) = slow(velGrid2indexF(min_im1_nzm1, min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc(18) = slow(velGrid2indexF(min_im1_nzm1, j1, k1, nzm1, nzm1_nxm1))
      sloc( 2) = slow(velGrid2indexF(i1,           max_jm1_1,    min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 4) = slow(velGrid2indexF(i1,           min_jm1_nxm1, min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 7) = slow(velGrid2indexF(max_im1_1,    j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(14) = slow(velGrid2indexF(i1, j1, min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 8) = slow(velGrid2indexF(min_im1_nzm1, j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      slowLoc(1) = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc(2) = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc(3) = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc(4) = MIN(sloc(13), sloc(14))
      slowLoc(5) = MIN(sloc(15), sloc(16))
      slowLoc(6) = MIN(sloc(17), sloc(18))
      slowLoc(7) = sloc(19)
      RETURN
      END SUBROUTINE
!>    @brief Fetches the slowness corresponding to the stencil at grid
!>           point (i,j,k) for the 3'th sweep.  Note, that this is machine
!>           generated code.
!>
!>    @param[in] i           iz'th grid point.  This is Fortran indexed.
!>    @param[in] j           ix'th grid point.  This is Fortran indexed.
!>    @param[in] k           iy'th grid point.  This is Fortran indexed.
!>    @param[in] nzm1        nz - 1.
!>    @param[in] nzm1_nxm1   (nz - 1)*(nx - 1).
!>    @param[in] slow        This is the slowness (s/m) in each cell.  This is a
!>                           a vector of dimension [(nz-1) x (nx-1) x (ny-1)].
!>
!>    @param[out] slowLoc    This is the slowness stencil for the (iz,ix,iy)'th
!>                           grid point for use in the local solver.  It is a vector
!>                           dimension [7].
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep3Slowness64fF(i, j, k, nz, nx, ny, &
                                                       nzm1, nzm1_nxm1, slow, slowLoc) &
                 BIND(C, NAME='fteik_prefetchSweep3Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
      REAL(C_DOUBLE) sloc(19)
      !DIR ATTRIBUTES ALIGN: 64:: sloc
      INTEGER(C_INT) i1, j1, k1,                            &
                     max_im1_1, max_jm1_1, max_km1_1,       &
                     min_im1_nzm1, min_jm1_nxm1, min_km1_nym1
      INTEGER(C_INT), PARAMETER :: sgnvz =  1
      INTEGER(C_INT), PARAMETER :: sgnvx =  1
      INTEGER(C_INT), PARAMETER :: sgnvy =  0
      i1 = i - sgnvz
      j1 = j - sgnvx
      k1 = k - sgnvy
      ! Mins and maxes; i think this is superfluous
      max_im1_1    = MAX(i - 1, 1)
      max_jm1_1    = MAX(j - 1, 1)
      max_km1_1    = MAX(k - 1, 1)
      min_im1_nzm1 = MIN(i, nz - 1)
      min_jm1_nxm1 = MIN(j, nx - 1)
      min_km1_nym1 = MIN(k, ny - 1)
      sloc( 1) = slow(velGrid2indexF(i1,           max_jm1_1,    max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 5) = slow(velGrid2indexF(max_im1_1,    j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc(13) = slow(velGrid2indexF(i1, j1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 6) = slow(velGrid2indexF(min_im1_nzm1, j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 3) = slow(velGrid2indexF(i1,           min_jm1_nxm1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 2) = slow(velGrid2indexF(i1,           max_jm1_1,    min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 7) = slow(velGrid2indexF(max_im1_1,    j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 9) = slow(velGrid2indexF(max_im1_1,    max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc(14) = slow(velGrid2indexF(i1, j1, min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(15) = slow(velGrid2indexF(i1, max_jm1_1,    k1, nzm1, nzm1_nxm1))
      sloc(17) = slow(velGrid2indexF(max_im1_1,    j1, k1, nzm1, nzm1_nxm1))
      sloc(19) = slow(velGrid2indexF(i1, j1, k1, nzm1, nzm1_nxm1))
      sloc( 8) = slow(velGrid2indexF(min_im1_nzm1, j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(11) = slow(velGrid2indexF(min_im1_nzm1, max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc(18) = slow(velGrid2indexF(min_im1_nzm1, j1, k1, nzm1, nzm1_nxm1))
      sloc( 4) = slow(velGrid2indexF(i1,           min_jm1_nxm1, min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(10) = slow(velGrid2indexF(max_im1_1,    min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc(16) = slow(velGrid2indexF(i1, min_jm1_nxm1, k1, nzm1, nzm1_nxm1))
      sloc(12) = slow(velGrid2indexF(min_im1_nzm1, min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      slowLoc(1) = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc(2) = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc(3) = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc(4) = MIN(sloc(13), sloc(14))
      slowLoc(5) = MIN(sloc(15), sloc(16))
      slowLoc(6) = MIN(sloc(17), sloc(18))
      slowLoc(7) = sloc(19)
      RETURN
      END SUBROUTINE
!>    @brief Fetches the slowness corresponding to the stencil at grid
!>           point (i,j,k) for the 4'th sweep.  Note, that this is machine
!>           generated code.
!>
!>    @param[in] i           iz'th grid point.  This is Fortran indexed.
!>    @param[in] j           ix'th grid point.  This is Fortran indexed.
!>    @param[in] k           iy'th grid point.  This is Fortran indexed.
!>    @param[in] nzm1        nz - 1.
!>    @param[in] nzm1_nxm1   (nz - 1)*(nx - 1).
!>    @param[in] slow        This is the slowness (s/m) in each cell.  This is a
!>                           a vector of dimension [(nz-1) x (nx-1) x (ny-1)].
!>
!>    @param[out] slowLoc    This is the slowness stencil for the (iz,ix,iy)'th
!>                           grid point for use in the local solver.  It is a vector
!>                           dimension [7].
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep4Slowness64fF(i, j, k, nz, nx, ny, &
                                                       nzm1, nzm1_nxm1, slow, slowLoc) &
                 BIND(C, NAME='fteik_prefetchSweep4Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
      REAL(C_DOUBLE) sloc(19)
      !DIR ATTRIBUTES ALIGN: 64:: sloc
      INTEGER(C_INT) i1, j1, k1,                            &
                     max_im1_1, max_jm1_1, max_km1_1,       &
                     min_im1_nzm1, min_jm1_nxm1, min_km1_nym1
      INTEGER(C_INT), PARAMETER :: sgnvz =  1
      INTEGER(C_INT), PARAMETER :: sgnvx =  0
      INTEGER(C_INT), PARAMETER :: sgnvy =  0
      i1 = i - sgnvz
      j1 = j - sgnvx
      k1 = k - sgnvy
      ! Mins and maxes; i think this is superfluous
      max_im1_1    = MAX(i - 1, 1)
      max_jm1_1    = MAX(j - 1, 1)
      max_km1_1    = MAX(k - 1, 1)
      min_im1_nzm1 = MIN(i, nz - 1)
      min_jm1_nxm1 = MIN(j, nx - 1)
      min_km1_nym1 = MIN(k, ny - 1)
      sloc( 1) = slow(velGrid2indexF(i1,           max_jm1_1,    max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 3) = slow(velGrid2indexF(i1,           min_jm1_nxm1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 5) = slow(velGrid2indexF(max_im1_1,    j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc(13) = slow(velGrid2indexF(i1, j1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 6) = slow(velGrid2indexF(min_im1_nzm1, j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 2) = slow(velGrid2indexF(i1,           max_jm1_1,    min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 9) = slow(velGrid2indexF(max_im1_1,    max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc(15) = slow(velGrid2indexF(i1, max_jm1_1,    k1, nzm1, nzm1_nxm1))
      sloc(11) = slow(velGrid2indexF(min_im1_nzm1, max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc( 4) = slow(velGrid2indexF(i1,           min_jm1_nxm1, min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 7) = slow(velGrid2indexF(max_im1_1,    j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(10) = slow(velGrid2indexF(max_im1_1,    min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc(14) = slow(velGrid2indexF(i1, j1, min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(16) = slow(velGrid2indexF(i1, min_jm1_nxm1, k1, nzm1, nzm1_nxm1))
      sloc(17) = slow(velGrid2indexF(max_im1_1,    j1, k1, nzm1, nzm1_nxm1))
      sloc(19) = slow(velGrid2indexF(i1, j1, k1, nzm1, nzm1_nxm1))
      sloc( 8) = slow(velGrid2indexF(min_im1_nzm1, j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(12) = slow(velGrid2indexF(min_im1_nzm1, min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc(18) = slow(velGrid2indexF(min_im1_nzm1, j1, k1, nzm1, nzm1_nxm1))
      slowLoc(1) = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc(2) = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc(3) = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc(4) = MIN(sloc(13), sloc(14))
      slowLoc(5) = MIN(sloc(15), sloc(16))
      slowLoc(6) = MIN(sloc(17), sloc(18))
      slowLoc(7) = sloc(19)
      RETURN
      END SUBROUTINE
!>    @brief Fetches the slowness corresponding to the stencil at grid
!>           point (i,j,k) for the 5'th sweep.  Note, that this is machine
!>           generated code.
!>
!>    @param[in] i           iz'th grid point.  This is Fortran indexed.
!>    @param[in] j           ix'th grid point.  This is Fortran indexed.
!>    @param[in] k           iy'th grid point.  This is Fortran indexed.
!>    @param[in] nzm1        nz - 1.
!>    @param[in] nzm1_nxm1   (nz - 1)*(nx - 1).
!>    @param[in] slow        This is the slowness (s/m) in each cell.  This is a
!>                           a vector of dimension [(nz-1) x (nx-1) x (ny-1)].
!>
!>    @param[out] slowLoc    This is the slowness stencil for the (iz,ix,iy)'th
!>                           grid point for use in the local solver.  It is a vector
!>                           dimension [7].
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep5Slowness64fF(i, j, k, nz, nx, ny, &
                                                       nzm1, nzm1_nxm1, slow, slowLoc) &
                 BIND(C, NAME='fteik_prefetchSweep5Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
      REAL(C_DOUBLE) sloc(19)
      !DIR ATTRIBUTES ALIGN: 64:: sloc
      INTEGER(C_INT) i1, j1, k1,                            &
                     max_im1_1, max_jm1_1, max_km1_1,       &
                     min_im1_nzm1, min_jm1_nxm1, min_km1_nym1
      INTEGER(C_INT), PARAMETER :: sgnvz =  0
      INTEGER(C_INT), PARAMETER :: sgnvx =  1
      INTEGER(C_INT), PARAMETER :: sgnvy =  1
      i1 = i - sgnvz
      j1 = j - sgnvx
      k1 = k - sgnvy
      ! Mins and maxes; i think this is superfluous
      max_im1_1    = MAX(i - 1, 1)
      max_jm1_1    = MAX(j - 1, 1)
      max_km1_1    = MAX(k - 1, 1)
      min_im1_nzm1 = MIN(i, nz - 1)
      min_jm1_nxm1 = MIN(j, nx - 1)
      min_km1_nym1 = MIN(k, ny - 1)
      sloc( 5) = slow(velGrid2indexF(max_im1_1,    j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 9) = slow(velGrid2indexF(max_im1_1,    max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc(17) = slow(velGrid2indexF(max_im1_1,    j1, k1, nzm1, nzm1_nxm1))
      sloc( 1) = slow(velGrid2indexF(i1,           max_jm1_1,    max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 6) = slow(velGrid2indexF(min_im1_nzm1, j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc(11) = slow(velGrid2indexF(min_im1_nzm1, max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc(13) = slow(velGrid2indexF(i1, j1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc(15) = slow(velGrid2indexF(i1, max_jm1_1,    k1, nzm1, nzm1_nxm1))
      sloc(18) = slow(velGrid2indexF(min_im1_nzm1, j1, k1, nzm1, nzm1_nxm1))
      sloc(19) = slow(velGrid2indexF(i1, j1, k1, nzm1, nzm1_nxm1))
      sloc(10) = slow(velGrid2indexF(max_im1_1,    min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc( 3) = slow(velGrid2indexF(i1,           min_jm1_nxm1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc(12) = slow(velGrid2indexF(min_im1_nzm1, min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc(16) = slow(velGrid2indexF(i1, min_jm1_nxm1, k1, nzm1, nzm1_nxm1))
      sloc( 7) = slow(velGrid2indexF(max_im1_1,    j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 2) = slow(velGrid2indexF(i1,           max_jm1_1,    min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 8) = slow(velGrid2indexF(min_im1_nzm1, j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(14) = slow(velGrid2indexF(i1, j1, min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 4) = slow(velGrid2indexF(i1,           min_jm1_nxm1, min_km1_nym1, nzm1, nzm1_nxm1))
      slowLoc(1) = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc(2) = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc(3) = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc(4) = MIN(sloc(13), sloc(14))
      slowLoc(5) = MIN(sloc(15), sloc(16))
      slowLoc(6) = MIN(sloc(17), sloc(18))
      slowLoc(7) = sloc(19)
      RETURN
      END SUBROUTINE
!>    @brief Fetches the slowness corresponding to the stencil at grid
!>           point (i,j,k) for the 6'th sweep.  Note, that this is machine
!>           generated code.
!>
!>    @param[in] i           iz'th grid point.  This is Fortran indexed.
!>    @param[in] j           ix'th grid point.  This is Fortran indexed.
!>    @param[in] k           iy'th grid point.  This is Fortran indexed.
!>    @param[in] nzm1        nz - 1.
!>    @param[in] nzm1_nxm1   (nz - 1)*(nx - 1).
!>    @param[in] slow        This is the slowness (s/m) in each cell.  This is a
!>                           a vector of dimension [(nz-1) x (nx-1) x (ny-1)].
!>
!>    @param[out] slowLoc    This is the slowness stencil for the (iz,ix,iy)'th
!>                           grid point for use in the local solver.  It is a vector
!>                           dimension [7].
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep6Slowness64fF(i, j, k, nz, nx, ny, &
                                                       nzm1, nzm1_nxm1, slow, slowLoc) &
                 BIND(C, NAME='fteik_prefetchSweep6Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
      REAL(C_DOUBLE) sloc(19)
      !DIR ATTRIBUTES ALIGN: 64:: sloc
      INTEGER(C_INT) i1, j1, k1,                            &
                     max_im1_1, max_jm1_1, max_km1_1,       &
                     min_im1_nzm1, min_jm1_nxm1, min_km1_nym1
      INTEGER(C_INT), PARAMETER :: sgnvz =  0
      INTEGER(C_INT), PARAMETER :: sgnvx =  0
      INTEGER(C_INT), PARAMETER :: sgnvy =  1
      i1 = i - sgnvz
      j1 = j - sgnvx
      k1 = k - sgnvy
      ! Mins and maxes; i think this is superfluous
      max_im1_1    = MAX(i - 1, 1)
      max_jm1_1    = MAX(j - 1, 1)
      max_km1_1    = MAX(k - 1, 1)
      min_im1_nzm1 = MIN(i, nz - 1)
      min_jm1_nxm1 = MIN(j, nx - 1)
      min_km1_nym1 = MIN(k, ny - 1)
      sloc( 9) = slow(velGrid2indexF(max_im1_1,    max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc( 1) = slow(velGrid2indexF(i1,           max_jm1_1,    max_km1_1,    nzm1, nzm1_nxm1))
      sloc(11) = slow(velGrid2indexF(min_im1_nzm1, max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc(15) = slow(velGrid2indexF(i1, max_jm1_1,    k1, nzm1, nzm1_nxm1))
      sloc( 5) = slow(velGrid2indexF(max_im1_1,    j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc(10) = slow(velGrid2indexF(max_im1_1,    min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc(17) = slow(velGrid2indexF(max_im1_1,    j1, k1, nzm1, nzm1_nxm1))
      sloc( 3) = slow(velGrid2indexF(i1,           min_jm1_nxm1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 6) = slow(velGrid2indexF(min_im1_nzm1, j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc(12) = slow(velGrid2indexF(min_im1_nzm1, min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc(13) = slow(velGrid2indexF(i1, j1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc(16) = slow(velGrid2indexF(i1, min_jm1_nxm1, k1, nzm1, nzm1_nxm1))
      sloc(18) = slow(velGrid2indexF(min_im1_nzm1, j1, k1, nzm1, nzm1_nxm1))
      sloc(19) = slow(velGrid2indexF(i1, j1, k1, nzm1, nzm1_nxm1))
      sloc( 2) = slow(velGrid2indexF(i1,           max_jm1_1,    min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 7) = slow(velGrid2indexF(max_im1_1,    j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 4) = slow(velGrid2indexF(i1,           min_jm1_nxm1, min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 8) = slow(velGrid2indexF(min_im1_nzm1, j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(14) = slow(velGrid2indexF(i1, j1, min_km1_nym1, nzm1, nzm1_nxm1))
      slowLoc(1) = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc(2) = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc(3) = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc(4) = MIN(sloc(13), sloc(14))
      slowLoc(5) = MIN(sloc(15), sloc(16))
      slowLoc(6) = MIN(sloc(17), sloc(18))
      slowLoc(7) = sloc(19)
      RETURN
      END SUBROUTINE
!>    @brief Fetches the slowness corresponding to the stencil at grid
!>           point (i,j,k) for the 7'th sweep.  Note, that this is machine
!>           generated code.
!>
!>    @param[in] i           iz'th grid point.  This is Fortran indexed.
!>    @param[in] j           ix'th grid point.  This is Fortran indexed.
!>    @param[in] k           iy'th grid point.  This is Fortran indexed.
!>    @param[in] nzm1        nz - 1.
!>    @param[in] nzm1_nxm1   (nz - 1)*(nx - 1).
!>    @param[in] slow        This is the slowness (s/m) in each cell.  This is a
!>                           a vector of dimension [(nz-1) x (nx-1) x (ny-1)].
!>
!>    @param[out] slowLoc    This is the slowness stencil for the (iz,ix,iy)'th
!>                           grid point for use in the local solver.  It is a vector
!>                           dimension [7].
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep7Slowness64fF(i, j, k, nz, nx, ny, &
                                                       nzm1, nzm1_nxm1, slow, slowLoc) &
                 BIND(C, NAME='fteik_prefetchSweep7Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
      REAL(C_DOUBLE) sloc(19)
      !DIR ATTRIBUTES ALIGN: 64:: sloc
      INTEGER(C_INT) i1, j1, k1,                            &
                     max_im1_1, max_jm1_1, max_km1_1,       &
                     min_im1_nzm1, min_jm1_nxm1, min_km1_nym1
      INTEGER(C_INT), PARAMETER :: sgnvz =  0
      INTEGER(C_INT), PARAMETER :: sgnvx =  1
      INTEGER(C_INT), PARAMETER :: sgnvy =  0
      i1 = i - sgnvz
      j1 = j - sgnvx
      k1 = k - sgnvy
      ! Mins and maxes; i think this is superfluous
      max_im1_1    = MAX(i - 1, 1)
      max_jm1_1    = MAX(j - 1, 1)
      max_km1_1    = MAX(k - 1, 1)
      min_im1_nzm1 = MIN(i, nz - 1)
      min_jm1_nxm1 = MIN(j, nx - 1)
      min_km1_nym1 = MIN(k, ny - 1)
      sloc( 5) = slow(velGrid2indexF(max_im1_1,    j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 1) = slow(velGrid2indexF(i1,           max_jm1_1,    max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 6) = slow(velGrid2indexF(min_im1_nzm1, j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc(13) = slow(velGrid2indexF(i1, j1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 3) = slow(velGrid2indexF(i1,           min_jm1_nxm1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 7) = slow(velGrid2indexF(max_im1_1,    j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 9) = slow(velGrid2indexF(max_im1_1,    max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc(17) = slow(velGrid2indexF(max_im1_1,    j1, k1, nzm1, nzm1_nxm1))
      sloc( 2) = slow(velGrid2indexF(i1,           max_jm1_1,    min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 8) = slow(velGrid2indexF(min_im1_nzm1, j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(11) = slow(velGrid2indexF(min_im1_nzm1, max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc(14) = slow(velGrid2indexF(i1, j1, min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(15) = slow(velGrid2indexF(i1, max_jm1_1,    k1, nzm1, nzm1_nxm1))
      sloc(18) = slow(velGrid2indexF(min_im1_nzm1, j1, k1, nzm1, nzm1_nxm1))
      sloc(19) = slow(velGrid2indexF(i1, j1, k1, nzm1, nzm1_nxm1))
      sloc(10) = slow(velGrid2indexF(max_im1_1,    min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc( 4) = slow(velGrid2indexF(i1,           min_jm1_nxm1, min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(12) = slow(velGrid2indexF(min_im1_nzm1, min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc(16) = slow(velGrid2indexF(i1, min_jm1_nxm1, k1, nzm1, nzm1_nxm1))
      slowLoc(1) = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc(2) = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc(3) = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc(4) = MIN(sloc(13), sloc(14))
      slowLoc(5) = MIN(sloc(15), sloc(16))
      slowLoc(6) = MIN(sloc(17), sloc(18))
      slowLoc(7) = sloc(19)
      RETURN
      END SUBROUTINE
!>    @brief Fetches the slowness corresponding to the stencil at grid
!>           point (i,j,k) for the 8'th sweep.  Note, that this is machine
!>           generated code.
!>
!>    @param[in] i           iz'th grid point.  This is Fortran indexed.
!>    @param[in] j           ix'th grid point.  This is Fortran indexed.
!>    @param[in] k           iy'th grid point.  This is Fortran indexed.
!>    @param[in] nzm1        nz - 1.
!>    @param[in] nzm1_nxm1   (nz - 1)*(nx - 1).
!>    @param[in] slow        This is the slowness (s/m) in each cell.  This is a
!>                           a vector of dimension [(nz-1) x (nx-1) x (ny-1)].
!>
!>    @param[out] slowLoc    This is the slowness stencil for the (iz,ix,iy)'th
!>                           grid point for use in the local solver.  It is a vector
!>                           dimension [7].
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_prefetchSweep8Slowness64fF(i, j, k, nz, nx, ny, &
                                                       nzm1, nzm1_nxm1, slow, slowLoc) &
                 BIND(C, NAME='fteik_prefetchSweep8Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
      REAL(C_DOUBLE) sloc(19)
      !DIR ATTRIBUTES ALIGN: 64:: sloc
      INTEGER(C_INT) i1, j1, k1,                            &
                     max_im1_1, max_jm1_1, max_km1_1,       &
                     min_im1_nzm1, min_jm1_nxm1, min_km1_nym1
      INTEGER(C_INT), PARAMETER :: sgnvz =  0
      INTEGER(C_INT), PARAMETER :: sgnvx =  0
      INTEGER(C_INT), PARAMETER :: sgnvy =  0
      i1 = i - sgnvz
      j1 = j - sgnvx
      k1 = k - sgnvy
      ! Mins and maxes; i think this is superfluous
      max_im1_1    = MAX(i - 1, 1)
      max_jm1_1    = MAX(j - 1, 1)
      max_km1_1    = MAX(k - 1, 1)
      min_im1_nzm1 = MIN(i, nz - 1)
      min_jm1_nxm1 = MIN(j, nx - 1)
      min_km1_nym1 = MIN(k, ny - 1)
      sloc( 1) = slow(velGrid2indexF(i1,           max_jm1_1,    max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 5) = slow(velGrid2indexF(max_im1_1,    j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 3) = slow(velGrid2indexF(i1,           min_jm1_nxm1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 6) = slow(velGrid2indexF(min_im1_nzm1, j1,           max_km1_1,    nzm1, nzm1_nxm1))
      sloc(13) = slow(velGrid2indexF(i1, j1, max_km1_1,    nzm1, nzm1_nxm1))
      sloc( 9) = slow(velGrid2indexF(max_im1_1,    max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc( 2) = slow(velGrid2indexF(i1,           max_jm1_1,    min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(11) = slow(velGrid2indexF(min_im1_nzm1, max_jm1_1,    k1,           nzm1, nzm1_nxm1))
      sloc(15) = slow(velGrid2indexF(i1, max_jm1_1,    k1, nzm1, nzm1_nxm1))
      sloc( 7) = slow(velGrid2indexF(max_im1_1,    j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(10) = slow(velGrid2indexF(max_im1_1,    min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc(17) = slow(velGrid2indexF(max_im1_1,    j1, k1, nzm1, nzm1_nxm1))
      sloc( 4) = slow(velGrid2indexF(i1,           min_jm1_nxm1, min_km1_nym1, nzm1, nzm1_nxm1))
      sloc( 8) = slow(velGrid2indexF(min_im1_nzm1, j1,           min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(12) = slow(velGrid2indexF(min_im1_nzm1, min_jm1_nxm1, k1,           nzm1, nzm1_nxm1))
      sloc(14) = slow(velGrid2indexF(i1, j1, min_km1_nym1, nzm1, nzm1_nxm1))
      sloc(16) = slow(velGrid2indexF(i1, min_jm1_nxm1, k1, nzm1, nzm1_nxm1))
      sloc(18) = slow(velGrid2indexF(min_im1_nzm1, j1, k1, nzm1, nzm1_nxm1))
      sloc(19) = slow(velGrid2indexF(i1, j1, k1, nzm1, nzm1_nxm1))
      slowLoc(1) = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc(2) = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc(3) = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc(4) = MIN(sloc(13), sloc(14))
      slowLoc(5) = MIN(sloc(15), sloc(16))
      slowLoc(6) = MIN(sloc(17), sloc(18))
      slowLoc(7) = sloc(19)
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep1LS64fF(linitk, ttimes, ierr) &
                 BIND(C, NAME='fteik_evaluateSweep1LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep1Slowness64fF, & 
                                 fteik_prefetchSweep1TravelTimes64fF
      USE FTEIK_UTILS64F, ONLY : levelPtr, lupd1, lupdInit1, ijkv1, slow, &
                                 dxi, dyi, dzi, &
                                 nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_UTILS64F, ONLY : chunkSize
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, k2, knode, level, loop, node, node1, node2, n2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:)
      INTEGER(C_INT), PARAMETER :: sweep = 1
      INTEGER(C_INT), PARAMETER :: sgntz = 1
      INTEGER(C_INT), PARAMETER :: sgntx = 1
      INTEGER(C_INT), PARAMETER :: sgnty = 1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
REAL(C_DOUBLE) sw(8), tw(8)
!DIR$ ATTRIBUTES align:64 :: sw, tw
      ierr = 0
      ! Some derivative items
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      slowWork(:) = 0.d0
      ttWork(:) = 0.d0
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv1, lupd1, l1, l2, k2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(n2, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, node, node1, node2, tt1, tw, sw) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=1,k2 !,chunkSize
               node1 = l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + n2 - 1)
!              DO node=node1,node2 !loop=1,n2
!                 loop = node - node1 + 1 !node = knode + loop - 1
!                 IF (lupd1(node)) THEN
!                    i    = ijkv1(4*(node-1)+1)
!                    j    = ijkv1(4*(node-1)+2)
!                    k    = ijkv1(4*(node-1)+3)
!                    !DIR$ FORCEINLINE
!                    CALL fteik_prefetchSweep1Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
!                                                          slow, slowWork(8*(loop-1)+1))
!                 ENDIF
!              ENDDO
!              DO node=node1,node2 !loop=1,n2
!                 loop = node - node1 + 1 !node = knode + loop - 1
!                 IF (lupd1(node)) THEN
!                    i    = ijkv1(4*(node-1)+1)
!                    j    = ijkv1(4*(node-1)+2)
!                    k    = ijkv1(4*(node-1)+3)
!                    !DIR$ FORCEINLINE
!                    CALL fteik_prefetchSweep1TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
!                                                             ttimes, ttWork(8*(loop-1)+1))
!                 ENDIF
!              ENDDO
               !DIR$ VECTOR ALWAYS
               DO node=node1,node2 !loop=1,n2
                  !loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd1(node)) THEN
                     i    = ijkv1(4*(node-1)+1)
                     j    = ijkv1(4*(node-1)+2)
                     k    = ijkv1(4*(node-1)+3)
                     indx = ijkv1(4*(node-1)+4)
                     CALL fteik_prefetchSweep1Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, sw) !slowWork(1))
                     CALL fteik_prefetchSweep1TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, tw) !ttWork(1))
                     ttimes(indx) = fteik_localSolver64fF(tw, sw, & !ttWork(1), slowWork(1), &
                                                          .FALSE., i, j, k, &
                                                          sgntz, sgntx, sgnty, sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
!                    ttimes(indx) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
!                                                 slowWork(8*(loop-1)+1),  &
!                                                 .FALSE.,                        &
!                                                 i, j, k,                       &
!                                                 sgntz, sgntx, sgnty,           &
!                                                 sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
!              !$OMP SIMD
!              DO node=node1,node2 !loop=1,n2
!                 loop = node - node1 + 1 !node = knode + loop - 1
!                 IF (lupd1(node)) THEN
!                    indx = ijkv1(4*(node-1)+4)
!                    ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
!                 ENDIF
!              ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 1       CONTINUE
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv1, lupdInit1, l1, l2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, n2, node, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=l1,l2!,chunkSize
               n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit1(node)) THEN
                     i    = ijkv1(4*(node-1)+1)
                     j    = ijkv1(4*(node-1)+2)
                     k    = ijkv1(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep1Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit1(node)) THEN
                     i    = ijkv1(4*(node-1)+1)
                     j    = ijkv1(4*(node-1)+2)
                     k    = ijkv1(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep1TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit1(node)) THEN
                     i    = ijkv1(4*(node-1)+1)
                     j    = ijkv1(4*(node-1)+2)
                     k    = ijkv1(4*(node-1)+3)
                     indx = ijkv1(4*(node-1)+4)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .TRUE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit1(node)) THEN
                     indx = ijkv1(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 11      CONTINUE
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep2LS64fF(linitk, ttimes, ierr) &
                 BIND(C, NAME='fteik_evaluateSweep2LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep2Slowness64fF, & 
                                 fteik_prefetchSweep2TravelTimes64fF
      USE FTEIK_UTILS64F, ONLY : levelPtr, lupd2, lupdInit2, ijkv2, slow, &
                                 dxi, dyi, dzi, &
                                 nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_UTILS64F, ONLY : chunkSize
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, k2, knode, level, loop, node, node1, node2, n2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:)
      INTEGER(C_INT), PARAMETER :: sweep = 2
      INTEGER(C_INT), PARAMETER :: sgntz = 1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty = 1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      ierr = 0
      ! Some derivative items
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      slowWork(:) = 0.d0
      ttWork(:) = 0.d0
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv2, lupd2, l1, l2, k2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(n2, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, node, node1, node2, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=1,k2 !,chunkSize
               node1 = l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd2(node)) THEN
                     i    = ijkv2(4*(node-1)+1)
                     j    = ijkv2(4*(node-1)+2)
                     k    = ijkv2(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep2Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd2(node)) THEN
                     i    = ijkv2(4*(node-1)+1)
                     j    = ijkv2(4*(node-1)+2)
                     k    = ijkv2(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep2TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd2(node)) THEN
                     i    = ijkv2(4*(node-1)+1)
                     j    = ijkv2(4*(node-1)+2)
                     k    = ijkv2(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .FALSE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd2(node)) THEN
                     indx = ijkv2(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 1       CONTINUE
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv2, lupdInit2, l1, l2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, n2, node, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=l1,l2!,chunkSize
               n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit2(node)) THEN
                     i    = ijkv2(4*(node-1)+1)
                     j    = ijkv2(4*(node-1)+2)
                     k    = ijkv2(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep2Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit2(node)) THEN
                     i    = ijkv2(4*(node-1)+1)
                     j    = ijkv2(4*(node-1)+2)
                     k    = ijkv2(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep2TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit2(node)) THEN
                     i    = ijkv2(4*(node-1)+1)
                     j    = ijkv2(4*(node-1)+2)
                     k    = ijkv2(4*(node-1)+3)
                     indx = ijkv2(4*(node-1)+4)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .TRUE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit2(node)) THEN
                     indx = ijkv2(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 11      CONTINUE
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep3LS64fF(linitk, ttimes, ierr) &
                 BIND(C, NAME='fteik_evaluateSweep3LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep3Slowness64fF, & 
                                 fteik_prefetchSweep3TravelTimes64fF
      USE FTEIK_UTILS64F, ONLY : levelPtr, lupd3, lupdInit3, ijkv3, slow, &
                                 dxi, dyi, dzi, &
                                 nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_UTILS64F, ONLY : chunkSize
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, k2, knode, level, loop, node, node1, node2, n2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:)
      INTEGER(C_INT), PARAMETER :: sweep = 3
      INTEGER(C_INT), PARAMETER :: sgntz = 1
      INTEGER(C_INT), PARAMETER :: sgntx = 1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      ierr = 0
      ! Some derivative items
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      slowWork(:) = 0.d0
      ttWork(:) = 0.d0
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv3, lupd3, l1, l2, k2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(n2, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, node, node1, node2, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=1,k2 !,chunkSize
               node1 = l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd3(node)) THEN
                     i    = ijkv3(4*(node-1)+1)
                     j    = ijkv3(4*(node-1)+2)
                     k    = ijkv3(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep3Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd3(node)) THEN
                     i    = ijkv3(4*(node-1)+1)
                     j    = ijkv3(4*(node-1)+2)
                     k    = ijkv3(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep3TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd3(node)) THEN
                     i    = ijkv3(4*(node-1)+1)
                     j    = ijkv3(4*(node-1)+2)
                     k    = ijkv3(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .FALSE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd3(node)) THEN
                     indx = ijkv3(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 1       CONTINUE
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv3, lupdInit3, l1, l2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, n2, node, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=l1,l2!,chunkSize
               n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit3(node)) THEN
                     i    = ijkv3(4*(node-1)+1)
                     j    = ijkv3(4*(node-1)+2)
                     k    = ijkv3(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep3Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit3(node)) THEN
                     i    = ijkv3(4*(node-1)+1)
                     j    = ijkv3(4*(node-1)+2)
                     k    = ijkv3(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep3TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit3(node)) THEN
                     i    = ijkv3(4*(node-1)+1)
                     j    = ijkv3(4*(node-1)+2)
                     k    = ijkv3(4*(node-1)+3)
                     indx = ijkv3(4*(node-1)+4)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .TRUE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit3(node)) THEN
                     indx = ijkv3(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 11      CONTINUE
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep4LS64fF(linitk, ttimes, ierr) &
                 BIND(C, NAME='fteik_evaluateSweep4LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep4Slowness64fF, & 
                                 fteik_prefetchSweep4TravelTimes64fF
      USE FTEIK_UTILS64F, ONLY : levelPtr, lupd4, lupdInit4, ijkv4, slow, &
                                 dxi, dyi, dzi, &
                                 nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_UTILS64F, ONLY : chunkSize
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, k2, knode, level, loop, node, node1, node2, n2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:)
      INTEGER(C_INT), PARAMETER :: sweep = 4
      INTEGER(C_INT), PARAMETER :: sgntz = 1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      ierr = 0
      ! Some derivative items
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      slowWork(:) = 0.d0
      ttWork(:) = 0.d0
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv4, lupd4, l1, l2, k2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(n2, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, node, node1, node2, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=1,k2 !,chunkSize
               node1 = l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd4(node)) THEN
                     i    = ijkv4(4*(node-1)+1)
                     j    = ijkv4(4*(node-1)+2)
                     k    = ijkv4(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep4Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd4(node)) THEN
                     i    = ijkv4(4*(node-1)+1)
                     j    = ijkv4(4*(node-1)+2)
                     k    = ijkv4(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep4TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd4(node)) THEN
                     i    = ijkv4(4*(node-1)+1)
                     j    = ijkv4(4*(node-1)+2)
                     k    = ijkv4(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .FALSE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd4(node)) THEN
                     indx = ijkv4(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 1       CONTINUE
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv4, lupdInit4, l1, l2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, n2, node, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=l1,l2!,chunkSize
               n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit4(node)) THEN
                     i    = ijkv4(4*(node-1)+1)
                     j    = ijkv4(4*(node-1)+2)
                     k    = ijkv4(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep4Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit4(node)) THEN
                     i    = ijkv4(4*(node-1)+1)
                     j    = ijkv4(4*(node-1)+2)
                     k    = ijkv4(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep4TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit4(node)) THEN
                     i    = ijkv4(4*(node-1)+1)
                     j    = ijkv4(4*(node-1)+2)
                     k    = ijkv4(4*(node-1)+3)
                     indx = ijkv4(4*(node-1)+4)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .TRUE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit4(node)) THEN
                     indx = ijkv4(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 11      CONTINUE
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep5LS64fF(linitk, ttimes, ierr) &
                 BIND(C, NAME='fteik_evaluateSweep5LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep5Slowness64fF, & 
                                 fteik_prefetchSweep5TravelTimes64fF
      USE FTEIK_UTILS64F, ONLY : levelPtr, lupd5, lupdInit5, ijkv5, slow, &
                                 dxi, dyi, dzi, &
                                 nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_UTILS64F, ONLY : chunkSize
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, k2, knode, level, loop, node, node1, node2, n2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:)
      INTEGER(C_INT), PARAMETER :: sweep = 5
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx = 1
      INTEGER(C_INT), PARAMETER :: sgnty = 1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      ierr = 0
      ! Some derivative items
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      slowWork(:) = 0.d0
      ttWork(:) = 0.d0
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv5, lupd5, l1, l2, k2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(n2, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, node, node1, node2, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=1,k2 !,chunkSize
               node1 = l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd5(node)) THEN
                     i    = ijkv5(4*(node-1)+1)
                     j    = ijkv5(4*(node-1)+2)
                     k    = ijkv5(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep5Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd5(node)) THEN
                     i    = ijkv5(4*(node-1)+1)
                     j    = ijkv5(4*(node-1)+2)
                     k    = ijkv5(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep5TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd5(node)) THEN
                     i    = ijkv5(4*(node-1)+1)
                     j    = ijkv5(4*(node-1)+2)
                     k    = ijkv5(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .FALSE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd5(node)) THEN
                     indx = ijkv5(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 1       CONTINUE
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv5, lupdInit5, l1, l2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, n2, node, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=l1,l2!,chunkSize
               n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit5(node)) THEN
                     i    = ijkv5(4*(node-1)+1)
                     j    = ijkv5(4*(node-1)+2)
                     k    = ijkv5(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep5Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit5(node)) THEN
                     i    = ijkv5(4*(node-1)+1)
                     j    = ijkv5(4*(node-1)+2)
                     k    = ijkv5(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep5TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit5(node)) THEN
                     i    = ijkv5(4*(node-1)+1)
                     j    = ijkv5(4*(node-1)+2)
                     k    = ijkv5(4*(node-1)+3)
                     indx = ijkv5(4*(node-1)+4)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .TRUE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit5(node)) THEN
                     indx = ijkv5(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 11      CONTINUE
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep6LS64fF(linitk, ttimes, ierr) &
                 BIND(C, NAME='fteik_evaluateSweep6LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep6Slowness64fF, & 
                                 fteik_prefetchSweep6TravelTimes64fF
      USE FTEIK_UTILS64F, ONLY : levelPtr, lupd6, lupdInit6, ijkv6, slow, &
                                 dxi, dyi, dzi, &
                                 nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_UTILS64F, ONLY : chunkSize
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, k2, knode, level, loop, node, node1, node2, n2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:)
      INTEGER(C_INT), PARAMETER :: sweep = 6
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty = 1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      ierr = 0
      ! Some derivative items
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      slowWork(:) = 0.d0
      ttWork(:) = 0.d0
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv6, lupd6, l1, l2, k2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(n2, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, node, node1, node2, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=1,k2 !,chunkSize
               node1 = l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd6(node)) THEN
                     i    = ijkv6(4*(node-1)+1)
                     j    = ijkv6(4*(node-1)+2)
                     k    = ijkv6(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep6Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd6(node)) THEN
                     i    = ijkv6(4*(node-1)+1)
                     j    = ijkv6(4*(node-1)+2)
                     k    = ijkv6(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep6TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd6(node)) THEN
                     i    = ijkv6(4*(node-1)+1)
                     j    = ijkv6(4*(node-1)+2)
                     k    = ijkv6(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .FALSE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd6(node)) THEN
                     indx = ijkv6(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 1       CONTINUE
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv6, lupdInit6, l1, l2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, n2, node, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=l1,l2!,chunkSize
               n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit6(node)) THEN
                     i    = ijkv6(4*(node-1)+1)
                     j    = ijkv6(4*(node-1)+2)
                     k    = ijkv6(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep6Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit6(node)) THEN
                     i    = ijkv6(4*(node-1)+1)
                     j    = ijkv6(4*(node-1)+2)
                     k    = ijkv6(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep6TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit6(node)) THEN
                     i    = ijkv6(4*(node-1)+1)
                     j    = ijkv6(4*(node-1)+2)
                     k    = ijkv6(4*(node-1)+3)
                     indx = ijkv6(4*(node-1)+4)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .TRUE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit6(node)) THEN
                     indx = ijkv6(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 11      CONTINUE
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep7LS64fF(linitk, ttimes, ierr) &
                 BIND(C, NAME='fteik_evaluateSweep7LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep7Slowness64fF, & 
                                 fteik_prefetchSweep7TravelTimes64fF
      USE FTEIK_UTILS64F, ONLY : levelPtr, lupd7, lupdInit7, ijkv7, slow, &
                                 dxi, dyi, dzi, &
                                 nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_UTILS64F, ONLY : chunkSize
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, k2, knode, level, loop, node, node1, node2, n2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:)
      INTEGER(C_INT), PARAMETER :: sweep = 7
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx = 1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      ierr = 0
      ! Some derivative items
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      slowWork(:) = 0.d0
      ttWork(:) = 0.d0
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv7, lupd7, l1, l2, k2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(n2, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, node, node1, node2, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=1,k2 !,chunkSize
               node1 = l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd7(node)) THEN
                     i    = ijkv7(4*(node-1)+1)
                     j    = ijkv7(4*(node-1)+2)
                     k    = ijkv7(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep7Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd7(node)) THEN
                     i    = ijkv7(4*(node-1)+1)
                     j    = ijkv7(4*(node-1)+2)
                     k    = ijkv7(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep7TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd7(node)) THEN
                     i    = ijkv7(4*(node-1)+1)
                     j    = ijkv7(4*(node-1)+2)
                     k    = ijkv7(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .FALSE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd7(node)) THEN
                     indx = ijkv7(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 1       CONTINUE
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv7, lupdInit7, l1, l2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, n2, node, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=l1,l2!,chunkSize
               n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit7(node)) THEN
                     i    = ijkv7(4*(node-1)+1)
                     j    = ijkv7(4*(node-1)+2)
                     k    = ijkv7(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep7Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit7(node)) THEN
                     i    = ijkv7(4*(node-1)+1)
                     j    = ijkv7(4*(node-1)+2)
                     k    = ijkv7(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep7TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit7(node)) THEN
                     i    = ijkv7(4*(node-1)+1)
                     j    = ijkv7(4*(node-1)+2)
                     k    = ijkv7(4*(node-1)+3)
                     indx = ijkv7(4*(node-1)+4)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .TRUE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit7(node)) THEN
                     indx = ijkv7(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 11      CONTINUE
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep8LS64fF(linitk, ttimes, ierr) &
                 BIND(C, NAME='fteik_evaluateSweep8LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep8Slowness64fF, & 
                                 fteik_prefetchSweep8TravelTimes64fF
      USE FTEIK_UTILS64F, ONLY : levelPtr, lupd8, lupdInit8, ijkv8, slow, &
                                 dxi, dyi, dzi, &
                                 nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_UTILS64F, ONLY : chunkSize
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, k2, knode, level, loop, node, node1, node2, n2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:)
      INTEGER(C_INT), PARAMETER :: sweep = 8
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      ierr = 0
      ! Some derivative items
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      slowWork(:) = 0.d0
      ttWork(:) = 0.d0
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv8, lupd8, l1, l2, k2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(n2, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, node, node1, node2, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=1,k2 !,chunkSize
               node1 = l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd8(node)) THEN
                     i    = ijkv8(4*(node-1)+1)
                     j    = ijkv8(4*(node-1)+2)
                     k    = ijkv8(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep8Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd8(node)) THEN
                     i    = ijkv8(4*(node-1)+1)
                     j    = ijkv8(4*(node-1)+2)
                     k    = ijkv8(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep8TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd8(node)) THEN
                     i    = ijkv8(4*(node-1)+1)
                     j    = ijkv8(4*(node-1)+2)
                     k    = ijkv8(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .FALSE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd8(node)) THEN
                     indx = ijkv8(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 1       CONTINUE
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP SHARED(ijkv8, lupdInit8, l1, l2, sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            !$OMP SHARED(nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
            !$OMP SHARED(slow, ttimes) &
            !$OMP PRIVATE(i, indx, j, k, knode, loop, n2, node, tt1) &
            !$OMP FIRSTPRIVATE(ttWork, slowWork)
            DO knode=l1,l2!,chunkSize
               n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit8(node)) THEN
                     i    = ijkv8(4*(node-1)+1)
                     j    = ijkv8(4*(node-1)+2)
                     k    = ijkv8(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep8Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, slowWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit8(node)) THEN
                     i    = ijkv8(4*(node-1)+1)
                     j    = ijkv8(4*(node-1)+2)
                     k    = ijkv8(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep8TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                              ttimes, ttWork(8*(loop-1)+1))
                  ENDIF
               ENDDO
               !DIR$ IVDEP
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit8(node)) THEN
                     i    = ijkv8(4*(node-1)+1)
                     j    = ijkv8(4*(node-1)+2)
                     k    = ijkv8(4*(node-1)+3)
                     indx = ijkv8(4*(node-1)+4)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                                                  slowWork(8*(loop-1)+1),  &
                                                  .TRUE.,                        &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO loop=1,n2
                  node = knode + loop - 1
                  IF (lupdInit8(node)) THEN
                     indx = ijkv8(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
 11      CONTINUE
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      RETURN
      END SUBROUTINE
