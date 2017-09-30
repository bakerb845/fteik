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
      PURE SUBROUTINE fteik_prefetchSweep1TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                         ttimes, tt1, tt2, tt3, tt4, &
                                                         tt5, tt6, tt7, tt8)         &
      BIND(C, NAME='fteik_prefetchSweep1TravelTimes64fF')
      !$OMP DECLARE SIMD(fteik_prefetchSweep1TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx
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
      PURE SUBROUTINE fteik_prefetchSweep2TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                         ttimes, tt1, tt2, tt3, tt4, &
                                                         tt5, tt6, tt7, tt8)         &
      BIND(C, NAME='fteik_prefetchSweep2TravelTimes64fF')
      !$OMP DECLARE SIMD(fteik_prefetchSweep2TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
      INTEGER(C_INT), PARAMETER :: sgntz =  1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty =  1
      tt6 = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      tt3 = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      tt7 = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      tt5 = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      tt1 = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      tt8 = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      tt4 = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      tt2 = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
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
      PURE SUBROUTINE fteik_prefetchSweep3TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                         ttimes, tt1, tt2, tt3, tt4, &
                                                         tt5, tt6, tt7, tt8)         &
      BIND(C, NAME='fteik_prefetchSweep3TravelTimes64fF')
      !$OMP DECLARE SIMD(fteik_prefetchSweep3TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
      INTEGER(C_INT), PARAMETER :: sgntz =  1
      INTEGER(C_INT), PARAMETER :: sgntx =  1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      tt4 = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      tt2 = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      tt1 = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      tt8 = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      tt7 = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      tt5 = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      tt6 = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      tt3 = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
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
      PURE SUBROUTINE fteik_prefetchSweep4TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                         ttimes, tt1, tt2, tt3, tt4, &
                                                         tt5, tt6, tt7, tt8)         &
      BIND(C, NAME='fteik_prefetchSweep4TravelTimes64fF')
      !$OMP DECLARE SIMD(fteik_prefetchSweep4TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
      INTEGER(C_INT), PARAMETER :: sgntz =  1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      tt1 = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      tt8 = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      tt4 = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      tt2 = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      tt6 = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      tt3 = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      tt7 = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      tt5 = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
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
      PURE SUBROUTINE fteik_prefetchSweep5TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                         ttimes, tt1, tt2, tt3, tt4, &
                                                         tt5, tt6, tt7, tt8)         &
      BIND(C, NAME='fteik_prefetchSweep5TravelTimes64fF')
      !$OMP DECLARE SIMD(fteik_prefetchSweep5TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx =  1
      INTEGER(C_INT), PARAMETER :: sgnty =  1
      tt5 = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      tt7 = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      tt3 = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      tt6 = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      tt2 = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      tt4 = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      tt8 = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      tt1 = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
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
      PURE SUBROUTINE fteik_prefetchSweep6TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                         ttimes, tt1, tt2, tt3, tt4, &
                                                         tt5, tt6, tt7, tt8)         &
      BIND(C, NAME='fteik_prefetchSweep6TravelTimes64fF')
      !$OMP DECLARE SIMD(fteik_prefetchSweep6TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty =  1
      tt3 = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      tt6 = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      tt5 = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      tt7 = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      tt8 = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      tt1 = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      tt2 = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      tt4 = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
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
      PURE SUBROUTINE fteik_prefetchSweep7TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                         ttimes, tt1, tt2, tt3, tt4, &
                                                         tt5, tt6, tt7, tt8)         &
      BIND(C, NAME='fteik_prefetchSweep7TravelTimes64fF')
      !$OMP DECLARE SIMD(fteik_prefetchSweep7TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx =  1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      tt2 = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      tt4 = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      tt8 = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      tt1 = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      tt5 = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      tt7 = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
      tt3 = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      tt6 = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
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
      PURE SUBROUTINE fteik_prefetchSweep8TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                         ttimes, tt1, tt2, tt3, tt4, &
                                                         tt5, tt6, tt7, tt8)         &
      BIND(C, NAME='fteik_prefetchSweep8TravelTimes64fF')
      !$OMP DECLARE SIMD(fteik_prefetchSweep8TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) 
      REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      tt8 = ttimes(grid2indexF(i,       j,       k, nz, nzx))
      tt1 = ttimes(grid2indexF(i-sgntz, j,       k, nz, nzx))
      tt2 = ttimes(grid2indexF(i,       j-sgntx, k, nz, nzx))
      tt4 = ttimes(grid2indexF(i-sgntz, j-sgntx, k, nz, nzx))
      tt3 = ttimes(grid2indexF(i,       j,       k-sgnty, nz, nzx))
      tt6 = ttimes(grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx))
      tt5 = ttimes(grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx))
      tt7 = ttimes(grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx))
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
                                                       nzm1, nzm1_nxm1, slow, & !slowLoc) &
                                                       slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                                       slowLoc5, slowLoc6, slowLoc7) &
      BIND(C, NAME='fteik_prefetchSweep1Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      !REAL(C_DOUBLE), INTENT(OUT) :: !slowLoc(7)
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                     slowLoc5, slowLoc6, slowLoc7
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
      slowLoc1 = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc2 = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc3 = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc4 = MIN(sloc(13), sloc(14))
      slowLoc5 = MIN(sloc(15), sloc(16))
      slowLoc6 = MIN(sloc(17), sloc(18))
      slowLoc7 = sloc(19)
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
                                                       nzm1, nzm1_nxm1, slow, & !slowLoc) &
                                                       slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                                       slowLoc5, slowLoc6, slowLoc7) &
      BIND(C, NAME='fteik_prefetchSweep2Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      !REAL(C_DOUBLE), INTENT(OUT) :: !slowLoc(7)
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                     slowLoc5, slowLoc6, slowLoc7
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
      slowLoc1 = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc2 = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc3 = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc4 = MIN(sloc(13), sloc(14))
      slowLoc5 = MIN(sloc(15), sloc(16))
      slowLoc6 = MIN(sloc(17), sloc(18))
      slowLoc7 = sloc(19)
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
                                                       nzm1, nzm1_nxm1, slow, & !slowLoc) &
                                                       slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                                       slowLoc5, slowLoc6, slowLoc7) &
      BIND(C, NAME='fteik_prefetchSweep3Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      !REAL(C_DOUBLE), INTENT(OUT) :: !slowLoc(7)
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                     slowLoc5, slowLoc6, slowLoc7
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
      slowLoc1 = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc2 = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc3 = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc4 = MIN(sloc(13), sloc(14))
      slowLoc5 = MIN(sloc(15), sloc(16))
      slowLoc6 = MIN(sloc(17), sloc(18))
      slowLoc7 = sloc(19)
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
                                                       nzm1, nzm1_nxm1, slow, & !slowLoc) &
                                                       slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                                       slowLoc5, slowLoc6, slowLoc7) &
      BIND(C, NAME='fteik_prefetchSweep4Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      !REAL(C_DOUBLE), INTENT(OUT) :: !slowLoc(7)
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                     slowLoc5, slowLoc6, slowLoc7
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
      slowLoc1 = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc2 = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc3 = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc4 = MIN(sloc(13), sloc(14))
      slowLoc5 = MIN(sloc(15), sloc(16))
      slowLoc6 = MIN(sloc(17), sloc(18))
      slowLoc7 = sloc(19)
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
                                                       nzm1, nzm1_nxm1, slow, & !slowLoc) &
                                                       slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                                       slowLoc5, slowLoc6, slowLoc7) &
      BIND(C, NAME='fteik_prefetchSweep5Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      !REAL(C_DOUBLE), INTENT(OUT) :: !slowLoc(7)
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                     slowLoc5, slowLoc6, slowLoc7
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
      slowLoc1 = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc2 = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc3 = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc4 = MIN(sloc(13), sloc(14))
      slowLoc5 = MIN(sloc(15), sloc(16))
      slowLoc6 = MIN(sloc(17), sloc(18))
      slowLoc7 = sloc(19)
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
                                                       nzm1, nzm1_nxm1, slow, & !slowLoc) &
                                                       slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                                       slowLoc5, slowLoc6, slowLoc7) &
      BIND(C, NAME='fteik_prefetchSweep6Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      !REAL(C_DOUBLE), INTENT(OUT) :: !slowLoc(7)
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                     slowLoc5, slowLoc6, slowLoc7
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
      slowLoc1 = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc2 = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc3 = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc4 = MIN(sloc(13), sloc(14))
      slowLoc5 = MIN(sloc(15), sloc(16))
      slowLoc6 = MIN(sloc(17), sloc(18))
      slowLoc7 = sloc(19)
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
                                                       nzm1, nzm1_nxm1, slow, & !slowLoc) &
                                                       slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                                       slowLoc5, slowLoc6, slowLoc7) &
      BIND(C, NAME='fteik_prefetchSweep7Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      !REAL(C_DOUBLE), INTENT(OUT) :: !slowLoc(7)
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                     slowLoc5, slowLoc6, slowLoc7
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
      slowLoc1 = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc2 = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc3 = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc4 = MIN(sloc(13), sloc(14))
      slowLoc5 = MIN(sloc(15), sloc(16))
      slowLoc6 = MIN(sloc(17), sloc(18))
      slowLoc7 = sloc(19)
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
                                                       nzm1, nzm1_nxm1, slow, & !slowLoc) &
                                                       slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                                       slowLoc5, slowLoc6, slowLoc7) &
      BIND(C, NAME='fteik_prefetchSweep8Slowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1
      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) 
      !REAL(C_DOUBLE), INTENT(OUT) :: !slowLoc(7)
      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc1, slowLoc2, slowLoc3, slowLoc4, &
                                     slowLoc5, slowLoc6, slowLoc7
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
      slowLoc1 = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowLoc2 = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowLoc3 = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowLoc4 = MIN(sloc(13), sloc(14))
      slowLoc5 = MIN(sloc(15), sloc(16))
      slowLoc6 = MIN(sloc(17), sloc(18))
      slowLoc7 = sloc(19)
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep1LS64fF(linitk, ttimes, ierr) &
      BIND(C, NAME='fteik_evaluateSweep1LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_LOCALSOLVER64F, ONLY : fteik_localSolver_noInit64fF, &
                                       fteik_localSolver_init64fF
      !USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF, fteik_localSolverExplicit64fF, &
      !                           fteik_localSolverNoInit64fF, fteik_localSolverInit64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep1Slowness64fF, & 
                                 fteik_prefetchSweep1TravelTimes64fF
      USE FTEIK_SOLVER64F, ONLY : levelPtr, lupd1, lupdInit1, ijkv1, nLevels
      !USE FTEIK_UTILS64F, ONLY : levelPtr, lupd1, lupdInit1, ijkv1, slow, &
      !                           nLevels
      USE FTEIK_MODEL64F, ONLY : slow, dx, dy, dz, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      !                           dx, dy, dz, &
      !                           nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, chunkSize, zero
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, kndx, knode, level, loop, node, node1, node2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:), &
                                     t1(:), t2(:), t3(:), t4(:), &
                                     t5(:), t6(:), t7(:), t8(:), &
                                     s1(:), s2(:), s3(:), s4(:), &
                                     s5(:), s6(:), s7(:)
      INTEGER(C_INT), PARAMETER :: sweep = 1
      INTEGER(C_INT), PARAMETER :: sgntz = 1
      INTEGER(C_INT), PARAMETER :: sgntx = 1
      INTEGER(C_INT), PARAMETER :: sgnty = 1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      !DIR$ ATTRIBUTES ALIGN:64 :: t1
      !DIR$ ATTRIBUTES ALIGN:64 :: t2
      !DIR$ ATTRIBUTES ALIGN:64 :: t3
      !DIR$ ATTRIBUTES ALIGN:64 :: t4
      !DIR$ ATTRIBUTES ALIGN:64 :: t5
      !DIR$ ATTRIBUTES ALIGN:64 :: t6
      !DIR$ ATTRIBUTES ALIGN:64 :: t7
      !DIR$ ATTRIBUTES ALIGN:64 :: t8
      !DIR$ ATTRIBUTES ALIGN:64 :: s1
      !DIR$ ATTRIBUTES ALIGN:64 :: s2
      !DIR$ ATTRIBUTES ALIGN:64 :: s3
      !DIR$ ATTRIBUTES ALIGN:64 :: s4
      !DIR$ ATTRIBUTES ALIGN:64 :: s5
      !DIR$ ATTRIBUTES ALIGN:64 :: s6
      !DIR$ ATTRIBUTES ALIGN:64 :: s7
      !DIR$ ATTRIBUTES ALIGN:64 :: tt1
      ierr = 0
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(dx, dy, dz, ijkv1, levelPtr, linitk, lupd1, lupdInit1) &
      !$OMP SHARED(nLevels, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
      !$OMP SHARED(slow, ttimes) & 
      !$OMP PRIVATE(sgnrx_dxi, sgnry_dyi, sgnrz_dzi) &
      !$OMP PRIVATE(i, level, indx, j, k, l1, l2, kndx, knode, loop, node, node1, node2, tt1) &
      !$OMP PRIVATE(t1, t2, t3, t4, t5, t6, t7, t8, ttWork, slowWork) &
      !$OMP PRIVATE(s1, s2, s3, s4, s5, s6, s7)
      ! Initialize
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      ALLOCATE(t1(chunkSize))
      ALLOCATE(t2(chunkSize))
      ALLOCATE(t3(chunkSize))
      ALLOCATE(t4(chunkSize))
      ALLOCATE(t5(chunkSize))
      ALLOCATE(t6(chunkSize))
      ALLOCATE(t7(chunkSize))
      ALLOCATE(t8(chunkSize))
      ALLOCATE(s1(chunkSize))
      ALLOCATE(s2(chunkSize))
      ALLOCATE(s3(chunkSize))
      ALLOCATE(s4(chunkSize))
      ALLOCATE(s5(chunkSize))
      ALLOCATE(s6(chunkSize))
      ALLOCATE(s7(chunkSize))
      slowWork(:) = zero
      ttWork(:) = zero
      tt1(:) = FTEIK_HUGE
      t1(:) = zero
      t2(:) = zero
      t3(:) = zero
      t4(:) = zero
      t5(:) = zero
      t6(:) = zero
      t7(:) = zero
      t8(:) = zero
      s1(:) = zero
      s2(:) = zero
      s3(:) = zero
      s4(:) = zero
      s5(:) = zero
      s6(:) = zero
      s7(:) = zero
      sgnrz_dzi = sgnrz/dz
      sgnrx_dxi = sgnrx/dx
      sgnry_dyi = sgnry/dy
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd1(node)) THEN
                     i    = ijkv1(4*(node-1)+1)
                     j    = ijkv1(4*(node-1)+2)
                     k    = ijkv1(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep1Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, & !slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd1(node)) THEN
                     i    = ijkv1(4*(node-1)+1)
                     j    = ijkv1(4*(node-1)+2)
                     k    = ijkv1(4*(node-1)+3)
                     kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep1TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                             ttimes, &
                                             t1(loop), t2(loop), t3(loop), t4(loop), &
                                             t5(loop), t6(loop), t7(loop), tt1(loop))
                                             !     ttWork(kndx+1), ttWork(kndx+2), &
                                             !     ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                             !     ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               loop = node2 - node1 + 1
               !DIR$ FORCEINLINE
               CALL fteik_localSolver_noInit64fF(loop,                           &
                                                t1, t2, t3, t4, t5, t6, t7, & !t8, &
                                                s1, s2, s3, s4, s5, s6, s7, tt1)
               !!DIR$ IVDEP !$OMP SIMD !DIR$ IVDEP
               !DO node=node1,node2 !loop=1,n2
               !   loop = node - node1 + 1 !node = knode + loop - 1
               !   IF (lupd1(node)) THEN
               !      !i    = ijkv1(4*(node-1)+1)
               !      !j    = ijkv1(4*(node-1)+2)
               !      !k    = ijkv1(4*(node-1)+3)
               !      kndx = 8*(loop - 1)
               !      !DIR$ FORCEINLINE
               !      tt1(loop) = fteik_localSolver_noInit64fF( &
               !                   t1(loop), t2(loop), t3(loop), t4(loop), &
               !                   t5(loop), t6(loop), t7(loop), t8(loop), &
               !                   s1(loop), s2(loop), s3(loop), s4(loop), &
               !                   s5(loop), s6(loop), s7(loop))
               !                   !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !                   !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !                   !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !                   !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7))
               !      !!DIR$ FORCEINLINE
               !      !tt1(loop) = fteik_localSolverExplicit64fF( &
               !      !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !      !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !      !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !      !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
               !      !             .FALSE.,                        &
               !      !              i, j, k,                       &
               !      !              sgntz, sgntx, sgnty,           &
               !      !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !      !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
               !      !                             slowWork(8*(loop-1)+1),  &
               !      !                             .FALSE.,                        &
               !      !                             i, j, k,                       &
               !      !                             sgntz, sgntx, sgnty,           &
               !      !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !   ENDIF
               !ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd1(node)) THEN
                     indx = ijkv1(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 1       CONTINUE ! loop on levels
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize !knode=l1,l2!,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               !n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit1(node)) THEN
                     i    = ijkv1(4*(node-1)+1)
                     j    = ijkv1(4*(node-1)+2)
                     k    = ijkv1(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep1Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, &!slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit1(node)) THEN
                     i    = ijkv1(4*(node-1)+1)
                     j    = ijkv1(4*(node-1)+2)
                     k    = ijkv1(4*(node-1)+3)
                     !kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep1TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                   ttimes, &
                                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                                   t5(loop), t6(loop), t7(loop), t8(loop))
                                                   !ttWork(kndx+1), ttWork(kndx+2), &
                                                   !ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                                   !ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               !!DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit1(node)) THEN
                     i    = ijkv1(4*(node-1)+1)
                     j    = ijkv1(4*(node-1)+2)
                     k    = ijkv1(4*(node-1)+3)
                     !indx = ijkv1(4*(node-1)+4)
                     !kndx = 8*(loop - 1)
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !              t1(loop), t2(loop), t3(loop), t4(loop), &
                     !              t5(loop), t6(loop), t7(loop), t8(loop), &
                     !              s1(loop), s2(loop), s3(loop), s4(loop), &
                     !              s5(loop), s6(loop), s7(loop), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver_init64fF( &
                                  !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                                  !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                                  !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                                  !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                   t5(loop), t6(loop), t7(loop), t8(loop), &
                                   s1(loop), s2(loop), s3(loop), s4(loop), &
                                   s5(loop), s6(loop), s7(loop), &
                                   i, j, k,                       &
                                   sgntz, sgntx, sgnty,           &
                                   sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !!DIR$ FORCEINLINE
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                     !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                     !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                     !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                     !                             slowWork(8*(loop-1)+1),  &
                     !                             .TRUE.,                        &
                     !                             i, j, k,                       &
                     !                             sgntz, sgntx, sgnty,           &
                     !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit1(node)) THEN
                     indx = ijkv1(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 11      CONTINUE ! loop on levels
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      DEALLOCATE(t1)
      DEALLOCATE(t2)
      DEALLOCATE(t3)
      DEALLOCATE(t4)
      DEALLOCATE(t5)
      DEALLOCATE(t6)
      DEALLOCATE(t7)
      DEALLOCATE(t8)
      DEALLOCATE(s1)
      DEALLOCATE(s2)
      DEALLOCATE(s3)
      DEALLOCATE(s4)
      DEALLOCATE(s5)
      DEALLOCATE(s6)
      DEALLOCATE(s7)
      !$OMP END PARALLEL
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep2LS64fF(linitk, ttimes, ierr) &
      BIND(C, NAME='fteik_evaluateSweep2LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_LOCALSOLVER64F, ONLY : fteik_localSolver_noInit64fF, &
                                       fteik_localSolver_init64fF
      !USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF, fteik_localSolverExplicit64fF, &
      !                           fteik_localSolverNoInit64fF, fteik_localSolverInit64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep2Slowness64fF, & 
                                 fteik_prefetchSweep2TravelTimes64fF
      USE FTEIK_SOLVER64F, ONLY : levelPtr, lupd2, lupdInit2, ijkv2, nLevels
      !USE FTEIK_UTILS64F, ONLY : levelPtr, lupd2, lupdInit2, ijkv2, slow, &
      !                           nLevels
      USE FTEIK_MODEL64F, ONLY : slow, dx, dy, dz, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      !                           dx, dy, dz, &
      !                           nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, chunkSize, zero
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, kndx, knode, level, loop, node, node1, node2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:), &
                                     t1(:), t2(:), t3(:), t4(:), &
                                     t5(:), t6(:), t7(:), t8(:), &
                                     s1(:), s2(:), s3(:), s4(:), &
                                     s5(:), s6(:), s7(:)
      INTEGER(C_INT), PARAMETER :: sweep = 2
      INTEGER(C_INT), PARAMETER :: sgntz = 1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty = 1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      !DIR$ ATTRIBUTES ALIGN:64 :: t1
      !DIR$ ATTRIBUTES ALIGN:64 :: t2
      !DIR$ ATTRIBUTES ALIGN:64 :: t3
      !DIR$ ATTRIBUTES ALIGN:64 :: t4
      !DIR$ ATTRIBUTES ALIGN:64 :: t5
      !DIR$ ATTRIBUTES ALIGN:64 :: t6
      !DIR$ ATTRIBUTES ALIGN:64 :: t7
      !DIR$ ATTRIBUTES ALIGN:64 :: t8
      !DIR$ ATTRIBUTES ALIGN:64 :: s1
      !DIR$ ATTRIBUTES ALIGN:64 :: s2
      !DIR$ ATTRIBUTES ALIGN:64 :: s3
      !DIR$ ATTRIBUTES ALIGN:64 :: s4
      !DIR$ ATTRIBUTES ALIGN:64 :: s5
      !DIR$ ATTRIBUTES ALIGN:64 :: s6
      !DIR$ ATTRIBUTES ALIGN:64 :: s7
      !DIR$ ATTRIBUTES ALIGN:64 :: tt1
      ierr = 0
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(dx, dy, dz, ijkv2, levelPtr, linitk, lupd2, lupdInit2) &
      !$OMP SHARED(nLevels, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
      !$OMP SHARED(slow, ttimes) & 
      !$OMP PRIVATE(sgnrx_dxi, sgnry_dyi, sgnrz_dzi) &
      !$OMP PRIVATE(i, level, indx, j, k, l1, l2, kndx, knode, loop, node, node1, node2, tt1) &
      !$OMP PRIVATE(t1, t2, t3, t4, t5, t6, t7, t8, ttWork, slowWork) &
      !$OMP PRIVATE(s1, s2, s3, s4, s5, s6, s7)
      ! Initialize
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      ALLOCATE(t1(chunkSize))
      ALLOCATE(t2(chunkSize))
      ALLOCATE(t3(chunkSize))
      ALLOCATE(t4(chunkSize))
      ALLOCATE(t5(chunkSize))
      ALLOCATE(t6(chunkSize))
      ALLOCATE(t7(chunkSize))
      ALLOCATE(t8(chunkSize))
      ALLOCATE(s1(chunkSize))
      ALLOCATE(s2(chunkSize))
      ALLOCATE(s3(chunkSize))
      ALLOCATE(s4(chunkSize))
      ALLOCATE(s5(chunkSize))
      ALLOCATE(s6(chunkSize))
      ALLOCATE(s7(chunkSize))
      slowWork(:) = zero
      ttWork(:) = zero
      tt1(:) = FTEIK_HUGE
      t1(:) = zero
      t2(:) = zero
      t3(:) = zero
      t4(:) = zero
      t5(:) = zero
      t6(:) = zero
      t7(:) = zero
      t8(:) = zero
      s1(:) = zero
      s2(:) = zero
      s3(:) = zero
      s4(:) = zero
      s5(:) = zero
      s6(:) = zero
      s7(:) = zero
      sgnrz_dzi = sgnrz/dz
      sgnrx_dxi = sgnrx/dx
      sgnry_dyi = sgnry/dy
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd2(node)) THEN
                     i    = ijkv2(4*(node-1)+1)
                     j    = ijkv2(4*(node-1)+2)
                     k    = ijkv2(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep2Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, & !slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd2(node)) THEN
                     i    = ijkv2(4*(node-1)+1)
                     j    = ijkv2(4*(node-1)+2)
                     k    = ijkv2(4*(node-1)+3)
                     kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep2TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                             ttimes, &
                                             t1(loop), t2(loop), t3(loop), t4(loop), &
                                             t5(loop), t6(loop), t7(loop), tt1(loop))
                                             !     ttWork(kndx+1), ttWork(kndx+2), &
                                             !     ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                             !     ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               loop = node2 - node1 + 1
               !DIR$ FORCEINLINE
               CALL fteik_localSolver_noInit64fF(loop,                           &
                                                t1, t2, t3, t4, t5, t6, t7, & !t8, &
                                                s1, s2, s3, s4, s5, s6, s7, tt1)
               !!DIR$ IVDEP !$OMP SIMD !DIR$ IVDEP
               !DO node=node1,node2 !loop=1,n2
               !   loop = node - node1 + 1 !node = knode + loop - 1
               !   IF (lupd2(node)) THEN
               !      !i    = ijkv2(4*(node-1)+1)
               !      !j    = ijkv2(4*(node-1)+2)
               !      !k    = ijkv2(4*(node-1)+3)
               !      kndx = 8*(loop - 1)
               !      !DIR$ FORCEINLINE
               !      tt1(loop) = fteik_localSolver_noInit64fF( &
               !                   t1(loop), t2(loop), t3(loop), t4(loop), &
               !                   t5(loop), t6(loop), t7(loop), t8(loop), &
               !                   s1(loop), s2(loop), s3(loop), s4(loop), &
               !                   s5(loop), s6(loop), s7(loop))
               !                   !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !                   !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !                   !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !                   !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7))
               !      !!DIR$ FORCEINLINE
               !      !tt1(loop) = fteik_localSolverExplicit64fF( &
               !      !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !      !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !      !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !      !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
               !      !             .FALSE.,                        &
               !      !              i, j, k,                       &
               !      !              sgntz, sgntx, sgnty,           &
               !      !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !      !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
               !      !                             slowWork(8*(loop-1)+1),  &
               !      !                             .FALSE.,                        &
               !      !                             i, j, k,                       &
               !      !                             sgntz, sgntx, sgnty,           &
               !      !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !   ENDIF
               !ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd2(node)) THEN
                     indx = ijkv2(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 1       CONTINUE ! loop on levels
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize !knode=l1,l2!,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               !n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit2(node)) THEN
                     i    = ijkv2(4*(node-1)+1)
                     j    = ijkv2(4*(node-1)+2)
                     k    = ijkv2(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep2Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, &!slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit2(node)) THEN
                     i    = ijkv2(4*(node-1)+1)
                     j    = ijkv2(4*(node-1)+2)
                     k    = ijkv2(4*(node-1)+3)
                     !kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep2TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                   ttimes, &
                                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                                   t5(loop), t6(loop), t7(loop), t8(loop))
                                                   !ttWork(kndx+1), ttWork(kndx+2), &
                                                   !ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                                   !ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               !!DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit2(node)) THEN
                     i    = ijkv2(4*(node-1)+1)
                     j    = ijkv2(4*(node-1)+2)
                     k    = ijkv2(4*(node-1)+3)
                     !indx = ijkv2(4*(node-1)+4)
                     !kndx = 8*(loop - 1)
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !              t1(loop), t2(loop), t3(loop), t4(loop), &
                     !              t5(loop), t6(loop), t7(loop), t8(loop), &
                     !              s1(loop), s2(loop), s3(loop), s4(loop), &
                     !              s5(loop), s6(loop), s7(loop), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver_init64fF( &
                                  !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                                  !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                                  !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                                  !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                   t5(loop), t6(loop), t7(loop), t8(loop), &
                                   s1(loop), s2(loop), s3(loop), s4(loop), &
                                   s5(loop), s6(loop), s7(loop), &
                                   i, j, k,                       &
                                   sgntz, sgntx, sgnty,           &
                                   sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !!DIR$ FORCEINLINE
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                     !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                     !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                     !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                     !                             slowWork(8*(loop-1)+1),  &
                     !                             .TRUE.,                        &
                     !                             i, j, k,                       &
                     !                             sgntz, sgntx, sgnty,           &
                     !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit2(node)) THEN
                     indx = ijkv2(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 11      CONTINUE ! loop on levels
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      DEALLOCATE(t1)
      DEALLOCATE(t2)
      DEALLOCATE(t3)
      DEALLOCATE(t4)
      DEALLOCATE(t5)
      DEALLOCATE(t6)
      DEALLOCATE(t7)
      DEALLOCATE(t8)
      DEALLOCATE(s1)
      DEALLOCATE(s2)
      DEALLOCATE(s3)
      DEALLOCATE(s4)
      DEALLOCATE(s5)
      DEALLOCATE(s6)
      DEALLOCATE(s7)
      !$OMP END PARALLEL
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep3LS64fF(linitk, ttimes, ierr) &
      BIND(C, NAME='fteik_evaluateSweep3LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_LOCALSOLVER64F, ONLY : fteik_localSolver_noInit64fF, &
                                       fteik_localSolver_init64fF
      !USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF, fteik_localSolverExplicit64fF, &
      !                           fteik_localSolverNoInit64fF, fteik_localSolverInit64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep3Slowness64fF, & 
                                 fteik_prefetchSweep3TravelTimes64fF
      USE FTEIK_SOLVER64F, ONLY : levelPtr, lupd3, lupdInit3, ijkv3, nLevels
      !USE FTEIK_UTILS64F, ONLY : levelPtr, lupd3, lupdInit3, ijkv3, slow, &
      !                           nLevels
      USE FTEIK_MODEL64F, ONLY : slow, dx, dy, dz, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      !                           dx, dy, dz, &
      !                           nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, chunkSize, zero
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, kndx, knode, level, loop, node, node1, node2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:), &
                                     t1(:), t2(:), t3(:), t4(:), &
                                     t5(:), t6(:), t7(:), t8(:), &
                                     s1(:), s2(:), s3(:), s4(:), &
                                     s5(:), s6(:), s7(:)
      INTEGER(C_INT), PARAMETER :: sweep = 3
      INTEGER(C_INT), PARAMETER :: sgntz = 1
      INTEGER(C_INT), PARAMETER :: sgntx = 1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      !DIR$ ATTRIBUTES ALIGN:64 :: t1
      !DIR$ ATTRIBUTES ALIGN:64 :: t2
      !DIR$ ATTRIBUTES ALIGN:64 :: t3
      !DIR$ ATTRIBUTES ALIGN:64 :: t4
      !DIR$ ATTRIBUTES ALIGN:64 :: t5
      !DIR$ ATTRIBUTES ALIGN:64 :: t6
      !DIR$ ATTRIBUTES ALIGN:64 :: t7
      !DIR$ ATTRIBUTES ALIGN:64 :: t8
      !DIR$ ATTRIBUTES ALIGN:64 :: s1
      !DIR$ ATTRIBUTES ALIGN:64 :: s2
      !DIR$ ATTRIBUTES ALIGN:64 :: s3
      !DIR$ ATTRIBUTES ALIGN:64 :: s4
      !DIR$ ATTRIBUTES ALIGN:64 :: s5
      !DIR$ ATTRIBUTES ALIGN:64 :: s6
      !DIR$ ATTRIBUTES ALIGN:64 :: s7
      !DIR$ ATTRIBUTES ALIGN:64 :: tt1
      ierr = 0
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(dx, dy, dz, ijkv3, levelPtr, linitk, lupd3, lupdInit3) &
      !$OMP SHARED(nLevels, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
      !$OMP SHARED(slow, ttimes) & 
      !$OMP PRIVATE(sgnrx_dxi, sgnry_dyi, sgnrz_dzi) &
      !$OMP PRIVATE(i, level, indx, j, k, l1, l2, kndx, knode, loop, node, node1, node2, tt1) &
      !$OMP PRIVATE(t1, t2, t3, t4, t5, t6, t7, t8, ttWork, slowWork) &
      !$OMP PRIVATE(s1, s2, s3, s4, s5, s6, s7)
      ! Initialize
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      ALLOCATE(t1(chunkSize))
      ALLOCATE(t2(chunkSize))
      ALLOCATE(t3(chunkSize))
      ALLOCATE(t4(chunkSize))
      ALLOCATE(t5(chunkSize))
      ALLOCATE(t6(chunkSize))
      ALLOCATE(t7(chunkSize))
      ALLOCATE(t8(chunkSize))
      ALLOCATE(s1(chunkSize))
      ALLOCATE(s2(chunkSize))
      ALLOCATE(s3(chunkSize))
      ALLOCATE(s4(chunkSize))
      ALLOCATE(s5(chunkSize))
      ALLOCATE(s6(chunkSize))
      ALLOCATE(s7(chunkSize))
      slowWork(:) = zero
      ttWork(:) = zero
      tt1(:) = FTEIK_HUGE
      t1(:) = zero
      t2(:) = zero
      t3(:) = zero
      t4(:) = zero
      t5(:) = zero
      t6(:) = zero
      t7(:) = zero
      t8(:) = zero
      s1(:) = zero
      s2(:) = zero
      s3(:) = zero
      s4(:) = zero
      s5(:) = zero
      s6(:) = zero
      s7(:) = zero
      sgnrz_dzi = sgnrz/dz
      sgnrx_dxi = sgnrx/dx
      sgnry_dyi = sgnry/dy
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd3(node)) THEN
                     i    = ijkv3(4*(node-1)+1)
                     j    = ijkv3(4*(node-1)+2)
                     k    = ijkv3(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep3Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, & !slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd3(node)) THEN
                     i    = ijkv3(4*(node-1)+1)
                     j    = ijkv3(4*(node-1)+2)
                     k    = ijkv3(4*(node-1)+3)
                     kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep3TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                             ttimes, &
                                             t1(loop), t2(loop), t3(loop), t4(loop), &
                                             t5(loop), t6(loop), t7(loop), tt1(loop))
                                             !     ttWork(kndx+1), ttWork(kndx+2), &
                                             !     ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                             !     ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               loop = node2 - node1 + 1
               !DIR$ FORCEINLINE
               CALL fteik_localSolver_noInit64fF(loop,                           &
                                                t1, t2, t3, t4, t5, t6, t7, & !t8, &
                                                s1, s2, s3, s4, s5, s6, s7, tt1)
               !!DIR$ IVDEP !$OMP SIMD !DIR$ IVDEP
               !DO node=node1,node2 !loop=1,n2
               !   loop = node - node1 + 1 !node = knode + loop - 1
               !   IF (lupd3(node)) THEN
               !      !i    = ijkv3(4*(node-1)+1)
               !      !j    = ijkv3(4*(node-1)+2)
               !      !k    = ijkv3(4*(node-1)+3)
               !      kndx = 8*(loop - 1)
               !      !DIR$ FORCEINLINE
               !      tt1(loop) = fteik_localSolver_noInit64fF( &
               !                   t1(loop), t2(loop), t3(loop), t4(loop), &
               !                   t5(loop), t6(loop), t7(loop), t8(loop), &
               !                   s1(loop), s2(loop), s3(loop), s4(loop), &
               !                   s5(loop), s6(loop), s7(loop))
               !                   !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !                   !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !                   !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !                   !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7))
               !      !!DIR$ FORCEINLINE
               !      !tt1(loop) = fteik_localSolverExplicit64fF( &
               !      !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !      !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !      !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !      !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
               !      !             .FALSE.,                        &
               !      !              i, j, k,                       &
               !      !              sgntz, sgntx, sgnty,           &
               !      !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !      !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
               !      !                             slowWork(8*(loop-1)+1),  &
               !      !                             .FALSE.,                        &
               !      !                             i, j, k,                       &
               !      !                             sgntz, sgntx, sgnty,           &
               !      !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !   ENDIF
               !ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd3(node)) THEN
                     indx = ijkv3(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 1       CONTINUE ! loop on levels
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize !knode=l1,l2!,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               !n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit3(node)) THEN
                     i    = ijkv3(4*(node-1)+1)
                     j    = ijkv3(4*(node-1)+2)
                     k    = ijkv3(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep3Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, &!slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit3(node)) THEN
                     i    = ijkv3(4*(node-1)+1)
                     j    = ijkv3(4*(node-1)+2)
                     k    = ijkv3(4*(node-1)+3)
                     !kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep3TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                   ttimes, &
                                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                                   t5(loop), t6(loop), t7(loop), t8(loop))
                                                   !ttWork(kndx+1), ttWork(kndx+2), &
                                                   !ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                                   !ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               !!DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit3(node)) THEN
                     i    = ijkv3(4*(node-1)+1)
                     j    = ijkv3(4*(node-1)+2)
                     k    = ijkv3(4*(node-1)+3)
                     !indx = ijkv3(4*(node-1)+4)
                     !kndx = 8*(loop - 1)
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !              t1(loop), t2(loop), t3(loop), t4(loop), &
                     !              t5(loop), t6(loop), t7(loop), t8(loop), &
                     !              s1(loop), s2(loop), s3(loop), s4(loop), &
                     !              s5(loop), s6(loop), s7(loop), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver_init64fF( &
                                  !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                                  !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                                  !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                                  !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                   t5(loop), t6(loop), t7(loop), t8(loop), &
                                   s1(loop), s2(loop), s3(loop), s4(loop), &
                                   s5(loop), s6(loop), s7(loop), &
                                   i, j, k,                       &
                                   sgntz, sgntx, sgnty,           &
                                   sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !!DIR$ FORCEINLINE
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                     !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                     !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                     !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                     !                             slowWork(8*(loop-1)+1),  &
                     !                             .TRUE.,                        &
                     !                             i, j, k,                       &
                     !                             sgntz, sgntx, sgnty,           &
                     !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit3(node)) THEN
                     indx = ijkv3(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 11      CONTINUE ! loop on levels
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      DEALLOCATE(t1)
      DEALLOCATE(t2)
      DEALLOCATE(t3)
      DEALLOCATE(t4)
      DEALLOCATE(t5)
      DEALLOCATE(t6)
      DEALLOCATE(t7)
      DEALLOCATE(t8)
      DEALLOCATE(s1)
      DEALLOCATE(s2)
      DEALLOCATE(s3)
      DEALLOCATE(s4)
      DEALLOCATE(s5)
      DEALLOCATE(s6)
      DEALLOCATE(s7)
      !$OMP END PARALLEL
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep4LS64fF(linitk, ttimes, ierr) &
      BIND(C, NAME='fteik_evaluateSweep4LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_LOCALSOLVER64F, ONLY : fteik_localSolver_noInit64fF, &
                                       fteik_localSolver_init64fF
      !USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF, fteik_localSolverExplicit64fF, &
      !                           fteik_localSolverNoInit64fF, fteik_localSolverInit64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep4Slowness64fF, & 
                                 fteik_prefetchSweep4TravelTimes64fF
      USE FTEIK_SOLVER64F, ONLY : levelPtr, lupd4, lupdInit4, ijkv4, nLevels
      !USE FTEIK_UTILS64F, ONLY : levelPtr, lupd4, lupdInit4, ijkv4, slow, &
      !                           nLevels
      USE FTEIK_MODEL64F, ONLY : slow, dx, dy, dz, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      !                           dx, dy, dz, &
      !                           nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, chunkSize, zero
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, kndx, knode, level, loop, node, node1, node2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:), &
                                     t1(:), t2(:), t3(:), t4(:), &
                                     t5(:), t6(:), t7(:), t8(:), &
                                     s1(:), s2(:), s3(:), s4(:), &
                                     s5(:), s6(:), s7(:)
      INTEGER(C_INT), PARAMETER :: sweep = 4
      INTEGER(C_INT), PARAMETER :: sgntz = 1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      !DIR$ ATTRIBUTES ALIGN:64 :: t1
      !DIR$ ATTRIBUTES ALIGN:64 :: t2
      !DIR$ ATTRIBUTES ALIGN:64 :: t3
      !DIR$ ATTRIBUTES ALIGN:64 :: t4
      !DIR$ ATTRIBUTES ALIGN:64 :: t5
      !DIR$ ATTRIBUTES ALIGN:64 :: t6
      !DIR$ ATTRIBUTES ALIGN:64 :: t7
      !DIR$ ATTRIBUTES ALIGN:64 :: t8
      !DIR$ ATTRIBUTES ALIGN:64 :: s1
      !DIR$ ATTRIBUTES ALIGN:64 :: s2
      !DIR$ ATTRIBUTES ALIGN:64 :: s3
      !DIR$ ATTRIBUTES ALIGN:64 :: s4
      !DIR$ ATTRIBUTES ALIGN:64 :: s5
      !DIR$ ATTRIBUTES ALIGN:64 :: s6
      !DIR$ ATTRIBUTES ALIGN:64 :: s7
      !DIR$ ATTRIBUTES ALIGN:64 :: tt1
      ierr = 0
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(dx, dy, dz, ijkv4, levelPtr, linitk, lupd4, lupdInit4) &
      !$OMP SHARED(nLevels, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
      !$OMP SHARED(slow, ttimes) & 
      !$OMP PRIVATE(sgnrx_dxi, sgnry_dyi, sgnrz_dzi) &
      !$OMP PRIVATE(i, level, indx, j, k, l1, l2, kndx, knode, loop, node, node1, node2, tt1) &
      !$OMP PRIVATE(t1, t2, t3, t4, t5, t6, t7, t8, ttWork, slowWork) &
      !$OMP PRIVATE(s1, s2, s3, s4, s5, s6, s7)
      ! Initialize
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      ALLOCATE(t1(chunkSize))
      ALLOCATE(t2(chunkSize))
      ALLOCATE(t3(chunkSize))
      ALLOCATE(t4(chunkSize))
      ALLOCATE(t5(chunkSize))
      ALLOCATE(t6(chunkSize))
      ALLOCATE(t7(chunkSize))
      ALLOCATE(t8(chunkSize))
      ALLOCATE(s1(chunkSize))
      ALLOCATE(s2(chunkSize))
      ALLOCATE(s3(chunkSize))
      ALLOCATE(s4(chunkSize))
      ALLOCATE(s5(chunkSize))
      ALLOCATE(s6(chunkSize))
      ALLOCATE(s7(chunkSize))
      slowWork(:) = zero
      ttWork(:) = zero
      tt1(:) = FTEIK_HUGE
      t1(:) = zero
      t2(:) = zero
      t3(:) = zero
      t4(:) = zero
      t5(:) = zero
      t6(:) = zero
      t7(:) = zero
      t8(:) = zero
      s1(:) = zero
      s2(:) = zero
      s3(:) = zero
      s4(:) = zero
      s5(:) = zero
      s6(:) = zero
      s7(:) = zero
      sgnrz_dzi = sgnrz/dz
      sgnrx_dxi = sgnrx/dx
      sgnry_dyi = sgnry/dy
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd4(node)) THEN
                     i    = ijkv4(4*(node-1)+1)
                     j    = ijkv4(4*(node-1)+2)
                     k    = ijkv4(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep4Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, & !slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd4(node)) THEN
                     i    = ijkv4(4*(node-1)+1)
                     j    = ijkv4(4*(node-1)+2)
                     k    = ijkv4(4*(node-1)+3)
                     kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep4TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                             ttimes, &
                                             t1(loop), t2(loop), t3(loop), t4(loop), &
                                             t5(loop), t6(loop), t7(loop), tt1(loop))
                                             !     ttWork(kndx+1), ttWork(kndx+2), &
                                             !     ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                             !     ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               loop = node2 - node1 + 1
               !DIR$ FORCEINLINE
               CALL fteik_localSolver_noInit64fF(loop,                           &
                                                t1, t2, t3, t4, t5, t6, t7, & !t8, &
                                                s1, s2, s3, s4, s5, s6, s7, tt1)
               !!DIR$ IVDEP !$OMP SIMD !DIR$ IVDEP
               !DO node=node1,node2 !loop=1,n2
               !   loop = node - node1 + 1 !node = knode + loop - 1
               !   IF (lupd4(node)) THEN
               !      !i    = ijkv4(4*(node-1)+1)
               !      !j    = ijkv4(4*(node-1)+2)
               !      !k    = ijkv4(4*(node-1)+3)
               !      kndx = 8*(loop - 1)
               !      !DIR$ FORCEINLINE
               !      tt1(loop) = fteik_localSolver_noInit64fF( &
               !                   t1(loop), t2(loop), t3(loop), t4(loop), &
               !                   t5(loop), t6(loop), t7(loop), t8(loop), &
               !                   s1(loop), s2(loop), s3(loop), s4(loop), &
               !                   s5(loop), s6(loop), s7(loop))
               !                   !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !                   !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !                   !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !                   !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7))
               !      !!DIR$ FORCEINLINE
               !      !tt1(loop) = fteik_localSolverExplicit64fF( &
               !      !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !      !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !      !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !      !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
               !      !             .FALSE.,                        &
               !      !              i, j, k,                       &
               !      !              sgntz, sgntx, sgnty,           &
               !      !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !      !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
               !      !                             slowWork(8*(loop-1)+1),  &
               !      !                             .FALSE.,                        &
               !      !                             i, j, k,                       &
               !      !                             sgntz, sgntx, sgnty,           &
               !      !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !   ENDIF
               !ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd4(node)) THEN
                     indx = ijkv4(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 1       CONTINUE ! loop on levels
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize !knode=l1,l2!,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               !n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit4(node)) THEN
                     i    = ijkv4(4*(node-1)+1)
                     j    = ijkv4(4*(node-1)+2)
                     k    = ijkv4(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep4Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, &!slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit4(node)) THEN
                     i    = ijkv4(4*(node-1)+1)
                     j    = ijkv4(4*(node-1)+2)
                     k    = ijkv4(4*(node-1)+3)
                     !kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep4TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                   ttimes, &
                                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                                   t5(loop), t6(loop), t7(loop), t8(loop))
                                                   !ttWork(kndx+1), ttWork(kndx+2), &
                                                   !ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                                   !ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               !!DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit4(node)) THEN
                     i    = ijkv4(4*(node-1)+1)
                     j    = ijkv4(4*(node-1)+2)
                     k    = ijkv4(4*(node-1)+3)
                     !indx = ijkv4(4*(node-1)+4)
                     !kndx = 8*(loop - 1)
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !              t1(loop), t2(loop), t3(loop), t4(loop), &
                     !              t5(loop), t6(loop), t7(loop), t8(loop), &
                     !              s1(loop), s2(loop), s3(loop), s4(loop), &
                     !              s5(loop), s6(loop), s7(loop), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver_init64fF( &
                                  !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                                  !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                                  !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                                  !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                   t5(loop), t6(loop), t7(loop), t8(loop), &
                                   s1(loop), s2(loop), s3(loop), s4(loop), &
                                   s5(loop), s6(loop), s7(loop), &
                                   i, j, k,                       &
                                   sgntz, sgntx, sgnty,           &
                                   sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !!DIR$ FORCEINLINE
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                     !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                     !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                     !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                     !                             slowWork(8*(loop-1)+1),  &
                     !                             .TRUE.,                        &
                     !                             i, j, k,                       &
                     !                             sgntz, sgntx, sgnty,           &
                     !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit4(node)) THEN
                     indx = ijkv4(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 11      CONTINUE ! loop on levels
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      DEALLOCATE(t1)
      DEALLOCATE(t2)
      DEALLOCATE(t3)
      DEALLOCATE(t4)
      DEALLOCATE(t5)
      DEALLOCATE(t6)
      DEALLOCATE(t7)
      DEALLOCATE(t8)
      DEALLOCATE(s1)
      DEALLOCATE(s2)
      DEALLOCATE(s3)
      DEALLOCATE(s4)
      DEALLOCATE(s5)
      DEALLOCATE(s6)
      DEALLOCATE(s7)
      !$OMP END PARALLEL
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep5LS64fF(linitk, ttimes, ierr) &
      BIND(C, NAME='fteik_evaluateSweep5LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_LOCALSOLVER64F, ONLY : fteik_localSolver_noInit64fF, &
                                       fteik_localSolver_init64fF
      !USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF, fteik_localSolverExplicit64fF, &
      !                           fteik_localSolverNoInit64fF, fteik_localSolverInit64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep5Slowness64fF, & 
                                 fteik_prefetchSweep5TravelTimes64fF
      USE FTEIK_SOLVER64F, ONLY : levelPtr, lupd5, lupdInit5, ijkv5, nLevels
      !USE FTEIK_UTILS64F, ONLY : levelPtr, lupd5, lupdInit5, ijkv5, slow, &
      !                           nLevels
      USE FTEIK_MODEL64F, ONLY : slow, dx, dy, dz, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      !                           dx, dy, dz, &
      !                           nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, chunkSize, zero
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, kndx, knode, level, loop, node, node1, node2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:), &
                                     t1(:), t2(:), t3(:), t4(:), &
                                     t5(:), t6(:), t7(:), t8(:), &
                                     s1(:), s2(:), s3(:), s4(:), &
                                     s5(:), s6(:), s7(:)
      INTEGER(C_INT), PARAMETER :: sweep = 5
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx = 1
      INTEGER(C_INT), PARAMETER :: sgnty = 1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      !DIR$ ATTRIBUTES ALIGN:64 :: t1
      !DIR$ ATTRIBUTES ALIGN:64 :: t2
      !DIR$ ATTRIBUTES ALIGN:64 :: t3
      !DIR$ ATTRIBUTES ALIGN:64 :: t4
      !DIR$ ATTRIBUTES ALIGN:64 :: t5
      !DIR$ ATTRIBUTES ALIGN:64 :: t6
      !DIR$ ATTRIBUTES ALIGN:64 :: t7
      !DIR$ ATTRIBUTES ALIGN:64 :: t8
      !DIR$ ATTRIBUTES ALIGN:64 :: s1
      !DIR$ ATTRIBUTES ALIGN:64 :: s2
      !DIR$ ATTRIBUTES ALIGN:64 :: s3
      !DIR$ ATTRIBUTES ALIGN:64 :: s4
      !DIR$ ATTRIBUTES ALIGN:64 :: s5
      !DIR$ ATTRIBUTES ALIGN:64 :: s6
      !DIR$ ATTRIBUTES ALIGN:64 :: s7
      !DIR$ ATTRIBUTES ALIGN:64 :: tt1
      ierr = 0
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(dx, dy, dz, ijkv5, levelPtr, linitk, lupd5, lupdInit5) &
      !$OMP SHARED(nLevels, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
      !$OMP SHARED(slow, ttimes) & 
      !$OMP PRIVATE(sgnrx_dxi, sgnry_dyi, sgnrz_dzi) &
      !$OMP PRIVATE(i, level, indx, j, k, l1, l2, kndx, knode, loop, node, node1, node2, tt1) &
      !$OMP PRIVATE(t1, t2, t3, t4, t5, t6, t7, t8, ttWork, slowWork) &
      !$OMP PRIVATE(s1, s2, s3, s4, s5, s6, s7)
      ! Initialize
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      ALLOCATE(t1(chunkSize))
      ALLOCATE(t2(chunkSize))
      ALLOCATE(t3(chunkSize))
      ALLOCATE(t4(chunkSize))
      ALLOCATE(t5(chunkSize))
      ALLOCATE(t6(chunkSize))
      ALLOCATE(t7(chunkSize))
      ALLOCATE(t8(chunkSize))
      ALLOCATE(s1(chunkSize))
      ALLOCATE(s2(chunkSize))
      ALLOCATE(s3(chunkSize))
      ALLOCATE(s4(chunkSize))
      ALLOCATE(s5(chunkSize))
      ALLOCATE(s6(chunkSize))
      ALLOCATE(s7(chunkSize))
      slowWork(:) = zero
      ttWork(:) = zero
      tt1(:) = FTEIK_HUGE
      t1(:) = zero
      t2(:) = zero
      t3(:) = zero
      t4(:) = zero
      t5(:) = zero
      t6(:) = zero
      t7(:) = zero
      t8(:) = zero
      s1(:) = zero
      s2(:) = zero
      s3(:) = zero
      s4(:) = zero
      s5(:) = zero
      s6(:) = zero
      s7(:) = zero
      sgnrz_dzi = sgnrz/dz
      sgnrx_dxi = sgnrx/dx
      sgnry_dyi = sgnry/dy
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd5(node)) THEN
                     i    = ijkv5(4*(node-1)+1)
                     j    = ijkv5(4*(node-1)+2)
                     k    = ijkv5(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep5Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, & !slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd5(node)) THEN
                     i    = ijkv5(4*(node-1)+1)
                     j    = ijkv5(4*(node-1)+2)
                     k    = ijkv5(4*(node-1)+3)
                     kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep5TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                             ttimes, &
                                             t1(loop), t2(loop), t3(loop), t4(loop), &
                                             t5(loop), t6(loop), t7(loop), tt1(loop))
                                             !     ttWork(kndx+1), ttWork(kndx+2), &
                                             !     ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                             !     ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               loop = node2 - node1 + 1
               !DIR$ FORCEINLINE
               CALL fteik_localSolver_noInit64fF(loop,                           &
                                                t1, t2, t3, t4, t5, t6, t7, & !t8, &
                                                s1, s2, s3, s4, s5, s6, s7, tt1)
               !!DIR$ IVDEP !$OMP SIMD !DIR$ IVDEP
               !DO node=node1,node2 !loop=1,n2
               !   loop = node - node1 + 1 !node = knode + loop - 1
               !   IF (lupd5(node)) THEN
               !      !i    = ijkv5(4*(node-1)+1)
               !      !j    = ijkv5(4*(node-1)+2)
               !      !k    = ijkv5(4*(node-1)+3)
               !      kndx = 8*(loop - 1)
               !      !DIR$ FORCEINLINE
               !      tt1(loop) = fteik_localSolver_noInit64fF( &
               !                   t1(loop), t2(loop), t3(loop), t4(loop), &
               !                   t5(loop), t6(loop), t7(loop), t8(loop), &
               !                   s1(loop), s2(loop), s3(loop), s4(loop), &
               !                   s5(loop), s6(loop), s7(loop))
               !                   !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !                   !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !                   !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !                   !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7))
               !      !!DIR$ FORCEINLINE
               !      !tt1(loop) = fteik_localSolverExplicit64fF( &
               !      !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !      !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !      !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !      !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
               !      !             .FALSE.,                        &
               !      !              i, j, k,                       &
               !      !              sgntz, sgntx, sgnty,           &
               !      !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !      !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
               !      !                             slowWork(8*(loop-1)+1),  &
               !      !                             .FALSE.,                        &
               !      !                             i, j, k,                       &
               !      !                             sgntz, sgntx, sgnty,           &
               !      !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !   ENDIF
               !ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd5(node)) THEN
                     indx = ijkv5(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 1       CONTINUE ! loop on levels
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize !knode=l1,l2!,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               !n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit5(node)) THEN
                     i    = ijkv5(4*(node-1)+1)
                     j    = ijkv5(4*(node-1)+2)
                     k    = ijkv5(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep5Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, &!slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit5(node)) THEN
                     i    = ijkv5(4*(node-1)+1)
                     j    = ijkv5(4*(node-1)+2)
                     k    = ijkv5(4*(node-1)+3)
                     !kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep5TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                   ttimes, &
                                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                                   t5(loop), t6(loop), t7(loop), t8(loop))
                                                   !ttWork(kndx+1), ttWork(kndx+2), &
                                                   !ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                                   !ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               !!DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit5(node)) THEN
                     i    = ijkv5(4*(node-1)+1)
                     j    = ijkv5(4*(node-1)+2)
                     k    = ijkv5(4*(node-1)+3)
                     !indx = ijkv5(4*(node-1)+4)
                     !kndx = 8*(loop - 1)
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !              t1(loop), t2(loop), t3(loop), t4(loop), &
                     !              t5(loop), t6(loop), t7(loop), t8(loop), &
                     !              s1(loop), s2(loop), s3(loop), s4(loop), &
                     !              s5(loop), s6(loop), s7(loop), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver_init64fF( &
                                  !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                                  !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                                  !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                                  !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                   t5(loop), t6(loop), t7(loop), t8(loop), &
                                   s1(loop), s2(loop), s3(loop), s4(loop), &
                                   s5(loop), s6(loop), s7(loop), &
                                   i, j, k,                       &
                                   sgntz, sgntx, sgnty,           &
                                   sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !!DIR$ FORCEINLINE
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                     !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                     !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                     !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                     !                             slowWork(8*(loop-1)+1),  &
                     !                             .TRUE.,                        &
                     !                             i, j, k,                       &
                     !                             sgntz, sgntx, sgnty,           &
                     !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit5(node)) THEN
                     indx = ijkv5(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 11      CONTINUE ! loop on levels
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      DEALLOCATE(t1)
      DEALLOCATE(t2)
      DEALLOCATE(t3)
      DEALLOCATE(t4)
      DEALLOCATE(t5)
      DEALLOCATE(t6)
      DEALLOCATE(t7)
      DEALLOCATE(t8)
      DEALLOCATE(s1)
      DEALLOCATE(s2)
      DEALLOCATE(s3)
      DEALLOCATE(s4)
      DEALLOCATE(s5)
      DEALLOCATE(s6)
      DEALLOCATE(s7)
      !$OMP END PARALLEL
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep6LS64fF(linitk, ttimes, ierr) &
      BIND(C, NAME='fteik_evaluateSweep6LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_LOCALSOLVER64F, ONLY : fteik_localSolver_noInit64fF, &
                                       fteik_localSolver_init64fF
      !USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF, fteik_localSolverExplicit64fF, &
      !                           fteik_localSolverNoInit64fF, fteik_localSolverInit64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep6Slowness64fF, & 
                                 fteik_prefetchSweep6TravelTimes64fF
      USE FTEIK_SOLVER64F, ONLY : levelPtr, lupd6, lupdInit6, ijkv6, nLevels
      !USE FTEIK_UTILS64F, ONLY : levelPtr, lupd6, lupdInit6, ijkv6, slow, &
      !                           nLevels
      USE FTEIK_MODEL64F, ONLY : slow, dx, dy, dz, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      !                           dx, dy, dz, &
      !                           nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, chunkSize, zero
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, kndx, knode, level, loop, node, node1, node2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:), &
                                     t1(:), t2(:), t3(:), t4(:), &
                                     t5(:), t6(:), t7(:), t8(:), &
                                     s1(:), s2(:), s3(:), s4(:), &
                                     s5(:), s6(:), s7(:)
      INTEGER(C_INT), PARAMETER :: sweep = 6
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty = 1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      !DIR$ ATTRIBUTES ALIGN:64 :: t1
      !DIR$ ATTRIBUTES ALIGN:64 :: t2
      !DIR$ ATTRIBUTES ALIGN:64 :: t3
      !DIR$ ATTRIBUTES ALIGN:64 :: t4
      !DIR$ ATTRIBUTES ALIGN:64 :: t5
      !DIR$ ATTRIBUTES ALIGN:64 :: t6
      !DIR$ ATTRIBUTES ALIGN:64 :: t7
      !DIR$ ATTRIBUTES ALIGN:64 :: t8
      !DIR$ ATTRIBUTES ALIGN:64 :: s1
      !DIR$ ATTRIBUTES ALIGN:64 :: s2
      !DIR$ ATTRIBUTES ALIGN:64 :: s3
      !DIR$ ATTRIBUTES ALIGN:64 :: s4
      !DIR$ ATTRIBUTES ALIGN:64 :: s5
      !DIR$ ATTRIBUTES ALIGN:64 :: s6
      !DIR$ ATTRIBUTES ALIGN:64 :: s7
      !DIR$ ATTRIBUTES ALIGN:64 :: tt1
      ierr = 0
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(dx, dy, dz, ijkv6, levelPtr, linitk, lupd6, lupdInit6) &
      !$OMP SHARED(nLevels, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
      !$OMP SHARED(slow, ttimes) & 
      !$OMP PRIVATE(sgnrx_dxi, sgnry_dyi, sgnrz_dzi) &
      !$OMP PRIVATE(i, level, indx, j, k, l1, l2, kndx, knode, loop, node, node1, node2, tt1) &
      !$OMP PRIVATE(t1, t2, t3, t4, t5, t6, t7, t8, ttWork, slowWork) &
      !$OMP PRIVATE(s1, s2, s3, s4, s5, s6, s7)
      ! Initialize
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      ALLOCATE(t1(chunkSize))
      ALLOCATE(t2(chunkSize))
      ALLOCATE(t3(chunkSize))
      ALLOCATE(t4(chunkSize))
      ALLOCATE(t5(chunkSize))
      ALLOCATE(t6(chunkSize))
      ALLOCATE(t7(chunkSize))
      ALLOCATE(t8(chunkSize))
      ALLOCATE(s1(chunkSize))
      ALLOCATE(s2(chunkSize))
      ALLOCATE(s3(chunkSize))
      ALLOCATE(s4(chunkSize))
      ALLOCATE(s5(chunkSize))
      ALLOCATE(s6(chunkSize))
      ALLOCATE(s7(chunkSize))
      slowWork(:) = zero
      ttWork(:) = zero
      tt1(:) = FTEIK_HUGE
      t1(:) = zero
      t2(:) = zero
      t3(:) = zero
      t4(:) = zero
      t5(:) = zero
      t6(:) = zero
      t7(:) = zero
      t8(:) = zero
      s1(:) = zero
      s2(:) = zero
      s3(:) = zero
      s4(:) = zero
      s5(:) = zero
      s6(:) = zero
      s7(:) = zero
      sgnrz_dzi = sgnrz/dz
      sgnrx_dxi = sgnrx/dx
      sgnry_dyi = sgnry/dy
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd6(node)) THEN
                     i    = ijkv6(4*(node-1)+1)
                     j    = ijkv6(4*(node-1)+2)
                     k    = ijkv6(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep6Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, & !slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd6(node)) THEN
                     i    = ijkv6(4*(node-1)+1)
                     j    = ijkv6(4*(node-1)+2)
                     k    = ijkv6(4*(node-1)+3)
                     kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep6TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                             ttimes, &
                                             t1(loop), t2(loop), t3(loop), t4(loop), &
                                             t5(loop), t6(loop), t7(loop), tt1(loop))
                                             !     ttWork(kndx+1), ttWork(kndx+2), &
                                             !     ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                             !     ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               loop = node2 - node1 + 1
               !DIR$ FORCEINLINE
               CALL fteik_localSolver_noInit64fF(loop,                           &
                                                t1, t2, t3, t4, t5, t6, t7, & !t8, &
                                                s1, s2, s3, s4, s5, s6, s7, tt1)
               !!DIR$ IVDEP !$OMP SIMD !DIR$ IVDEP
               !DO node=node1,node2 !loop=1,n2
               !   loop = node - node1 + 1 !node = knode + loop - 1
               !   IF (lupd6(node)) THEN
               !      !i    = ijkv6(4*(node-1)+1)
               !      !j    = ijkv6(4*(node-1)+2)
               !      !k    = ijkv6(4*(node-1)+3)
               !      kndx = 8*(loop - 1)
               !      !DIR$ FORCEINLINE
               !      tt1(loop) = fteik_localSolver_noInit64fF( &
               !                   t1(loop), t2(loop), t3(loop), t4(loop), &
               !                   t5(loop), t6(loop), t7(loop), t8(loop), &
               !                   s1(loop), s2(loop), s3(loop), s4(loop), &
               !                   s5(loop), s6(loop), s7(loop))
               !                   !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !                   !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !                   !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !                   !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7))
               !      !!DIR$ FORCEINLINE
               !      !tt1(loop) = fteik_localSolverExplicit64fF( &
               !      !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !      !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !      !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !      !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
               !      !             .FALSE.,                        &
               !      !              i, j, k,                       &
               !      !              sgntz, sgntx, sgnty,           &
               !      !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !      !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
               !      !                             slowWork(8*(loop-1)+1),  &
               !      !                             .FALSE.,                        &
               !      !                             i, j, k,                       &
               !      !                             sgntz, sgntx, sgnty,           &
               !      !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !   ENDIF
               !ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd6(node)) THEN
                     indx = ijkv6(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 1       CONTINUE ! loop on levels
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize !knode=l1,l2!,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               !n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit6(node)) THEN
                     i    = ijkv6(4*(node-1)+1)
                     j    = ijkv6(4*(node-1)+2)
                     k    = ijkv6(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep6Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, &!slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit6(node)) THEN
                     i    = ijkv6(4*(node-1)+1)
                     j    = ijkv6(4*(node-1)+2)
                     k    = ijkv6(4*(node-1)+3)
                     !kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep6TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                   ttimes, &
                                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                                   t5(loop), t6(loop), t7(loop), t8(loop))
                                                   !ttWork(kndx+1), ttWork(kndx+2), &
                                                   !ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                                   !ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               !!DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit6(node)) THEN
                     i    = ijkv6(4*(node-1)+1)
                     j    = ijkv6(4*(node-1)+2)
                     k    = ijkv6(4*(node-1)+3)
                     !indx = ijkv6(4*(node-1)+4)
                     !kndx = 8*(loop - 1)
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !              t1(loop), t2(loop), t3(loop), t4(loop), &
                     !              t5(loop), t6(loop), t7(loop), t8(loop), &
                     !              s1(loop), s2(loop), s3(loop), s4(loop), &
                     !              s5(loop), s6(loop), s7(loop), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver_init64fF( &
                                  !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                                  !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                                  !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                                  !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                   t5(loop), t6(loop), t7(loop), t8(loop), &
                                   s1(loop), s2(loop), s3(loop), s4(loop), &
                                   s5(loop), s6(loop), s7(loop), &
                                   i, j, k,                       &
                                   sgntz, sgntx, sgnty,           &
                                   sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !!DIR$ FORCEINLINE
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                     !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                     !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                     !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                     !                             slowWork(8*(loop-1)+1),  &
                     !                             .TRUE.,                        &
                     !                             i, j, k,                       &
                     !                             sgntz, sgntx, sgnty,           &
                     !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit6(node)) THEN
                     indx = ijkv6(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 11      CONTINUE ! loop on levels
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      DEALLOCATE(t1)
      DEALLOCATE(t2)
      DEALLOCATE(t3)
      DEALLOCATE(t4)
      DEALLOCATE(t5)
      DEALLOCATE(t6)
      DEALLOCATE(t7)
      DEALLOCATE(t8)
      DEALLOCATE(s1)
      DEALLOCATE(s2)
      DEALLOCATE(s3)
      DEALLOCATE(s4)
      DEALLOCATE(s5)
      DEALLOCATE(s6)
      DEALLOCATE(s7)
      !$OMP END PARALLEL
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep7LS64fF(linitk, ttimes, ierr) &
      BIND(C, NAME='fteik_evaluateSweep7LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_LOCALSOLVER64F, ONLY : fteik_localSolver_noInit64fF, &
                                       fteik_localSolver_init64fF
      !USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF, fteik_localSolverExplicit64fF, &
      !                           fteik_localSolverNoInit64fF, fteik_localSolverInit64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep7Slowness64fF, & 
                                 fteik_prefetchSweep7TravelTimes64fF
      USE FTEIK_SOLVER64F, ONLY : levelPtr, lupd7, lupdInit7, ijkv7, nLevels
      !USE FTEIK_UTILS64F, ONLY : levelPtr, lupd7, lupdInit7, ijkv7, slow, &
      !                           nLevels
      USE FTEIK_MODEL64F, ONLY : slow, dx, dy, dz, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      !                           dx, dy, dz, &
      !                           nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, chunkSize, zero
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, kndx, knode, level, loop, node, node1, node2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:), &
                                     t1(:), t2(:), t3(:), t4(:), &
                                     t5(:), t6(:), t7(:), t8(:), &
                                     s1(:), s2(:), s3(:), s4(:), &
                                     s5(:), s6(:), s7(:)
      INTEGER(C_INT), PARAMETER :: sweep = 7
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx = 1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      !DIR$ ATTRIBUTES ALIGN:64 :: t1
      !DIR$ ATTRIBUTES ALIGN:64 :: t2
      !DIR$ ATTRIBUTES ALIGN:64 :: t3
      !DIR$ ATTRIBUTES ALIGN:64 :: t4
      !DIR$ ATTRIBUTES ALIGN:64 :: t5
      !DIR$ ATTRIBUTES ALIGN:64 :: t6
      !DIR$ ATTRIBUTES ALIGN:64 :: t7
      !DIR$ ATTRIBUTES ALIGN:64 :: t8
      !DIR$ ATTRIBUTES ALIGN:64 :: s1
      !DIR$ ATTRIBUTES ALIGN:64 :: s2
      !DIR$ ATTRIBUTES ALIGN:64 :: s3
      !DIR$ ATTRIBUTES ALIGN:64 :: s4
      !DIR$ ATTRIBUTES ALIGN:64 :: s5
      !DIR$ ATTRIBUTES ALIGN:64 :: s6
      !DIR$ ATTRIBUTES ALIGN:64 :: s7
      !DIR$ ATTRIBUTES ALIGN:64 :: tt1
      ierr = 0
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(dx, dy, dz, ijkv7, levelPtr, linitk, lupd7, lupdInit7) &
      !$OMP SHARED(nLevels, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
      !$OMP SHARED(slow, ttimes) & 
      !$OMP PRIVATE(sgnrx_dxi, sgnry_dyi, sgnrz_dzi) &
      !$OMP PRIVATE(i, level, indx, j, k, l1, l2, kndx, knode, loop, node, node1, node2, tt1) &
      !$OMP PRIVATE(t1, t2, t3, t4, t5, t6, t7, t8, ttWork, slowWork) &
      !$OMP PRIVATE(s1, s2, s3, s4, s5, s6, s7)
      ! Initialize
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      ALLOCATE(t1(chunkSize))
      ALLOCATE(t2(chunkSize))
      ALLOCATE(t3(chunkSize))
      ALLOCATE(t4(chunkSize))
      ALLOCATE(t5(chunkSize))
      ALLOCATE(t6(chunkSize))
      ALLOCATE(t7(chunkSize))
      ALLOCATE(t8(chunkSize))
      ALLOCATE(s1(chunkSize))
      ALLOCATE(s2(chunkSize))
      ALLOCATE(s3(chunkSize))
      ALLOCATE(s4(chunkSize))
      ALLOCATE(s5(chunkSize))
      ALLOCATE(s6(chunkSize))
      ALLOCATE(s7(chunkSize))
      slowWork(:) = zero
      ttWork(:) = zero
      tt1(:) = FTEIK_HUGE
      t1(:) = zero
      t2(:) = zero
      t3(:) = zero
      t4(:) = zero
      t5(:) = zero
      t6(:) = zero
      t7(:) = zero
      t8(:) = zero
      s1(:) = zero
      s2(:) = zero
      s3(:) = zero
      s4(:) = zero
      s5(:) = zero
      s6(:) = zero
      s7(:) = zero
      sgnrz_dzi = sgnrz/dz
      sgnrx_dxi = sgnrx/dx
      sgnry_dyi = sgnry/dy
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd7(node)) THEN
                     i    = ijkv7(4*(node-1)+1)
                     j    = ijkv7(4*(node-1)+2)
                     k    = ijkv7(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep7Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, & !slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd7(node)) THEN
                     i    = ijkv7(4*(node-1)+1)
                     j    = ijkv7(4*(node-1)+2)
                     k    = ijkv7(4*(node-1)+3)
                     kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep7TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                             ttimes, &
                                             t1(loop), t2(loop), t3(loop), t4(loop), &
                                             t5(loop), t6(loop), t7(loop), tt1(loop))
                                             !     ttWork(kndx+1), ttWork(kndx+2), &
                                             !     ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                             !     ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               loop = node2 - node1 + 1
               !DIR$ FORCEINLINE
               CALL fteik_localSolver_noInit64fF(loop,                           &
                                                t1, t2, t3, t4, t5, t6, t7, & !t8, &
                                                s1, s2, s3, s4, s5, s6, s7, tt1)
               !!DIR$ IVDEP !$OMP SIMD !DIR$ IVDEP
               !DO node=node1,node2 !loop=1,n2
               !   loop = node - node1 + 1 !node = knode + loop - 1
               !   IF (lupd7(node)) THEN
               !      !i    = ijkv7(4*(node-1)+1)
               !      !j    = ijkv7(4*(node-1)+2)
               !      !k    = ijkv7(4*(node-1)+3)
               !      kndx = 8*(loop - 1)
               !      !DIR$ FORCEINLINE
               !      tt1(loop) = fteik_localSolver_noInit64fF( &
               !                   t1(loop), t2(loop), t3(loop), t4(loop), &
               !                   t5(loop), t6(loop), t7(loop), t8(loop), &
               !                   s1(loop), s2(loop), s3(loop), s4(loop), &
               !                   s5(loop), s6(loop), s7(loop))
               !                   !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !                   !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !                   !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !                   !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7))
               !      !!DIR$ FORCEINLINE
               !      !tt1(loop) = fteik_localSolverExplicit64fF( &
               !      !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !      !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !      !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !      !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
               !      !             .FALSE.,                        &
               !      !              i, j, k,                       &
               !      !              sgntz, sgntx, sgnty,           &
               !      !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !      !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
               !      !                             slowWork(8*(loop-1)+1),  &
               !      !                             .FALSE.,                        &
               !      !                             i, j, k,                       &
               !      !                             sgntz, sgntx, sgnty,           &
               !      !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !   ENDIF
               !ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd7(node)) THEN
                     indx = ijkv7(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 1       CONTINUE ! loop on levels
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize !knode=l1,l2!,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               !n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit7(node)) THEN
                     i    = ijkv7(4*(node-1)+1)
                     j    = ijkv7(4*(node-1)+2)
                     k    = ijkv7(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep7Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, &!slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit7(node)) THEN
                     i    = ijkv7(4*(node-1)+1)
                     j    = ijkv7(4*(node-1)+2)
                     k    = ijkv7(4*(node-1)+3)
                     !kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep7TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                   ttimes, &
                                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                                   t5(loop), t6(loop), t7(loop), t8(loop))
                                                   !ttWork(kndx+1), ttWork(kndx+2), &
                                                   !ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                                   !ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               !!DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit7(node)) THEN
                     i    = ijkv7(4*(node-1)+1)
                     j    = ijkv7(4*(node-1)+2)
                     k    = ijkv7(4*(node-1)+3)
                     !indx = ijkv7(4*(node-1)+4)
                     !kndx = 8*(loop - 1)
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !              t1(loop), t2(loop), t3(loop), t4(loop), &
                     !              t5(loop), t6(loop), t7(loop), t8(loop), &
                     !              s1(loop), s2(loop), s3(loop), s4(loop), &
                     !              s5(loop), s6(loop), s7(loop), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver_init64fF( &
                                  !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                                  !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                                  !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                                  !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                   t5(loop), t6(loop), t7(loop), t8(loop), &
                                   s1(loop), s2(loop), s3(loop), s4(loop), &
                                   s5(loop), s6(loop), s7(loop), &
                                   i, j, k,                       &
                                   sgntz, sgntx, sgnty,           &
                                   sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !!DIR$ FORCEINLINE
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                     !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                     !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                     !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                     !                             slowWork(8*(loop-1)+1),  &
                     !                             .TRUE.,                        &
                     !                             i, j, k,                       &
                     !                             sgntz, sgntx, sgnty,           &
                     !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit7(node)) THEN
                     indx = ijkv7(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 11      CONTINUE ! loop on levels
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      DEALLOCATE(t1)
      DEALLOCATE(t2)
      DEALLOCATE(t3)
      DEALLOCATE(t4)
      DEALLOCATE(t5)
      DEALLOCATE(t6)
      DEALLOCATE(t7)
      DEALLOCATE(t8)
      DEALLOCATE(s1)
      DEALLOCATE(s2)
      DEALLOCATE(s3)
      DEALLOCATE(s4)
      DEALLOCATE(s5)
      DEALLOCATE(s6)
      DEALLOCATE(s7)
      !$OMP END PARALLEL
      RETURN
      END SUBROUTINE
      SUBROUTINE fteik_evaluateSweep8LS64fF(linitk, ttimes, ierr) &
      BIND(C, NAME='fteik_evaluateSweep8LS64fF')
      USE ISO_C_BINDING
      USE FTEIK_LOCALSOLVER64F, ONLY : fteik_localSolver_noInit64fF, &
                                       fteik_localSolver_init64fF
      !USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF, fteik_localSolverExplicit64fF, &
      !                           fteik_localSolverNoInit64fF, fteik_localSolverInit64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep8Slowness64fF, & 
                                 fteik_prefetchSweep8TravelTimes64fF
      USE FTEIK_SOLVER64F, ONLY : levelPtr, lupd8, lupdInit8, ijkv8, nLevels
      !USE FTEIK_UTILS64F, ONLY : levelPtr, lupd8, lupdInit8, ijkv8, slow, &
      !                           nLevels
      USE FTEIK_MODEL64F, ONLY : slow, dx, dy, dz, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      !                           dx, dy, dz, &
      !                           nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, chunkSize, zero
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT) i, indx, j, k, l1, l2, kndx, knode, level, loop, node, node1, node2
      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:), &
                                     t1(:), t2(:), t3(:), t4(:), &
                                     t5(:), t6(:), t7(:), t8(:), &
                                     s1(:), s2(:), s3(:), s4(:), &
                                     s5(:), s6(:), s7(:)
      INTEGER(C_INT), PARAMETER :: sweep = 8
      INTEGER(C_INT), PARAMETER :: sgntz = -1
      INTEGER(C_INT), PARAMETER :: sgntx = -1
      INTEGER(C_INT), PARAMETER :: sgnty = -1
      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)
      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)
      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)
      !DIR$ ATTRIBUTES ALIGN:64 :: t1
      !DIR$ ATTRIBUTES ALIGN:64 :: t2
      !DIR$ ATTRIBUTES ALIGN:64 :: t3
      !DIR$ ATTRIBUTES ALIGN:64 :: t4
      !DIR$ ATTRIBUTES ALIGN:64 :: t5
      !DIR$ ATTRIBUTES ALIGN:64 :: t6
      !DIR$ ATTRIBUTES ALIGN:64 :: t7
      !DIR$ ATTRIBUTES ALIGN:64 :: t8
      !DIR$ ATTRIBUTES ALIGN:64 :: s1
      !DIR$ ATTRIBUTES ALIGN:64 :: s2
      !DIR$ ATTRIBUTES ALIGN:64 :: s3
      !DIR$ ATTRIBUTES ALIGN:64 :: s4
      !DIR$ ATTRIBUTES ALIGN:64 :: s5
      !DIR$ ATTRIBUTES ALIGN:64 :: s6
      !DIR$ ATTRIBUTES ALIGN:64 :: s7
      !DIR$ ATTRIBUTES ALIGN:64 :: tt1
      ierr = 0
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(dx, dy, dz, ijkv8, levelPtr, linitk, lupd8, lupdInit8) &
      !$OMP SHARED(nLevels, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &
      !$OMP SHARED(slow, ttimes) & 
      !$OMP PRIVATE(sgnrx_dxi, sgnry_dyi, sgnrz_dzi) &
      !$OMP PRIVATE(i, level, indx, j, k, l1, l2, kndx, knode, loop, node, node1, node2, tt1) &
      !$OMP PRIVATE(t1, t2, t3, t4, t5, t6, t7, t8, ttWork, slowWork) &
      !$OMP PRIVATE(s1, s2, s3, s4, s5, s6, s7)
      ! Initialize
      ALLOCATE(slowWork(8*chunkSize))
      ALLOCATE(ttWork(8*chunkSize))
      ALLOCATE(tt1(chunkSize))
      ALLOCATE(t1(chunkSize))
      ALLOCATE(t2(chunkSize))
      ALLOCATE(t3(chunkSize))
      ALLOCATE(t4(chunkSize))
      ALLOCATE(t5(chunkSize))
      ALLOCATE(t6(chunkSize))
      ALLOCATE(t7(chunkSize))
      ALLOCATE(t8(chunkSize))
      ALLOCATE(s1(chunkSize))
      ALLOCATE(s2(chunkSize))
      ALLOCATE(s3(chunkSize))
      ALLOCATE(s4(chunkSize))
      ALLOCATE(s5(chunkSize))
      ALLOCATE(s6(chunkSize))
      ALLOCATE(s7(chunkSize))
      slowWork(:) = zero
      ttWork(:) = zero
      tt1(:) = FTEIK_HUGE
      t1(:) = zero
      t2(:) = zero
      t3(:) = zero
      t4(:) = zero
      t5(:) = zero
      t6(:) = zero
      t7(:) = zero
      t8(:) = zero
      s1(:) = zero
      s2(:) = zero
      s3(:) = zero
      s4(:) = zero
      s5(:) = zero
      s6(:) = zero
      s7(:) = zero
      sgnrz_dzi = sgnrz/dz
      sgnrx_dxi = sgnrx/dx
      sgnry_dyi = sgnry/dy
      IF (.NOT.linitk) THEN
         DO 1 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd8(node)) THEN
                     i    = ijkv8(4*(node-1)+1)
                     j    = ijkv8(4*(node-1)+2)
                     k    = ijkv8(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep8Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, & !slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd8(node)) THEN
                     i    = ijkv8(4*(node-1)+1)
                     j    = ijkv8(4*(node-1)+2)
                     k    = ijkv8(4*(node-1)+3)
                     kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep8TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                             ttimes, &
                                             t1(loop), t2(loop), t3(loop), t4(loop), &
                                             t5(loop), t6(loop), t7(loop), tt1(loop))
                                             !     ttWork(kndx+1), ttWork(kndx+2), &
                                             !     ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                             !     ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               loop = node2 - node1 + 1
               !DIR$ FORCEINLINE
               CALL fteik_localSolver_noInit64fF(loop,                           &
                                                t1, t2, t3, t4, t5, t6, t7, & !t8, &
                                                s1, s2, s3, s4, s5, s6, s7, tt1)
               !!DIR$ IVDEP !$OMP SIMD !DIR$ IVDEP
               !DO node=node1,node2 !loop=1,n2
               !   loop = node - node1 + 1 !node = knode + loop - 1
               !   IF (lupd8(node)) THEN
               !      !i    = ijkv8(4*(node-1)+1)
               !      !j    = ijkv8(4*(node-1)+2)
               !      !k    = ijkv8(4*(node-1)+3)
               !      kndx = 8*(loop - 1)
               !      !DIR$ FORCEINLINE
               !      tt1(loop) = fteik_localSolver_noInit64fF( &
               !                   t1(loop), t2(loop), t3(loop), t4(loop), &
               !                   t5(loop), t6(loop), t7(loop), t8(loop), &
               !                   s1(loop), s2(loop), s3(loop), s4(loop), &
               !                   s5(loop), s6(loop), s7(loop))
               !                   !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !                   !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !                   !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !                   !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7))
               !      !!DIR$ FORCEINLINE
               !      !tt1(loop) = fteik_localSolverExplicit64fF( &
               !      !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
               !      !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
               !      !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
               !      !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
               !      !             .FALSE.,                        &
               !      !              i, j, k,                       &
               !      !              sgntz, sgntx, sgnty,           &
               !      !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !      !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
               !      !                             slowWork(8*(loop-1)+1),  &
               !      !                             .FALSE.,                        &
               !      !                             i, j, k,                       &
               !      !                             sgntz, sgntx, sgnty,           &
               !      !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               !   ENDIF
               !ENDDO
               !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupd8(node)) THEN
                     indx = ijkv8(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 1       CONTINUE ! loop on levels
      ELSE
         DO 11 level=1,nLevels
            l1 = levelPtr(level)
            l2 = levelPtr(level+1) - 1
            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)
            !n2 = MIN(l2 - l1 + 1, chunkSize)
            !$OMP DO
            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize !knode=l1,l2!,chunkSize
               node1 = knode !l1 + (knode - 1)*chunkSize
               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)
               !n2 = 1 !MIN(l2-l1+1, chunkSize) !1
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit8(node)) THEN
                     i    = ijkv8(4*(node-1)+1)
                     j    = ijkv8(4*(node-1)+2)
                     k    = ijkv8(4*(node-1)+3)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep8Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                           slow, &!slowWork(8*(loop-1)+1))
                                                           s1(loop), s2(loop), s3(loop), s4(loop), &
                                                           s5(loop), s6(loop), s7(loop))
                  ENDIF
               ENDDO
               DO node=node1,node2 !1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit8(node)) THEN
                     i    = ijkv8(4*(node-1)+1)
                     j    = ijkv8(4*(node-1)+2)
                     k    = ijkv8(4*(node-1)+3)
                     !kndx = 8*(loop-1)
                     !DIR$ FORCEINLINE
                     CALL fteik_prefetchSweep8TravelTimes64fF(i, j, k, nz, nx, ny, nzx, &
                                                   ttimes, &
                                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                                   t5(loop), t6(loop), t7(loop), t8(loop))
                                                   !ttWork(kndx+1), ttWork(kndx+2), &
                                                   !ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &
                                                   !ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))
                  ENDIF
               ENDDO
               !!DIR$ IVDEP !$OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit8(node)) THEN
                     i    = ijkv8(4*(node-1)+1)
                     j    = ijkv8(4*(node-1)+2)
                     k    = ijkv8(4*(node-1)+3)
                     !indx = ijkv8(4*(node-1)+4)
                     !kndx = 8*(loop - 1)
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !              t1(loop), t2(loop), t3(loop), t4(loop), &
                     !              t5(loop), t6(loop), t7(loop), t8(loop), &
                     !              s1(loop), s2(loop), s3(loop), s4(loop), &
                     !              s5(loop), s6(loop), s7(loop), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !DIR$ FORCEINLINE
                     tt1(loop) = fteik_localSolver_init64fF( &
                                  !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                                  !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                                  !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                                  !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                                   t1(loop), t2(loop), t3(loop), t4(loop), &
                                   t5(loop), t6(loop), t7(loop), t8(loop), &
                                   s1(loop), s2(loop), s3(loop), s4(loop), &
                                   s5(loop), s6(loop), s7(loop), &
                                   i, j, k,                       &
                                   sgntz, sgntx, sgnty,           &
                                   sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !!DIR$ FORCEINLINE
                     !tt1(loop) = fteik_localSolverExplicit64fF( &
                     !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &
                     !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &
                     !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &
                     !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &
                     !             .TRUE.,                        &
                     !              i, j, k,                       &
                     !              sgntz, sgntx, sgnty,           &
                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &
                     !                             slowWork(8*(loop-1)+1),  &
                     !                             .TRUE.,                        &
                     !                             i, j, k,                       &
                     !                             sgntz, sgntx, sgnty,           &
                     !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                  ENDIF
               ENDDO
               !OMP SIMD
               DO node=node1,node2 !loop=1,n2
                  loop = node - node1 + 1 !node = knode + loop - 1
                  IF (lupdInit8(node)) THEN
                     indx = ijkv8(4*(node-1)+4)
                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))
                  ENDIF
               ENDDO
            ENDDO ! parallel loop on nodes in level
            !$OMP END DO
            !$OMP BARRIER
 11      CONTINUE ! loop on levels
      ENDIF
      DEALLOCATE(slowWork)
      DEALLOCATE(ttWork)
      DEALLOCATE(tt1)
      DEALLOCATE(t1)
      DEALLOCATE(t2)
      DEALLOCATE(t3)
      DEALLOCATE(t4)
      DEALLOCATE(t5)
      DEALLOCATE(t6)
      DEALLOCATE(t7)
      DEALLOCATE(t8)
      DEALLOCATE(s1)
      DEALLOCATE(s2)
      DEALLOCATE(s3)
      DEALLOCATE(s4)
      DEALLOCATE(s5)
      DEALLOCATE(s6)
      DEALLOCATE(s7)
      !$OMP END PARALLEL
      RETURN
      END SUBROUTINE
