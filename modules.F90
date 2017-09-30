!     MODULE FTEIK_CONSTANTS64F
!        USE ISO_C_BINDING
!        IMPLICIT NONE
!        REAL(C_DOUBLE), PARAMETER :: FTEIK_HUGE = 99999.d0
!        REAL(C_DOUBLE), PARAMETER :: zero = 0.d0
!        REAL(C_DOUBLE), PARAMETER :: one = 1.d0
!        REAL(C_DOUBLE), PARAMETER :: DBL_EPSILON = EPSILON(1.d0)
!        REAL(C_DOUBLE), PARAMETER :: perturbSource = 0.0001d0
!        INTEGER(C_INT), PARAMETER :: chunkSize = 16
!     END MODULE !FTEIK_CONSTANTS64F
!                                                                                        !
!========================================================================================!
!                                                                                        !
!     MODULE FTEIK_MODEL64F
!        USE ISO_C_BINDING
!        USE FTEIK_CONSTANTS64F, ONLY : zero 
!        IMPLICIT NONE
!        REAL(C_DOUBLE), ALLOCATABLE, SAVE :: slow(:) !> Slowness model (seconds/meters)
!        !DIR$ ATTRIBUTES ALIGN: 64 :: slow
!        REAL(C_DOUBLE), SAVE :: dx = zero     !> x origin (meters)
!        REAL(C_DOUBLE), SAVE :: dy = zero     !> y origin (meters)
!        REAL(C_DOUBLE), SAVE :: dz = zero     !> z origin (meters)
!        REAL(C_DOUBLE), SAVE :: x0 = zero     !> x origin (meters)
!        REAL(C_DOUBLE), SAVE :: y0 = zero     !> y origin (meters)
!        REAL(C_DOUBLE), SAVE :: z0 = zero     !> z origin (meters)
!        INTEGER(C_INT), SAVE :: ncell =-1     !> Number of cells in slowness model
!        INTEGER(C_INT), SAVE :: ngrd =-1      !> Number of grid points (= nz*nx*ny)
!        INTEGER(C_INT), SAVE :: nx =-1        !> Number of x grid points in travel time field
!        INTEGER(C_INT), SAVE :: ny =-1        !> Number of y grid points in travel time field
!        INTEGER(C_INT), SAVE :: nz =-1        !> Number of z grid points in travel time field
!        INTEGER(C_INT), SAVE:: nzx =-1       !> nz*nx
!        INTEGER(C_INT), SAVE :: nzm1_nxm1 =-1 !> (nz - 1)*(nx - 1)
!        INTEGER(C_INT), SAVE :: nzm1 =-1      !> (nz - 1)
!        LOGICAL(C_BOOL), SAVE :: lhaveModel = .FALSE. !> If true slowness model was set

!        INTERFACE

!           SUBROUTINE fteik_model_finalizeF()             &
!                      BIND(C, NAME='fteik_model_finalizeF')
!           USE ISO_C_BINDING
!           END SUBROUTINE

!           SUBROUTINE fteik_model_getVelocityModel64fF(nwork, nv, vel, ierr) &
!                      BIND(C, NAME='fteik_model_getVelocityModel64fF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE 
!           INTEGER(C_INT), VALUE, INTENT(IN) :: nwork
!           REAL(C_DOUBLE), INTENT(OUT) :: vel(nwork)
!           INTEGER(C_INT), INTENT(OUT) :: nv, ierr
!           END SUBROUTINE

!           SUBROUTINE fteik_model_intializeGeometryF(nz, nx, ny,    &
!                                                     dz, dx, dy,    &
!                                                     z0, x0, y0,    &
!                                                     ierr)          &
!                      BIND(C, NAME='fteik_model_initializeGeometryF')
!           USE ISO_C_BINDING
!           INTEGER(C_INT), VALUE, INTENT(IN) :: nz, nx, ny
!           REAL(C_DOUBLE), VALUE, INTENT(IN) :: dz, dx, dy, z0, x0, y0
!           END SUBROUTINE

!           SUBROUTINE fteik_model_setGridSizeF(nzIn, nxIn, nyIn, ierr) &
!                      BIND(C, NAME='fteik_model_setGridSizeF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), VALUE, INTENT(IN) :: nzIn, nxIn, nyIn
!           INTEGER(C_INT), INTENT(OUT) :: ierr
!           END SUBROUTINE

!           SUBROUTINE fteik_model_setOriginF(z0In, x0In, y0In) &
!                      BIND(C, NAME='fteik_model_setOriginF')
!           USE ISO_C_BINDING
!           REAL(C_DOUBLE), VALUE, INTENT(IN) :: x0In, y0In, z0In 
!           END SUBROUTINE

!           SUBROUTINE fteik_model_setGridSpacingF(dzIn, dxIn, dyIn, ierr) &
!                      BIND(C, NAME='fteik_model_setGridSpacingF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           REAL(C_DOUBLE), VALUE, INTENT(IN) :: dzIn, dxIn, dyIn
!           INTEGER(C_INT), INTENT(OUT) :: ierr
!           END SUBROUTINE

!           SUBROUTINE fteik_model_setVelocityModel64fF(nv, vel, ierr) &
!                      BIND(C, NAME='fteik_model_setVelocityModel64fF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), VALUE, INTENT(IN) :: nv
!           REAL(C_DOUBLE), INTENT(IN) :: vel(nv)
!           INTEGER(C_INT), INTENT(OUT) :: ierr
!           END SUBROUTINE

!           SUBROUTINE fteik_model_setVelocityModel32fF(nv, vel4, ierr) &
!                      BIND(C, NAME='fteik_model_setVelocityModel32fF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), VALUE, INTENT(IN) :: nv
!           REAL(C_FLOAT), INTENT(IN) :: vel4(nv)
!           INTEGER(C_INT), INTENT(OUT) :: ierr
!           END SUBROUTINE

!           PURE INTEGER(C_INT) FUNCTION fteik_model_velGrid2indexF(i, j, k,        &
!                                                                nzm1, nzm1_nxm1)   &
!                            BIND(C, NAME='fteik_model_velGrid2indexF')             &
!                            RESULT(velGrid2IndexF)
!           !$OMP DECLARE SIMD(fteik_model_velGrid2indexF) UNIFORM(nzm1, nzm1_nxm1)
!           USE ISO_C_BINDING
!           IMPLICIT NONE 
!           INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nzm1, nzm1_nxm1
!           END FUNCTION

!           PURE SUBROUTINE fteik_model_index2gridF(igrd, i, j, k, ierr) &
!                           BIND(C, NAME='fteik_model_index2gridF')
!           !$OMP DECLARE SIMD(fteik_model_index2gridF) UNIFORM(ierr)
!           USE ISO_C_BINDING
!           INTEGER(C_INT), INTENT(IN), VALUE :: igrd
!           INTEGER(C_INT), INTENT(OUT) :: i, j, k, ierr
!           END SUBROUTINE

!           PURE INTEGER(C_INT) FUNCTION fteik_model_grid2indexF(i, j, k, nz, nzx) &
!                               BIND(C, NAME='fteik_model_grid2indexF')            &
!                               RESULT(grid2indexF)
!           !$OMP DECLARE SIMD(fteik_model_grid2indexF) UNIFORM(nz, nzx)
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nzx
!           END FUNCTION

!        END INTERFACE
!     END MODULE !FTEIK_MODEL64F
!                                                                                        !
!========================================================================================!
!                                                                                        !
!     MODULE FTEIK_RECEIVER64F
!        USE ISO_C_BINDING
!        ! receiver indices
!        REAL(C_DOUBLE), ALLOCATABLE, SAVE :: zdr(:), xdr(:), ydr(:)
!        INTEGER(C_INT), ALLOCATABLE, SAVE :: zri(:), xri(:), yri(:)
!        INTEGER(C_INT), SAVE :: nrec 
!        INTEGER(C_INT), PARAMETER :: INTERP_NEAREST = 0
!        INTEGER(C_INT), PARAMETER :: INTERP_LINEAR = 1
!        INTEGER(C_INT), PARAMETER :: INTERP_HIGHACCURACY = 2
!        INTERFACE
!           SUBROUTINE fteik_receiver_initialize64fF(nrecIn, z, x, y, ierr) &
!           BIND(C, NAME='fteik_receiver_initialize64fF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE 
!           INTEGER(C_INT), INTENT(IN), VALUE :: nrecIn
!           REAL(C_DOUBLE), INTENT(IN) :: z(nrecIn), x(nrecIn), y(nrecIn)
!           INTEGER(C_INT), INTENT(OUT) :: ierr 
!           END SUBROUTINE

!           SUBROUTINE fteik_receiver_getTravelTimes64fF(nrecIn, ttr, ierr) &
!           BIND(C, NAME='fteik_receiver_getTravelTimes64fF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE 
!           INTEGER(C_INT), VALUE, INTENT(IN) :: nrecIn
!           REAL(C_DOUBLE), INTENT(OUT) :: ttr(nrecIn)
!           INTEGER(C_INT), INTENT(OUT) :: ierr 
!           END SUBROUTINE

!           SUBROUTINE fteik_receiver_finalizeF() &
!           BIND(C, NAME='fteik_receiver_finalizeF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE 
!           END SUBROUTINE
!        END INTERFACE
!     END MODULE FTEIK_RECEIVER64F
!                                                                                        !
!========================================================================================!
!                                                                                        !
!     MODULE FTEIK_SOURCE64F
!        USE ISO_C_BINDING
!        IMPLICIT NONE
!        REAL(C_DOUBLE), ALLOCATABLE, SAVE :: zsv(:), xsv(:), ysv(:),    &
!                                             zsav(:), xsav(:), ysav(:), &
!                                             zsrc(:), xsrc(:), ysrc(:)
!        INTEGER(C_INT), ALLOCATABLE, SAVE :: zsiv(:), xsiv(:), ysiv(:)
!        !DIR$ ATTRIBUTES ALIGN: 64 :: zsv, xsv, ysv, zsav, xsav, ysav
!        !DIR$ ATTRIBUTES ALIGN: 64 :: zsrc, xsrc, ysrc, zsiv, xsiv, ysiv
!        INTEGER(C_INT), SAVE :: nsrc = 0
!        LOGICAL(C_BOOL), SAVE :: lhaveSource = .FALSE.

!        INTERFACE

!           SUBROUTINE fteik_source_finalizeF()             &
!                      BIND(C, NAME='fteik_source_finalizeF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           END SUBROUTINE

!           SUBROUTINE fteik_source_initialize64fF(nsrcIn, zsrcIn, xsrcIn, ysrcIn, ierr) &
!                      BIND(C, NAME='fteik_source_initialize64fF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), VALUE, INTENT(IN) :: nsrcIn
!           REAL(C_DOUBLE), INTENT(IN) :: zsrcIn(nsrcIn), xsrcIn(nsrcIn), ysrcIn(nsrcIn)
!           INTEGER(C_INT), INTENT(OUT) :: ierr
!           END SUBROUTINE

!           SUBROUTINE fteik_source_getSolverInfo64fF(isrc,          &
!                                                     zsi, xsi, ysi, &
!                                                     zsa, xsa, ysa, &
!                                                     szero, szero2, &
!                                                     ierr)          &
!                      BIND(C, NAME='fteik_source_getSolverInfo64fF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
!           REAL(C_DOUBLE), INTENT(OUT) :: zsa, xsa, ysa, szero, szero2
!           INTEGER(C_INT), INTENT(OUT) :: zsi, xsi, ysi, ierr
!           END SUBROUTINE

!           SUBROUTINE fteik_source_getSourceIndices64fF(isrc, zsa, xsa, ysa, ierr) &
!                      BIND(C, NAME='fteik_source_getSourceIndices64fF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
!           REAL(C_DOUBLE), INTENT(OUT) :: zsa, xsa, ysa
!           INTEGER(C_INT), INTENT(OUT) :: ierr
!           END SUBROUTINE

!           SUBROUTINE fteik_source_getSourceIndices32iF(isrc, zsi, xsi, ysi, ierr) &
!                      BIND(C, NAME='fteik_source_getSourceIndices32iF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
!           INTEGER(C_INT), INTENT(OUT) :: zsi, xsi, ysi, ierr
!           END SUBROUTINE

!           SUBROUTINE fteik_source_getSzero64fF(isrc, szero, szero2, ierr)   &
!                      BIND(C, NAME='fteik_source_getSzero64fF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
!           REAL(C_DOUBLE), INTENT(OUT) :: szero, szero2
!           INTEGER(C_INT), INTENT(OUT) :: ierr
!           END SUBROUTINE

!        END INTERFACE
!     END MODULE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!     MODULE FTEIK_SOLVER64F
!        USE FTEIK_CONSTANTS64F, ONLY : zero
!        USE ISO_C_BINDING
!        IMPLICIT NONE
!        !> Holds the travel-times.
!        REAL(C_DOUBLE), ALLOCATABLE, SAVE :: ttimes(:)
!        !DIR$ ATTRIBUTES ALIGN: 64 :: ttimes
!        !> Maps from level to start node.
!        INTEGER(C_INT), ALLOCATABLE, SAVE :: levelPtr(:)
!        !> Maps from node in level to (iz, ix, iy, and grid point).
!        INTEGER(C_INT), ALLOCATABLE, SAVE :: ijkv1(:), ijkv2(:), ijkv3(:), ijkv4(:), &
!                                             ijkv5(:), ijkv6(:), ijkv7(:), ijkv8(:)
!        !> True if the node'th node is to be updated in the given sweep.
!        LOGICAL(C_BOOL), ALLOCATABLE, SAVE :: lupd1(:), lupd2(:), lupd3(:), lupd4(:), &
!                                              lupd5(:), lupd6(:), lupd7(:), lupd8(:)
!        !> True if the node'th node is to be updated in the initalization phase.
!        LOGICAL(C_BOOL), ALLOCATABLE, SAVE :: lupdInit1(:), lupdInit2(:), &
!                                              lupdInit3(:), lupdInit4(:), &
!                                              lupdInit5(:), lupdInit6(:), &
!                                              lupdInit7(:), lupdInit8(:)
!        !> Defines transition from spherical to Cartesian during the intialization.
!        REAL(C_DOUBLE), SAVE :: epsS2C = zero
!        !> Number of Gauss-Sediel iterations.
!        INTEGER(C_INT), SAVE :: nsweep = 0 
!        !> Flag indicating whether or not traveltimes were computed.
!        LOGICAL(C_BOOL), SAVE :: lhaveTimes = .FALSE.
!        !> Number of leves in the level scheduling method.
!        INTEGER(C_INT), SAVE :: nLevels = 0
!        !> Largest level size in level scheduling method.
!        INTEGER(C_INT), SAVE :: maxLevelSize = 0

!        INTERFACE

!           SUBROUTINE fteik_solver_finalizeF()             &
!                      BIND(C, NAME='fteik_solver_finalizeF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           END SUBROUTINE

!           SUBROUTINE fteik_solver_initialize64fF(nzIn, nxIn, nyIn,      &
!                                                  z0In, x0In, y0In,      &
!                                                  dzIn, dxIn, dyIn,      &
!                                                  nsweepIn, epsIn, ierr) &
!                      BIND(C, NAME='fteik_solver_initialize64fF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           REAL(C_DOUBLE), VALUE, INTENT(IN) :: x0In, y0In, z0In
!           REAL(C_DOUBLE), VALUE, INTENT(IN) :: dxIn, dyIn, dzIn
!           REAL(C_DOUBLE), VALUE, INTENT(IN) :: epsIn
!           INTEGER(C_INT), VALUE, INTENT(IN) :: nxIn, nyIn, nzIn
!           INTEGER(C_INT), VALUE, INTENT(IN) :: nsweepIn
!           INTEGER(C_INT), INTENT(OUT) :: ierr
!           END SUBROUTINE

!           SUBROUTINE fteik_solver_setSphereToCartEpsilonF(epsIn, ierr)  &
!                      BIND(C, NAME='fteik_solver_setSphereToCartEpsilonF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           REAL(C_DOUBLE), VALUE, INTENT(IN) :: epsIn
!           INTEGER(C_INT), INTENT(OUT) :: ierr
!           END SUBROUTINE

!           SUBROUTINE fteik_solver_setNumberOfSweepsF(nsweepIn, ierr) &
!                      BIND(C, NAME='fteik_solver_setNumberOfSweepsF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), VALUE, INTENT(IN) :: nsweepIn
!           INTEGER(C_INT), INTENT(OUT) :: ierr
!           END SUBROUTINE

!           SUBROUTINE fteik_solver_setUpdateNodesF(sweep, nLevels, linitk, isrc, &
!                                                   levelPtr, ijkv, lupd, ierr)   &
!                      BIND(C, NAME='fteik_solver_setUpdateNodesF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), VALUE, INTENT(IN) :: isrc, sweep, nLevels
!           LOGICAL(C_BOOL), VALUE, INTENT(IN) :: linitk
!           INTEGER(C_INT), DIMENSION(:), INTENT(IN) :: ijkv, levelPtr
!           LOGICAL(C_BOOL), DIMENSION(:), INTENT(INOUT) :: lupd
!           INTEGER(C_INT), INTENT(OUT) :: ierr
!           END SUBROUTINE

!           SUBROUTINE fteik_solver_computeGraphF(ierr)           &
!                      BIND(C, NAME='fteik_solver_computeGraphF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE 
!           INTEGER(C_INT), INTENT(OUT) :: ierr 
!           END SUBROUTINE

!           SUBROUTINE fteik_solver_getSweepLimitsF(sweep, linitk,          &
!                                                   nz, nx, ny,             &
!                                                   zsi, xsi, ysi,          &
!                                                   z1, z2, x1, x2, y1, y2, &
!                                                   ierr)                   &
!                     BIND(C, NAME='fteik_solver_getSweepLimitsF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), VALUE, INTENT(IN) :: sweep, nx, ny, nz, xsi, ysi, zsi
!           LOGICAL(C_BOOL), VALUE, INTENT(IN) :: linitk
!           INTEGER(C_INT), INTENT(OUT) :: x1, x2, y1, y2, z1, z2, ierr
!           END SUBROUTINE

!           SUBROUTINE fteik_solver_meshConstants64fF() &
!                      BIND(C, NAME='fteik_solver_meshConstants64fF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           END SUBROUTINE

!        END INTERFACE 
!     END MODULE

!     MODULE FTEIK_LOCALSOLVER64F
!        USE ISO_C_BINDING
!        INTEGER(C_INT), SAVE :: zsi, xsi, ysi
!        REAL(C_DOUBLE), SAVE :: zsa, xsa, ysa 
!        REAL(C_DOUBLE), SAVE :: dz, dx, dy
!        REAL(C_DOUBLE), SAVE :: szero, szero2
!        REAL(C_DOUBLE), SAVE :: dxi, dyi, dzi
!        REAL(C_DOUBLE), SAVE :: dx2, dy2, dz2
!        REAL(C_DOUBLE), SAVE :: dx2i, dy2i, dz2i
!        REAL(C_DOUBLE), SAVE :: dsum, dsumi, epsSolver
!        REAL(C_DOUBLE), SAVE :: dz2dx2, dz2dy2, dx2dy2
!        REAL(C_DOUBLE), SAVE :: dz2i_dx2i, dz2i_dy2i, dx2i_dy2i
!        REAL(C_DOUBLE), SAVE :: dz2i_p_dx2i, dz2i_p_dy2i, dx2i_p_dy2i
!        REAL(C_DOUBLE), SAVE :: dz2i_p_dy2i_inv, dz2i_p_dx2i_inv, dx2i_p_dy2i_inv
!        INTERFACE

!           SUBROUTINE fteik_localSolver_initialize64fF(isrc, dest, ts8, ierr)   &
!           BIND(C, NAME='fteik_localSolver_initialize64fF')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
!           REAL(C_DOUBLE), INTENT(OUT) :: ts8(8)
!           INTEGER(C_INT), INTENT(OUT) :: dest(8), ierr
!           END SUBROUTINE

!           PURE REAL(C_DOUBLE)                                       &
!           FUNCTION fteik_localSolver_tAna64fF(i, j, k, dz, dx, dy,  &
!                                               zsa, xsa, ysa, szero) &
!           BIND(C, NAME='fteik_localSolver_tAna64fF')
!           !$OMP DECLARE SIMD(fteik_localSolver_tAna64fF) &
!           !$OMP UNIFORM(dz, dx, dy, zsa, xsa, ysa, szero)
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           REAL(C_DOUBLE), INTENT(IN), VALUE :: dz, dx, dy, szero, zsa, xsa, ysa 
!           INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k
!           END FUNCTION

!           PURE SUBROUTINE fteik_localSolver_tAnaD64fF(t_anad, tzc, txc, tyc, &
!                                                       i, j, k,               &
!                                                       dz, dx, dy,            &
!                                                       zsa, xsa, ysa, szero)  &
!                           BIND(C, NAME='fteik_localSolver_tAnaD64fF')
!           !$OMP DECLARE SIMD(fteik_localSolver_tAnaD64fF) &
!           !$OMP UNIFORM(dz, dx, dy, zsa, xsa, ysa, szero)
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           REAL(C_DOUBLE), INTENT(IN), VALUE :: dz, dx, dy, zsa, xsa, ysa, szero
!           INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k
!           REAL(C_DOUBLE), INTENT(OUT) :: t_anad, tzc, txc, tyc
!           END SUBROUTINE

!        END INTERFACE

!     END MODULE 
!                                                                                        !
!========================================================================================!
!                                                                                        !

      MODULE FTEIK_UTILS64F
         USE ISO_C_BINDING
         REAL(C_DOUBLE), ALLOCATABLE :: ttimes(:)
         REAL(C_DOUBLE), ALLOCATABLE :: slow(:)
         !!!DIR$ ATTRIBUTES ALIGN: 64 :: ttimes, slow
         REAL(C_DOUBLE) dx, dy, dz
         REAL(C_DOUBLE) xsa, ysa, zsa
         REAL(C_DOUBLE) x0, y0, z0 
         REAL(C_DOUBLE) epsS2C
         INTEGER(C_INT) ncell, ngrd, nx, ny, nz, nzx, nzm1, nzm1_nxm1
         INTEGER(C_INT), ALLOCATABLE :: levelPtr(:)
         INTEGER(C_INT), ALLOCATABLE :: ijkv1(:), ijkv2(:), ijkv3(:), ijkv4(:), &
                                        ijkv5(:), ijkv6(:), ijkv7(:), ijkv8(:)
         LOGICAL(C_BOOL), ALLOCATABLE :: lupd1(:), lupd2(:), lupd3(:), lupd4(:), &
                                         lupd5(:), lupd6(:), lupd7(:), lupd8(:), &
                                         lupdInit1(:), lupdInit2(:), &
                                         lupdInit3(:), lupdInit4(:), & 
                                         lupdInit5(:), lupdInit6(:), &
                                         lupdInit7(:), lupdInit8(:)
         INTEGER(C_INT) nLevels, maxLevelSize
 
         INTEGER(C_INT) nsweep ! number of sweeps over model
         ! initialization variables
         LOGICAL(C_BOOL) lhaveGrid, lhaveGridSpacing, lhaveSource, &
                         lhaveSlownessModel, lhaveTravelTimes
         ! receiver indices
         REAL(C_DOUBLE), ALLOCATABLE :: zdr(:), xdr(:), ydr(:)
         INTEGER(C_INT), ALLOCATABLE :: zri(:), xri(:), yri(:)
         INTEGER(C_INT) nrec
         ! source variables
         REAL(C_DOUBLE) szero
         INTEGER(C_INT) xsi, ysi, zsi
         ! private variables for solver to compute on startup
         REAL(C_DOUBLE), ALLOCATABLE :: tt1(:)
         !DIR$ ATTRIBUTES ALIGN: 64 :: tt1
         REAL(C_DOUBLE) szero2
         REAL(C_DOUBLE) dxi, dyi, dzi
         REAL(C_DOUBLE) dx2, dy2, dz2
         REAL(C_DOUBLE) dx2i, dy2i, dz2i
         REAL(C_DOUBLE) dsum, dsumi, epsSolver
         REAL(C_DOUBLE) dz2dx2, dz2dy2, dx2dy2
         REAL(C_DOUBLE) dz2i_dx2i, dz2i_dy2i, dx2i_dy2i
         REAL(C_DOUBLE) dz2i_p_dx2i, dz2i_p_dy2i, dx2i_p_dy2i 
         REAL(C_DOUBLE) dz2i_p_dy2i_inv, dz2i_p_dx2i_inv, dx2i_p_dy2i_inv
         ! permutations
         INTEGER(C_INT), PARAMETER :: ttPerm1(8) = [7, 5, 6, 3, 4, 2, 1, 8]
         INTEGER(C_INT), PARAMETER :: ttPerm2(8) = [6, 3, 7, 5, 1, 8, 4, 2]
         INTEGER(C_INT), PARAMETER :: ttPerm3(8) = [4, 2, 1, 8, 7, 5, 6, 3]
         INTEGER(C_INT), PARAMETER :: ttPerm4(8) = [1, 8, 4, 2, 6, 3, 7, 5]
         INTEGER(C_INT), PARAMETER :: ttPerm5(8) = [5, 7, 3, 6, 2, 4, 8, 1]
         INTEGER(C_INT), PARAMETER :: ttPerm6(8) = [3, 6, 5, 7, 8, 1, 2, 4]
         INTEGER(C_INT), PARAMETER :: ttPerm7(8) = [2, 4, 8, 1, 5, 7, 3, 6]
         INTEGER(C_INT), PARAMETER :: ttPerm8(8) = [8, 1, 2, 4, 3, 6, 5, 7]
         INTEGER(C_INT), PARAMETER :: v2l1(19) = [1,  5,  9, 13, 15, 17, 19,  6, 11, &
                                                 18,  3, 10, 16, 12,  2,  7, 14,  8,  4]
         INTEGER(C_INT), PARAMETER :: v2l2(19) = [1,  9, 15, 11,  3,  5, 10, 13, 16, &
                                                 17, 19,  6, 12, 18,  2,  4,  7, 14,  8]
         INTEGER(C_INT), PARAMETER :: v2l3(19) = [1,  5, 13,  6,  3,  2,  7,  9, 14, &
                                                 15, 17, 19,  8, 11, 18,  4, 10, 16, 12]
         INTEGER(C_INT), PARAMETER :: v2l4(19) = [1,  3,  5, 13,  6,  2,  9, 15, 11, &
                                                  4,  7, 10, 14, 16, 17, 19,  8, 12, 18]
         INTEGER(C_INT), PARAMETER :: v2l5(19) = [5,  9, 17,  1,  6, 11, 13, 15, 18, &
                                                 19, 10,  3, 12, 16,  7,  2,  8, 14,  4]
         INTEGER(C_INT), PARAMETER :: v2l6(19) = [9,  1, 11, 15,  5, 10, 17,  3,  6, &
                                                 12, 13, 16, 18, 19,  2,  7,  4,  8, 14]
         INTEGER(C_INT), PARAMETER :: v2l7(19) = [5,  1,  6, 13,  3,  7,  9, 17,  2, &
                                                  8, 11, 14, 15, 18, 19, 10,  4, 12, 16] 
         INTEGER(C_INT), PARAMETER :: v2l8(19) = [1,  5,  3,  6, 13,  9,  2, 11, 15, &
                                                  7, 10, 17,  4,  8, 12, 14, 16, 18, 19]
         INTEGER(C_INT), PARAMETER :: sgntzv(8) = [1, 1, 1, 1,-1,-1,-1,-1]
         INTEGER(C_INT), PARAMETER :: sgntxv(8) = [1,-1, 1,-1, 1,-1, 1,-1]
         INTEGER(C_INT), PARAMETER :: sgntyv(8) = [1, 1,-1,-1, 1, 1,-1,-1]
         INTEGER(C_INT), PARAMETER :: sgnvzv(8) = [1, 1, 1, 1, 0, 0, 0, 0]
         INTEGER(C_INT), PARAMETER :: sgnvxv(8) = [1, 0, 1, 0, 1, 0, 1, 0]
         INTEGER(C_INT), PARAMETER :: sgnvyv(8) = [1, 1, 0, 0, 1, 1, 0, 0]
         !DIR$ ATTRIBUTES ALIGN: 64:: ttPerm1, ttPerm2, ttPerm3, ttPerm4
         !DIR$ ATTRIBUTES ALIGN: 64:: ttPerm5, ttPerm6, ttPerm7, ttPerm8
         !DIR$ ATTRIBUTES ALIGN: 64:: v2l1, v2l2, v2l3, v2l4, v2l5, v2l6, v2l7, v2l8
         ! constants
         REAL(C_DOUBLE), PARAMETER :: FTEIK_HUGE = 99999.d0
         REAL(C_DOUBLE), PARAMETER :: zero = 0.d0
         REAL(C_DOUBLE), PARAMETER :: one = 1.d0
         REAL(C_DOUBLE), PARAMETER :: DBL_EPSILON = EPSILON(1.d0)
         INTEGER(C_INT), PARAMETER :: chunkSize = 16
         ! variables that must be saved
         SAVE :: ttimes, tt1
         SAVE :: ijkv1, ijkv2, ijkv3, ijvk4, ijkv5, ijvk5, ijkv6, ijkv7, ijkv8,     &
                 levelPtr, nLevels
         SAVE :: ludp1, lupd2, ludp3, lupd4, lupd5, lupd6, lupd7, lupd8
         SAVE :: lupdInit1, lupdInit2, lupdInit3, lupdInit4, &
                 lupdInit5, lupdInit6, lupdInit7, lupdInit8 
         SAVE :: dx, dy, dz, epsS2C, x0, y0, z0, xsa, ysa, zsa, xsi, ysi, zsi 
         SAVE :: ncell, ngrd, nx, ny, nzx, nzm1, nzm1_nxm1
         SAVE :: nsweep
         SAVE :: lhaveGrid, lhaveGridSpacing, lhaveSource,  &
                 lhaveSlownessModel, lhaveTravelTimes
         SAVE :: zdr, xdr, ydr, zri, xri, yri, nrec
         INTERFACE 

            SUBROUTINE fteik_meshConstants64fF() BIND(C, NAME='fteik_meshConstants64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            END


            PURE REAL(C_DOUBLE) FUNCTION fteik_tAna64fF(i, j, k, dz, dx, dy,   &
                                                        zsa, xsa, ysa, szero)  &
                                BIND(C, NAME='fteik_tAna64fF')
            !$OMP DECLARE SIMD(fteik_tAna64fF) UNIFORM(dz, dx, dy, zsa, xsa, ysa, szero) 
            USE ISO_C_BINDING
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN), VALUE :: dx, dy, dz, szero, xsa, ysa, zsa
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k
            END FUNCTION

            PURE SUBROUTINE fteik_tAnaD64fF(t_anad, tzc, txc, tyc, i, j, k,   &
                                            dz, dx, dy, zsa, xsa, ysa, szero) &
                            BIND(C, NAME='fteik_tAnaD64fF')
            !$OMP DECLARE SIMD(fteik_tAnaD64fF) UNIFORM(dz, dx, dy, zsa, xsa, ysa, szero)
            USE ISO_C_BINDING
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN), VALUE :: dz, dx, dy, zsa, xsa, ysa, szero
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k
            REAL(C_DOUBLE), INTENT(OUT) :: t_anad, tzc, txc, tyc
            END SUBROUTINE

!           REAL(C_DOUBLE) FUNCTION T_ANA64F(i, j, k) &
!                             BIND(C, NAME='t_ana64f')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), INTENT(IN) :: i, j, k
!           END FUNCTION

!           REAL(C_DOUBLE) FUNCTION T_ANAD64F(tzc, txc, tyc, i, j, k) &
!                                   BIND(C, NAME='t_anad64f')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           INTEGER(C_INT), INTENT(IN) :: i, j, k
!           REAL(C_DOUBLE), INTENT(OUT) :: tzc, txc, tyc
!           END FUNCTION

            SUBROUTINE fteik_finalizeF() BIND(C, NAME='fteik_finalizeF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            END SUBROUTINE

            SUBROUTINE fteik_initializeF(nzIn, nxIn, nyIn,      &
                                         z0In, x0In, y0In,      &
                                         dzIn, dxIn, dyIn,      &
                                         nsweepIn, epsIn, ierr) &
                       BIND(C, NAME='fteik_initializeF')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            REAL(C_DOUBLE), INTENT(IN) :: x0In, y0In, z0In
            REAL(C_DOUBLE), INTENT(IN) :: dxIn, dyIn, dzIn
            REAL(C_DOUBLE), INTENT(IN) :: epsIn
            INTEGER(C_INT), INTENT(IN) :: nxIn, nyIn, nzIn
            INTEGER(C_INT), INTENT(IN) :: nsweepIn
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            REAL(C_DOUBLE) &
            FUNCTION fteik_localSolver64fF(tt, slowLoc, linitk,             &
                                           i, j, k,                         &
                                           sgntz, sgntx, sgnty,             &
                                           sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
                                    RESULT(tupd) &
                                    BIND(C, NAME='fteik_localSolver64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: tt(8), slowLoc(7)
            REAL(C_DOUBLE), INTENT(IN) ::  sgnrz_dzi, sgnrx_dxi, sgnry_dyi
            INTEGER(C_INT), INTENT(IN) :: i, j, k, sgntx, sgnty, sgntz
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            END FUNCTION

            REAL(C_DOUBLE)                                                 &
            FUNCTION fteik_localSolverExplicit64fF(tv, te, tn, tev,                  &
                                                   ten, tnv, tnve, tt0,              &
                                                   slow1, slow2, slow3, slow4,       &
                                                   slow5, slow6, slow7,              &
                                                   linitk,                           &
                                                   i, j, k,                          &
                                                   sgntz, sgntx, sgnty,              &
                                                   sgnrz_dzi, sgnrx_dxi, sgnry_dyi)  &
            RESULT(tupd) BIND(C, NAME='fteik_localSolverExplicit64fF')
            !$OMP DECLARE SIMD(fteik_localSolverExplicit64fF) &
            !$OMP UNIFORM(linitk, sgntz, sgntx, sgnty, sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
            USE ISO_C_BINDING
            IMPLICIT NONE 
            REAL(C_DOUBLE), INTENT(IN), VALUE :: tv, te, tn, tev, ten, tnv, tnve, tt0
            REAL(C_DOUBLE), INTENT(IN), VALUE :: slow1, slow2, slow3, slow4, &
                                                 slow5, slow6, slow7
            REAL(C_DOUBLE), INTENT(IN), VALUE :: sgnrz_dzi, sgnrx_dxi, sgnry_dyi
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, sgntx, sgnty, sgntz
            LOGICAL(C_BOOL), INTENT(IN), VALUE :: linitk
            END FUNCTION

            !PURE REAL(C_DOUBLE)                                         &
            !FUNCTION fteik_localSolverNoInit64fF(tv, te, tn, tev,       &
            !                               ten, tnv, tnve, tt0,         &
            !                               slow1, slow2, slow3, slow4,  &
            !                               slow5, slow6, slow7) &
            !RESULT(tupd) BIND(C, NAME='fteik_localSolverNoInit64fF')
            !!$OMP DECLARE SIMD(fteik_localSolverNoInit64fF)
            !USE ISO_C_BINDING
            !IMPLICIT NONE
            !REAL(C_DOUBLE), INTENT(IN), VALUE :: tv, te, tn, tev, ten, tnv, tnve, tt0 
            !REAL(C_DOUBLE), INTENT(IN), VALUE :: slow1, slow2, slow3, slow4, &
            !                                     slow5, slow6, slow7
            !END FUNCTION

            PURE SUBROUTINE fteik_localSolverNoInit64fF(n, tv, te, tn, tev,         &
                                                        ten, tnv, tnve, tt0,        &
                                                        slow1, slow2, slow3, slow4, &
                                                        slow5, slow6, slow7,        &
                                                        tupd)                       &
            BIND(C, NAME='fteik_localSolverNoInit64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: n
            REAL(C_DOUBLE), INTENT(IN) :: tv(n), te(n), tn(n), tev(n),  &
                                          ten(n), tnv(n), tnve(n), tt0(n)
            REAL(C_DOUBLE), INTENT(IN) :: slow1(n), slow2(n), slow3(n), slow4(n), &
                                          slow5(n), slow6(n), slow7(n)
            REAL(C_DOUBLE), INTENT(OUT) :: tupd(n)
            !DIR$ ASSUME_ALIGNED tv: 64
            !DIR$ ASSUME_ALIGNED te: 64
            !DIR$ ASSUME_ALIGNED tn: 64
            !DIR$ ASSUME_ALIGNED tev: 64
            !DIR$ ASSUME_ALIGNED ten: 64
            !DIR$ ASSUME_ALIGNED tnv: 64
            !DIR$ ASSUME_ALIGNED tnve: 64
            !DIR$ ASSUME_ALIGNED tt0: 64
            !DIR$ ASSUME_ALIGNED slow1: 64
            !DIR$ ASSUME_ALIGNED slow2: 64
            !DIR$ ASSUME_ALIGNED slow3: 64
            !DIR$ ASSUME_ALIGNED slow4: 64
            !DIR$ ASSUME_ALIGNED slow5: 64
            !DIR$ ASSUME_ALIGNED slow6: 64
            !DIR$ ASSUME_ALIGNED slow7: 64
            !DIR$ ASSUME_ALIGNED tupd: 64
            END SUBROUTINE
 
            PURE REAL(C_DOUBLE)                                             &
            FUNCTION fteik_localSolverInit64fF(tv, te, tn, tev,             &
                                           ten, tnv, tnve, tt0,             &
                                           slow1, slow2, slow3, slow4,      &
                                           slow5, slow6, slow7,             &
                                           i, j, k,                         &
                                           sgntz, sgntx, sgnty,             &
                                           sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
            RESULT(tupd) BIND(C, NAME='fteik_localSolverInit64fF')
            !!$OMP DECLARE SIMD(fteik_localSolverInit64fF) &
            !!$OMP UNIFORM(sgntz, sgntx, sgnty, sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
            USE ISO_C_BINDING
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN), VALUE :: tv, te, tn, tev, ten, tnv, tnve, tt0
            REAL(C_DOUBLE), INTENT(IN), VALUE :: slow1, slow2, slow3, slow4, &
                                                 slow5, slow6, slow7
            REAL(C_DOUBLE), INTENT(IN), VALUE :: sgnrz_dzi, sgnrx_dxi, sgnry_dyi
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, sgntx, sgnty, sgntz
            END FUNCTION

            SUBROUTINE fteik_setGridSizeF(nzIn, nxIn, nyIn, ierr) &
                       BIND(C, NAME='fteik_setGridSizeF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: nzIn, nxIn, nyIn
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_setGridSpacingF(dzIn, dxIn, dyIn, ierr) &
                       BIND(C, NAME='fteik_setGridSpacingF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: dzIn, dxIn, dyIn
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_initializeTravelTimesNearSourceF(ttimes, ierr)      &
                       BIND(C, NAME='fteik_initializeTravelTimesNearSourceF')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_setSzeroF(ierr) BIND(C, NAME='fteik_setSzeroF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(OUT) :: ierr 
            END SUBROUTINE

            SUBROUTINE fteik_setModelOriginF(z0In, x0In, y0In) &
                       BIND(C, NAME='fteik_setModelOriginF')
            USE ISO_C_BINDING
            REAL(C_DOUBLE), INTENT(IN) :: x0In, y0In, z0In
            END SUBROUTINE

            SUBROUTINE fteik_setSourceLocationF(zs, xs, ys, ierr) &
                       BIND(C, NAME='fteik_setSourceLocationF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: xs, ys, zs
            INTEGER(C_INT), INTENT(OUT) :: ierr 
            REAL(C_DOUBLE) xmax, ymax, zmax, xsrc, ysrc, zsrc 
            END SUBROUTINE

            SUBROUTINE fteik_getSweepLimitsF(sweep, linitk,          &
                                             nz, nx, ny,             &
                                             zsi, xsi, ysi,          &
                                             z1, z2, x1, x2, y1, y2, &
                                             ierr)                   &
                       BIND(C, NAME='fteik_getSweepLimitsF')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            INTEGER(C_INT), INTENT(IN) :: sweep, nx, ny, nz, xsi, ysi, zsi
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            INTEGER(C_INT), INTENT(OUT) :: x1, x2, y1, y2, z1, z2, ierr
            END SUBROUTINE

            SUBROUTINE fteik_setUpdateNodesF(sweep, nLevels, linitk,     &
                                             levelPtr, ijkv, lupd, ierr) &
                       BIND(C, NAME='fteik_setUpdateNodesF')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            INTEGER(C_INT), INTENT(IN) :: sweep, nLevels
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            INTEGER(C_INT), DIMENSION(:), INTENT(IN) :: ijkv, levelPtr
            LOGICAL(C_BOOL), DIMENSION(:), INTENT(OUT) :: lupd 
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE fteik_setUpdateNodesF

            SUBROUTINE fteik_setSphericalToCartesianEpsilonF(epsIn, ierr) &
                       BIND(C, NAME='fteik_setSphericalToCartesianEpsilonF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: epsIn
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_setNumberOfSweepsF(nsweepIn, ierr)      &
                       BIND(C, NAME='fteik_setNumberOfSweepsF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: nsweepIn
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_setEpsSolverF(ierr) &
                       BIND(C, NAME='fteik_setEpsSolverF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_getTravelTimes64fF(ng, tt, ierr) &
                       BIND(C, NAME='fteik_getTravelTimes64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            INTEGER(C_INT), INTENT(IN), VALUE :: ng
            REAL(C_DOUBLE), INTENT(OUT) :: tt(ng)
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_setVelocityModel64fF(nv, vel, ierr) &
                       BIND(C, NAME='fteik_setVelocityModel64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: nv
            REAL(C_DOUBLE), INTENT(IN) :: vel(nv)
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_computeGraphF(ierr)           &
                       BIND(C, NAME='fteik_computeGraphF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_setVelocityModel32fF(nv, vel4, ierr) &
                       BIND(C, NAME='fteik_setVelocityModel32fF')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            INTEGER(C_INT), INTENT(IN) :: nv
            REAL(C_FLOAT), INTENT(IN) :: vel4(nv)
            INTEGER(C_INT), INTENT(OUT) :: ierr 
            END SUBROUTINE

            SUBROUTINE fteik_prefetchTravelTimes64fF(i, j, k, nz, nzx,    &
                                                     sgntz, sgntx, sgnty, &
                                                     perm, ttimes, ttLoc) &
                       BIND(C, NAME='fteik_prefetchTravelTimes64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            REAL(C_DOUBLE), DIMENSION(:), INTENT(IN) :: ttimes
            INTEGER(C_INT), DIMENSION(:), INTENT(IN) :: perm
            INTEGER(C_INT), INTENT(IN) ::  i, j, k, nz, nzx, sgntx, sgnty, sgntz
            REAL(C_DOUBLE), INTENT(OUT) :: ttLoc(8)
            END SUBROUTINE

            SUBROUTINE fteik_prefetchSlowness64fF(i, j, k, nzm1, nzm1_nxm1,     &
                                                  sgnvz, sgnvx, sgnvy,          &
                                                  slowPerm, slow, slowStencils) &
                       BIND(C, NAME='fteik_prefetchSlowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            INTEGER(C_INT), INTENT(IN) :: slowPerm(19), i, j, k, nzm1, nzm1_nxm1, &
                                          sgnvz, sgnvx, sgnvy
            REAL(C_DOUBLE), DIMENSION(:), INTENT(IN) :: slow 
            REAL(C_DOUBLE), INTENT(OUT) :: slowStencils(7)
            !DIR ATTRIBUTES ALIGN: 64:: slowStencils
            END SUBROUTINE

            SUBROUTINE fteik_getTravelTimePermF(sweep, perm, ierr) &
                       BIND(C, NAME='fteik_getTravelTimePermF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: sweep
            INTEGER(C_INT), INTENT(OUT) :: perm(8), ierr
            END SUBROUTINE

            SUBROUTINE fteik_getSlownessPermF(sweep, perm, ierr) &
                       BIND(C, NAME='fteik_getSlownessPermF')
            USE ISO_C_BINDING
            INTEGER(C_INT), INTENT(IN) :: sweep
            INTEGER(C_INT), INTENT(OUT) :: perm(19), ierr
            END SUBROUTINE

            SUBROUTINE fteik_getSweepSigns64fF(sweep,               &
                                               sgntz, sgntx, sgnty, &
                                               sgnvz, sgnvx, sgnvy, &
                                               sgnrz, sgnrx, sgnry) &
            BIND(C, NAME='fteik_getSweepSigns64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            INTEGER(C_INT), INTENT(IN) :: sweep
            INTEGER(C_INT), INTENT(OUT) :: sgntz, sgntx, sgnty, sgnvz, sgnvx, sgnvy
            REAL(C_DOUBLE), INTENT(OUT) :: sgnrz, sgnrx, sgnry
            END SUBROUTINE

            PURE INTEGER(C_INT) FUNCTION grid2indexF(i, j, k, nz, nzx)   &
                                BIND(C, NAME='fteik_grid2indexF')
            !$OMP DECLARE SIMD(grid2indexF) UNIFORM(nz, nzx)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nzx
            END FUNCTION

            PURE SUBROUTINE index2gridF(igrd, i, j, k, ierr) &
                            BIND(C, NAME='fteik_index2gridF')
            !$OMP DECLARE SIMD(index2gridF)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: igrd 
            INTEGER(C_INT), INTENT(OUT) :: i, j, k, ierr
            END SUBROUTINE

            PURE INTEGER(C_INT) FUNCTION velGrid2indexF(i, j, k, nzm1, nzm1_nxm1)  &
                                BIND(C, NAME='fteik_velGrid2indexF')
            !$OMP DECLARE SIMD(velGrid2indexF) UNIFORM(nzm1, nzm1_nxm1)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nzm1, nzm1_nxm1
            END FUNCTION

            SUBROUTINE fteik_evaluateLSSweep164fF(linitk,                     &
                                                  ngrd, ncell, nLevels,       &
                                                  levelPtr, ijkv, lupd, slow, &
                                                  ttimes, ierr)               &
                       BIND(C, NAME='fteik_evaluateLSSweep164fF')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            INTEGER(C_INT), INTENT(IN) :: ncell, ngrd, nLevels
            LOGICAL(C_BOOL), INTENT(IN) :: lupd(ngrd)
            INTEGER(C_INT), INTENT(IN) :: levelPtr(nLevels+1), ijkv(4*ngrd)
            REAL(C_DOUBLE), INTENT(IN) :: slow(ncell)
            REAL(C_DOUBLE), INTENT(INOUT) :: ttimes(ngrd)
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE
            SUBROUTINE fteik_evaluateLSSweep264fF(linitk,                     &
                                                  ngrd, ncell, nLevels,       &
                                                  levelPtr, ijkv, lupd, slow, &
                                                  ttimes, ierr)               &
                       BIND(C, NAME='fteik_evaluateLSSweep264fF')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            INTEGER(C_INT), INTENT(IN) :: ncell, ngrd, nLevels
            LOGICAL(C_BOOL), INTENT(IN) :: lupd(ngrd)
            INTEGER(C_INT), INTENT(IN) :: levelPtr(nLevels+1), ijkv(4*ngrd)
            REAL(C_DOUBLE), INTENT(IN) :: slow(ncell)
            REAL(C_DOUBLE), INTENT(INOUT) :: ttimes(ngrd)
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_evaluateLSSweep64fF(linitk, sweep,              &
                                                 ngrd, ncell, nLevels,       &
                                                 levelPtr, ijkv, lupd, slow, &
                                                 ttimes, ierr)               &
                       BIND(C, NAME='fteik_evaluateLSSweep64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            INTEGER(C_INT), INTENT(IN) :: ncell, ngrd, nLevels, sweep
            LOGICAL(C_BOOL), INTENT(IN) :: lupd(ngrd)
            INTEGER(C_INT), INTENT(IN) :: levelPtr(nLevels+1), ijkv(4*ngrd)
            REAL(C_DOUBLE), INTENT(IN) :: slow(ncell)
            REAL(C_DOUBLE), INTENT(INOUT) :: ttimes(ngrd)
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_getTravelTimesAtReceivers64fF(nrecIn, trec, ierr) &
                       BIND(C, NAME='fteik_getTravelTimesAtReceivers64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: nrecIn
            REAL(C_DOUBLE), INTENT(OUT) :: trec(nrecIn)
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_setReceivers64fF(nrecIn, z, x, y, ierr) &
                       BIND(C, NAME='fteik_setReceivers64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            INTEGER(C_INT), INTENT(IN), VALUE :: nrecIn
            REAL(C_DOUBLE), INTENT(IN) :: z(nrecIn), x(nrecIn), y(nrecIn)
            INTEGER(C_INT), INTENT(OUT) :: ierr 
            END SUBROUTINE

         END INTERFACE
      END MODULE FTEIK_UTILS64F

!                                                                                        !
!========================================================================================!
!                                                                                        !
      MODULE FTEIK_AUTOCODE
         INTERFACE
            PURE SUBROUTINE fteik_prefetchSweep1TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                               !ttimes, ttloc) &
                                                               ttimes, tt1, tt2, tt3, tt4, &
                                                               tt5, tt6, tt7, tt8)         &
            BIND(C, NAME='fteik_prefetchSweep1TravelTimes64fF')
            !$OMP DECLARE SIMD(fteik_prefetchSweep1TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx
            REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny)
            REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
            !REAL(C_DOUBLE), INTENT(OUT) :: ttloc(8)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep2TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                               !ttimes, ttloc) &
                                                               ttimes, tt1, tt2, tt3, tt4, &
                                                               tt5, tt6, tt7, tt8)         &
            BIND(C, NAME='fteik_prefetchSweep2TravelTimes64fF')
            !$OMP DECLARE SIMD(fteik_prefetchSweep2TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx 
            REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny)
            REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
            !REAL(C_DOUBLE), INTENT(OUT) :: ttloc(8)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep3TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                               !ttimes, ttloc) &
                                                               ttimes, tt1, tt2, tt3, tt4, &
                                                               tt5, tt6, tt7, tt8)         &
            BIND(C, NAME='fteik_prefetchSweep3TravelTimes64fF')
            !$OMP DECLARE SIMD(fteik_prefetchSweep3TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx 
            REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny)
            REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
            !REAL(C_DOUBLE), INTENT(OUT) :: ttloc(8)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep4TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                               !ttimes, ttloc) &
                                                               ttimes, tt1, tt2, tt3, tt4, &
                                                               tt5, tt6, tt7, tt8)         &
            BIND(C, NAME='fteik_prefetchSweep4TravelTimes64fF')
            !$OMP DECLARE SIMD(fteik_prefetchSweep4TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx 
            REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny)
            REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
            !REAL(C_DOUBLE), INTENT(OUT) :: ttloc(8)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep5TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                               !ttimes, ttloc) &
                                                               ttimes, tt1, tt2, tt3, tt4, &
                                                               tt5, tt6, tt7, tt8)         &
            BIND(C, NAME='fteik_prefetchSweep5TravelTimes64fF')
            !$OMP DECLARE SIMD(fteik_prefetchSweep5TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx 
            REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny)
            REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
            !REAL(C_DOUBLE), INTENT(OUT) :: ttloc(8)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep6TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                               !ttimes, ttloc) &
                                                               ttimes, tt1, tt2, tt3, tt4, &
                                                               tt5, tt6, tt7, tt8)         &
            BIND(C, NAME='fteik_prefetchSweep6TravelTimes64fF')
            !$OMP DECLARE SIMD(fteik_prefetchSweep6TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx
            REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny)
            REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
            !REAL(C_DOUBLE), INTENT(OUT) :: ttloc(8)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep7TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                               !ttimes, ttloc) &
                                                               ttimes, tt1, tt2, tt3, tt4, &
                                                               tt5, tt6, tt7, tt8)         &
            BIND(C, NAME='fteik_prefetchSweep7TravelTimes64fF')
            !$OMP DECLARE SIMD(fteik_prefetchSweep7TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx
            REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny)
            REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
            !REAL(C_DOUBLE), INTENT(OUT) :: ttloc(8)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep8TravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &
                                                               !ttimes, ttloc) &
                                                               ttimes, tt1, tt2, tt3, tt4, &
                                                               tt5, tt6, tt7, tt8)         &
            BIND(C, NAME='fteik_prefetchSweep8TravelTimes64fF')
            !$OMP DECLARE SIMD(fteik_prefetchSweep8TravelTimes64fF) UNIFORM(nz, nx, ny, nzx)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx
            REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny)
            REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8
            !REAL(C_DOUBLE), INTENT(OUT) :: ttloc(8)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep1Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, & !slowLoc) &
                                                             sl1, sl2, sl3, sl4, &
                                                             sl5, sl6, sl7) &
                       BIND(C, NAME='fteik_prefetchSweep1Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            !REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            REAL(C_DOUBLE), INTENT(OUT) :: sl1, sl2, sl3, sl4, sl5, sl6, sl7
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep2Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, & !slowLoc) &
                                                             sl1, sl2, sl3, sl4, &
                                                             sl5, sl6, sl7) &
                       BIND(C, NAME='fteik_prefetchSweep2Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            !REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            REAL(C_DOUBLE), INTENT(OUT) :: sl1, sl2, sl3, sl4, sl5, sl6, sl7
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep3Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, & !slowLoc) &
                                                             sl1, sl2, sl3, sl4, &
                                                             sl5, sl6, sl7) &
                       BIND(C, NAME='fteik_prefetchSweep3Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            !REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            REAL(C_DOUBLE), INTENT(OUT) :: sl1, sl2, sl3, sl4, sl5, sl6, sl7
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep4Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, & !slowLoc) &
                                                             sl1, sl2, sl3, sl4, &
                                                             sl5, sl6, sl7) &
                       BIND(C, NAME='fteik_prefetchSweep4Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            !REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            REAL(C_DOUBLE), INTENT(OUT) :: sl1, sl2, sl3, sl4, sl5, sl6, sl7
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep5Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, & !slowLoc) &
                                                             sl1, sl2, sl3, sl4, &
                                                             sl5, sl6, sl7) &
                       BIND(C, NAME='fteik_prefetchSweep5Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            !REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            REAL(C_DOUBLE), INTENT(OUT) :: sl1, sl2, sl3, sl4, sl5, sl6, sl7
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep6Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, & !slowLoc) &
                                                             sl1, sl2, sl3, sl4, &
                                                             sl5, sl6, sl7) &
                       BIND(C, NAME='fteik_prefetchSweep6Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            !REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            REAL(C_DOUBLE), INTENT(OUT) :: sl1, sl2, sl3, sl4, sl5, sl6, sl7
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep7Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, & !slowLoc) &
                                                             sl1, sl2, sl3, sl4, &
                                                             sl5, sl6, sl7) &
                       BIND(C, NAME='fteik_prefetchSweep7Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            !REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            REAL(C_DOUBLE), INTENT(OUT) :: sl1, sl2, sl3, sl4, sl5, sl6, sl7
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep8Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, & !slowLoc) &
                                                             sl1, sl2, sl3, sl4, &
                                                             sl5, sl6, sl7) &
                       BIND(C, NAME='fteik_prefetchSweep8Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            !REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            REAL(C_DOUBLE), INTENT(OUT) :: sl1, sl2, sl3, sl4, sl5, sl6, sl7
            END SUBROUTINE

            SUBROUTINE fteik_evaluateSweep1LS64fF(linitk, ttimes, ierr) &
                       BIND(C, NAME='fteik_evaluateSweep1LS64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_evaluateSweep2LS64fF(linitk, ttimes, ierr) &
                       BIND(C, NAME='fteik_evaluateSweep2LS64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_evaluateSweep3LS64fF(linitk, ttimes, ierr) &
                       BIND(C, NAME='fteik_evaluateSweep3LS64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_evaluateSweep4LS64fF(linitk, ttimes, ierr) &
                       BIND(C, NAME='fteik_evaluateSweep4LS64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_evaluateSweep5LS64fF(linitk, ttimes, ierr) &
                       BIND(C, NAME='fteik_evaluateSweep5LS64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_evaluateSweep6LS64fF(linitk, ttimes, ierr) &
                       BIND(C, NAME='fteik_evaluateSweep6LS64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_evaluateSweep7LS64fF(linitk, ttimes, ierr) &
                       BIND(C, NAME='fteik_evaluateSweep7LS64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE fteik_evaluateSweep8LS64fF(linitk, ttimes, ierr) &
                       BIND(C, NAME='fteik_evaluateSweep8LS64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            LOGICAL(C_BOOL), INTENT(IN) :: linitk
            REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE

         END INTERFACE
      END MODULE

      MODULE FTEIK_H5IO
         USE ISO_C_BINDING
         !INTEGER(C_INT64_T) h5fl
         CHARACTER(C_CHAR) :: fname(4096)
         LOGICAL(C_BOOL) linitH5FL
         SAVE linitH5FL, fname

         INTERFACE

            INTEGER(C_SIZE_T) FUNCTION strlenF(fileName) &
            BIND(C, NAME='strlen')
            USE ISO_C_BINDING
            CHARACTER(C_CHAR), INTENT(IN) :: fileName(*)
            END FUNCTION

            INTEGER(C_INT) FUNCTION fteik_h5io_initializeF(fileName) &
            BIND(C, NAME='fteik_h5io_initializeF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(C_CHAR), INTENT(IN) :: fileName(*)
            END FUNCTION

            ! internal functions

            SUBROUTINE fteik_h5io_ijkv2levelsF(nz, nx, ny, ngrd, nLevels, &
                                               levelPtr, ijkv, levels)    &
            BIND(C, NAME='fteik_h5io_ijkv2levelsF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: ngrd, nLevels , nz, nx, ny
            INTEGER(C_INT), INTENT(IN) :: levelPtr(nLevels+1), ijkv(4*ngrd)
            INTEGER(C_INT16_T), INTENT(OUT) :: levels(ngrd)
            INTEGER(C_INT) level, indx, ix, iy, iz, jz, node
            END SUBROUTINE 

            INTEGER(C_INT) &
            FUNCTION fteik_h5io_writeLevelScheduleF(h5fl, sweep,   &
                                                    nz, nx, ny,    &
                                                    levelSchedule) &
            BIND(C, NAME='fteik_h5io_writeLevelScheduleF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT64_T), INTENT(IN), VALUE :: h5fl
            INTEGER(C_INT), INTENT(IN), VALUE :: nz, nx, ny, sweep
            INTEGER(C_INT16_T), INTENT(IN) :: levelSchedule(nz*nx*ny)
            END FUNCTION

            INTEGER(C_INT) &
            FUNCTION fteik_h5io_writeTravelTimes32fF(h5fl, ttName, &
                                                     nz, nx, ny,   &
                                                     tt)           &
            BIND(C, NAME='fteik_h5io_writeTravelTimes32fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT64_T), INTENT(IN), VALUE :: h5fl
            CHARACTER(C_CHAR), INTENT(IN) :: ttName(*)
            INTEGER(C_INT), INTENT(IN), VALUE :: nx, ny, nz
            REAL(C_FLOAT), INTENT(IN) :: tt(nx*ny*nz)
            END FUNCTION

            INTEGER(C_INT) &
            FUNCTION fteik_h5io_writeVelocityModel16iF(h5fl, velName, &
                                                       nz, nx, ny,    &
                                                       vel)           &
            BIND(C, NAME='fteik_h5io_writeVelocityModel16iF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT64_T), INTENT(IN), VALUE :: h5fl
            CHARACTER(C_CHAR), INTENT(IN) :: velName(*)
            INTEGER(C_INT), INTENT(IN), VALUE :: nx, ny, nz
            INTEGER(C_INT16_T), INTENT(IN) :: vel((nx-1)*(ny-1)*(nz-1))
            END FUNCTION

            INTEGER(C_INT) FUNCTION fteik_h5io_closeFileF(h5fl) &
            BIND(C, NAME='h5io_closeFileF')
            USE ISO_C_BINDING 
            IMPLICIT NONE
            INTEGER(C_INT64_T), INTENT(IN), VALUE :: h5fl
            END FUNCTION

            INTEGER(C_INT64_T) FUNCTION fteik_h5io_createFileF(fileName) &
            BIND(C, NAME='h5io_createFileF')
            USE ISO_C_BINDING 
            IMPLICIT NONE
            CHARACTER(C_CHAR), INTENT(IN) :: fileName(*)
            END FUNCTION

            INTEGER(C_INT64_T) FUNCTION fteik_h5io_openFileReadWriteF(fileName) &
            BIND(C, NAME='h5io_openFileReadWriteF')
            USE ISO_C_BINDING 
            IMPLICIT NONE
            CHARACTER(C_CHAR), INTENT(IN) :: fileName(*)
            END FUNCTION

            INTEGER(C_INT) FUNCTION fteik_h5io_initialize(fileName) &
            BIND(C, NAME='fteik_h5io_initialize')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(C_CHAR), INTENT(IN) :: fileName(*)
            END FUNCTION

            INTEGER(C_INT) FUNCTION h5io_writeArray64fF(fid, dataName, n, x) &
            BIND(C, NAME='h5io_writeArray64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT64_T), INTENT(IN), VALUE :: fid
            INTEGER(C_INT), INTENT(IN), VALUE :: n
            CHARACTER(C_CHAR), INTENT(IN) :: dataName(*)
            REAL(C_DOUBLE), INTENT(IN) :: x(n)
            END FUNCTION

            INTEGER(C_INT) FUNCTION h5io_writeArray32iF(fid, dataName, n, x) &
            BIND(C, NAME='h5io_writeArray32iF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT64_T), INTENT(IN), VALUE :: fid 
            INTEGER(C_INT), INTENT(IN), VALUE :: n
            CHARACTER(C_CHAR), INTENT(IN) :: dataName(*)
            INTEGER(C_INT), INTENT(IN) :: x(n)
            END FUNCTION

            INTEGER(C_INT) FUNCTION h5io_writeArray16iF(fid, dataName, n, x) &
            BIND(C, NAME='h5io_writeArray16iF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT64_T), INTENT(IN), VALUE :: fid
            INTEGER(C_INT16_T), INTENT(IN), VALUE :: n
            CHARACTER(C_CHAR), INTENT(IN) :: dataName(*)
            REAL(C_INT), INTENT(IN) :: x(n)
            END FUNCTION

         END INTERFACE
      END MODULE 
