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
         ! source variables
         REAL(C_DOUBLE) szero
         INTEGER(C_INT) xsi, ysi, zsi
         ! private variables for solver to compute on startup
         REAL(C_DOUBLE), ALLOCATABLE :: tt1(:), tt2(:), tt3(:), tt4(:), &
                                        tt5(:), tt6(:), tt7(:), tt8(:)
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
!        INTEGER(C_INT), PARAMETER :: l2l1(19) = [1, 15, 11, 19,  2,  8, 16, 18,  3, &
!                                                12,  9, 14,  4, 17,  5, 13,  6, 10,  7]
!        INTEGER(C_INT), PARAMETER :: l2l2(19) = [1, 15,  5, 16,  6, 12, 17, 19,  2, &
!                                                 7,  4, 13,  8, 18,  3,  9, 10, 14, 11]
!        INTEGER(C_INT), PARAMETER :: l2l3(19) = [1,  6,  5, 16,  2,  4,  7, 13,  8, &
!                                                17, 14, 19,  3,  9, 10, 18, 11, 15, 12]
!        INTEGER(C_INT), PARAMETER :: l2l4(19) = [1,  6,  2, 10,  3,  5, 11, 17,  7, &
!                                                12,  9, 18,  4, 13,  8, 14, 15, 19, 16]
!        INTEGER(C_INT), PARAMETER :: l2l5(19) = [4, 16, 12, 19,  1,  5, 15, 17,  2, &
!                                                11,  6, 13,  7, 18,  8, 14,  3,  9, 10]
!        INTEGER(C_INT), PARAMETER :: l2l6(19) = [2, 15,  8, 17,  5,  9, 16, 18,  1, &
!                                                 6,  3, 10, 11, 19,  4, 12,  7, 13, 14]
!        INTEGER(C_INT), PARAMETER :: l2l7(19) = [2,  9,  5, 17,  1,  3,  6, 10,  7, &
!                                                16, 11, 18,  4, 12, 13, 19,  8, 14, 15]
!        INTEGER(C_INT), PARAMETER :: l2l8(19) = [1,  7,  3, 13,  2,  4, 10, 14,  6, &
!                                                11,  8, 15,  5, 16,  9, 17, 12, 18, 19]
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
         INTEGER(C_INT), PARAMETER :: chunkSize = 32
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

            PURE REAL(C_DOUBLE) &
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

            PURE REAL(C_DOUBLE)                                                 &
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
            INTEGER(C_INT), INTENT(IN) :: ng
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
         END INTERFACE
      END MODULE

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
                                                             slow, slowLoc) &
                       BIND(C, NAME='fteik_prefetchSweep1Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep2Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, slowLoc) &
                       BIND(C, NAME='fteik_prefetchSweep2Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep3Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, slowLoc) &
                       BIND(C, NAME='fteik_prefetchSweep3Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep4Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, slowLoc) &
                       BIND(C, NAME='fteik_prefetchSweep4Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep5Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, slowLoc) &
                       BIND(C, NAME='fteik_prefetchSweep5Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep6Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, slowLoc) &
                       BIND(C, NAME='fteik_prefetchSweep6Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep7Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, slowLoc) &
                       BIND(C, NAME='fteik_prefetchSweep7Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
            END SUBROUTINE

            PURE SUBROUTINE fteik_prefetchSweep8Slowness64fF(i, j, k, nz, nx, ny, &
                                                             nzm1, nzm1_nxm1, &
                                                             slow, slowLoc) &
                       BIND(C, NAME='fteik_prefetchSweep8Slowness64fF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: i, j, k, nzm1, nzm1_nxm1, nz, nx, ny
            REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1))
            REAL(C_DOUBLE), INTENT(OUT) :: slowLoc(7)
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
