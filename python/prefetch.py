#!/usr/bin/python
from numpy import array

def getTTsign(sweep):
   ttSignz = array([1, 1, 1, 1,-1,-1,-1,-1])
   ttSignx = array([1,-1, 1,-1, 1,-1, 1,-1])
   ttSigny = array([1, 1,-1,-1, 1, 1,-1,-1])
   return ttSignz[sweep-1], ttSignx[sweep-1], ttSigny[sweep-1]

def getTTPerm(sweep):
   ttPerm1 = array([7, 5, 6, 3, 4, 2, 1, 8]) - 1
   ttPerm2 = array([6, 3, 7, 5, 1, 8, 4, 2]) - 1
   ttPerm3 = array([4, 2, 1, 8, 7, 5, 6, 3]) - 1
   ttPerm4 = array([1, 8, 4, 2, 6, 3, 7, 5]) - 1
   ttPerm5 = array([5, 7, 3, 6, 2, 4, 8, 1]) - 1
   ttPerm6 = array([3, 6, 5, 7, 8, 1, 2, 4]) - 1
   ttPerm7 = array([2, 4, 8, 1, 5, 7, 3, 6]) - 1
   ttPerm8 = array([8, 1, 2, 4, 3, 6, 5, 7]) - 1
   if (sweep == 1):
      ttPerm = ttPerm1
   elif (sweep == 2):
      ttPerm = ttPerm2
   elif (sweep == 3):
      ttPerm = ttPerm3
   elif (sweep == 4):
      ttPerm = ttPerm4
   elif (sweep == 5):
      ttPerm = ttPerm5
   elif (sweep == 6):
      ttPerm = ttPerm6
   elif (sweep == 7):
      ttPerm = ttPerm7
   elif (sweep == 8):
      ttPerm = ttPerm8
   return ttPerm


def makeSSub(sweep):
   v2l1 = array([1,  5,  9, 13, 15, 17, 19,  6, 11, 18,  3, 10, 16, 12,  2,  7, 14,  8,  4]) - 1
   v2l2 = array([1,  9, 15, 11,  3,  5, 10, 13, 16, 17, 19,  6, 12, 18,  2,  4,  7, 14,  8]) - 1
   v2l3 = array([1,  5, 13,  6,  3,  2,  7,  9, 14, 15, 17, 19,  8, 11, 18,  4, 10, 16, 12]) - 1
   v2l4 = array([1,  3,  5, 13,  6,  2,  9, 15, 11,  4,  7, 10, 14, 16, 17, 19,  8, 12, 18]) - 1
   v2l5 = array([5,  9, 17,  1,  6, 11, 13, 15, 18, 19, 10,  3, 12, 16,  7,  2,  8, 14,  4]) - 1
   v2l6 = array([9,  1, 11, 15,  5, 10, 17,  3,  6, 12, 13, 16, 18, 19,  2,  7,  4,  8, 14]) - 1
   v2l7 = array([5,  1,  6, 13,  3,  7,  9, 17,  2,  8, 11, 14, 15, 18, 19, 10,  4, 12, 16]) - 1
   v2l8 = array([1,  5,  3,  6, 13,  9,  2, 11, 15,  7, 10, 17,  4,  8, 12, 14, 16, 18, 19]) - 1
   vSignz = array([1, 1, 1, 1, 0, 0, 0, 0]) 
   vSignx = array([1, 0, 1, 0, 1, 0, 1, 0]) 
   vSigny = array([1, 1, 0, 0, 1, 1, 0, 0])
   v2fs = []
   v2fs.append('velGrid2indexF(i1,           max_jm1_1,    max_km1_1,    nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(i1,           max_jm1_1,    min_km1_nym1, nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(i1,           min_jm1_nxm1, max_km1_1,    nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(i1,           min_jm1_nxm1, min_km1_nym1, nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(max_im1_1,    j1,           max_km1_1,    nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(min_im1_nzm1, j1,           max_km1_1,    nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(max_im1_1,    j1,           min_km1_nym1, nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(min_im1_nzm1, j1,           min_km1_nym1, nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(max_im1_1,    max_jm1_1,    k1,           nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(max_im1_1,    min_jm1_nxm1, k1,           nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(min_im1_nzm1, max_jm1_1,    k1,           nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(min_im1_nzm1, min_jm1_nxm1, k1,           nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(i1, j1, max_km1_1,    nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(i1, j1, min_km1_nym1, nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(i1, max_jm1_1,    k1, nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(i1, min_jm1_nxm1, k1, nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(max_im1_1,    j1, k1, nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(min_im1_nzm1, j1, k1, nzm1, nzm1_nxm1)')
   v2fs.append('velGrid2indexF(i1, j1, k1, nzm1, nzm1_nxm1)')
   header = '!>    @brief Fetches the slowness corresponding to the stencil at grid\n' \
          + "!>           point (i,j,k) for the %d'th sweep.  Note, that this is machine\n"%(sweep) \
          + '!>           generated code.\n' \
          + '!>\n' \
          + "!>    @param[in] i           iz'th grid point.  This is Fortran indexed.\n" \
          + "!>    @param[in] j           ix'th grid point.  This is Fortran indexed.\n" \
          + "!>    @param[in] k           iy'th grid point.  This is Fortran indexed.\n" \
          + "!>    @param[in] nzm1        nz - 1.\n" \
          + "!>    @param[in] nzm1_nxm1   (nz - 1)*(nx - 1).\n" \
          + "!>    @param[in] slow        This is the slowness (s/m) in each cell.  This is a\n" \
          + "!>                           a vector of dimension [(nz-1) x (nx-1) x (ny-1)].\n" \
          + '!>\n' \
          + "!>    @param[out] slowLoc    This is the slowness stencil for the (iz,ix,iy)'th\n" \
          + "!>                           grid point for use in the local solver.  It is a vector\n" \
          + "!>                           dimension [7].\n" \
          + "!>\n" \
          + "!>    @author Ben Baker\n" \
          + "!>\n" \
          + "!>    @copyright MIT\n" \
          + "!>\n"
   cline = '      PURE SUBROUTINE fteik_prefetchSweep%dSlowness64fF(i, j, k, nz, nx, ny, &\n'%(sweep) \
         + '                                                       nzm1, nzm1_nxm1, slow, & !slowLoc) &\n' \
         + '                                                       slowLoc1, slowLoc2, slowLoc3, slowLoc4, &\n' \
         + '                                                       slowLoc5, slowLoc6, slowLoc7) &\n' \
         + "      BIND(C, NAME='fteik_prefetchSweep%dSlowness64fF')\n"%(sweep)
   cline = cline \
         + '      USE ISO_C_BINDING\n' \
         + '      USE FTEIK_UTILS64F, ONLY : velgrid2indexF\n' \
         + '      IMPLICIT NONE\n' \
         + '      INTEGER(C_INT), INTENT(IN) :: i, j, k, nx, ny, nz, nzm1, nzm1_nxm1\n' \
         + '      REAL(C_DOUBLE), INTENT(IN) :: slow((nz-1)*(nx-1)*(ny-1)) \n' \
         + '      !REAL(C_DOUBLE), INTENT(OUT) :: !slowLoc(7)\n' \
         + '      REAL(C_DOUBLE), INTENT(OUT) :: slowLoc1, slowLoc2, slowLoc3, slowLoc4, &\n' \
         + '                                     slowLoc5, slowLoc6, slowLoc7\n' \
         + '      REAL(C_DOUBLE) sloc(19)\n' \
         + '      !DIR ATTRIBUTES ALIGN: 64:: sloc\n' \
         + '      INTEGER(C_INT) i1, j1, k1,                            &\n' \
         + '                     max_im1_1, max_jm1_1, max_km1_1,       &\n' \
         + '                     min_im1_nzm1, min_jm1_nxm1, min_km1_nym1\n'
   cline = cline + '      INTEGER(C_INT), PARAMETER :: sgnvz = %2d\n'%(vSignz[sweep-1])
   cline = cline + '      INTEGER(C_INT), PARAMETER :: sgnvx = %2d\n'%(vSignx[sweep-1])
   cline = cline + '      INTEGER(C_INT), PARAMETER :: sgnvy = %2d\n'%(vSigny[sweep-1])
   cline = cline + '      i1 = i - sgnvz\n' 
   cline = cline + '      j1 = j - sgnvx\n'
   cline = cline + '      k1 = k - sgnvy\n'
   cline = cline + '      ! Mins and maxes; i think this is superfluous\n'
   cline = cline + '      max_im1_1    = MAX(i - 1, 1)\n'
   cline = cline + '      max_jm1_1    = MAX(j - 1, 1)\n'
   cline = cline + '      max_km1_1    = MAX(k - 1, 1)\n'
   cline = cline + '      min_im1_nzm1 = MIN(i, nz - 1)\n'
   cline = cline + '      min_jm1_nxm1 = MIN(j, nx - 1)\n'
   cline = cline + '      min_km1_nym1 = MIN(k, ny - 1)\n'
   if (sweep == 1): 
      vPerm = v2l1
   elif (sweep == 2): 
      vPerm = v2l2
   elif (sweep == 3): 
      vPerm = v2l3
   elif (sweep == 4): 
      vPerm = v2l4
   elif (sweep == 5): 
      vPerm = v2l5
   elif (sweep == 6): 
      vPerm = v2l6
   elif (sweep == 7): 
      vPerm = v2l7
   elif (sweep == 8): 
      vPerm = v2l8
   for i in xrange(19):
      cline = cline + '      sloc(%2d) = slow(%s)\n'%(vPerm[i]+1, v2fs[vPerm[i]])
   cline = cline \
         + '      slowLoc1 = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))\n' \
         + '      slowLoc2 = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))\n' \
         + '      slowLoc3 = MIN(sloc(9), sloc(10), sloc(11), sloc(12))\n' \
         + '      slowLoc4 = MIN(sloc(13), sloc(14))\n' \
         + '      slowLoc5 = MIN(sloc(15), sloc(16))\n' \
         + '      slowLoc6 = MIN(sloc(17), sloc(18))\n' \
         + '      slowLoc7 = sloc(19)\n' \
         + '      RETURN\n' \
         + '      END SUBROUTINE\n' 
   program = header + cline
   return program
#############################################################################################

def makeTTSub(sweep):
   g2fs = []
   g2fs.append('grid2indexF(i-sgntz, j,       k, nz, nzx)')
   g2fs.append('grid2indexF(i,       j-sgntx, k, nz, nzx)')
   g2fs.append('grid2indexF(i,       j,       k-sgnty, nz, nzx)')
   g2fs.append('grid2indexF(i-sgntz, j-sgntx, k, nz, nzx)')
   g2fs.append('grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx)')
   g2fs.append('grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx)')
   g2fs.append('grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx)')
   g2fs.append('grid2indexF(i,       j,       k, nz, nzx)')
   header = "!>    @brief  Extracts the travel-times for the %d'th sweep.\n"%(sweep) \
          + "!>            Note, that this is machine generated code.\n" \
          + "!>\n" \
          + "!>    @param[in] i        iz'th grid point.  This is Fortran indexed.\n" \
          + "!>    @param[in] j        ix'th grid point.  This is Fortran indexed.\n" \
          + "!>    @param[in] k        iy'th grid point.  This is Fortran indexed.\n" \
          + "!>    @param[in] ttimes   The travel times at all points in the grid.\n" \
          + "!>                        this is a vector of dimension [nz x nx x ny].\n" \
          + "!>\n" \
          + "!>    @param[out] tt1     Travel time at node (i-sgntz,j,k) i.e, tv.\n" \
          + "!>    @param[out] tt2     Travel time at node (i,j-sgntx,k) i.e. te.\n" \
          + "!>    @param[out] tt3     Travel time at node (i,j,k-sgnty) i.e. tn.\n" \
          + "!>    @param[out] tt4     Travel time at node (i-sgntz,j-sgntx,k) i.e. tev.\n" \
          + "!>    @param[out] tt5     Travel time at node (i,j-sgntx,k-sgnty) i.e. ten.\n" \
          + "!>    @param[out] tt6     Travel time at node (i-sgntz,j,k-sgnty) i.e. tnv.\n" \
          + "!>    @param[out] tt7     Travel time at node (i-sgntz,j-sgntx,k-sgnty) i.e. tnve.\n" \
          + "!>    @param[out] tt8     travel time at (i,j,k).\n" \
          + "!>\n" \
          + "!>    @author Ben Baker\n" \
          + "!>\n" \
          + "!>    @copyright MIT\n" \
          + "!>\n"
   cline = "      PURE SUBROUTINE fteik_prefetchSweep%dTravelTimes64fF(i, j, k, nz, nx, ny, nzx,  &\n"%(sweep) \
         + "                                                         ttimes, tt1, tt2, tt3, tt4, &\n" \
         + "                                                         tt5, tt6, tt7, tt8)         &\n" \
         + "      BIND(C, NAME='fteik_prefetchSweep%dTravelTimes64fF')\n"%(sweep) \
         + "      !!$OMP DECLARE SIMD(fteik_prefetchSweep%dTravelTimes64fF) UNIFORM(nz, nx, ny, nzx)\n"%(sweep)
   cline = cline + '      USE ISO_C_BINDING\n'
   cline = cline + '      USE FTEIK_UTILS64F, ONLY : grid2indexF\n'
   cline = cline + '      IMPLICIT NONE\n'
   cline = cline + '      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nx, ny, nzx\n'
   cline = cline + '      REAL(C_DOUBLE), INTENT(IN) :: ttimes(nz*nx*ny) \n'
   cline = cline + '      REAL(C_DOUBLE), INTENT(OUT) :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8\n'
   #cline = cline + '      REAL(C_DOUBLE), INTENT(OUT) :: ttLoc(8)\n'
   ttSignz, ttSignx, ttSigny = getTTsign(sweep)
   cline = cline + '      INTEGER(C_INT), PARAMETER :: sgntz = %2d\n'%(ttSignz)
   cline = cline + '      INTEGER(C_INT), PARAMETER :: sgntx = %2d\n'%(ttSignx)
   cline = cline + '      INTEGER(C_INT), PARAMETER :: sgnty = %2d\n'%(ttSigny) 
   ttPerm = getTTPerm(sweep)
   #for i in xrange(len(ttPerm)):
   #   #if (ttPerm[i] + 1 == 8):
   #   #   continue
   #   cline = cline + '      ttLoc(%d) = ttimes(%s)\n'%(ttPerm[i]+1, g2fs[ttPerm[i]])
   for i in xrange(len(ttPerm)):
       #if (ttPerm[i] + 1 == 8):
       #   continue
       cline = cline + '      tt%d = ttimes(%s)\n'%(ttPerm[i]+1, g2fs[ttPerm[i]]) 
   cline = cline + '      RETURN\n'
   cline = cline + '      END SUBROUTINE\n'
   program = header + cline
   return program
####################################################################################################
def makeLS(sweep):
   ttPerm = getTTPerm(sweep)
   ttSignZ, ttSignX, ttSignY = getTTsign(sweep)
   cline = '      SUBROUTINE fteik_evaluateSweep%dLS64fF(linitk, ttimes, ierr) &\n'%(sweep) \
         + "      BIND(C, NAME='fteik_evaluateSweep%dLS64fF')\n"%(sweep) \
         + '      USE ISO_C_BINDING\n' \
         + '      USE FTEIK_LOCALSOLVER3D64F, ONLY : fteik_localSolver_noInit64fF, &\n' \
         + '                                         fteik_localSolver_init64fF\n' \
         + '      !USE FTEIK_UTILS64F, ONLY : fteik_localSolver64fF, fteik_localSolverExplicit64fF, &\n' \
         + '      !                           fteik_localSolverNoInit64fF, fteik_localSolverInit64fF\n' \
         + '      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep%dSlowness64fF, & \n'%(sweep) \
         + '                                 fteik_prefetchSweep%dTravelTimes64fF\n'%(sweep) \
         + '      USE FTEIK3D_SOLVER64F, ONLY : levelPtr, lupd%d, lupdInit%d, ijkv%d, nLevels\n'%(sweep, sweep, sweep) \
         + '      !USE FTEIK_UTILS64F, ONLY : levelPtr, lupd%d, lupdInit%d, ijkv%d, slow, &\n'%(sweep, sweep, sweep) \
         + '      !                           nLevels\n' \
         + '      USE FTEIK_MODEL64F, ONLY : slow, dx, dy, dz, nx, ny, nz, nzx, nzm1, nzm1_nxm1\n' \
         + '      !                           dx, dy, dz, &\n' \
         + '      !                           nLevels, nx, ny, nz, nzx, nzm1, nzm1_nxm1\n' \
         + '      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, chunkSize, zero\n' \
         + '      IMPLICIT NONE\n' \
         + '      LOGICAL(C_BOOL), INTENT(IN) :: linitk\n' \
         + '      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes\n' \
         + '      INTEGER(C_INT), INTENT(OUT) :: ierr\n' \
         + '      REAL(C_DOUBLE) sgnrz_dzi, sgnrx_dxi, sgnry_dyi\n' \
         + '      INTEGER(C_INT) i, indx, j, k, l1, l2, kndx, knode, level, loop, node, node1, node2\n' \
         + '      REAL(C_DOUBLE), ALLOCATABLE :: slowWork(:), ttWork(:), tt1(:), &\n' \
         + '                                     t1(:), t2(:), t3(:), t4(:), &\n' \
         + '                                     t5(:), t6(:), t7(:), t8(:), &\n' \
         + '                                     s1(:), s2(:), s3(:), s4(:), &\n' \
         + '                                     s5(:), s6(:), s7(:)\n' \
         + '      INTEGER(C_INT), PARAMETER :: sweep = %d\n'%(sweep) \
         + '      INTEGER(C_INT), PARAMETER :: sgntz = %d\n'%(ttSignZ) \
         + '      INTEGER(C_INT), PARAMETER :: sgntx = %d\n'%(ttSignX) \
         + '      INTEGER(C_INT), PARAMETER :: sgnty = %d\n'%(ttSignY) \
         + '      REAL(C_DOUBLE), PARAMETER :: sgnrz = DBLE(sgntz)\n' \
         + '      REAL(C_DOUBLE), PARAMETER :: sgnrx = DBLE(sgntx)\n' \
         + '      REAL(C_DOUBLE), PARAMETER :: sgnry = DBLE(sgnty)\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: t1\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: t2\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: t3\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: t4\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: t5\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: t6\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: t7\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: t8\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: s1\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: s2\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: s3\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: s4\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: s5\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: s6\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: s7\n' \
         + '      !DIR$ ATTRIBUTES ALIGN:64 :: tt1\n' \
         + '      ierr = 0\n' \
         + '      !$OMP PARALLEL DEFAULT(none) &\n' \
         + '      !$OMP SHARED(dx, dy, dz, ijkv%d, levelPtr, linitk, lupd%d, lupdInit%d) &\n'%(sweep,sweep,sweep) \
         + '      !$OMP SHARED(nLevels, nx, ny, nz, nzm1, nzm1_nxm1, nzx) &\n' \
         + '      !$OMP SHARED(slow, ttimes) & \n' \
         + '      !$OMP PRIVATE(sgnrx_dxi, sgnry_dyi, sgnrz_dzi) &\n' \
         + '      !$OMP PRIVATE(i, level, indx, j, k, l1, l2, kndx, knode, loop, node, node1, node2, tt1) &\n' \
         + '      !$OMP PRIVATE(t1, t2, t3, t4, t5, t6, t7, t8, ttWork, slowWork) &\n' \
         + '      !$OMP PRIVATE(s1, s2, s3, s4, s5, s6, s7)\n' \
         + '      ! Initialize\n' \
         + '      ALLOCATE(slowWork(8*chunkSize))\n' \
         + '      ALLOCATE(ttWork(8*chunkSize))\n' \
         + '      ALLOCATE(tt1(chunkSize))\n' \
         + '      ALLOCATE(t1(chunkSize))\n' \
         + '      ALLOCATE(t2(chunkSize))\n' \
         + '      ALLOCATE(t3(chunkSize))\n' \
         + '      ALLOCATE(t4(chunkSize))\n' \
         + '      ALLOCATE(t5(chunkSize))\n' \
         + '      ALLOCATE(t6(chunkSize))\n' \
         + '      ALLOCATE(t7(chunkSize))\n' \
         + '      ALLOCATE(t8(chunkSize))\n' \
         + '      ALLOCATE(s1(chunkSize))\n' \
         + '      ALLOCATE(s2(chunkSize))\n' \
         + '      ALLOCATE(s3(chunkSize))\n' \
         + '      ALLOCATE(s4(chunkSize))\n' \
         + '      ALLOCATE(s5(chunkSize))\n' \
         + '      ALLOCATE(s6(chunkSize))\n' \
         + '      ALLOCATE(s7(chunkSize))\n' \
         + '      slowWork(:) = zero\n' \
         + '      ttWork(:) = zero\n' \
         + '      tt1(:) = FTEIK_HUGE\n' \
         + '      t1(:) = zero\n' \
         + '      t2(:) = zero\n' \
         + '      t3(:) = zero\n' \
         + '      t4(:) = zero\n' \
         + '      t5(:) = zero\n' \
         + '      t6(:) = zero\n' \
         + '      t7(:) = zero\n' \
         + '      t8(:) = zero\n' \
         + '      s1(:) = zero\n' \
         + '      s2(:) = zero\n' \
         + '      s3(:) = zero\n' \
         + '      s4(:) = zero\n' \
         + '      s5(:) = zero\n' \
         + '      s6(:) = zero\n' \
         + '      s7(:) = zero\n' \
         + '      sgnrz_dzi = sgnrz/dz\n' \
         + '      sgnrx_dxi = sgnrx/dx\n' \
         + '      sgnry_dyi = sgnry/dy\n' \
         + '      IF (.NOT.linitk) THEN\n' \
         + '         DO 1 level=1,nLevels\n' \
         + '            l1 = levelPtr(level)\n' \
         + '            l2 = levelPtr(level+1) - 1\n' \
         + '            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)\n' \
         + '            !n2 = MIN(l2 - l1 + 1, chunkSize)\n' \
         + '            !$OMP DO\n' \
         + '            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize\n' \
         + '               node1 = knode !l1 + (knode - 1)*chunkSize\n' \
         + '               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)\n' \
         + '               DO node=node1,node2 !loop=1,n2\n' \
         + '                  loop = node - node1 + 1 !node = knode + loop - 1\n' \
         + '                  IF (lupd%d(node)) THEN\n'%(sweep) \
         + '                     i    = ijkv%d(4*(node-1)+1)\n'%(sweep) \
         + '                     j    = ijkv%d(4*(node-1)+2)\n'%(sweep) \
         + '                     k    = ijkv%d(4*(node-1)+3)\n'%(sweep) \
         + '                     !DIR$ FORCEINLINE\n' \
         + '                     CALL fteik_prefetchSweep%dSlowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &\n'%(sweep) \
         + '                                                           slow, & !slowWork(8*(loop-1)+1))\n' \
         + '                                                           s1(loop), s2(loop), s3(loop), s4(loop), &\n' \
         + '                                                           s5(loop), s6(loop), s7(loop))\n' \
         + '                  ENDIF\n' \
         + '               ENDDO\n' \
         + '               DO node=node1,node2 !loop=1,n2\n' \
         + '                  loop = node - node1 + 1 !node = knode + loop - 1\n' \
         + '                  IF (lupd%d(node)) THEN\n'%(sweep) \
         + '                     i    = ijkv%d(4*(node-1)+1)\n'%(sweep) \
         + '                     j    = ijkv%d(4*(node-1)+2)\n'%(sweep) \
         + '                     k    = ijkv%d(4*(node-1)+3)\n'%(sweep) \
         + '                     kndx = 8*(loop-1)\n' \
         + '                     !DIR$ FORCEINLINE\n' \
         + '                     CALL fteik_prefetchSweep%dTravelTimes64fF(i, j, k, nz, nx, ny, nzx, &\n'%(sweep) \
         + '                                             ttimes, &\n' \
         + '                                             t1(loop), t2(loop), t3(loop), t4(loop), &\n' \
         + '                                             t5(loop), t6(loop), t7(loop), tt1(loop))\n' \
         + '                                             !     ttWork(kndx+1), ttWork(kndx+2), &\n' \
         + '                                             !     ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &\n' \
         + '                                             !     ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))\n' \
         + '                  ENDIF\n' \
         + '               ENDDO\n' \
         + '               loop = node2 - node1 + 1\n' \
         + '               !DIR$ FORCEINLINE\n' \
         + '               CALL fteik_localSolver_noInit64fF(loop,                           &\n' \
         + '                                                t1, t2, t3, t4, t5, t6, t7, & !t8, &\n' \
         + '                                                s1, s2, s3, s4, s5, s6, s7, tt1)\n' \
         + '               !!DIR$ IVDEP !$OMP SIMD !DIR$ IVDEP\n' \
         + '               !DO node=node1,node2 !loop=1,n2\n' \
         + '               !   loop = node - node1 + 1 !node = knode + loop - 1\n' \
         + '               !   IF (lupd%d(node)) THEN\n'%(sweep) \
         + '               !      !i    = ijkv%d(4*(node-1)+1)\n'%(sweep) \
         + '               !      !j    = ijkv%d(4*(node-1)+2)\n'%(sweep) \
         + '               !      !k    = ijkv%d(4*(node-1)+3)\n'%(sweep) \
         + '               !      kndx = 8*(loop - 1)\n' \
         + '               !      !DIR$ FORCEINLINE\n' \
         + '               !      tt1(loop) = fteik_localSolver_noInit64fF( &\n' \
         + '               !                   t1(loop), t2(loop), t3(loop), t4(loop), &\n' \
         + '               !                   t5(loop), t6(loop), t7(loop), t8(loop), &\n' \
         + '               !                   s1(loop), s2(loop), s3(loop), s4(loop), &\n' \
         + '               !                   s5(loop), s6(loop), s7(loop))\n' \
         + '               !                   !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &\n' \
         + '               !                   !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &\n' \
         + '               !                   !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &\n' \
         + '               !                   !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7))\n' \
         + '               !      !!DIR$ FORCEINLINE\n' \
         + '               !      !tt1(loop) = fteik_localSolverExplicit64fF( &\n' \
         + '               !      !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &\n' \
         + '               !      !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &\n' \
         + '               !      !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &\n' \
         + '               !      !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &\n' \
         + '               !      !             .FALSE.,                        &\n' \
         + '               !      !              i, j, k,                       &\n' \
         + '               !      !              sgntz, sgntx, sgnty,           &\n' \
         + '               !      !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)\n' \
         + '               !      !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &\n' \
         + '               !      !                             slowWork(8*(loop-1)+1),  &\n' \
         + '               !      !                             .FALSE.,                        &\n' \
         + '               !      !                             i, j, k,                       &\n' \
         + '               !      !                             sgntz, sgntx, sgnty,           &\n' \
         + '               !      !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)\n' \
         + '               !   ENDIF\n' \
         + '               !ENDDO\n' \
         + '               !$OMP SIMD\n' \
         + '               DO node=node1,node2 !loop=1,n2\n' \
         + '                  loop = node - node1 + 1 !node = knode + loop - 1\n' \
         + '                  IF (lupd%d(node)) THEN\n'%(sweep) \
         + '                     indx = ijkv%d(4*(node-1)+4)\n'%(sweep) \
         + '                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))\n' \
         + '                  ENDIF\n' \
         + '               ENDDO\n' \
         + '            ENDDO ! parallel loop on nodes in level\n' \
         + '            !$OMP END DO\n' \
         + '            !$OMP BARRIER\n' \
         + ' 1       CONTINUE ! loop on levels\n' \
         + '      ELSE\n' \
         + '         DO 11 level=1,nLevels\n' \
         + '            l1 = levelPtr(level)\n' \
         + '            l2 = levelPtr(level+1) - 1\n' \
         + '            !k2 = INT(DBLE(l2 - l1 + 1)/DBLE(chunkSize) + 1.d0)\n' \
         + '            !n2 = MIN(l2 - l1 + 1, chunkSize)\n' \
         + '            !$OMP DO\n' \
         + '            DO knode=l1,l2,chunkSize !1,k2 !,chunkSize !knode=l1,l2!,chunkSize\n' \
         + '               node1 = knode !l1 + (knode - 1)*chunkSize\n' \
         + '               node2 = MIN(l2, node1 + chunkSize - 1) !MIN(l2, node1 + n2 - 1)\n' \
         + '               !n2 = 1 !MIN(l2-l1+1, chunkSize) !1\n' \
         + '               DO node=node1,node2 !loop=1,n2\n' \
         + '                  loop = node - node1 + 1 !node = knode + loop - 1\n' \
         + '                  IF (lupdInit%d(node)) THEN\n'%(sweep) \
         + '                     i    = ijkv%d(4*(node-1)+1)\n'%(sweep) \
         + '                     j    = ijkv%d(4*(node-1)+2)\n'%(sweep) \
         + '                     k    = ijkv%d(4*(node-1)+3)\n'%(sweep) \
         + '                     !DIR$ FORCEINLINE\n' \
         + '                     CALL fteik_prefetchSweep%dSlowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &\n'%(sweep) \
         + '                                                           slow, &!slowWork(8*(loop-1)+1))\n' \
         + '                                                           s1(loop), s2(loop), s3(loop), s4(loop), &\n' \
         + '                                                           s5(loop), s6(loop), s7(loop))\n' \
         + '                  ENDIF\n' \
         + '               ENDDO\n' \
         + '               DO node=node1,node2 !1,n2\n' \
         + '                  loop = node - node1 + 1 !node = knode + loop - 1\n' \
         + '                  IF (lupdInit%d(node)) THEN\n'%(sweep) \
         + '                     i    = ijkv%d(4*(node-1)+1)\n'%(sweep) \
         + '                     j    = ijkv%d(4*(node-1)+2)\n'%(sweep) \
         + '                     k    = ijkv%d(4*(node-1)+3)\n'%(sweep) \
         + '                     !kndx = 8*(loop-1)\n' \
         + '                     !DIR$ FORCEINLINE\n' \
         + '                     CALL fteik_prefetchSweep%dTravelTimes64fF(i, j, k, nz, nx, ny, nzx, &\n'%(sweep) \
         + '                                                   ttimes, &\n' \
         + '                                                   t1(loop), t2(loop), t3(loop), t4(loop), &\n' \
         + '                                                   t5(loop), t6(loop), t7(loop), t8(loop))\n' \
         + '                                                   !ttWork(kndx+1), ttWork(kndx+2), &\n' \
         + '                                                   !ttWork(kndx+3), ttWork(kndx+4), ttWork(kndx+5), &\n' \
         + '                                                   !ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8))\n' \
         + '                  ENDIF\n' \
         + '               ENDDO\n' \
         + '               !!DIR$ IVDEP !$OMP SIMD\n' \
         + '               DO node=node1,node2 !loop=1,n2\n' \
         + '                  loop = node - node1 + 1 !node = knode + loop - 1\n' \
         + '                  IF (lupdInit%d(node)) THEN\n'%(sweep) \
         + '                     i    = ijkv%d(4*(node-1)+1)\n'%(sweep) \
         + '                     j    = ijkv%d(4*(node-1)+2)\n'%(sweep) \
         + '                     k    = ijkv%d(4*(node-1)+3)\n'%(sweep) \
         + '                     !indx = ijkv%d(4*(node-1)+4)\n'%(sweep) \
         + '                     !kndx = 8*(loop - 1)\n' \
         + '                     !tt1(loop) = fteik_localSolverExplicit64fF( &\n' \
         + '                     !              t1(loop), t2(loop), t3(loop), t4(loop), &\n' \
         + '                     !              t5(loop), t6(loop), t7(loop), t8(loop), &\n' \
         + '                     !              s1(loop), s2(loop), s3(loop), s4(loop), &\n' \
         + '                     !              s5(loop), s6(loop), s7(loop), &\n' \
         + '                     !             .TRUE.,                        &\n' \
         + '                     !              i, j, k,                       &\n' \
         + '                     !              sgntz, sgntx, sgnty,           &\n' \
         + '                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)\n' \
         + '                     !DIR$ FORCEINLINE\n' \
         + '                     tt1(loop) = fteik_localSolver_init64fF( &\n' \
         + '                                  !ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &\n' \
         + '                                  !ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &\n' \
         + '                                  !slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &\n' \
         + '                                  !slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &\n' \
         + '                                   t1(loop), t2(loop), t3(loop), t4(loop), &\n' \
         + '                                   t5(loop), t6(loop), t7(loop), t8(loop), &\n' \
         + '                                   s1(loop), s2(loop), s3(loop), s4(loop), &\n' \
         + '                                   s5(loop), s6(loop), s7(loop), &\n' \
         + '                                   i, j, k,                       &\n' \
         + '                                   sgntz, sgntx, sgnty,           &\n' \
         + '                                   sgnrz_dzi, sgnrx_dxi, sgnry_dyi)\n' \
         + '                     !!DIR$ FORCEINLINE\n' \
         + '                     !tt1(loop) = fteik_localSolverExplicit64fF( &\n' \
         + '                     !             ttWork(kndx+1), ttWork(kndx+2), ttWork(kndx+3), ttWork(kndx+4), &\n' \
         + '                     !             ttWork(kndx+5), ttWork(kndx+6), ttWork(kndx+7), ttWork(kndx+8), &\n' \
         + '                     !             slowWork(kndx+1), slowWork(kndx+2), slowWork(kndx+3), slowWork(kndx+4), &\n' \
         + '                     !             slowWork(kndx+5), slowWork(kndx+6), slowWork(kndx+7), &\n' \
         + '                     !             .TRUE.,                        &\n' \
         + '                     !              i, j, k,                       &\n' \
         + '                     !              sgntz, sgntx, sgnty,           &\n' \
         + '                     !              sgnrz_dzi, sgnrx_dxi, sgnry_dyi)\n' \
         + '                     !tt1(loop) = fteik_localSolver64fF(ttWork(8*(loop-1)+1), &\n' \
         + '                     !                             slowWork(8*(loop-1)+1),  &\n' \
         + '                     !                             .TRUE.,                        &\n' \
         + '                     !                             i, j, k,                       &\n' \
         + '                     !                             sgntz, sgntx, sgnty,           &\n' \
         + '                     !                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi)\n' \
         + '                  ENDIF\n' \
         + '               ENDDO\n' \
         + '               !OMP SIMD\n' \
         + '               DO node=node1,node2 !loop=1,n2\n' \
         + '                  loop = node - node1 + 1 !node = knode + loop - 1\n' \
         + '                  IF (lupdInit%d(node)) THEN\n'%(sweep) \
         + '                     indx = ijkv%d(4*(node-1)+4)\n'%(sweep) \
         + '                     ttimes(indx) = tt1(loop) !DMIN1(ttimes(indx), tt1(loop))\n' \
         + '                  ENDIF\n' \
         + '               ENDDO\n' \
         + '            ENDDO ! parallel loop on nodes in level\n' \
         + '            !$OMP END DO\n' \
         + '            !$OMP BARRIER\n' \
         + ' 11      CONTINUE ! loop on levels\n' \
         + '      ENDIF\n' \
         + '      DEALLOCATE(slowWork)\n' \
         + '      DEALLOCATE(ttWork)\n' \
         + '      DEALLOCATE(tt1)\n' \
         + '      DEALLOCATE(t1)\n' \
         + '      DEALLOCATE(t2)\n' \
         + '      DEALLOCATE(t3)\n' \
         + '      DEALLOCATE(t4)\n' \
         + '      DEALLOCATE(t5)\n' \
         + '      DEALLOCATE(t6)\n' \
         + '      DEALLOCATE(t7)\n' \
         + '      DEALLOCATE(t8)\n' \
         + '      DEALLOCATE(s1)\n' \
         + '      DEALLOCATE(s2)\n' \
         + '      DEALLOCATE(s3)\n' \
         + '      DEALLOCATE(s4)\n' \
         + '      DEALLOCATE(s5)\n' \
         + '      DEALLOCATE(s6)\n' \
         + '      DEALLOCATE(s7)\n' \
         + '      !$OMP END PARALLEL\n' \
         + '      RETURN\n' \
         + '      END SUBROUTINE\n'

   return cline

program = '' 
for k in xrange(1,9):
   program = program + makeTTSub(k)
for k in xrange(1,9):
   program = program + makeSSub(k)
#print program
for k in xrange(1,9):
   program = program + makeLS(k)
#print program
of = open('autoCode.f90', 'w')
of.write(program)
of.close()
