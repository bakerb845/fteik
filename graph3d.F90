MODULE FTEIK_GRAPH3D
   USE ISO_C_BINDING
   IMPLICIT NONE
   !> Maps from start index of level'th level to corresponding node.
   INTEGER(C_INT), PROTECTED, SAVE, ALLOCATABLE :: levelPtr(:)
   !> Maps from the n'th node to the (iz, ix, iy, indx) for the first sweep.  This is
   !> an array of dimension [4 x ngrd] with leading dimension 4.
   INTEGER(C_INT), PROTECTED, SAVE, ALLOCATABLE :: ijkv1(:)
   !> Maps from the n'th node to the (iz, ix, iy, indx) for the second sweep.  This is
   !> an array of dimension [4 x ngrd] with leading dimension 4.
   INTEGER(C_INT), PROTECTED, SAVE, ALLOCATABLE :: ijkv2(:)
   !> Maps from the n'th node to the (iz, ix, iy, indx) for the third sweep.  This is
   !> an array of dimension [4 x ngrd] with leading dimension 4.
   INTEGER(C_INT), PROTECTED, SAVE, ALLOCATABLE :: ijkv3(:)
   !> Maps from the n'th node to the (iz, ix, iy, indx) for the fourth sweep.  This is
   !> an array of dimension [4 x ngrd] with leading dimension 4.
   INTEGER(C_INT), PROTECTED, SAVE, ALLOCATABLE :: ijkv4(:)
   !> Maps from the n'th node to the (iz, ix, iy, indx) for the fifth sweep.  This is
   !> an array of dimension [4 x ngrd] with leading dimension 4.
   INTEGER(C_INT), PROTECTED, SAVE, ALLOCATABLE :: ijkv5(:)
   !> Maps from the n'th node to the (iz, ix, iy, indx) for the first sweep.  This is
   !> an array of dimension [4 x ngrd] with leading dimension 6.
   INTEGER(C_INT), PROTECTED, SAVE, ALLOCATABLE :: ijkv6(:)
   !> Maps from the n'th node to the (iz, ix, iy, indx) for the first sweep.  This is
   !> an array of dimension [4 x ngrd] with leading dimension 7.
   INTEGER(C_INT), PROTECTED, SAVE, ALLOCATABLE :: ijkv7(:)
   !> Maps from the n'th node to the (iz, ix, iy, indx) for the first sweep.  This is
   !> an array of dimension [4 x ngrd] with leading dimension 8.
   INTEGER(C_INT), PROTECTED, SAVE, ALLOCATABLE :: ijkv8(:)
   !> Contains the level of the node'th node.  This is an array of dimension [ngrd].
   INTEGER(C_INT), PROTECTED, SAVE, ALLOCATABLE :: node2level(:)

   INTEGER(C_INT), PROTECTED, SAVE :: nLevels = 0 !< Number of levels in solver.
   INTEGER(C_INT), PROTECTED, SAVE :: maxLevelSize = 0 !< Max level size.
 
   INTEGER(C_INT), PRIVATE, SAVE :: nz = 0 !< Number of z grid points.
   INTEGER(C_INT), PRIVATE, SAVE :: nx = 0 !< Number of x grid points.
   INTEGER(C_INT), PRIVATE, SAVE :: ny = 0 !< Number of y grid points.
   INTEGER(C_INT), PRIVATE, SAVE :: ngrd = 0 !< Number of grid points.
   INTEGER(C_INT32_T), PRIVATE, PARAMETER :: SORT_ASCENDING = 0

   INTERFACE
      INTEGER(C_INT32_T) FUNCTION sorting_sort32i_finter(n, a, order) &
      BIND(C, NAME='sorting_sort32i_finter')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT32_T), VALUE, INTENT(IN) :: n, order
      INTEGER(C_INT32_T), INTENT(INOUT) :: a(n)
      END FUNCTION
   END INTERFACE
   CONTAINS
!----------------------------------------------------------------------------------------!
!                                   Begin the Code                                       !
!----------------------------------------------------------------------------------------!
!
!>    @brief Initializes the level set graph structure.
!>
!>    @param[in] nzIn    Number of grid points in z.  This must be greater than 2.
!>    @param[in] nxIn    Number of grid points in x.  This must be greater than 2.
!>    @param[in] nyIn    Number of grid points in y.  This must be greater than 2.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_graph3d_initializeF(nzIn, nxIn, nyIn, ierr) &
      BIND(C, NAME='fteik_graph3d_initializeF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nzIn, nxIn, nyIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      IF (nzIn < 3 .OR. nxIn < 3 .OR. nyIn < 3) THEN
         WRITE(*,*) 'fteik_graph3d_initializeF: Insufficient number of grid points'
         ierr =-1
         RETURN
      ENDIF
      nz = nzIn
      nx = nxIn
      ny = nyIn
      ngrd = nz*nx*ny
      RETURN
      END SUBROUTINE  
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the level set pointer and node ordering for each sweep.
!>
!>    @param[out] ierr    0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_graph3d_makeLevelStructuresF(ierr) &
      BIND(C, NAME='fteik_graph3d_makeLevelStructuresF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) ierrLoc, ijk, ix, iy, iz, i1, i2, j1, j2, k1, k2, level, node, np, &
                     sweep
      INTEGER, ALLOCATABLE :: linit(:), work(:)
      ierr = 0
      IF (nz < 1) THEN
         WRITE(*,*) 'fteik_graph3d_makeLevelStructuresF: Graph not initalized'
         ierr = 1
         RETURN
      ENDIF
      ! Compute some sizes and set space
      nLevels = nx + ny + nz - 2 ! Number of levels in 3D
      ALLOCATE(ijkv1(4*ngrd))
      ALLOCATE(ijkv2(4*ngrd))
      ALLOCATE(ijkv3(4*ngrd))
      ALLOCATE(ijkv4(4*ngrd))
      ALLOCATE(ijkv5(4*ngrd))
      ALLOCATE(ijkv6(4*ngrd))
      ALLOCATE(ijkv7(4*ngrd))
      ALLOCATE(ijkv8(4*ngrd))
      ALLOCATE(node2level(ngrd))
      ALLOCATE(linit(ngrd))
      ALLOCATE(levelPtr(nLevels+1))
      ALLOCATE(work(ngrd))
      node2level(:) = 0
      ijkv1(:) = 0
      linit(:) = 0
      levelPtr(:) = 0
      ! Loop on the levels
      levelPtr(1) = 1
      DO level=2,nx+ny+nz-1
         k1 = MAX(0, level - 2 - (nx - 1) - (nz - 1))
         k1 = MIN(k1, ny - 1)
         k2 = MIN(ny - 1, level - 2)
         np = 0
         DO iy=k1,k2
            j1 = MAX(0, level - 2 - iy - (nz - 1))
            j1 = MIN(j1, nx - 1)
            j2 = MIN(nx - 1, level - 2 - iy)
            DO ix=j1,j2
               i1 = MAX(0, level - 2 - iy - ix)
               i1 = MIN(i1, nz - 1)
               i2 = MIN(nz - 1, level - 2 - iy - ix)
               DO iz=i1,i2
                  np = np + 1
                  ijk = iy*nx*nz + ix*nz + iz + 1 ! C to Fortran numbering
                  work(np) = ijk
                  node2level(ijk) = level - 1
                  linit(ijk) = linit(ijk) + 1 
               ENDDO ! Loop on z
            ENDDO ! Loop on x
         ENDDO ! Loop on y
         ! Update level size and pointer
         levelPtr(level) = levelPtr(level-1) + np
         maxLevelSize = MAX(maxLevelSize, np)
         ! Sort the nodes in this level
         IF (np > 1) THEN
            ierr = sorting_sort32i_finter(np, work, SORT_ASCENDING)
            IF (ierr /= 0) THEN
               WRITE(*,*) 'fteik_graph3d_makeLevelStructuresF: Sort failed'
               EXIT
            ENDIF
         ENDIF
         ! Insert the nodes into ijkv
         ierr = 0
         !$OMP SIMD REDUCTION(+:ierr)
         DO np=levelPtr(level-1),levelPtr(level)-1
            node = np + 1 - levelPtr(level-1)
            CALL fteik_graph3d_index2gridF(work(node), iz, ix, iy, ierrLoc)
            IF (ierrLoc /= 0) ierr = 1
            ijkv1(4*(np-1)+1) = iz
            ijkv1(4*(np-1)+2) = ix
            ijkv1(4*(np-1)+3) = iy
            ijkv1(4*(np-1)+4) = work(node)
         ENDDO
         IF (ierr /= 0) THEN
            WRITE(*,*) 'fteik_graph3d_index2gridF: Failed to convert grid to (iz,ix,iy)'
            EXIT
         ENDIF 
      ENDDO
      IF (MINVAL(linit) < 1 .OR. MAXVAL(linit) > 1) THEN
         WRITE(*,*) 'fteik_graph3d_makeLevelStructuresF: Algorithm failure', &
         MINVAL(linit), MAXVAL(linit)
         ierr = 1
      ENDIF
      IF (MINVAL(ijkv1) < 1 .OR. MAXVAL(ijkv1) > ngrd) THEN
         WRITE(*,*) 'fteik_graph3d_makeLevelStructuresF: Failed to make ijkv1'
         ierr = 1
      ENDIF
      IF (MINVAL(node2level) < 1 .OR. MAXVAL(node2level) > nLevels) THEN
         WRITE(*,*) 'fteik_graph3d_makeLevelStructuresF: Failed to map node2leve'
         ierr = 1
      ENDIF
      IF (ALLOCATED(linit)) DEALLOCATE(linit)
      IF (ALLOCATED(work))  DEALLOCATE(work)
      ! Generate the IJKV orderings
      !$OMP PARALLEL DO SHARED(ijkv2, ijkv3, ijkv4, ijkv5, ijkv6, ijkv7, ijkv8) &
      !$OMP PRIVATE(ierr, sweep) &
      !$OMP DEFAULT(NONE)
      DO sweep=2,8
         IF (sweep == 2) CALL fteik_graph3d_generateIJKV(2, ijkv2, ierr)
         IF (sweep == 3) CALL fteik_graph3d_generateIJKV(3, ijkv3, ierr)
         IF (sweep == 4) CALL fteik_graph3d_generateIJKV(4, ijkv4, ierr)
         IF (sweep == 5) CALL fteik_graph3d_generateIJKV(5, ijkv5, ierr)
         IF (sweep == 6) CALL fteik_graph3d_generateIJKV(6, ijkv6, ierr)
         IF (sweep == 7) CALL fteik_graph3d_generateIJKV(7, ijkv7, ierr)
         IF (sweep == 8) CALL fteik_graph3d_generateIJKV(8, ijkv8, ierr)
      ENDDO
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Generates the IJKV orderings for the given sweep.
!>
!>    @param[in] sweep     Sweep number.  This must be in the range [1,8].
!>    @param[out] ijkvOut  IJKV ordering for the given sweep.   This is an array
!>                         of dimension [4 x ngrd] with leading dimension 4.
!>
!>    @param[out] ierr     0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_graph3d_generateIJKV(sweep, ijkvOut, ierr) &
      BIND(C, NAME='fteik_graph3d_generateIJKV') 
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: sweep
      INTEGER(C_INT), INTENT(OUT) :: ijkvOut(4*ngrd)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT32_T), ALLOCATABLE :: work(:)
      INTEGER ierrLoc, il, indx, iz, ix, iy, j, jz, jx, jy, node, nsort
      ierr = 0
      IF (ngrd < 1 .OR. .NOT.ALLOCATED(levelPtr)) THEN
         WRITE(*,*) 'fteik_graph3d_generateIJKV: Invalid initialization'
         ierr = 1
         RETURN
      ENDIF
      IF (sweep < 1 .OR. sweep > 8) THEN
         WRITE(*,*) 'fteik_graph3d_generateIJKV: Invalid sweep:', sweep
         ierr = 1
         RETURN
      ENDIF
      ALLOCATE(work(maxLevelSize))
      ijkvOut(:) = 0
      ierr = 0
      DO il=1,nLevels
         ! Permute the nodes in this level
         j = 0
         DO node=levelPtr(il),levelPtr(il+1)-1
            iz = ijkv1(4*(node-1)+1)
            ix = ijkv1(4*(node-1)+2)
            iy = ijkv1(4*(node-1)+3)
            indx = ijkv1(4*(node-1)+4)
            IF (node2level(indx) == il) THEN
               IF (sweep == 1) THEN
                  jz = iz
                  jx = ix
                  jy = iy
               ELSEIF (sweep == 2) THEN
                  jz = iz
                  jx = nx + 1 - ix
                  jy = iy
               ELSEIF (sweep == 3) THEN
                  jz = iz
                  jx = ix
                  jy = ny + 1 - iy
               ELSEIF (sweep == 4) THEN
                  jz = iz
                  jx = nx + 1 - ix
                  jy = ny + 1 - iy
               ELSEIF (sweep == 5) THEN
                  jz = nz + 1 - iz
                  jx = ix
                  jy = iy
               ELSEIF (sweep == 6) THEN
                  jz = nz + 1 - iz
                  jx = nx + 1 - ix
                  jy = iy
               ELSEIF (sweep == 7) THEN
                  jz = nz + 1 - iz
                  jx = ix
                  jy = ny + 1 - iy
               ELSE !if (sweep == 8)
                  jz = nz + 1 - iz
                  jx = nx + 1 - ix
                  jy = ny + 1 - iy
               ENDIF 
               j = j + 1
               work(j) =  (jy - 1)*nx*nz + (jx - 1)*nz + jz
            ENDIF ! end check on level match
         ENDDO ! loop on grid points
         ! Sort the nodes in this level
         nsort = j
         IF (levelPtr(il+1) - levelPtr(il) /= nsort) THEN
            WRITE(*,*) 'fteik_graph3d_generateIJKV: levelPtr not set correctly', &
                       nsort, levelPtr(il+1) - levelPtr(il)
            ierr = 1
            EXIT
         ENDIF
         IF (nsort > 1) THEN
            ierr = sorting_sort32i_finter(nsort, work, SORT_ASCENDING)
            IF (ierr /= 0) THEN
               WRITE(*,*) 'fteik_graph3d_generateIJKV: Failed to sort work'
               EXIT
            ENDIF
         ENDIF
         ! Put the sorted nodes into the ijkv array
         ierr = 0
         DO j=1,nsort
            indx = levelPtr(il) - 1 + j
            CALL fteik_graph3d_index2gridF(work(j), iz, ix, iy, ierrLoc)
            IF (ierrLoc /= 0) ierr = 1 
            ijkvOut(4*(indx-1)+1) = iz
            ijkvOut(4*(indx-1)+2) = ix
            ijkvOut(4*(indx-1)+3) = iy
            ijkvOut(4*(indx-1)+4) = work(j)
         ENDDO
      ENDDO ! loop on levels
      IF (ALLOCATED(work)) DEALLOCATE(work)
      IF (MINVAL(ijkvOut) < 1 .OR. MAXVAL(ijkvOut) > ngrd) THEN
         WRITE(*,*) 'fteik_graph3d_generateIJKV: Failed to make ijkv for sweep', sweep
         ierr = 1
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Frees/releases memory for the 3D graph.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_graph3d_finalizeF( ) &
      BIND(C, NAME='fteik_graph3d_finalizeF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      IF (ngrd < 1) THEN
         WRITE(*,*) 'fteik_graph3d_finalizeF: Graph never intialized'
         RETURN
      ENDIF
      nz = 0
      nx = 0
      ny = 0
      ngrd = 0
      nLevels = 0
      maxLevelSize = 0
      IF (ALLOCATED(node2level)) DEALLOCATE(node2level)
      IF (ALLOCATED(levelPtr))   DEALLOCATE(levelPtr)
      IF (ALLOCATED(ijkv1))      DEALLOCATE(ijkv1)
      IF (ALLOCATED(ijkv2))      DEALLOCATE(ijkv2)
      IF (ALLOCATED(ijkv3))      DEALLOCATE(ijkv3)
      IF (ALLOCATED(ijkv4))      DEALLOCATE(ijkv4)
      IF (ALLOCATED(ijkv5))      DEALLOCATE(ijkv5)
      IF (ALLOCATED(ijkv6))      DEALLOCATE(ijkv6)
      IF (ALLOCATED(ijkv7))      DEALLOCATE(ijkv7)
      IF (ALLOCATED(ijkv8))      DEALLOCATE(ijkv8)
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Converts from the grid index in the travel time array to the (i, j, k)'th
!>           grid indices.
!>
!>    @param[in] igrd   Fortran indxed grid point in travel time mesh to convert.
!>
!>    @param[out] i     iz'th grid point.  This is Fortran indexed.
!>    @param[out] j     ix'th grid point.  This is Fortran indexed.
!>    @param[out] k     iy'th grid point.  This is Fortran indexed.
!>    @param[out] ierr  0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      PURE SUBROUTINE fteik_graph3d_index2gridF(igrd, i, j, k, ierr) &
      BIND(C, NAME='fteik_graph3d_index2gridF')
      !$OMP DECLARE SIMD(fteik_graph3d_index2gridF) UNIFORM(ierr)
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(IN), VALUE :: igrd
      INTEGER(C_INT), INTENT(OUT) :: i, j, k, ierr
      INTEGER(C_INT) igrd1, nzx
      igrd1 = igrd - 1 ! F to C
      nzx = nz*nx
      k = (igrd1)/nzx
      j = (igrd1 - k*nzx)/nz
      i =  igrd1 - k*nzx - j*nz
      k = k + 1 !C to F
      j = j + 1 !C to F
      i = i + 1 !C to F
      ierr = 0 
      IF (i < 1 .OR. i > nz) ierr = ierr + 1 
      IF (j < 1 .OR. j > nx) ierr = ierr + 1 
      IF (k < 1 .OR. k > ny) ierr = ierr + 1 
      RETURN
      END
END MODULE
