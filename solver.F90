!>    @brief This program initializes the solver.  It will also perform the graph
!>           reordering.
!>
!>    @param[in] nzIn      Number of z grid points in model.  This must be at least 3.
!>    @param[in] nxIn      Number of x grid points in model.  This must be at least 3.
!>    @param[in] nyIn      Number of y grid points in model.  This must be at least 3.
!>    @param[in] z0In      z0 model origin (meters).
!>    @param[in] x0In      x0 model origin (meters).
!>    @param[in] y0In      y0 model origin (meters).
!>    @param[in] dzIn      Grid spacing (meters) in z.  This must be positive.
!>    @param[in] dxIn      Grid spacing (meters) in x.  This must be positive.
!>    @param[in] dyIn      Grid spacing (meters) in y.  This must be positive.
!>    @param[in] nsweepIn  This is the number of sweeps or iterations of the 
!>                         Gauss-Seidel method to be performed.  One is generally
!>                         sufficient.
!>    @param[in] epsIn     This is the radius, in number of grid points, around source
!>                         where the spherical approximation finite difference stencils
!>                         will be used.  This cannot be negative. 
!>
!>    @param[out] ierr     0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solver_initialize64fF(nzIn, nxIn, nyIn,      &
                                             z0In, x0In, y0In,      &
                                             dzIn, dxIn, dyIn,      &
                                             nsweepIn, epsIn, ierr) &
                 BIND(C, NAME='fteik_solver_initialize64fF')
      USE ISO_C_BINDING
      USE FTEIK_MODEL64F, ONLY : fteik_model_intializeGeometryF
      USE FTEIK_SOLVER64F, ONLY : fteik_solver_setSphereToCartEpsilonF, &
                                  fteik_solver_setNumberOfSweepsF, &
                                  fteik_solver_computeGraphF, &
                                  fteik_solver_finalizeF
      IMPLICIT NONE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x0In, y0In, z0In
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dxIn, dyIn, dzIn
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: epsIn
      INTEGER(C_INT), VALUE, INTENT(IN) :: nxIn, nyIn, nzIn
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsweepIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      WRITE(*,*) 'fteik_solver_initialize64fF: Initializing...'
      CALL fteik_solver_finalizeF()
      CALL fteik_model_intializeGeometryF(nzIn, nxIn, nyIn, &
                                          dzIn, dxIn, dyIn, &
                                          z0In, x0In, y0In, &
                                          ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver_initialize64fF: Error setting grid geometry'
         RETURN
      ENDIF
      CALL fteik_solver_setNumberOfSweepsF(nsweepIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver_initialize64fF: Failed to set max number of sweeps'
         RETURN
      ENDIF
      CALL fteik_solver_setSphereToCartEpsilonF(epsIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver_initialize64fF: Failed to set epsilon'
         RETURN
      ENDIF
      CALL fteik_solver_computeGraphF(ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver_initialize64fF: Failed to compute graph'
         RETURN
      ENDIF
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the graph and pointers to be used by the level schedule solver.
!>
!>    @result 0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solver_computeGraphF(ierr)          &
                 BIND(C, NAME='fteik_solver_computeGraphF') 
      USE ISO_C_BINDING
      USE FTEIK_MODEL64F, ONLY : ngrd, nx, ny, nz
      USE FTEIK_SOLVER64F, ONLY : nLevels, maxLevelSize
      USE FTEIK_SOLVER64F, ONLY : ijkv1, ijkv2, ijkv3, ijkv4, ijkv5, ijkv6, ijkv7, ijkv8
      USE FTEIK_SOLVER64F, ONLY : lupd1, lupd2, lupd3, lupd4, lupd5, lupd6, lupd7, lupd8
      USE FTEIK_SOLVER64F, ONLY : levelPtr
!     USE FTEIK_SOLVER64F, ONLY : tt1 
      USE FTEIK_SOLVER64F, ONLY : fteik_solver_setUpdateNodesF
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE FTEIK_GRAPH
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: ierr
      TYPE(GRAPH_TYPE) graph
      INTEGER(C_INT) ierrs(8), i
      ierr = 0
      IF (nx < 1 .OR. ny < 1 .OR. nz < 1 .OR. ngrd < 1) THEN
         WRITE(*,*) 'fteik_solver_computeGraph: Error grid not yet set'
         ierr = 1
         RETURN
      ENDIF
      CALL INIT(graph, nz, nx, ny, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver_computeGraphF: Error initializing graph'
         ierr = 1
         RETURN
      ENDIF
      nLevels = NUMBER_OF_LEVELS(graph)
      IF (.NOT.ALLOCATED(ijkv1)) ALLOCATE(ijkv1(4*ngrd))
      IF (.NOT.ALLOCATED(ijkv2)) ALLOCATE(ijkv2(4*ngrd))
      IF (.NOT.ALLOCATED(ijkv3)) ALLOCATE(ijkv3(4*ngrd))
      IF (.NOT.ALLOCATED(ijkv4)) ALLOCATE(ijkv4(4*ngrd))
      IF (.NOT.ALLOCATED(ijkv5)) ALLOCATE(ijkv5(4*ngrd))
      IF (.NOT.ALLOCATED(ijkv6)) ALLOCATE(ijkv6(4*ngrd))
      IF (.NOT.ALLOCATED(ijkv7)) ALLOCATE(ijkv7(4*ngrd))
      IF (.NOT.ALLOCATED(ijkv8)) ALLOCATE(ijkv8(4*ngrd))
      ierr = IJKV(graph, 1, ijkv1)
      ierr = IJKV(graph, 2, ijkv2) + ierr
      ierr = IJKV(graph, 3, ijkv3) + ierr
      ierr = IJKV(graph, 4, ijkv4) + ierr
      ierr = IJKV(graph, 5, ijkv5) + ierr
      ierr = IJKV(graph, 6, ijkv6) + ierr
      ierr = IJKV(graph, 7, ijkv7) + ierr
      ierr = IJKV(graph, 8, ijkv8) + ierr
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver_computeGraphF: Errors encountered getting ijkv'
         ierr = 1
         RETURN
      ENDIF
      ierr = LEVEL_PTR(graph, 1, levelPtr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver_computeGraphF: Error encountered getting levelPtr'
         ierr = 1
         RETURN
      ENDIF
      ! finished with the graph
      CALL DELETE(graph)
      maxLevelSize = 0
      DO 1 i=1,nLevels
         maxLevelSize = MAX(maxLevelSize, levelPtr(i+1) - levelPtr(i))
    1 CONTINUE
!     ALLOCATE(tt1(maxLevelSize))
!     tt1(:) = zero
      ! Set the update grid.  This will mask nodes on the boundary.
      IF (.NOT.ALLOCATED(lupd1)) ALLOCATE(lupd1(ngrd))
      IF (.NOT.ALLOCATED(lupd2)) ALLOCATE(lupd2(ngrd))
      IF (.NOT.ALLOCATED(lupd3)) ALLOCATE(lupd3(ngrd))
      IF (.NOT.ALLOCATED(lupd4)) ALLOCATE(lupd4(ngrd))
      IF (.NOT.ALLOCATED(lupd5)) ALLOCATE(lupd5(ngrd))
      IF (.NOT.ALLOCATED(lupd6)) ALLOCATE(lupd6(ngrd))
      IF (.NOT.ALLOCATED(lupd7)) ALLOCATE(lupd7(ngrd))
      IF (.NOT.ALLOCATED(lupd8)) ALLOCATE(lupd8(ngrd))
      lupd1(:) = .FALSE.
      lupd2(:) = .FALSE.
      lupd3(:) = .FALSE.
      lupd4(:) = .FALSE.
      lupd5(:) = .FALSE.
      lupd6(:) = .FALSE.
      lupd7(:) = .FALSE.
      lupd8(:) = .FALSE.
      CALL fteik_solver_setUpdateNodesF(1, nlevels, .FALSE., 0, levelPtr, &
                                        ijkv1, lupd1, ierrs(1))
      CALL fteik_solver_setUpdateNodesF(2, nlevels, .FALSE., 0, levelPtr, &
                                        ijkv2, lupd2, ierrs(2))
      CALL fteik_solver_setUpdateNodesF(3, nlevels, .FALSE., 0, levelPtr, &
                                        ijkv3, lupd3, ierrs(3))
      CALL fteik_solver_setUpdateNodesF(4, nlevels, .FALSE., 0, levelPtr, &
                                        ijkv4, lupd4, ierrs(4))
      CALL fteik_solver_setUpdateNodesF(5, nlevels, .FALSE., 0, levelPtr, &
                                        ijkv5, lupd5, ierrs(5))
      CALL fteik_solver_setUpdateNodesF(6, nlevels, .FALSE., 0, levelPtr, &
                                        ijkv6, lupd6, ierrs(6))
      CALL fteik_solver_setUpdateNodesF(7, nlevels, .FALSE., 0, levelPtr, &
                                        ijkv7, lupd7, ierrs(7))
      CALL fteik_solver_setUpdateNodesF(8, nlevels, .FALSE., 0, levelPtr, &
                                        ijkv8, lupd8, ierrs(8))
      IF (MAXVAL(ABS(ierrs)) /= 0) THEN
         WRITE(*,*) 'fteik_solver_computeGraphF: Error setting update nodes'
         ierr = 1
      ENDIF
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility routine for setting the nodes to update in the sweep.
!>
!>    @param[in] sweep      Sweep number.  This must be in the range [1,8].
!>    @param[in] nLevels    Number of levels in level set method.
!>    @param[in] linitk     If true then this sets the initialization nodes for the
!>                          given source.
!>                          Otherwise, this sets the initialization nodes for the
!>                          general solver.
!>    @param[in] isrc       Source number.
!>    @param[in] levelPtr   Maps from level'th level to start node in ijkv.
!>                          This is a vector of dimension [nLevels+1].
!>    @param[in] ijkv       Maps from in'th node in to the (iz,ix,iy,igrd) indices.
!>                          This is a vector of length [4*ngrd] where ngrd is the
!>                          number of grid points.
!>
!>    @param[in,out] lupd   On input this contains the nodes to update.  \n
!>                          On output this contains a new list of nodes to update based
!>                          on the sweep and potentially the source information.
!>
!>    @param[out] ierr      0 indicates success. 
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>    
      SUBROUTINE fteik_solver_setUpdateNodesF(sweep, nLevels, linitk, isrc, &
                                              levelPtr, ijkv, lupd, ierr)   &
                 BIND(C, NAME='fteik_solver_setUpdateNodesF')
      USE ISO_C_BINDING
      USE FTEIK_MODEL64F, ONLY : nz, nx, ny
      USE FTEIK_SOURCE64F, ONLY : fteik_source_getSourceIndices32iF
      USE FTEIK_SOLVER64F, ONLY : fteik_solver_getSweepLimitsF
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc, sweep, nLevels
      LOGICAL(C_BOOL), VALUE, INTENT(IN) :: linitk
      INTEGER(C_INT), DIMENSION(:), INTENT(IN) :: ijkv, levelPtr
      LOGICAL(C_BOOL), DIMENSION(:), INTENT(INOUT) :: lupd
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) i, ix, iy, iz, node, maxx, maxy, maxz, minx, miny, minz, &
                     x1, x2, xsi, y1, y2, ysi, z1, z2, zsi
      IF (linitk) THEN
         CALL fteik_source_getSourceIndices32iF(isrc, zsi, xsi, ysi, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'fteik_solver_setUpdateNodesF: Error with source node init'
            ierr = 1
            RETURN
         ENDIF
      ENDIF
      CALL fteik_solver_getSweepLimitsF(sweep, linitk,          &
                                        nz, nx, ny,             &
                                        zsi, xsi, ysi,          &
                                        z1, z2, x1, x2, y1, y2, &
                                        ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver_setUpdateNodesF: Invalid sweep', sweep
         ierr = 1
         RETURN
      ENDIF
      minz = MIN(z1, z2)
      maxz = MAX(z1, z2)
      minx = MIN(x1, x2)
      maxx = MAX(x1, x2)
      miny = MIN(y1, y2)
      maxy = MAX(y1, y2)
      DO 1 i=1,nLevels
         DO 2 node=levelPtr(i),levelPtr(i+1)-1 
            iz   = ijkv(4*(node - 1) + 1)
            ix   = ijkv(4*(node - 1) + 2)
            iy   = ijkv(4*(node - 1) + 3)
            IF (iz >= minz .AND. iz <= maxz .AND. &
                ix >= minx .AND. ix <= maxx .AND. &
                iy >= miny .AND. iy <= maxy) THEN
               lupd(node) = .TRUE.
            ENDIF 
    2    CONTINUE
    1 CONTINUE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Returns the sweep limits the given sweep.
!>
!>    @param[in] sweep    Sweep number.  This must be in the range [1,8].
!>    @param[in] linitk   If true then the sweep ranges are limited by the source
!>                        position.  This is useful on start-up when causality prevents
!>                        a chunk of nodes from being able to be updated.
!>                        Otherwise, these are the start and stop indices for the
!>                        entire grid.
!>    @param[in] nz       Number of z grid points in mesh.
!>    @param[in] nx       Number of x grid points in mesh.
!>    @param[in] ny       Number of y grid points in mesh.
!>    @param[in] zsi      z source grid point (Fortran numbered).
!>    @param[in] xsi      x source grid point (Fortran numbered).
!>    @param[in] ysi      y source grid point (Fortran numbered).
!>
!>    @param[out] z1      Start index in grid for loop on z (Fortran numbered).
!>    @param[out] z2      End index in grid for loop on z (Fortran numbered).
!>    @param[out] x1      Start index in grid for loop on x (Fortran numbered).
!>    @param[out] x2      End index in grid for loop on x (Fortran numbered).
!>    @param[out] y1      Start index in grid for loop on y (Fortran numbered).
!>    @param[out] y2      End index in grid for loop on y (Fortran numbered).
!>    @param[out] ierr    0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solver_getSweepLimitsF(sweep, linitk,          &
                                              nz, nx, ny,             &
                                              zsi, xsi, ysi,          &
                                              z1, z2, x1, x2, y1, y2, &
                                              ierr)                   &
                 BIND(C, NAME='fteik_solver_getSweepLimitsF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: sweep, nx, ny, nz, xsi, ysi, zsi
      LOGICAL(C_BOOL), VALUE, INTENT(IN) :: linitk
      INTEGER(C_INT), INTENT(OUT) :: x1, x2, y1, y2, z1, z2, ierr
      ierr = 0
      IF (sweep == 1) THEN
         y1 = 2;  x1 = 2;  z1 = 2;
         y2 = ny; x2 = nx; z2 = nz;
         IF (linitk) THEN
            y1 = MAX(2, ysi); x1 = MAX(2, xsi); z1 = MAX(2, zsi);
            y2 = ny;          x2 = nx;          z2 = nz;
         ENDIF
      ELSE IF (sweep == 2) THEN
         y1 = 2;  x1 = nx - 1; z1 = 2; 
         y2 = ny; x2 = 1;      z2 = nz;
         IF (linitk) THEN
            y1 = MAX(2, ysi); x1 = xsi + 1;     z1 = MAX(2, zsi);
            y2 = ny;          x2 = 1;           z2 = nz;
         ENDIF
      ELSE IF (sweep == 3) THEN
         y1 = ny - 1; x1 = 2;  z1 = 2;
         y2 = 1;      x2 = nx; z2 = nz;
         IF (linitk) THEN
            y1 = ysi + 1;     x1 = MAX(2, xsi); z1 = MAX(2, zsi);
            y2 = 1;           x2 = nx;          z2 = nz;
         ENDIF
      ELSE IF (sweep == 4) THEN
         y1 = ny - 1; x1 = nx - 1; z1 = 2;
         y2 = 1;      x2 = 1;      z2 = nz;
         IF (linitk) THEN
            y1 = ysi + 1;     x1 = xsi + 1;     z1 = MAX(2, zsi);
            y2 = 1;           x2 = 1;           z2 = nz;
         ENDIF
      ELSE IF (sweep == 5) THEN
         y1 = 2;  x1 = 2;  z1 = nz - 1
         y2 = ny; x2 = nx; z2 = 1;
         IF (linitk) THEN 
            y1 = MAX(2, ysi); x1 = MAX(2, xsi); z1 = zsi + 1
            y2 = ny;          x2 = nx;          z2 = 1;
         ENDIF
      ELSE IF (sweep == 6) THEN
         y1 = 2;  x1 = nx - 1; z1 = nz - 1;
         y2 = ny; x2 = 1;      z2 = 1;
         IF (linitk) THEN 
            y1 = MAX(2, ysi); x1 = xsi + 1;     z1 = zsi + 1
            y2 = ny;          x2 = 1;           z2 = 1;
         ENDIF
      ELSE IF (sweep == 7) THEN
         y1 = ny - 1; x1 = 2;  z1 = nz - 1;
         y2 = 1;      x2 = nx; z2 = 1;
         IF (linitk) THEN 
            y1 = ysi + 1;     x1 = MAX(2, xsi); z1 = zsi + 1
            y2 = 1;           x2 = nx;          z2 = 1;
         ENDIF
      ELSE IF (sweep == 8) THEN
         y1 = ny - 1; x1 = nx - 1; z1 = nz - 1; 
         y2 = 1;      x2 = 1;      z2 = 1;
         IF (linitk) THEN 
            y1 = ysi + 1;     x1 = xsi + 1;     z1 = zsi + 1
            y2 = 1;           x2 = 1;           z2 = 1;
         ENDIF
      ELSE
         WRITE(*,*) 'fteik_solver_getSweepLimitsF: Invalid sweep', sweep
         ierr = 1
      ENDIF 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases all memory associated with the solver and resets all variables.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solver_finalizeF()             &
                 BIND(C, NAME='fteik_solver_finalizeF')
      USE FTEIK_SOLVER64F, ONLY : ttimes,  &
                                  ijkv1, ijkv2, ijkv3, ijkv4, &
                                  ijkv5, ijkv6, ijkv7, ijkv8, &
                                  lupd1, lupd2, lupd3, lupd4, &
                                  lupd5, lupd6, lupd7, lupd8, &
                                  epsS2C, nsweep, lhaveTimes
      USE FTEIK_RECEIVER64F, ONLY : fteik_receiver_finalizeF
      USE FTEIK_SOURCE64F, ONLY : fteik_source_finalizeF
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      lhaveTimes = .FALSE.
      epsS2C = zero
      nsweep = 0
      IF (ALLOCATED(ttimes)) DEALLOCATE(ttimes)
      IF (ALLOCATED(lupd1)) DEALLOCATE(lupd1)
      IF (ALLOCATED(lupd2)) DEALLOCATE(lupd2)
      IF (ALLOCATED(lupd3)) DEALLOCATE(lupd3)
      IF (ALLOCATED(lupd4)) DEALLOCATE(lupd4)
      IF (ALLOCATED(lupd5)) DEALLOCATE(lupd5)
      IF (ALLOCATED(lupd6)) DEALLOCATE(lupd6)
      IF (ALLOCATED(lupd7)) DEALLOCATE(lupd7)
      IF (ALLOCATED(lupd8)) DEALLOCATE(lupd8)
      IF (ALLOCATED(ijkv1)) DEALLOCATE(ijkv1)
      IF (ALLOCATED(ijkv2)) DEALLOCATE(ijkv2)
      IF (ALLOCATED(ijkv3)) DEALLOCATE(ijkv3)
      IF (ALLOCATED(ijkv4)) DEALLOCATE(ijkv4)
      IF (ALLOCATED(ijkv5)) DEALLOCATE(ijkv5)
      IF (ALLOCATED(ijkv6)) DEALLOCATE(ijkv6)
      IF (ALLOCATED(ijkv7)) DEALLOCATE(ijkv7)
      IF (ALLOCATED(ijkv8)) DEALLOCATE(ijkv8)
      CALL fteik_receiver_finalizeF()
      CALL fteiK_source_finalizeF()
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the spherical approximation transition tolerance around the source.
!>
!>    @param[in] epsIn   Radius in number of grid points around source where the
!>                       spherical approximation finite difference stencils will 
!>                       be used.  This is non-negative.
!>    @param[out] ierr   0 indicates success.
!>
!>    @author Mark Noble, Alexandrine Gesret, and Ben Baker
!>
!>    @version 2
!>
!>    @date July 2017
!>
!>    @copyright CeCILL-3
!>
!>    This is now a subroutine and incorporated with the fteik Fortran modules.
!>
      SUBROUTINE fteik_solver_setSphereToCartEpsilonF(epsIn, ierr) &
                 BIND(C, NAME='fteik_solver_setSphereToCartEpsilonF')
      USE ISO_C_BINDING
      USE FTEIK_MODEL64F, ONLY : nz, nx, ny
      USE FTEIK_SOLVER64F, ONLY : epsS2C
      USE FTEIK_CONSTANTS64F, ONLY : zero
      IMPLICIT NONE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: epsIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      epsS2C = zero
      IF (nx < 1 .OR. ny < 1 .OR. nz < 1) THEN
         WRITE(*,*) 'fteik_solver_setSphereToCartEpsilonF: Grid not yet set'
         ierr = 1
         RETURN
      ENDIF
      IF (INT(epsIn) > nz .OR. INT(epsIn) > nx .OR. INT(epsIn) > ny) THEN
         IF (INT(epsIn) > nz) THEN
            WRITE(*,*) 'fteik_solver_setSphereToCartEpsilonF: eps bigger than nz', &
                       INT(epsIn), nz
         ENDIF
         IF (INT(epsIn) > nx) THEN
            WRITE(*,*) 'fteik_solver_setSphereToCartEpsilonF: eps bigger than nx', &
                       INT(epsIn), nx
         ENDIF
         IF (INT(epsIn) > ny) THEN
            WRITE(*,*) 'fteik_solver_setSphereToCartEpsilonF: eps bigger than ny', &
                       INT(epsIn), ny
         ENDIF
         ierr = 1
         RETURN
      ENDIF
      epsS2C = epsIn
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the number of Gauss-Seidel iterations  in the fast sweeping method.
!>
!>    @param[in] nsweepIn    Number of sweeps.  This cannot be negative and 1 is
!>                           usually sufficient.
!>
!>    @param[out] ierr       0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT 
!>
      SUBROUTINE fteik_solver_setNumberOfSweepsF(nsweepIn, ierr) &
                 BIND(C, NAME='fteik_solver_setNumberOfSweepsF')
      USE FTEIK_SOLVER64F, ONLY : nsweep
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsweepIn 
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      nsweep = 0
      IF (nsweepIn < 0) THEN
         WRITE(*,*) 'fteik_solver_setNumberOfSweepsF: nsweep must be positive', nsweep
         ierr = 1
         RETURN
      ENDIF
      nsweep = nsweepIn
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !

