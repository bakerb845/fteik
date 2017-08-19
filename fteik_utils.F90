!                                                                                        !
!========================================================================================!
!                                                                                        !
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
      SUBROUTINE fteik_initializeF(nzIn, nxIn, nyIn,      &
                                   z0In, x0In, y0In,      &
                                   dzIn, dxIn, dyIn,      &
                                   nsweepIn, epsIn, ierr) &
                 BIND(C, NAME='fteik_initializeF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_setGridSizeF, fteik_setGridSpacingF, &
                                 fteik_setModelOriginF, fteik_setNumberOfSweepsF, &
                                 fteik_setSphericalToCartesianEpsilonF, &
                                 fteik_computeGraphF, fteik_setReceivers64fF
      USE FTEIK_UTILS64F, ONLY : lhaveGrid, lhaveGridSpacing, lhaveSource, &
                                 lhaveSlownessModel, lhaveTravelTimes
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN) :: x0In, y0In, z0In
      REAL(C_DOUBLE), INTENT(IN) :: dxIn, dyIn, dzIn
      REAL(C_DOUBLE), INTENT(IN) :: epsIn
      INTEGER(C_INT), INTENT(IN) :: nxIn, nyIn, nzIn
      INTEGER(C_INT), INTENT(IN) :: nsweepIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) xt(1), yt(1), zt(1)
      lhaveGrid = .FALSE.
      lhaveGridSpacing = .FALSE.
      lhaveSource = .FALSE.
      lhaveSlownessModel = .FALSE.
      lhaveTravelTimes = .FALSE.
      CALL fteik_setGridSizeF(nzIn, nxIn, nyIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_initializeF: Error setting the grid size'
         RETURN
      ENDIF
      CALL fteik_setGridSpacingF(dzIn, dxIn, dyIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_initializeF: Error setting grid spacing'
         RETURN
      ENDIF
      CALL fteik_setNumberOfSweepsF(nsweepIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_initializeF: Failed to set max number of sweeps'
         RETURN
      ENDIF
      CALL fteik_setSphericalToCartesianEpsilonF(epsIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_initializeF: Failed to set epsilon'
         RETURN
      ENDIF
      CALL fteik_setModelOriginF(z0In, x0In, y0In)
      CALL fteik_computeGraphF(ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_initializeF: Failed to compute graph'
         RETURN
      ENDIF
      CALL fteik_setReceivers64fF(-1, xt, yt, zt, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_initializeF: Strange error in zero-ing out receivers'
         ierr = 0
      ENDIF
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE fteik_computeGraphF(ierr)          &
                 BIND(C, NAME='fteik_computeGraphF') 
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : lhaveGrid, nx, ny, nz
      USE FTEIK_UTILS64F, ONLY : nLevels, maxLevelSize
      USE FTEIK_UTILS64F, ONLY : ijkv1, ijkv2, ijkv3, ijkv4, ijkv5, ijkv6, ijkv7, ijkv8
      USE FTEIK_UTILS64F, ONLY : lupd1, lupd2, lupd3, lupd4, lupd5, lupd6, lupd7, lupd8
      USE FTEIK_UTILS64F, ONLY : levelPtr
      USE FTEIK_UTILS64F, ONLY : tt1 
      USE FTEIK_UTILS64F, ONLY : zero
      USE FTEIK_UTILS64F, ONLY : fteik_setUpdateNodesF
      USE FTEIK_GRAPH
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: ierr
      TYPE(GRAPH_TYPE) graph
      INTEGER(C_INT) ierrs(8), i
      ierr = 0
      IF (.NOT.lhaveGrid) THEN
         WRITE(*,*) 'fteik_computeGraph: Error grid not yet set'
         ierr = 1
         RETURN
      ENDIF
      CALL INIT(graph, nz, nx, ny, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_computeGraphF: Error initializing graph'
         ierr = 1
         RETURN
      ENDIF
      nLevels = NUMBER_OF_LEVELS(graph)
      ierr = IJKV(graph, 1, ijkv1)
      ierr = IJKV(graph, 2, ijkv2) + ierr
      ierr = IJKV(graph, 3, ijkv3) + ierr
      ierr = IJKV(graph, 4, ijkv4) + ierr
      ierr = IJKV(graph, 5, ijkv5) + ierr
      ierr = IJKV(graph, 6, ijkv6) + ierr
      ierr = IJKV(graph, 7, ijkv7) + ierr
      ierr = IJKV(graph, 8, ijkv8) + ierr
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_computeGraphF: Errors encountered getting ijkv'
         ierr = 1
         RETURN
      ENDIF
      ierr = LEVEL_PTR(graph, 1, levelPtr)
      ! TODO all levelPtrs better be the same; should only copy 1
      !ierr = LEVEL_PTR(graph, 1, levelPtr1)
      !ierr = LEVEL_PTR(graph, 2, levelPtr2) + ierr
      !ierr = LEVEL_PTR(graph, 3, levelPtr3) + ierr
      !ierr = LEVEL_PTR(graph, 4, levelPtr4) + ierr
      !ierr = LEVEL_PTR(graph, 5, levelPtr5) + ierr
      !ierr = LEVEL_PTR(graph, 6, levelPtr6) + ierr
      !ierr = LEVEL_PTR(graph, 7, levelPtr7) + ierr
      !ierr = LEVEL_PTR(graph, 8, levelPtr8) + ierr
      !print *, 'bigloc:', maxval(levelPtr4 - levelPtr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_computeGraphF: Error encountered getting levelPtr'
         ierr = 1
         RETURN
      ENDIF
      ! finished with the graph
      CALL DELETE(graph)
      maxLevelSize = 0
      DO 1 i=1,nLevels
         maxLevelSize = MAX(maxLevelSize, levelPtr(i+1) - levelPtr(i))
    1 CONTINUE
print *, maxLevelSize
      ALLOCATE(tt1(maxLevelSize))
      tt1(:) = zero
      ! Set the update grid.  This will mask nodes on the boundary.
      CALL fteik_setUpdateNodesF(1, nlevels, .FALSE., levelPtr, ijkv1, lupd1, ierrs(1))
      CALL fteik_setUpdateNodesF(2, nlevels, .FALSE., levelPtr, ijkv2, lupd2, ierrs(2))
      CALL fteik_setUpdateNodesF(3, nlevels, .FALSE., levelPtr, ijkv3, lupd3, ierrs(3))
      CALL fteik_setUpdateNodesF(4, nlevels, .FALSE., levelPtr, ijkv4, lupd4, ierrs(4))
      CALL fteik_setUpdateNodesF(5, nlevels, .FALSE., levelPtr, ijkv5, lupd5, ierrs(5))
      CALL fteik_setUpdateNodesF(6, nlevels, .FALSE., levelPtr, ijkv6, lupd6, ierrs(6))
      CALL fteik_setUpdateNodesF(7, nlevels, .FALSE., levelPtr, ijkv7, lupd7, ierrs(7))
      CALL fteik_setUpdateNodesF(8, nlevels, .FALSE., levelPtr, ijkv8, lupd8, ierrs(8))
      IF (MAXVAL(ABS(ierrs)) /= 0) THEN
         WRITE(*,*) 'fteik_computeGraphF: Error setting update nodes'
         ierr = 1
      ENDIF
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>
!>    @brief Gets the travel-times at the receiver locations.
!>
!>    @param[int] nrecIn    Length of array trec.  This must be at least nrec. 
!>
!>    @param[out] trec      Travel-time from source to receiver (seconds).  This
!>                          has dimension [nrecIn].  If nrecIn > nrec then the
!>                          unused elements of trec will be set to FTEIK_HUGE.
!>
!>    @param[out] ierr      0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_getTravelTimesAtReceivers64fF(nrecIn, trec, ierr) &
                 BIND(C, NAME='fteik_getTravelTimesAtReceivers64fF')
      USE FTEIK_UTILS64F, ONLY : ttimes, xri, yri, zri, xdr, ydr, zdr, &
                                 nz, nzx, nrec, lhaveTravelTimes
      USE FTEIK_UTILS64F, ONLY : one, FTEIK_HUGE
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: nrecIn
      REAL(C_DOUBLE), INTENT(OUT) :: trec(nrecIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) v000, v001, v010, v011, v100, v101, v110, v111
      REAL(C_DOUBLE) c0, c1, c00, c10, c01, c11
      INTEGER(C_INT) i
      ! Some error checks
      ierr = 0
      IF (nrec < 1) THEN
         WRITE(*,*) 'fteik_getTravelTimesAtReceivers64fF: recv indices not yet set'
         ierr = 1
      ENDIF
      IF (nrecIn < nrec) THEN
         WRITE(*,*) 'fteik_getTravelTimesAtRecievers64fF: nrecIn < nrec!'
         ierr = 1
         RETURN
      ENDIF
      IF (nrecIn > nrec) trec(nrec+1:nrecIn) = FTEIK_HUGE
      IF (.NOT. lhaveTravelTimes) THEN
         WRITE(*,*) 'fteik_getTravelTimesAtReceivers64fF: ttimes not yet computed'
         ierr = 1
         RETURN
      ENDIF
      ! Compute travel-times
      DO 1 i=1,nrec
         ! Try to make the memory accesses as sequential as possible
         v000 = ttimes(grid2indexF(zri(i),   xri(i),   yri(i),   nz, nzx))
         v100 = ttimes(grid2indexF(zri(i)+1, xri(i),   yri(i),   nz, nzx))
         v010 = ttimes(grid2indexF(zri(i),   xri(i)+1, yri(i),   nz, nzx))
         v110 = ttimes(grid2indexF(zri(i)+1, xri(i)+1, yri(i),   nz, nzx))
         v001 = ttimes(grid2indexF(zri(i),   xri(i),   yri(i)+1, nz, nzx))
         v101 = ttimes(grid2indexF(zri(i)+1, xri(i),   yri(i)+1, nz, nzx))
         v011 = ttimes(grid2indexF(zri(i),   xri(i)+1, yri(i)+1, nz, nzx))
         v111 = ttimes(grid2indexF(zri(i)+1, xri(i)+1, yri(i)+1, nz, nzx)) 
         ! Interpolate in z
         c00 = V000*(one - zdr(i)) + V100*zdr(i);
         c10 = V010*(one - zdr(i)) + V110*zdr(i);
         c01 = V001*(one - zdr(i)) + V101*zdr(i);
         c11 = V011*(one - zdr(i)) + V111*zdr(i);
         ! Inteprolate in x
         c0 = c00*(one - xdr(i)) + c10*xdr(i);
         c1 = c01*(one - xdr(i)) + c11*xdr(i);
         ! Interpolate in y
         trec(i) = c0*(one - ydr(i)) + c1*ydr(i);
    1 CONTINUE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the receiver indices in the grid.  This will make for a quick 
!>           linear interpolation.
!>           https://github.com/bottero/IMCMCrun/blob/master/src/functions.cpp
!>
!>    @param[in] nrecIn   Number of receivers.  If this is -1 then this function will
!>                        simply deallocate/reset the receiver information.
!>    @param[in] z        Receiver positions in z (meters).  This has dimension [nrecIn].
!>    @param[in] x        Receiver positions in x (meters).  This has dimension [nrecIn].
!>    @param[in] y        Receiver positions in y (meters).  This has dimension [nrecIn].
!>
!>    @param[out] ierr    0 indicate success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_setReceivers64fF(nrecIn, z, x, y, ierr) &
                 BIND(C, NAME='fteik_setReceivers64fF')
      USE FTEIK_UTILS64F, ONLY : dx, dy, dz, nx, ny, nz, x0, y0, z0, lhaveGrid
      USE FTEIK_UTILS64F, ONLY : nrec, xri, yri, zri, xdr, ydr, zdr
      USE FTEIK_UTILS64F, ONLY : one
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: nrecIn
      REAL(C_DOUBLE), INTENT(IN) :: z(nrecIn), x(nrecIn), y(nrecIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) xr, yr, zr
      INTEGER(C_INT) irec, ix, iy, iz
      ! Initialize 
      ierr = 0 
      nrec = 0
      IF (ALLOCATED(zri)) DEALLOCATE(zri)
      IF (ALLOCATED(xri)) DEALLOCATE(xri)
      IF (ALLOCATED(yri)) DEALLOCATE(yri)
      IF (nrecIn ==-1) RETURN
      ! Check that the model has been initialized
      IF (.NOT.lhaveGrid) THEN 
         WRITE(*,*) 'fteik_setReceivers64fF: Grid not yet set'
         ierr = 1
         RETURN
      ENDIF
      ! Verify the inputs
      IF (nrecIn < 1) THEN
         WRITE(*,*) 'fteik_setReceivers64fF: No receivers'
         ierr = 1
         RETURN
      ENDIF
      ! Initialize
      nrec = nrecIn
      ALLOCATE(zri(nrec))
      ALLOCATE(xri(nrec))
      ALLOCATE(yri(nrec))
      ALLOCATE(zdr(nrec))
      ALLOCATE(xdr(nrec))
      ALLOCATE(ydr(nrec))
      ! Collocate the receivers
      DO 1 irec=1,nrecIn
         IF (z(irec) < z0 .OR. z(irec) > z0 + dz*DBLE(nz - 1)) THEN
            WRITE(*,*) 'fteik_setReceivers: Warning z position is out of model', z(irec)
            ierr = ierr + 1
         ENDIF
         IF (x(irec) < x0 .OR. x(irec) > x0 + dx*DBLE(nx - 1)) THEN
            WRITE(*,*) 'fteik_setReceivers: Warning x position is out of bounds', x(irec)
            ierr = ierr + 1
         ENDIF
         IF (y(irec) < y0 .OR. y(irec) > y0 + dy*DBLE(ny - 1)) THEN
            WRITE(*,*) 'fteik_setReceivers: Warning y position is out of bounds', y(irec)
            ierr = ierr + 1
         ENDIF
         ! Get the grid point position
         zr = (z(irec) - z0)/dz
         xr = (x(irec) - x0)/dx
         yr = (y(irec) - y0)/dy
         iz = INT(zr) + 1
         ix = INT(xr) + 1
         iy = INT(yr) + 1
         ! Interpolation goes from (iz,ix,iy) to (iz+1,ix+1,iy+1).  This is an issue when
         ! the receiver is at the edge.  In this case; shift it to model bounds - 1. 
         ! Additionally, being extra generous, I'll force the receiver into the model
         ! if the receiver is out of bounds.
         zri(irec) = MAX(1, MIN(nz-1, iz))
         xri(irec) = MAX(1, MIN(nx-1, ix))
         yri(irec) = MAX(1, MIN(ny-1, iy))
         zdr(irec) = (z(irec) - (z0 + dz*DBLE(zri(irec)-1)))/dz
         xdr(irec) = (x(irec) - (x0 + dx*DBLE(xri(irec)-1)))/dx
         ydr(irec) = (y(irec) - (y0 + dy*DBLE(yri(irec)-1)))/dy
    1 CONTINUE
      RETURN
      END
! int kz=iz+1, kx=ix+1, ky=iy+1;
! double dz1=Z-(double)(iz+1),dx1=X-(double)(ix+1),dy1=Y-(double)(iy+1);
! double dz2=1.0-dz1,dx2=1.0-dx1,dy2=1.0-dy1;
! double t = dz2 * dx2 * dy2 * tt3d->get(iz,ix,iy)
!       +  dz1 * dx2 * dy2 * tt3d->get(kz,ix,iy)
!       +  dz2 * dx1 * dy2 * tt3d->get(iz,kx,iy)
!       +  dz2 * dx2 * dy1 * tt3d->get(iz,ix,ky)
!       +  dz2 * dx1 * dy1 * tt3d->get(iz,kx,ky)
!       +  dz1 * dx2 * dy1 * tt3d->get(kz,ix,ky)
!       +  dz1 * dx1 * dy2 * tt3d->get(kz,kx,iy)
!       +  dz1 * dx1 * dy1 * tt3d->get(kz,kx,ky);
!     RETURN
!     END
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
!>    @param[in] levelPtr   Maps from level'th level to start node in ijkv.
!>                          This is a vector of dimension [nLevels+1].
!>    @param[in] ijkv       Maps from in'th node in to the (iz,ix,iy,igrd) indices.
!>                          This is a vector of length [4*ngrd] where ngrd is the
!>                          number of grid points.
!>
!>    @param[out] lupd      Determines if the in'th node is to be updated for the 
!>                          given sweep.  This is a vector of length [ngrd].
!>    @param[out] ierr      0 indicates success. 
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>    
      SUBROUTINE fteik_setUpdateNodesF(sweep, nLevels, linitk,     &
                                       levelPtr, ijkv, lupd, ierr) &
                 BIND(C, NAME='fteik_setUpdateNodesF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : nx, ny, nz, xsi, ysi, zsi, lhaveSource
      USE FTEIK_UTILS64F, ONLY : fteik_getSweepLimitsF, index2gridF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: sweep, nLevels
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      INTEGER(C_INT), DIMENSION(:), INTENT(IN) :: ijkv, levelPtr
      LOGICAL(C_BOOL), DIMENSION(:), INTENT(OUT) :: lupd
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) i, igrd, ix, iy, iz, node, maxx, maxy, maxz, minx, miny, minz, &
                     x1, x2, y1, y2, z1, z2
      lupd(:) = .FALSE.
      IF (linitk) THEN
         IF (.NOT.lhaveSource) THEN
            WRITE(*,*) 'fteik_setUpdateNodesF: Error source not initialized'
            ierr = 1
            RETURN
         ENDIF
      ENDIF
      CALL fteik_getSweepLimitsF(sweep, linitk,          &
                                 nz, nx, ny,             &
                                 zsi, xsi, ysi,          &
                                 z1, z2, x1, x2, y1, y2, &
                                 ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_setUpdateNodesF: Invalid sweep', sweep
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
            igrd = ijkv(4*(node - 1) + 4)
            !CALL index2gridF(igrd, iz, ix, iy, ierr)
            IF (ierr == 0) THEN
               IF (iz >= minz .AND. iz <= maxz .AND. &
                   ix >= minx .AND. ix <= maxx .AND. &
                   iy >= miny .AND. iy <= maxy) THEN
                  lupd(node) = .TRUE.
               ENDIF 
            ENDIF
    2    CONTINUE
    1 CONTINUE
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
      SUBROUTINE fteik_setSphericalToCartesianEpsilonF(epsIn, ierr) &
                 BIND(C, NAME='fteik_setSphericalToCartesianEpsilonF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : epsS2C, nx, ny, nz, lhaveGrid
      USE FTEIK_UTILS64F, ONLY : zero
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN) :: epsIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      epsS2C = zero
      IF (.NOT.lhaveGrid) THEN
         WRITE(*,*) 'fteik_setSphericalToCartesianEpsilonF: grid not yet set'
         ierr = 1
         RETURN
      ENDIF
      IF (INT(epsIn) > nz .OR. INT(epsIn) > nx .OR. INT(epsIn) > ny) THEN
         IF (INT(epsIn) > nz) THEN
            WRITE(*,*) 'fteik_setSphericalToCartesianEpsilonF: eps bigger than nz', &
                       INT(epsIn), nz
         ENDIF
         IF (INT(epsIn) > nx) THEN
            WRITE(*,*) 'fteik_setSphericalToCartesianEpsilonF: eps bigger than nx', &
                       INT(epsIn), nx
         ENDIF
         IF (INT(epsIn) > ny) THEN
            WRITE(*,*) 'fteik_setSphericalToCartesianEpsilonF: eps bigger than ny', &
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
!>    @brief Frees the memory on the solver and resets variables to 0.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_finalizeF() BIND(C, NAME='fteik_finalizeF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : slow, ttimes
      USE FTEIK_UTILS64F, ONLY : lupd1, lupd2, lupd3, lupd4, lupd5, lupd6, lupd7, lupd8
      USE FTEIK_UTILS64F, ONLY : lupdInit1, lupdInit2, lupdInit3, lupdInit4, &
                                 lupdInit5, lupdInit6, lupdInit7, lupdInit8
      USE FTEIK_UTILS64F, ONLY : xsi, ysi, zsi, xsa, ysa, zsa
      USE FTEIK_UTILS64F, ONLY : epsS2C
      USE FTEIK_UTILS64F, ONLY : dx, dy, dz, x0, y0, z0
      USE FTEIK_UTILS64F, ONLY : ngrd, ncell, nx, ny, nz, nzm1, nzm1_nxm1, nzx
      USE FTEIK_UTILS64F, ONLY : lhaveGrid, lhaveGridSpacing, nx, ny, nz
      USE FTEIK_UTILS64F, ONLY : nLevels, maxLevelSize, nsweep
      USE FTEIK_UTILS64F, ONLY : ijkv1, ijkv2, ijkv3, ijkv4, ijkv5, ijkv6, ijkv7, ijkv8
      USE FTEIK_UTILS64F, ONLY : levelPtr
      USE FTEIK_UTILS64F, ONLY : tt1
      USE FTEIK_UTILS64F, ONLY : lhaveGrid, lhaveSource, lhaveSlownessModel
      USE FTEIK_UTILS64F, ONLY : zero
      USE FTEIK_UTILS64F, ONLY : fteik_setReceivers64fF
      IMPLICIT NONE
      REAL(C_DOUBLE) xt(1), yt(1), zt(1)
      INTEGER(C_INT) ierr
      IF (ALLOCATED(slow))      DEALLOCATE(slow)
      IF (ALLOCATED(ttimes))    DEALLOCATE(ttimes)
      IF (ALLOCATED(lupd1))     DEALLOCATE(lupd1)
      IF (ALLOCATED(lupd2))     DEALLOCATE(lupd2)
      IF (ALLOCATED(lupd3))     DEALLOCATE(lupd3)
      IF (ALLOCATED(lupd4))     DEALLOCATE(lupd4)
      IF (ALLOCATED(lupd5))     DEALLOCATE(lupd5)
      IF (ALLOCATED(lupd6))     DEALLOCATE(lupd6)
      IF (ALLOCATED(lupd7))     DEALLOCATE(lupd7)
      IF (ALLOCATED(lupd8))     DEALLOCATE(lupd8)
      IF (ALLOCATED(lupdInit1)) DEALLOCATE(lupdInit1)
      IF (ALLOCATED(lupdInit2)) DEALLOCATE(lupdInit2)
      IF (ALLOCATED(lupdInit3)) DEALLOCATE(lupdInit3)
      IF (ALLOCATED(lupdInit4)) DEALLOCATE(lupdInit4)
      IF (ALLOCATED(lupdInit5)) DEALLOCATE(lupdInit5)
      IF (ALLOCATED(lupdInit6)) DEALLOCATE(lupdInit6)
      IF (ALLOCATED(lupdInit7)) DEALLOCATE(lupdInit7)
      IF (ALLOCATED(lupdInit8)) DEALLOCATE(lupdInit8)
      IF (ALLOCATED(levelPtr))  DEALLOCATE(levelPtr)
      IF (ALLOCATED(ijkv1))     DEALLOCATE(ijkv1)
      IF (ALLOCATED(ijkv2))     DEALLOCATE(ijkv2)
      IF (ALLOCATED(ijkv3))     DEALLOCATE(ijkv3)
      IF (ALLOCATED(ijkv4))     DEALLOCATE(ijkv4)
      IF (ALLOCATED(ijkv5))     DEALLOCATE(ijkv5)
      IF (ALLOCATED(ijkv6))     DEALLOCATE(ijkv6)  
      IF (ALLOCATED(ijkv7))     DEALLOCATE(ijkv7)
      IF (ALLOCATED(ijkv8))     DEALLOCATE(ijkv8)
      IF (ALLOCATED(tt1))       DEALLOCATE(tt1)
      CALL fteik_setReceivers64fF(-1, xt, yt, zt, ierr)
      dz = zero; dx = zero; dy = zero
      zsa = zero; xsa = zero; ysa = zero
      z0 = zero; x0 = zero; y0 = zero
      zsi = 0; xsi = 0; ysi = 0
      epsS2C = zero
      nz = 0; nx = 0; ny = 0
      ngrd = 0; ncell = 0; nzx = 0; nzm1 = 0; nzm1_nxm1 = 0
      nLevels = 0; maxLevelSize = 0; nsweep = 0
      lhaveGrid = .FALSE.
      lhaveGridSpacing = .FALSE.
      lhaveSource = .FALSE.
      lhaveSlownessModel = .FALSE.
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the number of refinement sweeps (iterations) in the fast
!>           sweeping method.
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
      SUBROUTINE fteik_setNumberOfSweepsF(nsweepIn, ierr) &
                 BIND(C, NAME='fteik_setNumberOfSweepsF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : nsweep
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: nsweepIn 
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      nsweep = 0
      IF (nsweepIn < 0) THEN
         WRITE(*,*) 'fteik_setNumberOfSweepsF: nsweep must be positive', nsweep
         ierr = 1
         RETURN
      ENDIF
      nsweep = nsweepIn
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Convenience utility to get the internal grid sizes.
!>
!>    @param[out] nzOut     Number of z grid points in travel time field
!>    @param[out] nxOut     Number of x grid points in travel time field.
!>    @param[out] nyOUt     Number of y grid points in travel time field.
!>    @param[out] ngrdOut   Number of grid points in travel time field (=nz*nx*ny).
!>    @param[out] ncellOut  Number of cells in velocity field (=(nz-1)*(nx-1)*(ny-1)).
!>    @param[out] ierr      0 indicates success.
!> 
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_getGridSizeF(nzOut, nxOut, nyOut, ngrdOut, ncellOut, ierr) &
                 BIND(C, NAME='fteik_getGridSizeF')
      USE FTEIK_UTILS64F, ONLY : nz, nx, ny, ngrd, ncell, lhaveGrid
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: nzOut, nxOut, nyOut, ngrdOut, ncellOut, ierr
      ierr = 0
      nzOut = 0
      nxOut = 0
      nyOut = 0
      ngrdOut = 0
      ncellOut = 0
      IF (.NOT.lhaveGrid) THEN
         WRITE(*,*) 'fteik_getGridSize: Grid not yet set'
         ierr = 1
         RETURN
      ENDIF
      nzOut = nz
      nxOut = nx
      nyOut = ny
      ngrdOut = ngrd
      ncellOut = ncell
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Returns the number of levels in the level-scheduling method.
!>
!>    @param[out] nLevelsOut   On successful exit this is the number of levels. \n
!>                             Otherwise, it is 0.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_getNumberOfLevelsF(nLevelsOut) &
                 BIND(C, NAME='fteik_getNumberOfLevelsF')
      USE FTEIK_UTILS64F, ONLY : nLevels, lhaveGrid
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(OUT) :: nLevelsOut
      nLevelsOut = 0
      IF (.NOT.lhaveGrid) THEN
         WRITE(*,*) 'fteik_getNumberOfLevelsF: Graph not initialized'
         RETURN
      ENDIF
      nLevelsOut = nLevels
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Returns the graph pointers fo rthe level-scheduling method.
!>
!>    @param[in] sweep         Sweep number.  This must be in the range of [1,8].
!>    @param[in] nLevelsIn     Workspace-1 of levelPtrOut.  This must be at least nLevels.
!>    @param[in] ngrd4         Workspace of ijkvOut.  This must be at least 4*ngrd.
!>    @param[out] levelPtrOut  Maps from the level'th level to the start of ijkvOut.
!>                             This is an array of dimension [nLevelsIn+1].
!>    @param[out] ijkvOut      Returns the (iz, ix, iy, global node number) of the 
!>                             node'th point in the level scheduling method.  This is
!>                             an array of dimension [ngrd4].
!>    @param[out] ierr         0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_getGraphPointersF(sweep, nLevelsIn, ngrd4,     &
                                         levelPtrOut, ijkvOut, ierr)  &
                 BIND(C, NAME='fteik_getGraphPointersF')
      USE FTEIK_UTILS64F, ONLY : ijkv1, ijkv2, ijkv3, ijkv4, &
                                 ijkv5, ijkv6, ijkv7, ijkv8, &
                                 levelPtr, ngrd, nLevels, lhaveGrid
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(IN), VALUE :: sweep, nLevelsIn
      INTEGER(C_INT), INTENT(OUT) :: levelPtrOut(nLevelsIn+1), ijkvOut(ngrd4), ierr
      ierr = 0
      IF (.NOT.lhaveGrid) THEN
         WRITE(*,*) 'fteik_getGraphPointersF: Graph not initialized'
         ierr = 1
         RETURN
      ENDIF
      IF (ngrd4 < 4*ngrd) THEN
         WRITE(*,*) 'fteik_getGraphPointersF: Insufficient space for ijkv' 
         ierr = 1
         RETURN
      ENDIF
      IF (nLevelsIn < nLevels) THEN
         WRITE(*,*) 'fteik_getGraphPointersF: Insufficient space for levelPtrOut'
         ierr = 1
         RETURN
      ENDIF
      levelPtrOut(1:nLevels+1) = levelPtr(1:nLevels+1)
      IF (sweep == 1) THEN
         ijkvOut(1:4*ngrd) = ijkv1(1:4*ngrd) 
      ELSEIF (sweep == 2) THEN
         ijkvOut(1:4*ngrd) = ijkv2(1:4*ngrd)
      ELSEIF (sweep == 3) THEN
         ijkvOut(1:4*ngrd) = ijkv3(1:4*ngrd)
      ELSEIF (sweep == 4) THEN
         ijkvOut(1:4*ngrd) = ijkv4(1:4*ngrd)
      ELSEIF (sweep == 5) THEN
         ijkvOut(1:4*ngrd) = ijkv5(1:4*ngrd)
      ELSEIF (sweep == 6) THEN
         ijkvOut(1:4*ngrd) = ijkv6(1:4*ngrd)
      ELSEIF (sweep == 7) THEN
         ijkvOut(1:4*ngrd) = ijkv7(1:4*ngrd)
      ELSEIF (sweep == 8) THEN
         ijkvOut(1:4*ngrd) = ijkv8(1:4*ngrd)
      ELSE
         WRITE(*,*) 'fteik_getGraphPointersF: Invalid sweep', sweep
         ierr = 1
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the model origin.
!>
!>    @param[in] z0In    z origin (meters).
!>    @param[in] x0In    x origin (meters).
!>    @param[in] y0In    y origin (meters).
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_setModelOriginF(z0In, x0In, y0In) &
                 BIND(C, NAME='fteik_setModelOriginF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : x0, y0, z0
      REAL(C_DOUBLE), INTENT(IN) :: x0In, y0In, z0In
      z0 = z0In
      x0 = x0In
      y0 = y0In
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief This sets the number of grid points in each dimension of the travel-time
!>           field.  This function must be called to properly start the solver.
!>
!>    @param[in] nzIn     Number of z grid points in the travel-time field.
!>                        This must be at least 3.
!>    @param[in] nxIn     Number of x grid points in the travel-time field.
!>                        This must be at least 3.
!>    @param[in] nyIn     Number of y grid points in the travel-time field.
!>                        This must be at least 3.
!>
!>    @param[out] ierr    0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_setGridSizeF(nzIn, nxIn, nyIn, ierr) &
                 BIND(C, NAME='fteik_setGridSizeF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : lupd1, lupd2, lupd3, lupd4, &
                                 lupd5, lupd6, lupd7, lupd8, &
                                 slow, ttimes
      USE FTEIK_UTILS64F, ONLY : lupdInit1, lupdInit2, lupdInit3, lupdInit4, &
                                 lupdInit5, lupdInit6, lupdInit7, lupdInit8
      USE FTEIK_UTILS64F, ONLY : nz, nx, ny, ncell, ngrd, nzx, nzm1, nzm1_nxm1
      USE FTEIK_UTILS64F, ONLY : lhaveGrid
      USE FTEIK_UTILS64F, ONLY : FTEIK_HUGE, zero
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: nzIn, nxIn, nyIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      nz = 0
      nx = 0
      ny = 0
      ngrd = 0
      ncell = 0
      lhaveGrid = .FALSE.
      IF (nzIn < 3 .OR. nxIn < 3 .OR. nyIn < 3) THEN 
         IF (nzIn < 3) WRITE(*,*) 'fteik_setGridSize: ERROR nz is too small', nzIn
         IF (nxIn < 3) WRITE(*,*) 'fteik_setGridSize: ERROR nx is too small', nxIn
         IF (nyIn < 3) WRITE(*,*) 'fteik_setGridSize: ERROR ny is too small', nyIn
         ierr = 1
         RETURN
      ENDIF
      nz = nzIn
      nx = nxIn
      ny = nyIn
      ncell = (nz - 1)*(nx - 1)*(ny - 1)
      ngrd = nz*nx*ny 
      nzx = nz*nx
      nzm1 = nz - 1
      nzm1_nxm1 = (nz - 1)*(nx - 1)
      ALLOCATE(ttimes(ngrd))
      ALLOCATE(slow(ncell))
      ALLOCATE(lupd1(ngrd))
      ALLOCATE(lupd2(ngrd))
      ALLOCATE(lupd3(ngrd))
      ALLOCATE(lupd4(ngrd))
      ALLOCATE(lupd5(ngrd))
      ALLOCATE(lupd6(ngrd))
      ALLOCATE(lupd7(ngrd))
      ALLOCATE(lupd8(ngrd))
      ALLOCATE(lupdInit1(ngrd))
      ALLOCATE(lupdInit2(ngrd))
      ALLOCATE(lupdInit3(ngrd))
      ALLOCATE(lupdInit4(ngrd))
      ALLOCATE(lupdInit5(ngrd))
      ALLOCATE(lupdInit6(ngrd))
      ALLOCATE(lupdInit7(ngrd))
      ALLOCATE(lupdInit8(ngrd))
      ttimes(:) = FTEIK_HUGE
      slow(:) = zero 
      lupd1(:) = .FALSE.
      lupd2(:) = .FALSE.
      lupd3(:) = .FALSE.
      lupd4(:) = .FALSE.
      lupd5(:) = .FALSE.
      lupd6(:) = .FALSE.
      lupd7(:) = .FALSE.
      lupd8(:) = .FALSE.
      lupdInit1(:) = .FALSE.
      lupdInit2(:) = .FALSE.
      lupdInit3(:) = .FALSE.
      lupdInit4(:) = .FALSE.
      lupdInit5(:) = .FALSE.
      lupdInit6(:) = .FALSE.
      lupdInit7(:) = .FALSE.
      lupdInit8(:) = .FALSE.
      lhaveGrid = .TRUE.
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the grid spacing of the solver.  This must be called to properly
!>           start the solver.
!>
!>    @param[in] dzIn      Mesh spacing in z (meters).  This must be positive.
!>    @param[in] dxIn      Mesh spacing in x (meters).  This must be positive.
!>    @param[in] dyIn      Mesh spacing in y (meters).  This must be positive.
!>
!>    @param[out] ierr     0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_setGridSpacingF(dzIn, dxIn, dyIn, ierr) &
                 BIND(C, NAME='fteik_setGridSpacingF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : lhaveGridSpacing, dx, dy, dz, zero
      USE FTEIK_UTILS64F, ONLY : fteik_meshConstants64fF
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN) :: dzIn, dxIn, dyIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      lhaveGridSpacing = .FALSE.
      dz = zero
      dx = zero
      dy = zero
      IF (dzIn <= zero .OR. dxIn <= zero .OR. dyIn <= zero) THEN
         IF (dzIn <= zero) WRITE(*,*) 'fteik_setGridSpacing: dz is too small', dzIn
         IF (dxIn <= zero) WRITE(*,*) 'fteik_setGridSpacing: dx is too small', dxIn
         IF (dyIn <= zero) WRITE(*,*) 'fteik_setGridSpacing: dy is too small', dyIn
         ierr = 1
         RETURN
      ENDIF
      dz = dzIn
      dx = dxIn
      dy = dyIn
      lhaveGridSpacing = .TRUE.
      CALL fteik_meshConstants64fF()
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the source location in the model.
!>
!>    @param[in] zs     z source location (meters).
!>    @param[in] xs     x source location (meters).
!>    @param[in] ys     y source location (meters).
!>
!>    @param[out] ierr  0 indicates success.
!>
!>    @author Mark Noble, Alexandrine Gesret, and Ben Baker
!>
!>    @copyright CeCILL-3
      SUBROUTINE fteik_setSourceLocationF(zs, xs, ys, ierr) &
                 BIND(C, NAME='fteik_setSourceLocationF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : dx, dy, dz, nx, ny, nz, &
                                 xsi, ysi, zsi, xsa, ysa, zsa, &
                                 x0, y0, z0, lhaveSource, & 
                                 ijkv1, ijkv2, ijkv3, ijkv4, &
                                 ijkv5, ijkv6, ijkv7, ijkv8, &
                                 lupdInit1, lupdInit2, lupdInit3, lupdInit4, &
                                 lupdInit5, lupdInit6, lupdInit7, lupdInit8, &
                                 levelPtr, nLevels
      USE FTEIk_UTILS64F, ONLY : fteik_setUpdateNodesF
      USE FTEIK_UTILS64F, ONLY : DBL_EPSILON, one, zero
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN) :: xs, ys, zs
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) xmax, ymax, zmax, xsrc, ysrc, zsrc
      INTEGER(C_INT) ierrs(8)
      REAL(C_DOUBLE), PARAMETER :: perturbSource = 0.0001d0
      lhaveSource = .FALSE.
      ierr = 0
      zsi = 0
      xsi = 0
      ysi = 0
      zsa =-one
      xsa =-one
      ysa =-one
      ! Create a relative position in the grid
      zsrc = zs - z0
      xsrc = xs - x0
      ysrc = ys - y0
      zmax = DBLE(nz - 1)*dz
      xmax = DBLE(nx - 1)*dx
      ymax = DBLE(ny - 1)*dy
      IF (zsrc < zero .OR. zsrc > zmax .OR. &
          xsrc < zero .OR. xsrc > xmax .OR. &
          ysrc < zero .OR. ysrc > ymax) THEN
         IF (zsrc < zero .OR. zsrc > zmax) THEN
            WRITE(*,*) 'fteik_setSourceLocationF: zsrc is out of bounds', zsrc, zero, zmax
         ENDIF 
         IF (xsrc < zero .OR. xsrc > xmax) THEN
            WRITE(*,*) 'fteik_setSourceLocationF: xsrc is out of bounds', xsrc, zero, xmax
         ENDIF
         IF (ysrc < zero .OR. ysrc > ymax) THEN
            WRITE(*,*) 'fteik_setSourceLocationF: ysrc is out of bounds', ysrc, zero, ymax
         ENDIF 
         ierr = 1
         RETURN
      ENDIF
      ! Convert source position to grid position (Fortran indexing)
      zsa = (zsrc/dz) + one 
      xsa = (xsrc/dx) + one
      ysa = (ysrc/dy) + one
      ! If a source falls on a grid point then this will push it into a cell
      IF (ABS(zsa - one) < DBL_EPSILON) zsa = zsa + perturbSource
      IF (zsa >= DBLE(nz))              zsa = zsa - perturbSource
      IF (ABS(xsa - one) < DBL_EPSILON) xsa = xsa + perturbSource
      IF (xsa >= DBLE(nx))              xsa = xsa - perturbSource
      IF (ABS(ysa - one) < DBL_EPSILON) ysa = ysa + perturbSource
      IF (ysa >= DBLE(ny))              ysa = ysa - perturbSource
      zsi = INT(zsa)
      xsi = INT(xsa)
      ysi = INT(ysa)
      ! verify bounds make sense
      IF (zsi < 1 .OR. zsi > nz .OR. &
          xsi < 1 .OR. xsi > nx .OR. &
          ysi < 1 .OR. ysi > ny) THEN
         WRITE(*,*) 'fteik_setSourceLocationF: Point (',zsi,xsi,ysi, ') ouf of bounds'
         ierr = 1
         RETURN
      ENDIF
      IF (zsi > nz - 1 .OR. xsi > nx - 1 .OR. ysi > ny - 1) THEN
         WRITE(*,*) 'fteik_setSourceLocationF: Warning solver may segfault'
      ENDIF
      ! set the update grid
      lhaveSource = .TRUE.
      CALL fteik_setUpdateNodesF(1, nlevels, .TRUE., levelPtr, ijkv1, lupdInit1, ierrs(1))
      CALL fteik_setUpdateNodesF(2, nlevels, .TRUE., levelPtr, ijkv2, lupdInit2, ierrs(2))
      CALL fteik_setUpdateNodesF(3, nlevels, .TRUE., levelPtr, ijkv3, lupdInit3, ierrs(3))
      CALL fteik_setUpdateNodesF(4, nlevels, .TRUE., levelPtr, ijkv4, lupdInit4, ierrs(4))
      CALL fteik_setUpdateNodesF(5, nlevels, .TRUE., levelPtr, ijkv5, lupdInit5, ierrs(5))
      CALL fteik_setUpdateNodesF(6, nlevels, .TRUE., levelPtr, ijkv6, lupdInit6, ierrs(6))
      CALL fteik_setUpdateNodesF(7, nlevels, .TRUE., levelPtr, ijkv7, lupdInit7, ierrs(7))
      CALL fteik_setUpdateNodesF(8, nlevels, .TRUE., levelPtr, ijkv8, lupdInit8, ierrs(8))
      IF (MAXVAL(ABS(ierrs)) /= 0) THEN 
         WRITE(*,*) 'fteik_setSourceLocationF: Error setting update nodes'
         ierr = 1
         lhaveSource = .FALSE.
         RETURN
      ENDIF
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the travel times analytically near the source.
!>
!>    @param[in,out] ttimes    On input the travel times have been initialized to a
!>                             large number.
!>                             On output the travel times near the source have been
!>                             computed analytically.  This is a vector of dimension
!>                             [nz*nx*ny].
!>    @param[out] ierr         0 indicates success.
!>
!>    @author Mark Noble, Alexandrine Gesret, and Ben Baker
!>
!>    @version 2
!>
!>    @date July 2017
!>
!>    This now fills the travel-times in a cache friendlier way and now requires
!>    the source with fteik_setSourceLocationF prior to calling.
!> 
      SUBROUTINE fteik_initializeTravelTimesNearSourceF(ttimes, ierr)      &
                 BIND(C, NAME='fteik_initializeTravelTimesNearSourceF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_tAna64fF, grid2indexF
      USE FTEIK_UTILS64F, ONLY : dx, dy, dz, szero, xsa, ysa, zsa, &
                                 xsi, ysi, zsi, nz, nzx, lhaveSource, lhaveSlownessModel
      IMPLICIT NONE
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: ttimes
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) ts8(8)
      INTEGER(C_INT) dest(8), i
      !dir$ attributes align:64 :: ts8, dest
      ierr = 0
      IF (.NOT.lhaveSource) THEN
         WRITE(*,*) 'fteik_initializeTravelTimesNearSourceF: Source not yet set'
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.lhaveSlownessModel) THEN
         WRITE(*,*) 'fteik_initializeTravelTimesNearSourceF: Slowness model not yet set'
         ierr = 1
         RETURN
      ENDIF
      ts8(1) = fteik_tAna64fF(zsi,   xsi,   ysi,              &
                              dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(2) = fteik_tAna64fF(zsi+1, xsi,   ysi,              &
                              dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(3) = fteik_tAna64fF(zsi,   xsi+1, ysi,              &
                              dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(4) = fteik_tAna64fF(zsi+1, xsi+1, ysi,              &
                              dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(5) = fteik_tAna64fF(zsi,   xsi,   ysi+1,            &
                              dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(6) = fteik_tAna64fF(zsi+1, xsi,   ysi+1,            &
                              dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(7) = fteik_tAna64fF(zsi,   xsi+1, ysi+1,            &
                              dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(8) = fteik_tAna64fF(zsi+1, xsi+1, ysi+1,            &
                              dz, dx, dy, zsa, xsa, ysa, szero)
      dest(1) = grid2indexF(zsi,   xsi  , ysi, nz, nzx)
      dest(2) = grid2indexF(zsi+1, xsi  , ysi, nz, nzx)
      dest(3) = grid2indexF(zsi  , xsi+1, ysi, nz, nzx)
      dest(4) = grid2indexF(zsi+1, xsi+1, ysi, nz, nzx)
      dest(5) = grid2indexF(zsi,   xsi  , ysi+1, nz, nzx)
      dest(6) = grid2indexF(zsi+1, xsi  , ysi+1, nz, nzx)
      dest(7) = grid2indexF(zsi,   xsi+1, ysi+1, nz, nzx)
      dest(8) = grid2indexF(zsi+1, xsi+1, ysi+1, nz, nzx)
      DO 1 i=1,8
         ttimes(dest(i)) = ts8(i)
         !print *, ts8(i)
    1 CONTINUE 
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the slowness near the source.
!>
!>    @result 0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_setSzeroF(ierr) BIND(C, NAME='fteik_setSzeroF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : slow, zsi, xsi, ysi, &
                                 lhaveSource, lhaveSlownessModel, &
                                 nzm1, nzm1_nxm1
      USE FTEIK_UTILS64F, ONlY : szero, szero2
      USE FTEIK_UTILS64F, ONLY : zero
      USE FTEIK_UTILS64F, ONLY : velGrid2indexF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) indx
      szero = zero
      ierr = 1
      IF (.NOT.lhaveSource) THEN
         WRITE(*,*) 'fteik_setSzeroF: Source not yet initialized'
         RETURN
      ENDIF
      IF (.NOT.lhaveSlownessModel) THEN
         WRITE(*,*) 'fteik_setSzeroF: Slowness model not yet set'
         RETURN
      ENDIF
      ierr = 0
      indx = velGrid2indexF(zsi, xsi, ysi, nzm1, nzm1_nxm1)
      szero = slow(indx) 
      szero2 = szero*szero
print *, indx, szero
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Converts the spherical to epsilon transition (epsS2C) in grid points
!>           to a transition in seconds.  This is for internal solver use only.
!>
!>    @param[out] ierr   0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_setEpsSolverF(ierr) &
                 BIND(C, NAME='fteik_setEpsSolverF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : dx, dy, dz, epsS2C, epsSolver, lhaveSource, &
                                 lhaveSlownessModel, szero
      USE FTEIK_UTILS64F, ONLY : zero
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 1
      epsSolver = zero
      IF (.NOT.lhaveSource) THEN 
         WRITE(*,*) 'fteik_setEpsSolver: Source not yet initialized'
         RETURN
      ENDIF
      IF (.NOT.lhaveSlownessModel) THEN 
         WRITE(*,*) 'fteik_setSzeroF: Slowness model not yet set'
         RETURN
      ENDIF
      ierr = 0
      epsSolver = epsS2C*MIN(dz, dx, dy)*szero
      !print *, 'epsSolver, szero:', epsSolver, szero
      RETURN
      END
!                           
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief This is a debugging solver for testing certain subroutines and not
!>           appropriate for general use.
!>
!>    @param[out] ierr   0 indicates success.
!>
!>    @author Ben Baker
!>
      SUBROUTINE fteik_solveEikonalDebugF(ierr) &
                 BIND(C, NAME='fteik_solveEikonalDebugF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_initializeTravelTimesNearSourceF, &
                                 fteik_setSzeroF, fteik_setEpsSolverF, &
                                 fteik_getSweepSigns64fF, fteik_getTravelTimePermF, &
                                 fteik_getSlownessPermF, fteik_getSweepLimitsF, &
                                 fteik_prefetchSlowness64fF, fteik_localSolver64fF, &
                                 fteik_prefetchTravelTimes64fF, &
                                 fteik_meshConstants64fF, grid2indexF
      USE FTEIK_UTILS64F, ONLY : ttimes, slow, lhaveGrid, lhaveTravelTimes, nsweep, &
                                 nx, ny, nz, nzx, nzm1, nzm1_nxm1, xsi, ysi, zsi
      USE FTEIK_UTILS64F, ONLY : dxi, dyi, dzi
      USE FTEIK_UTILS64F, ONLY : FTEIK_HUGE
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) ttWork(8), slowWork(7), sgnrz, sgnrx, sgnry, &
                     sgnrz_dzi, sgnrx_dxi, sgnry_dyi, tupd
      INTEGER(C_INT) slowPerm(19), ttPerm(8), i, j, k, kndx, &
                     kiter, sgntz, sgntx, sgnty, sgnvz, sgnvx, sgnvy, sweep, &
                     x1, x2, y1, y2, z1, z2
      LOGICAL(C_BOOL) linitk
      !DIR$ ATTRIBUTES ALIGN: 64:: slowPerm, ttPerm, slowWork, ttWork
      ierr = 0
      lhaveTravelTimes = .FALSE.
      ! some easy checks
      IF (.NOT. lhaveGrid) THEN
         WRITE(*,*) 'fteik_solveEikonalDebugF: Grid not yet set'
         ierr = 1
         RETURN
      ENDIF
      ! don't let mesh constants go out of scope
      CALL fteik_meshConstants64fF()
      ! initialize travel times to a big number
      ttimes(:) = FTEIK_HUGE
      ! set the slowness around the source for the local solver
      CALL fteik_setSzeroF(ierr) 
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solveEikonalDebugF: Error initializing szero'
         RETURN
      ENDIF
      ! set the spherical to cartesian transition in time 
      CALL fteik_setEpsSolverF(ierr) 
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solveEikonalDebugF: Error setting epsSolver'
         RETURN
      ENDIF
      ! Initialize the points around the source
      ttimes(:) = FTEIK_HUGE
      CALL fteik_initializeTravelTimesNearSourceF(ttimes, ierr)
      IF (ierr /= 0) THEN 
         WRITE(*,*) 'fteik_solveEikonalDebugF: Error initializing travel-time'
         RETURN
      ENDIF
      ! iterative loop
      linitk = .TRUE.
      DO kiter=0,nsweep
         IF (kiter > 0) linitk = .FALSE.
         ! sweep directions
         DO sweep=1,8
            ! Get some stencil information for this sweep
            CALL fteik_getSweepSigns64fF(sweep,               &
                                         sgntz, sgntx, sgnty, &
                                         sgnvz, sgnvx, sgnvy, &
                                         sgnrz, sgnrx, sgnry)
            ! Some derivative items
            sgnrz_dzi = sgnrz*dzi
            sgnrx_dxi = sgnrx*dxi
            sgnry_dyi = sgnry*dyi
            ! Get the slowness and travel time permutations
            CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
            CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
            ! Get the loop limits
            CALL fteik_getSweepLimitsF(sweep, linitk,          &
                                       nz, nx, ny,             &
                                       zsi, xsi, ysi,          &
                                       z1, z2, x1, x2, y1, y2, &
                                       ierr)
            DO k=y1,y2,sgnty
               DO j=x1,x2,sgntx
                  DO i=z1,z2,sgntz
                     ! prefetch the slownesses and traveltimes
                     CALL fteik_prefetchSlowness64fF(i, j, k, nzm1, nzm1_nxm1,    &
                                                     sgnvz, sgnvx, sgnvy,    &
                                                     slowPerm, slow, slowWork)
                     CALL fteik_prefetchTravelTimes64fF(i, j, k, nz, nzx,     &
                                                        sgntz, sgntx, sgnty,  &
                                                        ttPerm, ttimes, ttWork)
!print *, i, j, k
!print *, sngl(ttWork)
                     ! apply the finite differencing
                     tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &
                                                  i, j, k,                       &
                                                  sgntz, sgntx, sgnty,           &
                                                  sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
                     ! update ttimes
                     kndx = grid2indexF(i, j, k, nz, nzx)
                     ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
                  ENDDO
               ENDDO
            ENDDO
!return
         ENDDO ! loop on sweep directions 
      ENDDO ! end iterative loop on sweeps 
print *, minval(ttimes), maxval(ttimes)
      lhaveTravelTimes = .TRUE.
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Solves the eikonal equation using the level set based fast sweeping method.
!>
!>    @result 0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solveEikonalF(ierr) BIND(C, NAME='fteik_solveEikonalF')
      USE ISO_C_BINDING 
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the velocity model on the solver from a float velocity model.
!>
!>    @param[in] nv     Number of cells.  This must be equal to [(nz-1)*(nx-1)*(ny-1)].
!>
!>    @param[in] vel4   Velocity model (s/m).  This is a vector of dimension
!>                      [(nz-1) x (nx-1) x (ny-1)] with leading dimension 1 of (nz-1)
!>                      and leading dimension 2 of (nx-1).
!>
!>    @param[out] ierr  0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_setVelocityModel32fF(nv, vel4, ierr) &
                 BIND(C, NAME='fteik_setVelocityModel32fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : slow, ncell, lhaveSlownessModel, one
      INTEGER(C_INT), INTENT(IN) :: nv
      REAL(C_FLOAT), INTENT(IN) :: vel4(nv)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) i
      ierr = 0
      lhaveSlownessModel = .FALSE.
      IF (nv /= ncell) THEN 
         WRITE(*,*) 'fteik_setVelocityModel32fF: ERROR - ncell /= nv', ncell, nv
         ierr = 1
         RETURN
      ENDIF
      IF (nv < 1) THEN 
         WRITE(*,*) 'fteik_setVelocityMode32flF: ERROR - no cells in vel', nv
         ierr = 1
         RETURN
      ENDIF
      IF (MINVAL(vel4) <= 0.0) THEN
         WRITE(*,*) 'fteik_setVelocityModel32fF: ERROR - all velocities must be positive'
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.ALLOCATED(slow)) ALLOCATE(slow(ncell))
      DO 1 i=1,ncell
         slow(i) = one/DBLE(vel4(i))
    1 CONTINUE
      lhaveSlownessModel = .TRUE.
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility routine for extracting the travel times from the solver.
!>
!>    @param[in] ng     Number of grid points in tt.  This must be nz*nx*ny.
!>
!>    @param[out] tt    Travel times at each grid point.  This is a vector of 
!>                      dimension [nz*nx*ny] where the first leading dimension is
!>                      nz and the second leading dimension is nx.
!>    @param[out] ierr  0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_getTravelTimes64fF(ng, tt, ierr) &
                 BIND(C, NAME='fteik_getTravelTimes64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : ttimes, ngrd, lhaveTravelTimes, zero
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: ng
      REAL(C_DOUBLE), INTENT(OUT) :: tt(ng)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      IF (ng /= ngrd) THEN
         WRITE(*,*) 'fteik_getTravelTimes64fF: ERROR - ng /= ngrd', ng, ngrd
         ierr = 1
         RETURN
      ENDIF
      IF (ng < 1) THEN
         WRITE(*,*) 'fteik_getTravelTiems64fF: ERROR - no grid points', ng
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT. lhaveTravelTimes) THEN
         WRITE(*,*) 'fteik_getTravelTimes64fF: ERROR - travel times not yet computed'
         tt(:) = zero
         ierr = 1
      ENDIF
      tt(1:ngrd) = ttimes(1:ngrd)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the velocity model on the solver.
!>
!>    @param[in] nv     Number of cells.  This must be equal to [(nz-1)*(nx-1)*(ny-1)].
!>
!>    @param[in] vel    Velocity model (s/m).  This is a vector of dimension
!>                      [(nz-1) x (nx-1) x (ny-1)] with leading dimension 1 of (nz-1)
!>                      and leading dimension 2 of (nx-1).
!>
!>    @param[out] ierr  0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_setVelocityModel64fF(nv, vel, ierr) &
                 BIND(C, NAME='fteik_setVelocityModel64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : slow, lhaveSlownessModel, ncell, one, zero
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: nv
      REAL(C_DOUBLE), INTENT(IN) :: vel(nv)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) i
      ierr = 0
      lhaveSlownessModel = .FALSE.
      IF (nv /= ncell) THEN
         WRITE(*,*) 'fteik_setVelocityModel64fF: ERROR - ncell /= nv', ncell, nv
         ierr = 1
         RETURN
      ENDIF
      IF (nv < 1) THEN
         WRITE(*,*) 'fteik_setVelocityModel64fF: ERROR - no cells in vel', nv
         ierr = 1
         RETURN
      ENDIF
      IF (MINVAL(vel) <= zero) THEN
         WRITE(*,*) 'fteik_setVelocityModel64fF: ERROR - all velocities must be positive'
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.ALLOCATED(slow)) ALLOCATE(slow(ncell))
      DO 1 i=1,ncell
         slow(i) = one/vel(i)
    1 CONTINUE 
      lhaveSlownessModel = .TRUE.
      RETURN
      END SUBROUTINE
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
         WRITE(*,*) 'fteiK_get_sweepLimitsF: Invalid sweep', sweep
         ierr = 1
      ENDIF 
      RETURN
      END
!========================================================================================!
!                                                                                        !
!>    @brief Gets the permutation for extracting the slowness for a given sweep.
!>
!>    @param[in] sweep     Sweep number.  This must be in the range [1,8].
!>
!>    @param[out] perm     This is permutation so that the slowness prefetch can
!>                         extract nodes in increasing grid point number as to 
!>                         improve cache-coherency.  This is a vector of dimension [19].
!>    @param[out] ierr     0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_getSlownessPermF(sweep, perm, ierr) &
                 BIND(C, NAME='fteik_getSlownessPermF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : v2l1, v2l2, v2l3, v2l4, v2l5, v2l6, v2l7, v2l8
      INTEGER(C_INT), INTENT(IN) :: sweep
      INTEGER(C_INT), INTENT(OUT) :: perm(19), ierr
      ierr = 0
      IF (sweep == 1) THEN
         perm(:) = v2l1(:)
      ELSE IF (sweep == 2) THEN 
         perm(:) = v2l2(:)
      ELSE IF (sweep == 3) THEN
         perm(:) = v2l3(:)
      ELSE IF (sweep == 4) THEN
         perm(:) = v2l4(:)
      ELSE IF (sweep == 5) THEN
         perm(:) = v2l5(:)
      ELSE IF (sweep == 6) THEN
         perm(:) = v2l6(:)
      ELSE IF (sweep == 7) THEN 
         perm(:) = v2l7(:)
      ELSE IF (sweep == 8) THEN 
         perm(:) = v2l8(:)
      ELSE 
         WRITE(*,*) 'fteik_getSlownessPermF: Invalid sweep number'
         ierr = 1
      ENDIF
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gets the permutation for extracting the travel times for a given sweep.
!>
!>    @param[in] sweep     Sweep number.  This must be in the range [1,8].
!>
!>    @param[out] perm     This is permutation so that the travel time prefetch can
!>                         extract nodes in increasing grid point number as to 
!>                         improve cache-coherency.  This is a vector of dimension [8].
!>    @param[out] ierr     0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_getTravelTimePermF(sweep, perm, ierr) &
                 BIND(C, NAME='fteik_getTravelTimePermF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : ttPerm1, ttPerm2, ttPerm3, ttPerm4, &
                                 ttPerm5, ttPerm6, ttPerm7, ttPerm8
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: sweep
      INTEGER(C_INT), INTENT(OUT) :: perm(8), ierr
      ierr = 0
      IF (sweep == 1) THEN
         perm(:) = ttPerm1(:)
      ELSE IF (sweep == 2) THEN
         perm(:) = ttPerm2(:)
      ELSE IF (sweep == 3) THEN
         perm(:) = ttPerm3(:)
      ELSE IF (sweep == 4) THEN
         perm(:) = ttPerm4(:)
      ELSE IF (sweep == 5) THEN
         perm(:) = ttPerm5(:)
      ELSE IF (sweep == 6) THEN
         perm(:) = ttPerm6(:)
      ELSE IF (sweep == 7) THEN
         perm(:) = ttPerm7(:)
      ELSE IF (sweep == 8) THEN
         perm(:) = ttPerm8(:)
      ELSE
         WRITE(*,*) 'fteik_getTravelTimePermF: Invalid sweep number'
         ierr = 1
      ENDIF
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief This applies the finite difference operators to a grid point. 
!>
!>    @param[in] tt         This is an array of length 8 holding the travel-times (s)
!>                          to be used by the FD stencils.  It should have been
!>                          filled in tandem with fteik_prefetchTravelTimes64fF
!>                          and fteik_getTravelTimePermF for the appropriate sweep.
!>    @param[in] slowLoc    This is an array of length 7 holding the slownesses (s/m) 
!>                          to be used by the FD stencils.  It should have been
!>                          filled in tandeom with fteik_prefetchSlowness64fF
!>                          and fteik_getSlownessPermF for the appropriate sweep.
!>    @param[in] linitk     If true then this flag indicates that this is the first
!>                          iteration of the Gauss-Seidel method. \n
!>                          Otherwise, this is a higher iteration.
!>    @param[in] i          Grid point in z.  This is Fortran indexed.
!>    @param[in] j          Grid point in x.  This is Fortran indexed.
!>    @param[in] k          Grid point in y.  This is Fortran indexed.
!>    @param[in] sgntz      Sign to add to i'th grid point for analytic computations.
!>    @param[in] sgntx      Sign to add to j'th grid point for analytic computations.
!>    @param[in] sgnty      Sign to add to k'th grid point for analytic computations.
!>    @param[in] sgnrz_dzi  sgnrz/dz
!>    @param[in] sgnrx_dxi  sgnrx/dx
!>    @param[in] sgnry_dyi  sgnry/dy
!>
!>    @result The minimal travel-time of all the finite difference stensils 
!>            at the (i,j,k)'th grid node.
!>
!>    @author Mark Noble, Alexandrine Gesret, and Ben Baker
!>
!>    @copyright Mines ParisTech, France distributed under CeCILL-C license.
!>
!>    @version 3
!>
!>    @date July 2017
!> 
!>    Note, this version has been substantially modified from the original
!>    Include_FTeik3d_2.0.f of Nobel and Gesret.  In particular: \n
!>     - Made this into a function instead of a segment of code to be included. \n
!>     - Reworked preamble to be consistent with Doxygen. \n
!>     - Changed BIG to FTEIK_HUGE which is defined in the module. \n
!>     - Changed things like **0.5 to SQRT to avoid the potential use of logarithms. \n
!>     - Tried to precompute as many terms as possible. \n
!>     - Operate on a local velocity/traveltime stencil to provide user opportunity
!>       to emphasize cache-coherency; particularly in level set sweeping.
!>             
      PURE REAL(C_DOUBLE)  &
      FUNCTION fteik_localSolver64fF(tt, slowLoc, linitk,             &
                                     i, j, k,                         &
                                     sgntz, sgntx, sgnty,             &
                                     sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
      RESULT(tupd) BIND(C, NAME='fteik_localSolver64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_localSolverExplicit64fF
      IMPLICIT NONE 
      REAL(C_DOUBLE), INTENT(IN) :: tt(8), slowLoc(7)
      REAL(C_DOUBLE), INTENT(IN) :: sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT), INTENT(IN) :: i, j, k, sgntx, sgnty, sgntz
      LOGICAL(C_BOOL), INTENT(IN) :: linitk
      tupd = fteik_localSolverExplicit64fF(tt(1), tt(2), tt(3), tt(4),         &
                                           tt(5), tt(6), tt(7), tt(8),         &
                                           slowLoc(1), slowLoc(2), slowLoc(3), slowLoc(4), &
                                           slowLoc(5), slowLoc(6), slowLoc(7),             &
                                           linitk,                                         &
                                           i, j, k,                          &
                                           sgntz, sgntx, sgnty,              &
                                           sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
      RETURN
      END FUNCTION



      PURE REAL(C_DOUBLE)                                                     &
      FUNCTION fteik_localSolverExplicit64fF(tv, te, tn, tev,                 &
                                             ten, tnv, tnve, tt0,             &
                                             slow1, slow2, slow3, slow4,      &
                                             slow5, slow6, slow7,             &
                                             linitk,                          &
                                             i, j, k,                         &
                                             sgntz, sgntx, sgnty,             &
                                             sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
      RESULT(tupd) BIND(C, NAME='fteik_localSolverExplicit64fF')
      !$OMP DECLARE SIMD(fteik_localSolverExplicit64fF) &
      !$OMP UNIFORM(linitk, sgntz, sgntx, sgnty, sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_tAnaD64fF, fteik_tAna64fF
      USE FTEIK_UTILS64F, ONLY : FTEIK_HUGE
      USE FTEIK_UTILS64F, ONLY : dx, dy, dz, xsa, ysa, zsa, &
                                 dx2i, dy2i, dz2i, &
                                 dsum, dsumi, epsSolver, &
                                 dz2dx2, dz2dy2, dx2dy2, &
                                 dz2i_dx2i, dz2i_dy2i, dx2i_dy2i, &
                                 dz2i_p_dx2i, dz2i_p_dy2i, dx2i_p_dy2i, &
                                 dz2i_p_dy2i_inv, dz2i_p_dx2i_inv, dx2i_p_dy2i_inv, &
                                 szero, szero2
      USE FTEIK_UTILS64F, ONLY : xsi, ysi, zsi
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN), VALUE :: tv, te, tn, tev, ten, tnv, tnve, tt0
      REAL(C_DOUBLE), INTENT(IN), VALUE :: slow1, slow2, slow3, slow4, &
                                           slow5, slow6, slow7
      REAL(C_DOUBLE), INTENT(IN), VALUE :: sgnrz_dzi, sgnrx_dxi, sgnry_dyi !, &
      !                                     dx, dy, dz, xsa, ysa, zsa, &
      !                                     dx2i, dy2i, dz2i, &
      !                                     dsum, dsumi, epsSolver, &
      !                                     dz2dx2, dz2dy2, dx2dy2, &
      !                                     dz2i_dx2i, dz2i_dy2i, dx2i_dy2i, &
      !                                     dz2i_p_dx2i, dz2i_p_dy2i, dx2i_p_dy2i, & 
      !                                     dz2i_p_dy2i_inv, dz2i_p_dx2i_inv,      &
      !                                     dx2i_p_dy2i_inv, szero, szero2
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, sgntx, sgnty, sgntz !, xsi, ysi, zsi
      LOGICAL(C_BOOL), INTENT(IN), VALUE :: linitk
      ! local variables
      REAL(C_DOUBLE) apoly, bpoly, cpoly, dpoly,                                &
                     four_sref2, sref, sref2,                                   &
                     t0c, t1, t1d, t1d_t2d_min, t2, t2d, t3, t3d,               &
                     ta, tb, tab, tab2, tac, tac2, tbc, tbc2, tc,               &
                     taue, tauev, tauen, taun, taunv, taunve, tauv,             &
                     tmin, txc, tyc, tzc
      LOGICAL(C_BOOL) ltmin_eps, ltmin_eps_notk
      ! Initialze the output to a large value
      tupd = FTEIK_HUGE
      ! Get the local times of surrounding points
      !tv   = tt(1) ! tt(i-sgntz,j,k)
      !te   = tt(2) ! tt(i,j-sgntx,k)
      !tn   = tt(3) ! tt(i,j,k-sgnty)
      !tev  = tt(4) ! tt(i-sgntz,j-sgntx,k)
      !ten  = tt(5) ! tt(i,j-sgntx,k-sgnty)
      !tnv  = tt(6) ! tt(i-sgntz,j,k-sgnty)
      !tnve = tt(7) ! tt(i-sgntz,j-sgntx,k-sgnty)
      !tt0  = tt(8) ! tt(i, j, k)
      ! Check time to see when to switch to plane approximation
      tmin = MIN(tv, te, tn)
      ! This gets computed over and over and over...
      ltmin_eps = .FALSE.
      ltmin_eps_notk = .FALSE.
      IF (tmin > epsSolver) ltmin_eps = .TRUE.
      IF (ltmin_eps .OR. .NOT.linitk) ltmin_eps_notk = .TRUE.
      ! Get analytical solution, if using pertubation
      t0c = 0.d0
      tzc = 0.d0
      txc = 0.d0
      tyc = 0.d0
      !tauv = FTEIK_HUGE
      !taue = FTEIK_HUGE
      !taun = FTEIK_HUGE
      !taun = FTEIK_HUGE
      !tauev = FTEIK_HUGE
      !tauen = FTEIK_HUGE
      !taunv = FTEIK_HUGE
      !taunve = FTEIK_HUGE
      IF (tmin <= epsSolver .OR. linitk) THEN
         !t0c = t_anad64f(tzc, txc, tyc, i, j, k)
         CALL fteik_tAnaD64fF(t0c, tzc, txc, tyc, i, j, k,    &
                              dz, dx, dy, zsa, xsa, ysa, szero)
         ! Convert times into pertubations
         tauv   = tv   - fteik_tAna64fF(i-sgntz, j,       k,            &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
         taue   = te   - fteik_tAna64fF(i,       j-sgntx, k,            &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
         taun   = tn   - fteik_tAna64fF(i,       j,       k-sgnty,      &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
         tauev  = tev  - fteik_tAna64fF(i-sgntz, j-sgntx, k,            &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
         tauen  = ten  - fteik_tAna64fF(i,       j-sgntx, k-sgnty,      &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
         taunv  = tnv  - fteik_tAna64fF(i-sgntz, j,       k-sgnty,      &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
         taunve = tnve - fteik_tAna64fF(i-sgntz, j-sgntx, k-sgnty,      &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
      ENDIF
      !--------------------1D operators, (refracted times),set times to BIG--------------!
      !t1d = FTEIK_HUGE ! Big
      !t1  = FTEIK_HUGE ! Big
      !t2  = FTEIK_HUGE ! Big
      !t3  = FTEIK_HUGE ! Big
      ! V plane
      t1 = tv + dz*slow1 !t1 = tv + dz*slowLoc(1)
      ! WE plane
      t2 = te + dx*slow2 !t2 = te + dx*slowLoc(2)
      ! NS plane
      t3 = tn + dy*slow3 !t3 = tn + dy*slowLoc(3)
      t1d = MIN(t1, t2, t3)
      !-------------------------------------2D operators---------------------------------!
      t2d = FTEIK_HUGE  !Big;
      t1  = FTEIK_HUGE  !Big;
      t2  = FTEIK_HUGE  !Big;
      t3  = FTEIK_HUGE  !Big;
      ! ZX (VE) plane
      sref = slow4 !sref = slowLoc(4)
      IF ( (tv < te + dx*sref) .AND. (te < tv + dz*sref) ) THEN
         sref2 = sref*sref
         IF (ltmin_eps_notk .OR. k /= ysi) THEN
            ta = tev + te - tv
            tb = tev - te + tv
            tab = ta - tb
            tab2 = tab*tab
            four_sref2 = 4.d0*sref2
            t1 = ( (tb*dz2i + ta*dx2i) &
                  + SQRT(four_sref2*dz2i_p_dx2i - dz2i_dx2i*tab2) )*dz2i_p_dx2i_inv
         ELSE
            ta = tauev + taue - tauv
            tb = tauev - taue + tauv
            apoly = dz2i_p_dx2i
            bpoly = 4.d0*(sgnrx_dxi*txc + sgnrz_dzi*tzc) - 2.d0*(ta*dx2i + tb*dz2i)
            cpoly = ((ta*ta)*dx2i) + ((tb*tb)*dz2i)            &
                  - 4.d0*(sgnrx_dxi*txc*ta + sgnrz_dzi*tzc*tb) &
                  + 4.d0*(szero2 - sref2 + tyc*tyc)
            dpoly = bpoly*bpoly - 4.d0*apoly*cpoly
            IF (dpoly >= 0.d0) t1 = (SQRT(dpoly) - bpoly)/(2.d0*apoly) + t0c
            IF (t1 < tv .OR. t1 < te) t1 = FTEIK_HUGE
         ENDIF 
      ENDIF
      ! ZY (VN) plane
      sref = slow5 !sref = slowLoc(5)
      IF ( (tv < tn + dy*sref) .AND. (tn < tv + dz*sref) ) THEN
         sref2 = sref*sref
         IF (ltmin_eps_notk .OR. j /= xsi) THEN
            ta = tv - tn + tnv
            tb = tn - tv + tnv
            tab = ta - tb
            tab2 = tab*tab
            four_sref2 = 4.d0*sref2
            t2 = ( (ta*dz2i + tb*dy2i) &
                  + SQRT(four_sref2*dz2i_p_dy2i - dz2i_dy2i*tab2) )*dz2i_p_dy2i_inv
         ELSE
            ta = tauv - taun + taunv ! Z
            tb = taun - tauv + taunv ! Y
            apoly = dz2i_p_dy2i
            bpoly = 4.d0*(sgnry_dyi*tyc + sgnrz_dzi*tzc) - 2.d0*(ta*dz2i + tb*dy2i)
            cpoly = ((ta*ta)*dz2i) + ((tb*tb)*dy2i)            &
                  - 4.d0*(sgnrz_dzi*tzc*ta + sgnry_dyi*tyc*tb) &
                  + 4.d0*(szero2 - sref2 + txc*txc)
            dpoly = bpoly*bpoly - 4.d0*apoly*cpoly
            IF (dpoly >= 0.d0) t2 = (SQRT(dpoly) - bpoly)/(2.d0*apoly) + t0c
            IF (t2 < tv .OR. t2 < tn) t2 = FTEIK_HUGE
         ENDIF
      ENDIF
      ! XY (EN) plane
      sref = slow6 !sref = slowLoc(6)
      IF ( (te < tn + dy*sref) .AND. (tn < te + dx*sref) ) THEN
         sref2 = sref*sref
         IF (ltmin_eps_notk .OR. i /= zsi) THEN
            ta = te - tn + ten
            tb = tn - te + ten
            tab = ta - tb
            tab2 = tab*tab
            four_sref2 = 4.d0*sref2
            t3 = ( (ta*dx2i + tb*dy2i) &
               + SQRT(four_sref2*dx2i_p_dy2i - dx2i_dy2i*tab2) )*dx2i_p_dy2i_inv
         ELSE
            ta = taue - taun + tauen ! X
            tb = taun - taue + tauen ! Y
            apoly = dx2i_p_dy2i
            bpoly = 4.d0*(sgnry_dyi*tyc + sgnrx_dxi*txc) - 2.d0*(ta*dx2i + tb*dy2i)
            cpoly = ((ta*ta)*dx2i)+((tb*tb)*dy2i)              &
                  - 4.d0*(sgnrx_dxi*txc*ta + sgnry_dyi*tyc*tb) &
                  + 4.d0*(szero2 - sref2 + tzc*tzc);
            dpoly = bpoly*bpoly - 4.d0*apoly*cpoly;
            IF (dpoly >= 0.d0) t3 = (SQRT(dpoly) - bpoly)/(2.d0*apoly) + t0c
            IF (t3 < te .OR. t3 < tn) t3 = FTEIK_HUGE
         ENDIF
      ENDIF
      t2d = MIN(t1, t2, t3)
      !-------------------------------------3D operators---------------------------------!
      t3d = FTEIK_HUGE
      sref = slow7 !sref = slowLoc(7)
      t1d_t2d_min = MIN(t1d, t2d)
      IF ( t1d_t2d_min > MAX(tv, te, tn) ) THEN
         sref2 = sref*sref
         IF (ltmin_eps_notk) THEN
            ta = te + 0.5d0*(-tn + ten - tv + tev) - tnv + tnve ! X
            tb = tv + 0.5d0*(-tn + tnv - te + tev) - ten + tnve ! Z
            tc = tn + 0.5d0*(-te + ten - tv + tnv) - tev + tnve ! Y
            tab = ta - tb
            tbc = tb - tc
            tac = ta - tc
            t2 = sref2*dsum*9.d0
            tab2 = tab*tab
            tbc2 = tbc*tbc
            tac2 = tac*tac
            t3 = dz2dx2*tab2 + dz2dy2*tbc2 + dx2dy2*tac2
            IF (t2 >= t3) THEN 
               t1 = tb*dz2i + ta*dx2i + tc*dy2i
               t3d = (t1 + SQRT(t2 - t3))*dsumi
            ENDIF
         ELSE
            ta = taue + 0.5d0*(-taun + tauen - tauv + tauev) - taunv + taunve ! X
            tb = tauv + 0.5d0*(-taun + taunv - taue + tauev) - tauen + taunve ! Z
            tc = taun + 0.5d0*(-taue + tauen - tauv + taunv) - tauev + taunve ! Y
            apoly = dx2i + dz2i + dy2i
            bpoly = 6.d0*(sgnrz_dzi*tzc + sgnry_dyi*tyc + sgnrx_dxi*txc) &
                  - 2.d0*(ta*dx2i + tb*dz2i + tc*dy2i)
            cpoly = ((ta*ta)*dx2i) + ((tb*tb)*dz2i) + ((tc*tc)*dy2i)               &
                  - 6.d0*( sgnrx_dxi*txc*ta + sgnrz_dzi*tzc*tb + sgnry_dyi*tyc*tc) &
                  + 9.d0*(szero2 - sref2)
            dpoly = bpoly*bpoly - 4.d0*apoly*cpoly
            IF (dpoly >= 0.d0) t3d = (SQRT(dpoly) - bpoly)/(2.d0*apoly) + t0c
            IF (t3d < te .OR. t3d < tn .OR. t3d < tv) t3d = FTEIK_HUGE
         ENDIF
      ENDIF
      tupd = MIN(tt0, t1d_t2d_min, t3d)
      !tupd = DMIN1(t1d_t2d_min, t3d)
      RETURN 
      END FUNCTION
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Fetches the travel times corresponding to the stencil at grid point
!>           (i, j, k) for use in the local solver.
!>
!>    @param[in] i          iz'th grid point.  This is Fortran indexed.
!>    @param[in] j          ix'th grid point.  This is Fortran indexed.
!>    @param[in] k          iy'th grid point.  This is Fortran indexed.
!>    @param[in] nz         Number of z grid points in travel time field.
!>    @param[in] nzx        Number of z and x grid points in travel time field (=nx*nz).
!>    @param[in] sgntz      The travel time stencil shift in z for the given sweep.
!>    @param[in] sgntx      The travel time stencil shift in x for the given sweep.
!>    @param[in] sgnty      The travel time stencil shift in y for the given sweep.
!>    @param[in] perm       This is the permuation for the given sweep so that
!>                          the travel times are extracted in increasing grid point
!>                          order.  This is a vector of length 8.
!>    @param[in] ttimes     These are the travel times at all points in the grid.
!>                          This is a vector of dimension [nz*nx*ny].
!>
!>    @param[out] ttLoc     This is the travel-time stencil for the (iz,ix,iy)'th grid
!>                          point appropriate for use in the local solver.  This is a
!>                          vector of dimension [8].
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_prefetchTravelTimes64fF(i, j, k, nz, nzx, &
                                            sgntz, sgntx, sgnty, &
                                            perm, ttimes, ttLoc) &
                 BIND(C, NAME='fteik_prefetchTravelTimes64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      IMPLICIT NONE
      REAL(C_DOUBLE), DIMENSION(:), INTENT(IN) :: ttimes
      INTEGER(C_INT), DIMENSION(:), INTENT(IN) :: perm
      INTEGER(C_INT), INTENT(IN) :: i, j, k, nz, nzx, sgntx, sgnty, sgntz
      REAL(C_DOUBLE), INTENT(OUT) :: ttLoc(8)
      INTEGER(C_INT) iloc(8), kiter
      !DIR$ ATTRIBUTES ALIGN: 64:: iloc
      iloc(1) = grid2indexF(i-sgntz, j,       k,       nz, nzx)
      iloc(2) = grid2indexF(i,       j-sgntx, k,       nz ,nzx)
      iloc(3) = grid2indexF(i,       j,       k-sgnty, nz, nzx)
      iloc(4) = grid2indexF(i-sgntz, j-sgntx, k,       nz, nzx)
      iloc(5) = grid2indexF(i,       j-sgntx, k-sgnty, nz, nzx)
      iloc(6) = grid2indexF(i-sgntz, j,       k-sgnty, nz, nzx)
      iloc(7) = grid2indexF(i-sgntz, j-sgntx, k-sgnty, nz, nzx)
      iloc(8) = grid2indexF(i,       j,       k,       nz, nzx)
      DO 1 kiter=1,8
         ttLoc(perm(kiter)) = ttimes(iloc(perm(kiter)))
    1 CONTINUE
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Fetches the slowness corresponding to the stencil at grid point (i,j,k)
!>           for use in the local solver.
!>
!>    @param[in] i              iz'th grid point.  This is Fortran indexed.
!>    @param[in] j              ix'th grid point.  This is Fortran indexed.
!>    @param[in] k              iy'th grid point.  This is Fortran indexed.
!>    @param[in] nzm1           Number of cells in velocity field in z (=nz - 1).
!>    @param[in] nzm1_nxm1      Number of cells in velocity in x and z (=(nz-1)*(nx-1)).
!>    @param[in] sgnvz          The velocity stencil shift in z for the given sweep.
!>    @param[in] sgnvx          The velocity stencil shift in x for the given sweep.
!>    @param[in] sgnvy          The velocity stencil shift in y for the given sweep.
!>    @param[in] slowPerm       This is the permuation for the given sweep so that
!>                              the slownesses are extracted in increasing grid point
!>                              order This is a vector of length 19.
!>    @param[in] slow           This is the slowness (s/m) in each cell.  This is a
!>                              a vector of dimension [(nz-1)*(nx-1)*(ny-1)].
!>
!>    @param[out] slowStencils  This is the slowness stencil for the (iz,ix,iy)'th grid
!>                              point appropriate for use in the local solver.  This
!>                              is a vector of dimension [7].
!>    @author Ben Baker
!>
!>    @copyright MIT
!
      SUBROUTINE fteik_prefetchSlowness64fF(i, j, k, nzm1, nzm1_nxm1,     &
                                            sgnvz, sgnvx, sgnvy,          &
                                            slowPerm, slow, slowStencils) &
                 BIND(C, NAME='fteik_prefetchSlowness64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : velgrid2indexF
      USE FTEIK_UTILS64F, ONLY : nx, ny, nz
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: slowPerm(19), i, j, k, nzm1, nzm1_nxm1, & 
                                    sgnvz, sgnvx, sgnvy
      REAL(C_DOUBLE), DIMENSION(:), INTENT(IN) :: slow
      REAL(C_DOUBLE), INTENT(OUT) :: slowStencils(7)
      !DIR ATTRIBUTES ALIGN: 64:: slowStencils
      REAL(C_DOUBLE) sloc(19)
      INTEGER(C_INT) indx(19), i1, it, j1, k1,              &
                     max_im1_1, max_jm1_1, max_km1_1,       &
                     min_im1_nzm1, min_jm1_nxm1, min_km1_nym1
      !DIR ATTRIBUTES ALIGN: 64:: indx
      !DIR ATTRIBUTES ALIGN: 64:: sloc
      i1 = i - sgnvz
      j1 = j - sgnvx
      k1 = k - sgnvy
      ! Mins and maxes
      max_im1_1    = MAX(i - 1, 1)
      max_jm1_1    = MAX(j - 1, 1)
      max_km1_1    = MAX(k - 1, 1)
      min_im1_nzm1 = MIN(i, nz - 1)
      min_jm1_nxm1 = MIN(j, nx - 1)
      min_km1_nym1 = MIN(k, ny - 1)
      ! 1D operators - V plane
      indx(1) = velGrid2indexF(i1, max_jm1_1,    max_km1_1,    nzm1, nzm1_nxm1)
      indx(2) = velGrid2indexF(i1, max_jm1_1,    min_km1_nym1, nzm1, nzm1_nxm1)
      indx(3) = velGrid2indexF(i1, min_jm1_nxm1, max_km1_1,    nzm1, nzm1_nxm1)
      indx(4) = velGrid2indexF(i1, min_jm1_nxm1, min_km1_nym1, nzm1, nzm1_nxm1)
      ! WE plane
      indx(5) = velGrid2indexF(max_im1_1,    j1, max_km1_1,    nzm1, nzm1_nxm1)
      indx(6) = velGrid2indexF(min_im1_nzm1, j1, max_km1_1,    nzm1, nzm1_nxm1)
      indx(7) = velGrid2indexF(max_im1_1,    j1, min_km1_nym1, nzm1, nzm1_nxm1)
      indx(8) = velGrid2indexF(min_im1_nzm1, j1, min_km1_nym1, nzm1, nzm1_nxm1)
      ! NS plane
      indx(9) = velGrid2indexF(max_im1_1,    max_jm1_1,    k1, nzm1, nzm1_nxm1)
      indx(10)= velGrid2indexF(max_im1_1,    min_jm1_nxm1, k1, nzm1, nzm1_nxm1)
      indx(11)= velGrid2indexF(min_im1_nzm1, max_jm1_1,    k1, nzm1, nzm1_nxm1)
      indx(12)= velGrid2indexF(min_im1_nzm1, min_jm1_nxm1, k1, nzm1, nzm1_nxm1)
      ! 2D operators - VE plane
      indx(13) = velGrid2indexF(i1, j1, max_km1_1,    nzm1, nzm1_nxm1)
      indx(14) = velGrid2indexF(i1, j1, min_km1_nym1, nzm1, nzm1_nxm1)
      ! ZY plane
      indx(15) = velGrid2indexF(i1, max_jm1_1,    k1, nzm1, nzm1_nxm1)
      indx(16) = velGrid2indexF(i1, min_jm1_nxm1, k1, nzm1, nzm1_nxm1)
      ! XY plane
      indx(17) = velGrid2indexF(max_im1_1,    j1, k1, nzm1, nzm1_nxm1)
      indx(18) = velGrid2indexF(min_im1_nzm1, j1, k1, nzm1, nzm1_nxm1)
      ! 3D Operator
      indx(19) = velGrid2indexF(i1, j1, k1, nzm1, nzm1_nxm1)
      DO 1 it=1,19
         sloc(slowPerm(it)) = slow(indx(slowPerm(it)))
    1 CONTINUE
      ! And get the slownesses for the stencils
      !slowStencils(1) = MIN(MIN(sloc(1), sloc(2)),  MIN(sloc(3), sloc(4)))
      !slowStencils(2) = MIN(MIN(sloc(5), sloc(6)),  MIN(sloc(7), sloc(8)))
      !slowStencils(3) = MIN(MIN(sloc(9), sloc(10)), MIN(sloc(11), sloc(12)))
      slowStencils(1) = MIN(sloc(1), sloc(2),  sloc(3),  sloc(4))
      slowStencils(2) = MIN(sloc(5), sloc(6),  sloc(7),  sloc(8))
      slowStencils(3) = MIN(sloc(9), sloc(10), sloc(11), sloc(12))
      slowStencils(4) = MIN(sloc(13), sloc(14))
      slowStencils(5) = MIN(sloc(15), sloc(16))
      slowStencils(6) = MIN(sloc(17), sloc(18))
      slowStencils(7) = sloc(19)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Converts the (i, j, k)'th grid indices to the corresonding grid index
!>           in the travel time array. 
!>
!>    @param[in] i    iz'th grid point.  This is Fortran indexed.
!>    @param[in] j    ix'th grid point.  This is Fortran indexed.
!>    @param[in] k    iy'th grid point.  This is Fortran indexed.
!>    @param[in] nz   Number of z grid points in grid.
!>    @param[in] nzx  Number of z and x grid points in grid (i.e.z nzx = nz x nx).
!>
!>    @result The grid point, g = (k-1)*nx*nz + (j-1)*nz + i, corresponiding to the
!>            (i, j, k)'th grid indices.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE INTEGER(C_INT) FUNCTION grid2indexF(i, j, k, nz, nzx) &
                          BIND(C, NAME='fteik_grid2indexF')
      !$OMP DECLARE SIMD(grid2indexF) UNIFORM(nz, nzx)
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nzx
      grid2indexF = (k - 1)*nzx + (j - 1)*nz + i
      RETURN
      END FUNCTION
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Converts from the grid index in the travel time array to the (i, j, k)'th
!>           grid indices.
!>    @param[in] igrd   Fortran indxed grid point in travel time mesh to convert.
!>
!>    @param[out] i     iz'th grid point.  This is Fortran indexed.
!>    @param[out] j     ix'th grid point.  This is Fortran indexed.
!>    @param[out] k     iy'th grid point.  This is Fortran indexed.
!>    @param[out] ierr  0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE index2gridF(igrd, i, j, k, ierr) &
                      BIND(C, NAME='fteik_index2gridF')
      !$OMP DECLARE SIMD(index2gridF)
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : nx, ny, nz, nz, nzx
      INTEGER(C_INT), INTENT(IN), VALUE :: igrd
      INTEGER(C_INT), INTENT(OUT) :: i, j, k, ierr
      INTEGER(C_INT) igrd1
      igrd1 = igrd - 1 ! F to C
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
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility for routine for getting the sweep direction signs and offsets
!>           for the stencils.
!>
!>    @param[in] sweep     Sweep number.  This must be on the interval [1,8].
!>
!>    @param[out] sgntz    Stencil direction in z for travel times for this sweep.
!>    @param[out] sgntx    Stencil direction in x for travel times for this sweep.
!>    @param[out] sgnty    Stencil direction in y for travel times for this sweep.
!>    @param[out] sgnvz    Stencil offset in z for velocity model for this sweep.
!>    @param[out] sgnvx    Stencil offset in x for velocity model for this sweep.
!>    @param[out] sgnvy    Stencil offset in y for velocity model for this sweep.
!>    @param[out] sgnrz    Double precision version of sgntz for local solver.
!>    @param[out] sgnrx    Double precision version of sgntx for local solver.
!>    @param[out] sgnry    Double precision version of sgnty for local solver. 
!>
!>    @author Ben Baker
!>
!>    @copyright MIT 
!>
      SUBROUTINE fteik_getSweepSigns64fF(sweep,               &
                                         sgntz, sgntx, sgnty, &
                                         sgnvz, sgnvx, sgnvy, &
                                         sgnrz, sgnrx, sgnry) &
      BIND(C, NAME='fteik_getSweepSigns64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : sgntzv, sgntxv, sgntyv, sgnvzv, sgnvxv, sgnvyv
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: sweep
      INTEGER(C_INT), INTENT(OUT) :: sgntz, sgntx, sgnty, sgnvz, sgnvx, sgnvy
      REAL(C_DOUBLE), INTENT(OUT) :: sgnrz, sgnrx, sgnry
      sgntz = sgntzv(sweep)
      sgntx = sgntxv(sweep)
      sgnty = sgntyv(sweep)

      sgnvz = sgnvzv(sweep)
      sgnvx = sgnvxv(sweep)
      sgnvy = sgnvyv(sweep)

      sgnrz = DBLE(sgntz)
      sgnrx = DBLE(sgntx)
      sgnry = DBLE(sgnty)
      RETURN
      END SUBROUTINE 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      PURE INTEGER(C_INT) FUNCTION velGrid2indexF(i, j, k, nzm1, nzm1_nxm1)   &
                          BIND(C, NAME='fteik_velGrid2indexF')
      !$OMP DECLARE SIMD(velGrid2indexF) UNIFORM(nzm1, nzm1_nxm1)
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nzm1, nzm1_nxm1
      velGrid2indexF = (k - 1)*nzm1_nxm1 + (j - 1)*nzm1 + i
      RETURN
      END FUNCTION

!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility routine for computing some grid spacing constants use by the 
!>           finite difference stencils in the local solver.
!>
!>    @author Mark Noble, Alexandrine Gesret, and Ben Baker
!>
!>    @version 2
!>
!>    @date July 2017
!>
!>    Note, this now updates private variables in the fteik module and computes a 
!>    few additional terms that are constantly used in the local solver.
!>
      SUBROUTINE fteik_meshConstants64fF() &
                 BIND(C, NAME='fteik_meshConstants64fF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : dx, dy, dz,                                      &
                                 dx2, dy2, dz2, dxi, dyi, dzi,                    &
                                 dx2i, dy2i, dz2i, dsum, dsumi,                   &
                                 dz2dx2, dz2dy2, dx2dy2,                          &
                                 dz2i_dx2i, dz2i_dy2i, dx2i_dy2i,                 &
                                 dz2i_p_dx2i, dz2i_p_dy2i, dx2i_p_dy2i,           &
                                 dz2i_p_dy2i_inv, dz2i_p_dx2i_inv, dx2i_p_dy2i_inv
      USE FTEIK_UTILS64F, ONLY : one
      IMPLICIT NONE
      dz2 = dz*dz
      dx2 = dx*dx
      dy2 = dy*dy
      dzi = one/dz
      dxi = one/dx
      dyi = one/dy
      dz2i = one/dz2
      dx2i = one/dx2
      dy2i = one/dy2
      dsum = dz2i + dx2i + dy2i
      dsumi = one/dsum
      dz2dx2 = one/(dz2*dx2)
      dz2dy2 = one/(dz2*dy2)
      dx2dy2 = one/(dx2*dy2)
      dz2i_dx2i = dz2i*dx2i
      dz2i_dy2i = dz2i*dy2i
      dx2i_dy2i = dx2i*dy2i
      dz2i_p_dy2i_inv = one/(dz2i + dy2i)
      dz2i_p_dx2i_inv = one/(dz2i + dx2i)
      dx2i_p_dy2i_inv = one/(dx2i + dy2i)
      dz2i_p_dx2i = dz2i + dx2i
      dz2i_p_dy2i = dz2i + dy2i
      dx2i_p_dy2i = dx2i + dy2i
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Compute analytic travel time at the point (i,j,,k) in a homogeneous model.
!>
!>    @param[in] i      iz'th grid point (Fortran numbered).
!>    @param[in] j      ix'th grid point (Fortran numbered).
!>    @param[in] k      iy'th grid point (Fortran numbered).
!>    @param[in] dz     Grid spacing (meters) in z.
!>    @param[in] dx     Grid spacing (meters) in x.
!>    @param[in] dy     Grid spacing (meters) in y.
!>    @param[in] zsa    Source offset (meters) in z.
!>    @param[in] xsa    Source offset (meters) in x.
!>    @param[in] ysa    Source offset (meters) in y.
!>    @param[in] szero  Slowness at source (s/m).  
!>
!>    @result The travel time from the source at (zsa,xsa,ysa) to the (i,j,k)'th grid
!>            point given a constant slowness around the source.
!>
!>    @author Mark Noble, Alexandrine, Gesret, and Ben Baker.
!>
!>    @version 2
!>
!>    @date July 2017
!>
!>    @copyright CeCILL-3
!>
!>
      PURE REAL(C_DOUBLE) FUNCTION fteik_tAna64fF(i, j, k, dz, dx, dy,   &
                                                  zsa, xsa, ysa, szero)  &
                          BIND(C, NAME='fteik_tAna64fF')
      !$OMP DECLARE SIMD(fteik_tAna64fF) UNIFORM(dz, dx, dy, zsa, xsa, ysa, szero)
      USE ISO_C_BINDING
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN), VALUE :: dz, dx, dy, szero, zsa, xsa, ysa
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k
      REAL(C_DOUBLE) d0, diffz, diffz2, diffx, diffx2, diffy, diffy2
      diffz = (DBLE(i) - zsa)*dz
      diffx = (DBLE(j) - xsa)*dx
      diffy = (DBLE(k) - ysa)*dy
      diffz2 = diffz*diffz
      diffx2 = diffx*diffx
      diffy2 = diffy*diffy
      d0 = diffz2 + diffx2 + diffy2
      fteik_tAna64fF = szero*SQRT(d0) !sqrt(diffz2 + diffx2 + diffy2);
      RETURN
      END FUNCTION
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Compute derivative of analytic travel time and derivative of times
!>           at point (i,j,k) in a homogeneous model.
!>
!>    @param[in] i      iz'th grid point (Fortran numbered).
!>    @param[in] j      ix'th grid point (Fortran numbered).
!>    @param[in] k      iy'th grid point (Fortran numbered).
!>    @param[in] dz     Grid spacing (meters) in z.
!>    @param[in] dx     Grid spacing (meters) in x.
!>    @param[in] dy     Grid spacing (meters) in y.
!>    @param[in] zsa    Source offset (meters) in z.
!>    @param[in] xsa    Source offset (meters) in x.
!>    @param[in] ysa    Source offset (meters) in y.
!>    @param[in] szero  Slowness at source (s/m).  
!>
!>    @param[out] t_anad   Derivative of analytic travel time at point (i,j,k).
!>    @param[out] tzc      Derivative of analytic travel time in z.
!>    @param[out] txc      Derivative of analytic travel time in x.
!>    @param[out] tyc      Derivative of analytic travel time in y.
!>
!>    @author Mark Noble, Alexandrine, Gesret, and Ben Baker.
!>
!>    @version 2
!>
!>    @date July 2017
!>
!>    @copyright CeCILL-3
!>
!>
      PURE SUBROUTINE fteik_tAnaD64fF(t_anad, tzc, txc, tyc, i, j, k,   &
                                      dz, dx, dy, zsa, xsa, ysa, szero) &
                              BIND(C, NAME='fteik_tAnaD64fF')
      !$OMP DECLARE SIMD(fteik_tAnaD64fF) UNIFORM(dz, dx, dy, zsa, xsa, ysa, szero)
      USE ISO_C_BINDING
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN), VALUE :: dz, dx, dy, zsa, xsa, ysa, szero
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k
      REAL(C_DOUBLE), INTENT(OUT) :: t_anad, tzc, txc, tyc
      REAL(C_DOUBLE) d0, d0i_szero, diffz, diffz2, diffx, diffx2, diffy, diffy2, sqrtd0
      REAL(C_DOUBLE), PARAMETER :: zero = 0.d0
      diffz = (DBLE(i) - zsa)*dz
      diffx = (DBLE(j) - xsa)*dx
      diffy = (DBLE(k) - ysa)*dy
      diffz2 = diffz*diffz
      diffx2 = diffx*diffx
      diffy2 = diffy*diffy
      d0 = diffz2 + diffx2 + diffy2

      d0i_szero = zero
      t_anad = zero
      IF (d0 > zero) THEN
         sqrtd0 = SQRT(d0)
         t_anad = szero*sqrtd0
         d0i_szero = szero/sqrtd0
      ENDIF
      tzc = d0i_szero*diffz
      txc = d0i_szero*diffx
      tyc = d0i_szero*diffy
      RETURN
      END SUBROUTINE 
