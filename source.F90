!>    @brief Deallocates memory on the source module and resets all variables.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!> 
      SUBROUTINE fteik_source_finalizeF()             &
                 BIND(C, NAME='fteik_source_finalizeF')
      USE FTEIK_SOURCE64F, ONLY : zsrc, xsrc, ysrc, zsiv, xsiv, ysiv, zsav, xsav, ysav, &
                                  nsrc, lhaveSource
      USE ISO_C_BINDING
      IMPLICIT NONE
      lhaveSource = .FALSE.
      nsrc = 0 
      IF (ALLOCATED(zsrc)) DEALLOCATE(zsrc)
      IF (ALLOCATED(xsrc)) DEALLOCATE(xsrc)
      IF (ALLOCATED(ysrc)) DEALLOCATE(ysrc)
      IF (ALLOCATED(zsiv)) DEALLOCATE(zsiv)
      IF (ALLOCATED(xsiv)) DEALLOCATE(xsiv)
      IF (ALLOCATED(ysiv)) DEALLOCATE(ysiv)
      IF (ALLOCATED(zsav)) DEALLOCATE(zsav)
      IF (ALLOCATED(xsav)) DEALLOCATE(xsav)
      IF (ALLOCATED(ysav)) DEALLOCATE(ysav)
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Initializes the source positions on the module.  The model must be set
!>           prior to calling this.
!>
!>    @param[in] nsrcIn     Number of sources.
!>    @param[in] zsrcIn     z locations of source in model (meters).  This is a vector
!>                          of dimension [nsrcIn].
!>    @param[in] xsrcIn     x locations of source in model (meters).  This is a vector
!>                          of dimension [nsrcIn].
!>    @param[in] ysrcIn     y locations of source in model (meters).  This is a vector
!>                          of dimension [nsrcIn].
!>
!>    @param[out] ierr      0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_source_initialize64fF(nsrcIn, zsrcIn, xsrcIn, ysrcIn, ierr) &
                 BIND(C, NAME='fteik_source_initialize64fF')
      USE FTEIK_SOURCE64F, ONLY : zsrc, xsrc, ysrc, &
                                  zsav, xsav, ysav, &
                                  zsiv, xsiv, ysiv, &
                                  nsrc, lhaveSource
      USE FTEIK_MODEL64F, ONLY : dz, dx, dy, z0, x0, y0, nz, nx, ny, lhaveModel
      USE FTEIK_CONSTANTS64F, ONLY : DBL_EPSILON, perturbSource, one, zero 
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsrcIn
      REAL(C_DOUBLE), INTENT(IN) :: zsrcIn(nsrcIn), xsrcIn(nsrcIn), ysrcIn(nsrcIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) xsa, ysa, zsa, xsrcLoc, ysrcLoc, zsrcLoc, xmax, ymax, zmax
      INTEGER(C_INT) isrc, xsi, ysi, zsi
      !----------------------------------------------------------------------------------!
      ierr = 0
      lhaveSource = .FALSE.
      nsrc = 0
      IF (ALLOCATED(zsrc)) DEALLOCATE(zsrc)
      IF (ALLOCATED(xsrc)) DEALLOCATE(xsrc)
      IF (ALLOCATED(ysrc)) DEALLOCATE(ysrc)
      IF (ALLOCATED(zsiv)) DEALLOCATE(zsiv)
      IF (ALLOCATED(xsiv)) DEALLOCATE(xsiv)
      IF (ALLOCATED(ysiv)) DEALLOCATE(ysiv)
      IF (ALLOCATED(zsav)) DEALLOCATE(zsav)
      IF (ALLOCATED(xsav)) DEALLOCATE(xsav)
      IF (ALLOCATED(ysav)) DEALLOCATE(ysav)
      IF (nsrcIn < 1) THEN
         WRITE(*,*) 'fteik_source_initialize64fF: No sources'
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.lhaveModel) THEN
         WRITE(*,*) 'fteik_source_initialize64fF: Model not yet set'
         ierr = 1
         RETURN
      ENDIF
      nsrc = nsrcIn
      ALLOCATE(zsrc(nsrc))
      ALLOCATE(xsrc(nsrc))
      ALLOCATE(ysrc(nsrc))
      ALLOCATE(zsiv(nsrc))
      ALLOCATE(xsiv(nsrc))
      ALLOCATE(ysiv(nsrc)) 
      ALLOCATE(zsav(nsrc))
      ALLOCATE(xsav(nsrc))
      ALLOCATE(ysav(nsrc))
      zsiv(:) = 0
      xsiv(:) = 0
      ysiv(:) = 0
      zsav(:) =-one
      xsav(:) =-one
      ysav(:) =-one
      ! Compute max of grid
      zmax = DBLE(nz - 1)*dz
      xmax = DBLE(nx - 1)*dx
      ymax = DBLE(ny - 1)*dy
      DO 1 isrc=1,nsrc
         ! Create a relative position in the grid
         zsrcLoc = zsrcIn(isrc) - z0
         xsrcLoc = xsrcIn(isrc) - x0
         ysrcLoc = ysrcIn(isrc) - y0
         ! Check the bounds
         IF (zsrcLoc < zero .OR. zsrcLoc > zmax .OR. &
             xsrcLoc < zero .OR. xsrcLoc > xmax .OR. &
             ysrcLoc < zero .OR. ysrcLoc > ymax) THEN
            IF (zsrcLoc < zero .OR. zsrcLoc > zmax) THEN
               WRITE(*,*) 'fteik_source_initialize64fF: zsrc is out of bounds', &
                          isrc, zsrcLoc, zero + z0, zmax + z0
            ENDIF
            IF (xsrcLoc < zero .OR. xsrcLoc > xmax) THEN
               WRITE(*,*) 'fteik_source_initialize64fF: xsrc is out of bounds', &
                          isrc, xsrcLoc, zero + x0, xmax + x0
            ENDIF
            IF (ysrcLoc < zero .OR. ysrcLoc > ymax) THEN
               WRITE(*,*) 'fteik_source_initialize64fF: ysrc is out of bounds', &
                          isrc, ysrcLoc, zero + y0, ymax + y0
            ENDIF
            ierr = 1
            RETURN
         ENDIF
         ! Convert source position to grid position (Fortran indexing)
         zsa = (zsrcLoc/dz) + one
         xsa = (xsrcLoc/dx) + one
         ysa = (ysrcLoc/dy) + one 
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
            WRITE(*,*)'fteik_source_initialize64fF: Point (',zsi,xsi,ysi,') out of bounds'
            ierr = 1 
            RETURN
         ENDIF
         IF (zsi > nz - 1 .OR. xsi > nx - 1 .OR. ysi > ny - 1) THEN
            WRITE(*,*)'fteik_source_initialize64fF: Warning solver may segfault', isrc
         ENDIF
         ! Set on the arrays
         zsrc(isrc) = zsa + z0 !zsrcIn(isrc)
         xsrc(isrc) = xsa + x0 !xsrcIn(isrc)
         ysrc(isrc) = ysa + z0 !ysrcIn(isrc)
         zsiv(isrc) = zsi
         xsiv(isrc) = xsi
         ysiv(isrc) = ysi
         zsav(isrc) = zsa
         xsav(isrc) = xsa
         ysav(isrc) = ysa
         WRITE(*,*) 'fteik_source_initilialize64fF: Source translation: (dz,dx,dy)=', &
                    ABS(zsrc(isrc) - zsrcIn(isrc)), ABS(xsrc(isrc) - xsrcIn(isrc)),   &
                    ABS(ysrc(isrc) - ysrcIn(isrc))
    1 CONTINUE
      lhaveSource = .TRUE.
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE fteik_source_getSolverInfo64fF(isrc,          &
                                                zsi, xsi, ysi, &
                                                zsa, xsa, ysa, &
                                                szero, szero2, &
                                                ierr)          &
                 BIND(C, NAME='fteik_source_getSolverInfo64fF')
      USE FTEIK_SOURCE64F, ONLY : fteik_source_getSourceIndices32iF, &
                                  fteik_source_getSourceIndices64fF, &
                                  fteik_source_getSzero64fF
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      REAL(C_DOUBLE), INTENT(OUT) :: zsa, xsa, ysa, szero, szero2
      INTEGER(C_INT), INTENT(OUT) :: zsi, xsi, ysi, ierr
      ierr = 0
      zsi = 0
      xsi = 0
      ysi = 0
      zsa = zero
      xsa = zero
      ysa = zero
      szero = zero
      szero2 = zero
      CALL fteik_source_getSourceIndices32iF(isrc, zsi, xsi, ysi, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_source_getSolverInfo64fF: Error getting source indices'
         RETURN
      ENDIF
      CALL fteik_source_getSourceIndices64fF(isrc, zsa, xsa, ysa, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_source_getSolverInfo64fF: Error getting source double indices'
         RETURN
      ENDIF 
      CALL fteik_source_getSzero64fF(isrc, szero, szero2, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_source_getSolverInfo64fF: Error getting szero'
         RETURN
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gets the source indices.
!>
!>    @param[in] isrc    Source number.
!>
!>    @param[out] zsa    z source grid point (double precision version).
!>    @param[out] xsa    x source grid point (double precision version).
!>    @param[out] ysa    y source grid point (double precision version).
!>    @param[out] ierr   0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_source_getSourceIndices64fF(isrc, zsa, xsa, ysa, ierr) &
                 BIND(C, NAME='fteik_source_getSourceIndices64fF')
      USE FTEIK_SOURCE64F, ONLY : zsav, xsav, ysav, nsrc, lhaveSource 
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      REAL(C_DOUBLE), INTENT(OUT) :: zsa, xsa, ysa
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0 
      zsa = zero 
      xsa = zero
      ysa = zero
      IF (.NOT.lhaveSource) THEN
         WRITE(*,*) 'fteik_source_getSourceIndices64fF: Source not initialized'
         ierr = 1 
         RETURN
      ENDIF
      IF (isrc < 1 .OR. isrc > nsrc) THEN
         WRITE(*,*) 'fteik_source_getSourceIndices64fF: Invalid source number', isrc, nsrc
         ierr = 1 
         RETURN
      ENDIF
      zsa = zsav(isrc)
      xsa = xsav(isrc)
      ysa = ysav(isrc)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gets the source indices.
!>
!>    @param[in] isrc    Source number.
!>
!>    @param[out] zsi    z source grid point.
!>    @param[out] xsi    x source grid point.
!>    @param[out] ysi    y source grid point.
!>    @param[out] ierr   0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_source_getSourceIndices32iF(isrc, zsi, xsi, ysi, ierr) &
                 BIND(C, NAME='fteik_source_getSourceIndices32iF')
      USE FTEIK_SOURCE64F, ONLY : zsiv, xsiv, ysiv, nsrc, lhaveSource 
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      INTEGER(C_INT), INTENT(OUT) :: zsi, xsi, ysi, ierr
      ierr = 0
      zsi = 0
      xsi = 0
      ysi = 0
      IF (.NOT.lhaveSource) THEN
         WRITE(*,*) 'fteik_source_getSourceIndices32iF: Source not initialized'
         ierr = 1
         RETURN
      ENDIF
      IF (isrc < 1 .OR. isrc > nsrc) THEN
         WRITE(*,*) 'fteik_source_getSourceIndices32iF: Invalid source number', isrc, nsrc
         ierr = 1
         RETURN
      ENDIF
      zsi = zsiv(isrc)
      xsi = xsiv(isrc)
      ysi = ysiv(isrc)
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gets the slowness at the source location.
!>
!>    @param[in] isrc        Source number.
!>
!>    @param[out] szero      Slowness, \f$ s_0 \f$, at source grid point.
!>    @param[out] szero2     \f$ s_{0}^2 \f$.
!>    @param[out] ierr       0 indicates success.
!>
!>    @result 0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_source_getSzero64fF(isrc, szero, szero2, ierr)   &
                 BIND(C, NAME='fteik_source_getSzero64fF')
      USE FTEIK_MODEL64F, ONLY : slow, nzm1, nzm1_nxm1, lhaveModel
      USE FTEIK_SOURCE64F, ONLY : zsiv, xsiv, ysiv, nsrc, lhaveSource
      USE FTEIK_CONSTANTS64F, ONLY : zero 
      USE FTEIK_UTILS64F, ONLY : velGrid2indexF
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      REAL(C_DOUBLE), INTENT(OUT) :: szero, szero2
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) indx
      szero = zero
      szero2 = zero
      ierr = 1
      IF (.NOT.lhaveSource) THEN
         WRITE(*,*) 'fteik_source_getSzero64fF: Source not yet initialized'
         RETURN
      ENDIF
      IF (.NOT.lhaveModel) THEN
         WRITE(*,*) 'fteik_source_getSzero64fF: Slowness model not yet set'
         RETURN
      ENDIF
      IF (isrc < 1 .OR. isrc > nsrc) THEN
         WRITE(*,*) 'fteik_source_getSzero64fF: Invalid source number', isrc, 1, nsrc
         RETURN
      ENDIF
      ierr = 0
      indx = velGrid2indexF(zsiv(isrc), xsiv(isrc), ysiv(isrc), nzm1, nzm1_nxm1)
      szero = slow(indx)
      szero2 = szero*szero
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
      SUBROUTINE fteik_source_setLocationF(zs, xs, ys, ierr) &
                 BIND(C, NAME='fteik_source_setLocationF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : dx, dy, dz, nx, ny, nz, &
                                 xsi, ysi, zsi, xsa, ysa, zsa, &
                                 x0, y0, z0, lhaveSource, &
                                 ijkv1, ijkv2, ijkv3, ijkv4, &
                                 ijkv5, ijkv6, ijkv7, ijkv8, &
                                 lupdInit1, lupdInit2, lupdInit3, lupdInit4, &
                                 lupdInit5, lupdInit6, lupdInit7, lupdInit8, &
                                 levelPtr, nLevels
      USE FTEIK_UTILS64F, ONLY : fteik_setUpdateNodesF
      USE FTEIK_UTILS64F, ONLY : DBL_EPSILON, one, zero
      IMPLICIT NONE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: xs, ys, zs
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
            WRITE(*,*) 'fteik_source_setLocationF: zsrc is out of bounds', &
                       zsrc, zero, zmax
         ENDIF
         IF (xsrc < zero .OR. xsrc > xmax) THEN
            WRITE(*,*) 'fteik_source_setLocationF: xsrc is out of bounds', &
                       xsrc, zero, xmax
         ENDIF
         IF (ysrc < zero .OR. ysrc > ymax) THEN
            WRITE(*,*) 'fteik_source_setLocationF: ysrc is out of bounds', &
                       ysrc, zero, ymax
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
         WRITE(*,*) 'fteik_source_setLocationF: Point (',zsi,xsi,ysi,') ouf of bounds'
         ierr = 1
         RETURN
      ENDIF
      IF (zsi > nz - 1 .OR. xsi > nx - 1 .OR. ysi > ny - 1) THEN
         WRITE(*,*) 'fteik_source_setLocationF: Warning solver may segfault'
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
         WRITE(*,*) 'fteik_source_setLocationF: Error setting update nodes'
         ierr = 1
         lhaveSource = .FALSE.
         RETURN
      ENDIF
      RETURN
      END SUBROUTINE

