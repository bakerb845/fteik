MODULE FTEIK_SOURCE64F
  USE ISO_C_BINDING
  IMPLICIT NONE 
  !> Double precision versions of the source grid points
  REAL(C_DOUBLE), PROTECTED, ALLOCATABLE, SAVE :: zsav(:), xsav(:), ysav(:)
  !> Source locations
  REAL(C_DOUBLE), PROTECTED, ALLOCATABLE, SAVE :: zsrc(:), xsrc(:), ysrc(:)
  !> Source grid points 
  INTEGER(C_INT), PROTECTED, ALLOCATABLE, SAVE :: zsiv(:), xsiv(:), ysiv(:)
  !DIR$ ATTRIBUTES ALIGN: 64 :: zsav, xsav, ysav
  !DIR$ ATTRIBUTES ALIGN: 64 :: zsrc, xsrc, ysrc
  !DIR$ ATTRIBUTES ALIGN: 64 :: zsiv, xsiv, ysiv
  !> Number of sources
  INTEGER(C_INT), PROTECTED, SAVE :: nsrc = 0
  !> Flag indicating the sources have been set -> this should be removed.
  LOGICAL(C_BOOL), PROTECTED, SAVE :: lhaveSource = .FALSE.
  CONTAINS
!----------------------------------------------------------------------------------------!
!                                   Begin the Code                                       !
!----------------------------------------------------------------------------------------!
!>
!>    @brief Deallocates memory on the source module and resets all variables.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!> 
      SUBROUTINE fteik_source_finalizeF()             &
                 BIND(C, NAME='fteik_source_finalizeF')
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
      USE FTEIK_MODEL64F, ONLY : dz, dx, dy, z0, x0, y0, nz, nx, ny
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
      CALL fteik_source_finalizeF()
      IF (nsrcIn < 1) THEN
         WRITE(*,*) 'fteik_source_initialize64fF: No sources'
         ierr = 1
         RETURN
      ENDIF
      IF (nz < 1 .OR. nx < 1 .OR. ny < 1 .OR.             &
          dz <= zero .OR. dx <= zero .OR. dy <= zero) THEN
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
                          isrc, zsrcLoc + z0, zero + z0, zmax + z0
            ENDIF
            IF (xsrcLoc < zero .OR. xsrcLoc > xmax) THEN
               WRITE(*,*) 'fteik_source_initialize64fF: xsrc is out of bounds', &
                          isrc, xsrcLoc + x0, zero + x0, xmax + x0
            ENDIF
            IF (ysrcLoc < zero .OR. ysrcLoc > ymax) THEN
               WRITE(*,*) 'fteik_source_initialize64fF: ysrc is out of bounds', &
                          isrc, ysrcLoc + y0, zero + y0, ymax + y0
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
         zsrc(isrc) = DBLE(zsi - 1)*dz + z0 !zsrcIn(isrc)
         xsrc(isrc) = DBLE(xsi - 1)*dx + x0 !xsrcIn(isrc)
         ysrc(isrc) = DBLE(ysi - 1)*dy + y0 !ysrcIn(isrc)
         zsiv(isrc) = zsi
         xsiv(isrc) = xsi
         ysiv(isrc) = ysi
         zsav(isrc) = zsa
         xsav(isrc) = xsa
         ysav(isrc) = ysa
print *, z0, x0, y0, zsi, xsi, ysi
print *, dz, dx, dy, xsa
         WRITE(*,900) zsrcIn(isrc), xsrcIn(isrc), ysrcIn(isrc)
         WRITE(*,901) zsrc(isrc), xsrc(isrc), ysrc(isrc)
         WRITE(*,902) ABS(zsrc(isrc) - zsrcIn(isrc)), ABS(xsrc(isrc) - xsrcIn(isrc)),   &
                      ABS(ysrc(isrc) - ysrcIn(isrc))
         WRITE(*,*)
    1 CONTINUE
      lhaveSource = .TRUE.
  900 FORMAT(' fteik_source_initialize64fF: Original source coordinates (z,x,y)=', &
             3F12.2, ' (m)') 
  901 FORMAT(' fteik_source_initialize64fF: Grid source coordinates (z,x,y)=', &
             3F12.2, ' (m)')
  902 FORMAT(' fteik_source_initialize64fF: Source translation: (dz,dx,dy)=', &
             3F12.2, ' (m)')
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility function to extract the number of sources.
!>
!>    @param[out] nsrcOut   Number of sources initialized on module.
!>    @param[out] ierr      0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_source_getNumberOfSources(nsrcOut, ierr) &
      BIND(C, NAME='fteik_source_getNumberOfSources')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: nsrcOut, ierr
      ierr = 0
      nsrcOut = 0
      IF (.NOT.lhaveSource) THEN
         WRITE(*,*) 'fteik_source_getNumberOfSources: Never initialized'
         ierr = 1
         RETURN
      ENDIF 
      nsrcOut = nsrc
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility routine which returns the source properties for the solver.
!>
!>    @param[in] isrc     Number number.  This is in range [1,nsrc].
!>
!>    @param[out] zsi     z source index in grid.
!>    @param[out] xsi     x source index in grid.
!>    @param[out] ysi     y source index in grid.
!>    @param[out] zsa     z source offset from origin in grid (meters).
!>    @param[out] xsa     x source offset from origin in grid (meters).
!>    @param[out] ysa     y source offset from origin in grid (meters).
!>    @param[out] szero   Slowness (s/m) at the source.
!>    @param[out] szero2  Slowness squared (s^2/m^2) at the source.
!>    @param[out] ierr    0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_source_getSolverInfo64fF(isrc,          &
                                                zsi, xsi, ysi, &
                                                zsa, xsa, ysa, &
                                                szero, szero2, &
                                                ierr)          &
                 BIND(C, NAME='fteik_source_getSolverInfo64fF')
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
      USE FTEIK_MODEL64F, ONLY : fteik_model_velGrid2indexF
      USE FTEIK_MODEL64F, ONLY : slow, nzm1, nzm1_nxm1, lhaveModel
      USE FTEIK_CONSTANTS64F, ONLY : zero 
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
      indx = fteik_model_velGrid2indexF(zsiv(isrc), xsiv(isrc), ysiv(isrc), &
                                        nzm1, nzm1_nxm1)
      szero = slow(indx)
      szero2 = szero*szero
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------------------------!
!                                      End the Code                                      !
!----------------------------------------------------------------------------------------!
END MODULE
