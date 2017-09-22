
      SUBROUTINE fteik_source_setSources(nsrcIn, zsrc, xsrc, ysrc)

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
      SUBROUTINE fteik_source_setSzeroF(ierr) BIND(C, NAME='fteik_source_setSzeroF')
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
         WRITE(*,*) 'fteik_source_setSzeroF: Source not yet initialized'
         RETURN
      ENDIF
      IF (.NOT.lhaveSlownessModel) THEN
         WRITE(*,*) 'fteik_source_setSzeroF: Slowness model not yet set'
         RETURN
      ENDIF
      ierr = 0
      indx = velGrid2indexF(zsi, xsi, ysi, nzm1, nzm1_nxm1)
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

