MODULE FTEIK_MODEL64F
  USE ISO_C_BINDING
  USE FTEIK_CONSTANTS64F, ONLY : zero
  IMPLICIT NONE
  REAL(C_DOUBLE), ALLOCATABLE, SAVE :: slow(:)  !< Slowness model (seconds/meter).
                                                !< This has dimension [ncell].
  !DIR$ ATTRIBUTES ALIGN: 64 :: slow
  REAL(C_DOUBLE), PROTECTED, SAVE :: dx = zero  !< x origin (meters).
  REAL(C_DOUBLE), PROTECTED, SAVE :: dy = zero  !< y origin (meters).
  REAL(C_DOUBLE), PROTECTED, SAVE :: dz = zero  !< z origin (meters).
  REAL(C_DOUBLE), PROTECTED, SAVE :: x0 = zero  !< x origin (meters).
  REAL(C_DOUBLE), PROTECTED, SAVE :: y0 = zero  !< y origin (meters).
  REAL(C_DOUBLE), PROTECTED, SAVE :: z0 = zero  !< z origin (meters).
  INTEGER(C_INT), PROTECTED, SAVE :: ncell =-1  !< Number of cells in slowness model.
  INTEGER(C_INT), PROTECTED, SAVE :: ngrd =-1   !< Number of grid points (= nz*nx*ny).
  INTEGER(C_INT), PROTECTED, SAVE :: nx =-1     !< Number of x grid points in 
                                                !< travel time field.
  INTEGER(C_INT), PROTECTED, SAVE :: ny =-1     !< Number of y grid points in 
                                                !< travel time field.
  INTEGER(C_INT), PROTECTED, SAVE :: nz =-1     !< Number of z grid points in 
                                                !< travel time field.
  INTEGER(C_INT), PROTECTED, SAVE:: nzx =-1     !< =nz*nx
  INTEGER(C_INT), PROTECTED, SAVE :: nzm1_nxm1 =-1 !< =(nz - 1)*(nx - 1)
  INTEGER(C_INT), PROTECTED, SAVE :: nzm1 =-1      !< =(nz - 1)
  LOGICAL(C_BOOL), PROTECTED, SAVE :: lhaveModel = .FALSE. !< If true slowness model 
                                                           !> was set.

  CONTAINS
!----------------------------------------------------------------------------------------!
!                                 Begin the Code                                         !
!----------------------------------------------------------------------------------------!
!>    @brief Convenience function to initialize the model geometry.
!>
!>    @param[in] nz       Number of z grid points in the travel-time field.
!>                        This must be at least 3.
!>    @param[in] nx       Number of x grid points in the travel-time field.
!>                        This must be at least 3.
!>    @param[in] ny       Number of y grid points in the travel-time field.
!>                        This must be at least 3.
!>    @param[in] dz       Mesh spacing in z (meters).  This must be positive.
!>    @param[in] dx       Mesh spacing in x (meters).  This must be positive.
!>    @param[in] dy       Mesh spacing in y (meters).  This must be positive.
!>    @param[in] z0       z origin (meters).
!>    @param[in] x0       x origin (meters).
!>    @param[in] y0       y origin (meters).
!>
!>    @param[out] ierr    0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_model_intializeGeometryF(nz, nx, ny,    &
                                                dz, dx, dy,    &
                                                z0, x0, y0,    &
                                                ierr)          &
                 BIND(C, NAME='fteik_model_initializeGeometryF')
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: nz, nx, ny
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dz, dx, dy, z0, x0, y0
      INTEGER(C_INT), INTENT(OUT) :: ierr
      CALL fteik_model_setGridSizeF(nz, nx, ny, ierr) 
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_model_initializeGeometryF: Error setting gridsize'
         RETURN
      ENDIF 
      CALL fteik_model_setGridSpacingF(dz, dx, dy, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_model_initializeGeomtryF: Error setting grid spacing'
         RETURN
      ENDIF
      CALL fteik_model_setOriginF(z0, x0, y0)
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
      SUBROUTINE fteik_model_setGridSpacingF(dzIn, dxIn, dyIn, ierr) &
                 BIND(C, NAME='fteik_model_setGridSpacingF')
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dzIn, dxIn, dyIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      dz = zero
      dx = zero 
      dy = zero
      IF (dzIn <= zero .OR. dxIn <= zero .OR. dyIn <= zero) THEN
         IF (dzIn <= zero) WRITE(*,*) 'fteik_model_setGridSpacing: dz is too small', dzIn
         IF (dxIn <= zero) WRITE(*,*) 'fteik_model_setGridSpacing: dx is too small', dxIn
         IF (dyIn <= zero) WRITE(*,*) 'fteik_model_setGridSpacing: dy is too small', dyIn
         ierr = 1
         RETURN
      ENDIF
      dz = dzIn
      dx = dxIn
      dy = dyIn
      RETURN
      END
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
      SUBROUTINE fteik_model_setGridSizeF(nzIn, nxIn, nyIn, ierr) &
                 BIND(C, NAME='fteik_model_setGridSizeF')
      USE ISO_C_BINDING
      USE FTEIK_CONSTANTS64F, ONLY : zero
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nzIn, nxIn, nyIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      nz = 0
      nx = 0
      ny = 0
      ngrd = 0
      ncell = 0
      IF (nzIn < 3 .OR. nxIn < 3 .OR. nyIn < 3) THEN
         IF (nzIn < 3) WRITE(*,*) 'fteik_model_setGridSize: ERROR nz is too small', nzIn
         IF (nxIn < 3) WRITE(*,*) 'fteik_model_setGridSize: ERROR nx is too small', nxIn
         IF (nyIn < 3) WRITE(*,*) 'fteik_model_setGridSize: ERROR ny is too small', nyIn
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
      IF (ALLOCATED(slow)) DEALLOCATE(slow)
      ALLOCATE(slow(ncell))
      slow(:) = zero
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
      SUBROUTINE fteik_model_setOriginF(z0In, x0In, y0In) &
                 BIND(C, NAME='fteik_model_setOriginF')
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x0In, y0In, z0In 
      z0 = z0In
      x0 = x0In
      y0 = y0In
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
      SUBROUTINE fteik_model_getGridSizeF(nzOut, nxOut, nyOut,     &
                                          ngrdOut, ncellOut, ierr) &
                 BIND(C, NAME='fteik_model_getGridSizeF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: nzOut, nxOut, nyOut, ngrdOut, ncellOut, ierr
      ierr = 0
      nzOut = 0
      nxOut = 0
      nyOut = 0
      ngrdOut = 0
      ncellOut = 0
      IF (nx < 1 .OR. ny < 1 .OR. nz < 1) THEN
         WRITE(*,*) 'fteik_model_getGridSize: Grid not yet set'
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
      SUBROUTINE fteik_model_setVelocityModel64fF(nv, vel, ierr) &
                 BIND(C, NAME='fteik_model_setVelocityModel64fF')
      USE FTEIK_CONSTANTS64F, ONLY : one, zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nv
      REAL(C_DOUBLE), INTENT(IN) :: vel(nv)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) i
      ierr = 0
      lhaveModel = .FALSE.
      IF (nv /= ncell) THEN
         WRITE(*,*) 'fteik_model_setVelocityModel64fF: ERROR - ncell /= nv', ncell, nv
         ierr = 1 
         RETURN
      ENDIF
      IF (nv < 1) THEN
         WRITE(*,*) 'fteik_model_setVelocityModel64fF: ERROR - no cells in vel', nv
         ierr = 1 
         RETURN
      ENDIF
      IF (MINVAL(vel) <= zero) THEN
         WRITE(*,*) 'fteik_model_setVelocityModel64fF: All velocities must be positive'
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.ALLOCATED(slow)) ALLOCATE(slow(ncell))
      DO 1 i=1,ncell
         slow(i) = one/vel(i)
    1 CONTINUE  
      lhaveModel = .TRUE.
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Finalizes the model module by releasing memory and setting variables to
!>           undefined types.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_model_finalizeF() &
                 BIND(C, NAME='fteik_model_finalizeF')
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      dx = zero
      dy = zero
      dz = zero
      x0 = zero
      y0 = zero
      z0 = zero
      ncell =-1
      ngrd =-1
      nx =-1
      ny =-1
      nz =-1
      nzx =-1
      nzm1_nxm1 =-1
      nzm1 =-1
      lhaveModel = .FALSE.
      IF (ALLOCATED(slow)) DEALLOCATE(slow)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Convenience function to copy the velocity model from the 
!>
!>    @param[in] nwork     If nwork =-1 then this is a space query for nv. \n
!>                         Otherwise, this is the max space allocated to vel.
!>
!>    @param[out] nv       Number of cells in the velocity model.
!>    @param[out] vel      Velocity model (m/s).  This is a vector of length [nwork].
!>                         If nwork > nv then the velocity model is in vel(1:nv).
!>    @param[out] ierr     0 indicates usccess. 
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_model_getVelocityModel64fF(nwork, nv, vel, ierr) &
                 BIND(C, NAME='fteik_model_getVelocityModel64fF')
      USE FTEIK_CONSTANTS64F, ONLY : one, zero
      USE ISO_C_BINDING
      IMPLICIT NONE 
      INTEGER(C_INT), VALUE, INTENT(IN) :: nwork
      REAL(C_DOUBLE), INTENT(OUT) :: vel(nwork)
      INTEGER(C_INT), INTENT(OUT) :: nv, ierr
      INTEGER(C_INT) i
      ierr = 0
      nv = 0
      IF (nwork ==-1) THEN
         nv = ncell
         RETURN
      ENDIF
      IF (nwork < ncell) THEN
         WRITE(*,*) 'fteik_model_getVelocityModel64fF: nv < ncell'
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.lhaveModel) THEN
         WRITE(*,*) 'fteik_model_getVelocityModel64fF: Velocity model not yet set'
         ierr = 1
         vel(:) = zero
         RETURN
      ENDIF 
      IF (nv > ncell) vel(ncell+1:nwork) = zero
      DO 1 i=1,ncell
         vel(i) = one/slow(i) 
    1 CONTINUE
      RETURN
      END
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
      SUBROUTINE fteik_model_setVelocityModel32fF(nv, vel4, ierr) &
                 BIND(C, NAME='fteik_model_setVelocityModel32fF')
      USE FTEIK_CONSTANTS64F, ONLY : one, zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nv
      REAL(C_FLOAT), INTENT(IN) :: vel4(nv)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) i
      ierr = 0
      lhaveModel = .FALSE.
      IF (nv /= ncell) THEN
         WRITE(*,*) 'fteik_model_setVelocityModel32fF: ERROR - ncell /= nv', ncell, nv
         ierr = 1
         RETURN
      ENDIF
      IF (nv < 1) THEN
         WRITE(*,*) 'fteik_model_setVelocityMode32flF: ERROR - no cells in vel', nv
         ierr = 1
         RETURN
      ENDIF
      IF (MINVAL(vel4) <= 0.0) THEN
         WRITE(*,*) 'fteik_model_setVelocityModel32fF: All velocities must be positive'
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.ALLOCATED(slow)) ALLOCATE(slow(ncell))
      DO 1 i=1,ncell
         slow(i) = one/DBLE(vel4(i))
    1 CONTINUE
      lhaveModel = .TRUE.
      RETURN
      END SUBROUTINE
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
      PURE INTEGER(C_INT)                                       &
      FUNCTION fteik_model_grid2indexF(i, j, k, nz, nzx)        &
      BIND(C, NAME='fteik_model_grid2indexF')                   &
      RESULT(grid2indexF)
      !$OMP DECLARE SIMD(fteik_model_grid2indexF) UNIFORM(nz, nzx)
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
      PURE SUBROUTINE fteik_model_index2gridF(igrd, i, j, k, ierr) &
      BIND(C, NAME='fteik_model_index2gridF')
      !$OMP DECLARE SIMD(fteik_model_index2gridF) UNIFORM(ierr)
      USE ISO_C_BINDING
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
!>    @brief Converts the (i, j, k)'th model index to the velocity cell.
!>
!>    @param[in] i           iz'th grid point.  This is in range [1, nz-1].
!>    @param[in] j           ix'th grid point.  This is in range [1, nx-1].
!>    @param[in] k           iy'th grid point.  This is in range [1, ny-1].
!>    @param[in] nzm1        nz - 1.
!>    @param[in] nzm1_nxm1   (nz - 1)*(nx - 1)
!>
!>    @result The velocity cell corresponding to the (i, j, k)'th model index.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      PURE INTEGER(C_INT)                                             &
      FUNCTION fteik_model_velGrid2indexF(i, j, k, nzm1, nzm1_nxm1)   &
      BIND(C, NAME='fteik_model_velGrid2indexF')                      &
      RESULT(velGrid2IndexF)
      !$OMP DECLARE SIMD(fteik_model_velGrid2indexF) UNIFORM(nzm1, nzm1_nxm1)
      USE ISO_C_BINDING
      IMPLICIT NONE 
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nzm1, nzm1_nxm1
      velGrid2indexF = (k - 1)*nzm1_nxm1 + (j - 1)*nzm1 + i
      RETURN
      END FUNCTION
!----------------------------------------------------------------------------------------!
!                                    End the Code                                        !
!----------------------------------------------------------------------------------------!
END MODULE
