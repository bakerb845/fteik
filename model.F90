!> @defgroup model Model
!> @ingroup solver2d
!> @ingroup solver3d
!> @brief Defines the velocity model and physical domain.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE FTEIK_MODEL64F
  USE ISO_FORTRAN_ENV
  USE ISO_C_BINDING
  USE FTEIK_CONSTANTS64F, ONLY : zero
  USE FTEIK_CONSTANTS64F, ONLY : FTEIK_NATURAL_ORDERING, &
                                 FTEIK_ZXY_ORDERING,     &    
                                 FTEIK_XYZ_ORDERING,     &    
                                 FTEIK_ZYX_ORDERING,     &
                                 FTEIK_ZX_ORDERING,      &
                                 FTEIK_XZ_ORDERING
  IMPLICIT NONE
  REAL(C_DOUBLE), PROTECTED, ALLOCATABLE, SAVE :: slow(:)  !< Slowness model in
                                                           !< (seconds/meter).
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
  INTEGER(C_INT), PROTECTED, SAVE :: verbose = 0   !< Controsl verbosity.
  LOGICAL(C_BOOL), PROTECTED, SAVE :: lhaveModel = .FALSE. !< If true slowness model 
                                                           !< was set.
  LOGICAL(C_BOOL), PROTECTED, SAVE :: lis3dModel !< If true then the model is 3D.
                                                 !< Otherwise, it is 2D. 
  ! Label the routines
  PUBLIC :: fteik_model_initializeGeometry
  PUBLIC :: fteik_model_setVelocityModel64f
  PUBLIC :: fteik_model_setVelocityModel32f
  PUBLIC :: fteik_model_getGridSize
  PUBLIC :: fteik_model_getVelocityModel64f
  PUBLIC :: fteik_model_setNodalVelocityModel64f
  PUBLIC :: fteik_model_setCellVelocityModel64f
  PUBLIC :: fteik_model_isModel3D
  PUBLIC :: fteik_model_free
  PUBLIC :: fteik_model_setGridSpacing
  PUBLIC :: fteik_model_setOrigin3D
  PUBLIC :: fteik_model_setOrigin2D
  PUBLIC :: fteik_model_setGridSize

  PUBLIC :: fteik_model_grid2indexF
  PUBLIC :: fteik_model_index2gridF
  PUBLIC :: fteik_model_velGrid2indexF
  CONTAINS
!----------------------------------------------------------------------------------------!
!                                 Begin the Code                                         !
!----------------------------------------------------------------------------------------!
!>    @brief Convenience function to initialize the model geometry.
!>
!>    @param[in] lis3d    If true then this is a 3D model.  \n
!>                        Otherwise, it is a 2D model and the y variables will be
!>                        ignored. 
!>    @param[in] nz       Number of z grid points in the travel-time field.
!>                        This must be at least 3.
!>    @param[in] nx       Number of x grid points in the travel-time field.
!>                        This must be at least 3.
!>    @param[in] ny       Number of y grid points in the travel-time field.
!>                        This must be at least 3.  This is only accessed for 3D models.
!>    @param[in] dz       Mesh spacing in z (meters).  This must be positive.
!>    @param[in] dx       Mesh spacing in x (meters).  This must be positive.
!>    @param[in] dy       Mesh spacing in y (meters).  This must be positive and is only
!>                        accessed for 3D models.
!>    @param[in] z0       z origin (meters).
!>    @param[in] x0       x origin (meters).
!>    @param[in] y0       y origin (meters).
!>
!>    @param[out] ierr    0 indicates success.
!>
!>    @ingroup model
      SUBROUTINE fteik_model_initializeGeometry(lis3d,         &
                                                nz, nx, ny,    &
                                                z0, x0, y0,    &
                                                dz, dx, dy,    &
                                                ierr)          &
                 BIND(C, NAME='fteik_model_initializeGeometry')
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: nz, nx, ny
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dz, dx, dy, z0, x0, y0
      LOGICAL(C_BOOl), VALUE, INTENT(IN) :: lis3d
      INTEGER(C_INT), INTENT(OUT) :: ierr
      lis3dModel = lis3d
      CALL fteik_model_setGridSize(lis3d, nz, nx, ny, ierr) 
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         RETURN
      ENDIF 
      CALL fteik_model_setGridSpacing(lis3d, dz, dx, dy, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,901)
         RETURN
      ENDIF
      IF (lis3d) THEN
         CALL fteik_model_setOrigin3D(z0, x0, y0)
      ELSE
         CALL fteik_model_setOrigin2D(z0, x0)
      ENDIF
  900 FORMAT('fteik_model_initializeGeometry: Error setting gridsize')
  901 FORMAT('fteik_model_initializeGeomtry: Error setting grid spacing')
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the grid spacing of the solver.  This must be called to properly
!>           start the solver.
!>
!>    @param[in] lis3d     If true then this is a 3D model. \n
!>                         Otherwise, it is a 2D model and y will not be accessed.
!>    @param[in] dzIn      Mesh spacing in z (meters).  This must be positive.
!>    @param[in] dxIn      Mesh spacing in x (meters).  This must be positive.
!>    @param[in] dyIn      Mesh spacing in y (meters).  This must be positive if 
!>                         the model is 3D.
!>
!>    @param[out] ierr     0 indicates success.
!>
!>    @ingroup model
      SUBROUTINE fteik_model_setGridSpacing(lis3d, dzIn, dxIn, dyIn, ierr) &
                 BIND(C, NAME='fteik_model_setGridSpacing')
      USE FTEIK_CONSTANTS64F, ONLY : zero, one
      USE ISO_C_BINDING
      IMPLICIT NONE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dzIn, dxIn, dyIn
      LOGICAL(C_BOOL), VALUE, INTENT(IN) :: lis3d
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      dz = zero
      dx = zero 
      dy = zero
      IF (dzIn <= zero .OR. dxIn <= zero .OR. (lis3d .AND. dyIn <= zero)) THEN
         IF (dzIn <= zero) WRITE(ERROR_UNIT,900) dzIn
         IF (dxIn <= zero) WRITE(ERROR_UNIT,901) dxIn
         IF (lis3d .AND. dyIn <= zero) WRITE(ERROR_UNIT,902)  dyIn
         ierr = 1
         RETURN
      ENDIF
      dz = dzIn
      dx = dxIn
      IF (lis3d) THEN
         dy = dyIn
      ELSE
         dy = one
      ENDIF
  900 FORMAT('fteik_model_setGridSpacing: dz is too small', E12.5)
  901 FORMAT('fteik_model_setGridSpacing: dy is too small', E12.5)
  902 FORMAT('fteik_model_setGridSpacing: dx is too small', E12.5)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief This sets the number of grid points in each dimension of the travel-time
!>           field.  This function must be called to properly start the solver.
!>
!>    @param[in] lis3d    If true then this is a 3D model. \n
!>                        Otherwise, it is a 2D model and y will not be accessed.
!>    @param[in] nzIn     Number of z grid points in the travel-time field.
!>                        This must be at least 3.
!>    @param[in] nxIn     Number of x grid points in the travel-time field.
!>                        This must be at least 3.
!>    @param[in] nyIn     Number of y grid points in the travel-time field.
!>                        This must be at least 3 if lis3d is true.
!>
!>    @param[out] ierr    0 indicates success.
!>
!>    @ingroup model
      SUBROUTINE fteik_model_setGridSize(lis3d, nzIn, nxIn, nyIn, ierr) &
                 BIND(C, NAME='fteik_model_setGridSize')
      USE ISO_C_BINDING
      USE FTEIK_CONSTANTS64F, ONLY : zero
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nzIn, nxIn, nyIn
      LOGICAL(C_BOOL), VALUE, INTENT(IN) :: lis3d
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      nz = 0
      nx = 0
      ny = 0
      ngrd = 0
      ncell = 0
      IF (nzIn < 3 .OR. nxIn < 3 .OR. (lis3d .AND. nyIn < 3)) THEN
         IF (nzIn < 3) WRITE(ERROR_UNIT,900) nzIn
         IF (nxIn < 3) WRITE(ERROR_UNIT,901) nxIn
         IF (lis3d .AND. nyIn < 3) WRITE(ERROR_UNIT,902) nyIn
         ierr = 1
         RETURN
      ENDIF
      nz = nzIn
      nx = nxIn
      ny = nyIn
      IF (.NOT. lis3d) ny = 1
      ncell = (nz - 1)*(nx - 1)*MAX(1, (ny - 1))
      ngrd = nz*nx*ny
      nzx = nz*nx
      nzm1 = nz - 1
      nzm1_nxm1 = (nz - 1)*(nx - 1)
      IF (ALLOCATED(slow)) DEALLOCATE(slow)
      ALLOCATE(slow(ncell))
      slow(:) = zero
  900 FORMAT('fteik_model_setGridSize: ERROR nz is too small', I0)
  901 FORMAT('fteik_model_setGridSize: ERROR nx is too small', I0)
  902 FORMAT('fteik_model_setGridSize: ERROR ny is too small', I0)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the model origin.
!>    @param[in] z0In    z origin (meters).
!>    @param[in] x0In    x origin (meters).
!>    @param[in] y0In    y origin (meters).
!>    @ingroup model
      SUBROUTINE fteik_model_setOrigin3D(z0In, x0In, y0In) &
                 BIND(C, NAME='fteik_model_setOrigin3D')
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
!>    @brief Sets the model origin.
!>    @param[in] z0In    z origin (meters).
!>    @param[in] x0In    x origin (meters).
!>    @ingroup model
      SUBROUTINE fteik_model_setOrigin2D(z0In, x0In) &
                 BIND(C, NAME='fteik_model_setOrigin2D')
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x0In, z0In
      z0 = z0In
      x0 = x0In
      y0 = zero 
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
!>    @ingroup model
      SUBROUTINE fteik_model_getGridSize(nzOut, nxOut, nyOut,     &
                                         ngrdOut, ncellOut, ierr) &
                 BIND(C, NAME='fteik_model_getGridSize')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: nzOut, nxOut, nyOut, ngrdOut, ncellOut, ierr
      ierr = 0
      nzOut = 0
      nxOut = 0
      nyOut = 0
      ngrdOut = 0
      ncellOut = 0
      IF (nx < 1 .OR. nz < 1) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      nzOut = nz
      nxOut = nx
      nyOut = ny
      ngrdOut = ngrd
      ncellOut = ncell
  900 FORMAT('fteik_model_getGridSize: Grid not yet set')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the nodal velocity model.
!>
!>    @param[in] ng     Number of grid points in velocity.  This should be equal to
!>                      [nz x nx x ny] for 3D or [nz x nx] for 2D.
!>    @param[in] order  Defines the column major order of vel. \n
!>                      If order == FTEIK_XYZ_ORDERING then vel has dimension 
!>                      [nz x ny x nx] where nz is leading dimension 1 and ny is
!>                      leading dimension 2. \n
!>                      If order == FTEIK_ZYX_ORDERING then vel has dimension
!>                      [nx x ny x nz] where nx is leading dimension 1 and ny is
!>                      leading dimension 2. \n
!>                      If order == FTEIK_ZXY_ORDERING or FTEIK_NATURAL_ORDERING
!>                      then vel has dimension [nz x nx x ny] where nz is leading
!>                      dimension 1 and nx is leading dimension 2. \n
!>                      If order == FTEIK_XZ_ORDERING then vel has dimension [nx x nz]
!>                      where nx is the leading dimension and the model is 2D. \n
!>                      If order == FTEIK_ZX_ORDERING or FTEIK_NATURAL_ORDERING
!>                      then vel has dimension [nz x nx] where nz is the leading
!>                      dimension and the model is 2D.
!>    @param[in] vel    Nodal velocity model whose leading dimension are given
!>                      by order.
!>
!>    @param[out] ierr  0 indicates success.
!>
!>    @ingroup model
      SUBROUTINE fteik_model_setNodalVelocityModel64f(ng, order, vel, ierr)  &
      BIND(C, NAME='fteik_model_setNodalVelocityModel64f')
      USE FTEIK_CONSTANTS64F, ONLY : one, zero
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: ng, order
      REAL(C_DOUBLE), INTENT(IN) :: vel(ng)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) vmin
      INTEGER(C_INT) indx(8), icell, ix, iy, iz
      ierr = 0
      lhaveModel = .FALSE.
      IF (ng < ngrd) THEN
         WRITE(ERROR_UNIT,*) &
         'fteik_model_setNodalVelocityModel64f: ng should be at least', ngrd 
         ierr = 1
         RETURN
      ENDIF
      IF (ng < 1) THEN
         WRITE(ERROR_UNIT,*) &
         'fteik_model_setNodalVelocityModel64f: No points in model', ng
         ierr = 1
         RETURN
      ENDIF
      vmin = MINVAL(vel(1:ng))
      IF (vmin <= zero) THEN
         WRITE(ERROR_UNIT,*) &
         'fteik_model_setNodalVelocityModel64f: vmin must be positive', vmin
         ierr = 1
         RETURN
      ENDIF
      ! Average the nodes to make the cells
      slow(:) = zero
      IF (lis3dModel) THEN
         IF (order == FTEIK_XYZ_ORDERING) THEN
            DO iy=1,ny-1
               DO ix=1,nx-1
                  DO iz=1,nz-1
                     indx(1) = (iz - 1)*nx*ny + (iy - 1)*nx + ix
                     indx(2) = (iz - 1)*nx*ny + (iy - 1)*nx + ix + 1
                     indx(3) = (iz - 1)*nx*ny + (iy    )*nx + ix
                     indx(4) = (iz - 1)*nx*ny + (iy    )*nx + ix + 1
                     indx(5) = (iz    )*nx*ny + (iy - 1)*nx + ix
                     indx(6) = (iz    )*nx*ny + (iy - 1)*nx + ix + 1
                     indx(7) = (iz    )*nx*ny + (iy    )*nx + ix
                     indx(8) = (iz    )*nx*ny + (iy    )*nx + ix + 1
                     icell = fteik_model_velGrid2indexF(iz, ix, iy, nzm1, nzm1_nxm1)
                     slow(icell) = SUM(vel(indx(1:8)))
                  ENDDO
               ENDDO
            ENDDO
         ELSEIF (order == FTEIK_ZYX_ORDERING) THEN
            DO iy=1,ny-1
               DO ix=1,nx-1
                  DO iz=1,nz-1
                     indx(1) = (ix - 1)*ny*nz + (iy - 1)*nz + iz
                     indx(2) = (ix - 1)*ny*nz + (iy - 1)*nz + iz + 1
                     indx(3) = (ix - 1)*ny*nz + (iy    )*nz + iz
                     indx(4) = (ix - 1)*ny*nz + (iy    )*nz + iz + 1
                     indx(5) = (ix    )*ny*nz + (iy - 1)*nz + iz
                     indx(6) = (ix    )*ny*nz + (iy - 1)*nz + iz + 1 
                     indx(7) = (ix    )*ny*nz + (iy    )*nz + iz
                     indx(8) = (ix    )*ny*nz + (iy    )*nz + iz + 1
                     icell = fteik_model_velGrid2indexF(iz, ix, iy, nzm1, nzm1_nxm1)
                     slow(icell) = SUM(vel(indx(1:8)))*0.125d0
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ((order /= FTEIK_NATURAL_ORDERING .AND. order /= FTEIK_ZXY_ORDERING).AND. &
                verbose > 0) THEN
               WRITE(OUTPUT_UNIT,*) &
               'fteik_model_setNodalVelocityModel64f: Defaulting ot zxy' 
            ENDIF 
            DO iy=1,ny-1
               DO ix=1,nx-1
                  DO iz=1,nz-1
                     indx(1) = fteik_model_grid2indexF(iz,   ix,   iy,   nz, nzx)
                     indx(2) = fteik_model_grid2indexF(iz+1, ix,   iy,   nz, nzx)
                     indx(3) = fteik_model_grid2indexF(iz,   ix+1, iy,   nz, nzx)
                     indx(4) = fteik_model_grid2indexF(iz+1, ix+1, iy,   nz, nzx)
                     indx(5) = fteik_model_grid2indexF(iz,   ix,   iy,   nz, nzx)
                     indx(6) = fteik_model_grid2indexF(iz+1, ix,   iy,   nz, nzx)
                     indx(7) = fteik_model_grid2indexF(iz,   ix+1, iy,   nz, nzx)
                     indx(8) = fteik_model_grid2indexF(iz+1, ix+1, iy,   nz, nzx)
                     icell = fteik_model_velGrid2indexF(iz, ix, iy, nzm1, nzm1_nxm1)
                     slow(icell) = SUM(vel(indx(1:8)))*0.125d0
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ! Same activity but for 2D
      ELSE
         IF (order == FTEIK_XZ_ORDERING) THEN
            DO ix=1,nx-1
               DO iz=1,nz-1
                  indx(1) = (iz - 1)*nx + ix
                  indx(2) = (iz - 1)*nx + ix + 1
                  indx(3) = (iz    )*nx + ix
                  indx(4) = (iz    )*nx + ix + 1
                  icell = (ix - 1)*(nz - 1) + iz
                  slow(icell) = SUM(vel(indx(1:4)))*0.25d0
               ENDDO
            ENDDO
         ELSE
            IF ((order /= FTEIK_NATURAL_ORDERING .AND. order /= FTEIK_ZXY_ORDERING).AND. &
                verbose > 0) THEN
               WRITE(OUTPUT_UNIT,*) &
               'fteik_model_setNodalVelocityModel64f: Defaulting to zx'
            ENDIF
            DO iz=1,nz-1
               DO ix=1,nx-1
                  indx(1) = (ix - 1)*nz + iz
                  indx(2) = (ix - 1)*nz + iz + 1
                  indx(3) = (ix    )*nz + iz
                  indx(4) = (ix    )*nz + iz + 1
                  icell = (ix - 1)*(nz - 1) + iz
                  slow(icell) = SUM(vel(indx(1:4)))*0.25d0
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      IF (MINVAL(slow(1:ncell)) <= zero) THEN
         WRITE(ERROR_UNIT,*) 'fteik_model_setNodalVelocityModel64f: Internal error'
         ierr = 1
         RETURN
      ENDIF
      slow(1:ncell) = one/slow(1:ncell)
      lhaveModel = .TRUE.
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the cell-based velocity model.
!>
!>    @param[in] ng     Number of grid points in velocity.  This should be equal to
!>                      [nz-1 x nx-1 x ny-1] for 3D or [nz-1 x nx-1] for 2D.
!>    @param[in] order  Defines the column major order of vel. \n
!>                      If order == FTEIK_XYZ_ORDERING then vel has dimension 
!>                      [nz-1 x ny-1 x nx-1] where nz-1 is leading dimension 1 and ny-1 is
!>                      leading dimension 2. \n
!>                      If order == FTEIK_ZYX_ORDERING then vel has dimension
!>                      [nx-1 x ny-1 x nz-1] where nx-1 is leading dimension 1 and ny-1 is
!>                      leading dimension 2. \n
!>                      If order == FTEIK_ZXY_ORDERING or FTEIK_NATURAL_ORDERING
!>                      then vel has dimension [nz-1 x nx-1 x ny-1] where nz-1 is leading
!>                      dimension 1 and nx-1 is leading dimension 2. \n
!>                      If order == FTEIK_XZ_ORDERING then vel has dimension [nx-1 x nz-1]
!>                      where nx-1 is the leading dimension and the model is 2D. \n
!>                      If order == FTEIK_ZX_ORDERING or FTEIK_NATURAL_ORDERING
!>                      then vel has dimension [nz-1 x nx-1] where nz-1 is the leading
!>                      dimension and the model is 2D.
!>    @param[in] vel    Cell-based velocity model whose leading dimension are given
!>                      by order.
!>
!>    @param[out] ierr  0 indicates success.
!>
!>    @ingroup model
      SUBROUTINE fteik_model_setCellVelocityModel64f(nc, order, vel, ierr)  &
      BIND(C, NAME='fteik_model_setCellVelocityModel64f')
      USE FTEIK_CONSTANTS64F, ONLY : one, zero
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: nc, order
      REAL(C_DOUBLE), INTENT(IN) :: vel(nc)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) vmin
      INTEGER(C_INT) icell, ix, iy, iz, jcell
      ierr = 0
      lhaveModel = .FALSE.
      IF (nc < ncell) THEN
         WRITE(ERROR_UNIT,*) &
         'fteik_model_setCellVelocityModel64f: nc should be at least', ncell
         ierr = 1
         RETURN
      ENDIF
      IF (nc < 1) THEN
         WRITE(ERROR_UNIT,*) &
         'fteik_model_setCellVelocityModel64f: No points in model', nc
         ierr = 1
         RETURN
      ENDIF
      vmin = MINVAL(vel(1:nc))
      IF (vmin <= zero) THEN
         WRITE(ERROR_UNIT,*) &
         'fteik_model_setCellVelocityModel64f: vmin must be positive', vmin
         ierr = 1
         RETURN
      ENDIF
      IF (lis3dModel) THEN
         IF (order == FTEIK_XYZ_ORDERING) THEN
            DO iy=1,ny-1
               DO ix=1,nx-1
                  DO iz=1,nz-1
                     jcell = (iz - 1)*(nx - 1)*(ny - 1) + (iy - 1)*(nx - 1) + ix
                     icell = fteik_model_velGrid2indexF(iz, ix, iy, nzm1, nzm1_nxm1)
                     slow(icell) = vel(jcell)
                  ENDDO
               ENDDO
            ENDDO
         ELSEIF (order == FTEIK_ZYX_ORDERING) THEN
            DO iy=1,ny-1
               DO ix=1,nx-1
                  DO iz=1,nz-1
                     jcell = (ix - 1)*(nz - 1)*(ny - 1) + (iy - 1)*(nz - 1) + iz
                     icell = fteik_model_velGrid2indexF(iz, ix, iy, nzm1, nzm1_nxm1)
                     slow(icell) = vel(jcell)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ((order /= FTEIK_NATURAL_ORDERING .AND. order /= FTEIK_ZXY_ORDERING).AND. &
                verbose > 0) THEN
               WRITE(OUTPUT_UNIT,*) &
               'fteik_model_setCellVelocityModel64f: Defaulting ot zxy'
            ENDIF
            slow(1:ncell) = one/vel(1:ncell)
            lhaveModel = .TRUE.
            RETURN
         ENDIF
      ! Same activity but for 2D
      ELSE
         IF (order == FTEIK_XZ_ORDERING) THEN
            DO ix=1,nx-1
               DO iz=1,nz-1
                  jcell = (iz - 1)*(nx - 1) + ix
                  icell = (ix - 1)*(nz - 1) + iz
                  slow(icell) = vel(jcell)
               ENDDO
            ENDDO
         ELSE
            IF ((order /= FTEIK_NATURAL_ORDERING .AND. order /= FTEIK_ZXY_ORDERING).AND. &
                verbose > 0) THEN
               WRITE(OUTPUT_UNIT,*) &
               'fteik_model_setCellVelocityModel64f: Defaulting ot zx' 
            ENDIF
            slow(1:ncell) = one/vel(1:ncell)
            lhaveModel = .TRUE.
            RETURN
         ENDIF
      ENDIF
      IF (MINVAL(slow(1:ncell)) <= zero) THEN
         WRITE(ERROR_UNIT,*) 'fteik_model_setCellVelocityModel64f: Internal error'
         ierr = 1 
         RETURN
      ENDIF
      slow(1:ncell) = one/slow(1:ncell)
      lhaveModel = .TRUE.
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
!>    @ingroup model
      SUBROUTINE fteik_model_setVelocityModel64f(nv, vel, ierr) &
                 BIND(C, NAME='fteik_model_setVelocityModel64f')
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
         WRITE(ERROR_UNIT,900) ncell, nv 
         ierr = 1 
         RETURN
      ENDIF
      IF (nv < 1) THEN
         WRITE(ERROR_UNIT,901) nv
         ierr = 1 
         RETURN
      ENDIF
      IF (MINVAL(vel) <= zero) THEN
         WRITE(ERROR_UNIT,902)
         ierr = 1
         RETURN
      ENDIF
      !IF (.NOT.ALLOCATED(slow)) ALLOCATE(slow(ncell))
      DO 1 i=1,ncell
         slow(i) = one/vel(i)
    1 CONTINUE  
      lhaveModel = .TRUE.
  900 FORMAT('fteik_model_setVelocityModel64f: ERROR - ncell /= nv', I0, I0)
  901 FORMAT('fteik_model_setVelocityModel64f: ERROR - no cells in vel', I0)
  902 FORMAT('fteik_model_setVelocityModel64f: All velocities must be positive')
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Finalizes the model module by releasing memory and setting variables to
!>           undefined types.
!>
!>    @ingroup model
      SUBROUTINE fteik_model_free() &
                 BIND(C, NAME='fteik_model_free')
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      lis3dModel = .TRUE.
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
!>    @ingroup model
      SUBROUTINE fteik_model_getVelocityModel64f(nwork, nv, vel, ierr) &
                 BIND(C, NAME='fteik_model_getVelocityModel64f')
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
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.lhaveModel) THEN
         WRITE(ERROR_UNIT,901)
         ierr = 1
         vel(:) = zero
         RETURN
      ENDIF 
      IF (nv > ncell) vel(ncell+1:nwork) = zero
      DO 1 i=1,ncell
         vel(i) = one/slow(i) 
    1 CONTINUE
  900 FORMAT('fteik_model_getVelocityModel64f: nv < ncell')
  901 FORMAT('fteik_model_getVelocityModel64f: Velocity model not yet set')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Determines if the model is 2D or 3D.
!>    @param[out] lis3d    If true then the model is 3D.
!>    @param[out] ierr     0 indicates usccess.
!>    @ingroup model
      SUBROUTINE fteik_model_isModel3D(lis3d, ierr) &
      BIND(C, NAME='fteik_model_isModel3D')
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(OUT) :: ierr
      LOGICAL(C_BOOL), INTENT(OUT) :: lis3d
      ierr = 0
      lis3d = lis3dModel
      IF (nx < 1 .OR. nz < 1) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
  900 FORMAT('fteik_model_isModel3D: Model not yet initialized')
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
!>    @ingroup model
      SUBROUTINE fteik_model_setVelocityModel32f(nv, vel4, ierr) &
                 BIND(C, NAME='fteik_model_setVelocityModel32f')
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
         WRITE(ERROR_UNIT,900) ncell, nv
         ierr = 1
         RETURN
      ENDIF
      IF (nv < 1) THEN
         WRITE(ERROR_UNIT,901) nv
         ierr = 1
         RETURN
      ENDIF
      IF (MINVAL(vel4) <= 0.0) THEN
         WRITE(ERROR_UNIT,902)
         ierr = 1
         RETURN
      ENDIF
      !IF (.NOT.ALLOCATED(slow)) ALLOCATE(slow(ncell))
      DO 1 i=1,ncell
         slow(i) = one/DBLE(vel4(i))
    1 CONTINUE
      lhaveModel = .TRUE.
  900 FORMAT('fteik_model_setVelocityModel64f: ERROR - ncell /= nv', I0, I0) 
  901 FORMAT('fteik_model_setVelocityModel64f: ERROR - no cells in vel', I0) 
  902 FORMAT('fteik_model_setVelocityModel64f: All velocities must be positive')
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
!>    @ingroup model
      INTEGER(C_INT) FUNCTION fteik_model_grid2indexF(i, j, k, nz, nzx) &
      BIND(C, NAME='fteik_model_grid2indexF')                           &
      RESULT(grid2indexF)
      !!$OMP DECLARE SIMD(fteik_model_grid2indexF) UNIFORM(nz, nzx)
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
!>    @ingroup model
      SUBROUTINE fteik_model_index2gridF(igrd, i, j, k, ierr) &
      BIND(C, NAME='fteik_model_index2gridF')
      !!$OMP DECLARE SIMD(fteik_model_index2gridF) UNIFORM(ierr)
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
!>    @ingroup model
      INTEGER(C_INT) FUNCTION fteik_model_velGrid2indexF(i, j, k, nzm1, nzm1_nxm1)   &
      BIND(C, NAME='fteik_model_velGrid2indexF')                                     &
      RESULT(velGrid2IndexF)
      !!$OMP DECLARE SIMD(fteik_model_velGrid2indexF) UNIFORM(nzm1, nzm1_nxm1)
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
