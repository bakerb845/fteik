!> @defgroup h5io HDF5 I/O
!> @ingroup solver2d
!> @ingroup solver3d
!> @brief HDF5-based file input/output.
!> @copyright Ben Baker distributed under the MIT license.
MODULE FTEIK_H5IO64F
    USE ISO_C_BINDING
    IMPLICIT NONE
    CHARACTER(C_CHAR), PRIVATE, SAVE :: fname(4096)
    LOGICAL(C_BOOL), PRIVATE, SAVE :: linitH5FL  


    INTERFACE

         INTEGER(C_SIZE_T) FUNCTION strlenF(fileName) &
         BIND(C, NAME='strlen')
         USE ISO_C_BINDING
         CHARACTER(C_CHAR), INTENT(IN) :: fileName(*)
         END FUNCTION

         INTEGER(C_INT) FUNCTION fteik_h5io_finalize() &
         BIND(C, NAME='fteik_h5io_finalize')
         USE ISO_C_BINDING
         END FUNCTION

         INTEGER(C_INT) &
         FUNCTION fteik_h5io_initialize(fileName,   &
                                        nz, nx, ny, &
                                        dz, dx, dy, &
                                        z0, x0, y0) &
         BIND(C, NAME='fteik_h5io_initialize')
         USE ISO_C_BINDING
         IMPLICIT NONE 
         CHARACTER(C_CHAR), INTENT(IN) :: fileName(*)
         REAL(C_DOUBLE), VALUE, INTENT(IN) :: dz, dx, dy, z0, x0, y0
         INTEGER(C_INT), VALUE, INTENT(IN) :: nz, nx, ny
         END FUNCTION

         INTEGER(C_INT64_T) FUNCTION fteik_h5io_openFileReadWriteF(fileName) &
         BIND(C, NAME='h5io_openFileReadWriteF')
         USE ISO_C_BINDING 
         IMPLICIT NONE 
         CHARACTER(C_CHAR), INTENT(IN) :: fileName(*)
         END FUNCTION

         INTEGER(C_INT) &
         FUNCTION fteik_h5io_writeLevelScheduleF(h5fl, sweep,   &
                                                 nz, nx, ny,    &
                                                 levelSchedule) &
         BIND(C, NAME='fteik_h5io_writeLevelScheduleF')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT64_T), INTENT(IN), VALUE :: h5fl
         INTEGER(C_INT), INTENT(IN), VALUE :: nz, nx, ny, sweep
         INTEGER(C_INT16_T), INTENT(IN) :: levelSchedule(nz*nx*ny)
         END FUNCTION

         INTEGER(C_INT) &
         FUNCTION fteik_h5io_writeTravelTimes32fF(h5fl, ttName, &
                                                  nz, nx, ny,   &
                                                  tt)           &
         BIND(C, NAME='fteik_h5io_writeTravelTimes32fF')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT64_T), INTENT(IN), VALUE :: h5fl
         CHARACTER(C_CHAR), INTENT(IN) :: ttName(*)
         INTEGER(C_INT), INTENT(IN), VALUE :: nx, ny, nz
         REAL(C_FLOAT), INTENT(IN) :: tt(nx*ny*nz)
         END FUNCTION

         INTEGER(C_INT) FUNCTION fteik_h5io_closeFileF(h5fl) &
         BIND(C, NAME='h5io_closeFileF')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT64_T), INTENT(IN), VALUE :: h5fl
         END FUNCTION

         INTEGER(C_INT) FUNCTION h5io_writeArray16iF(fid, dataName, n, x) &
         BIND(C, NAME='h5io_writeArray16iF')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT64_T), INTENT(IN), VALUE :: fid
         INTEGER(C_INT), INTENT(IN), VALUE :: n
         CHARACTER(C_CHAR), INTENT(IN) :: dataName(*)
         INTEGER(C_INT16_T), INTENT(IN) :: x(n)
         END FUNCTION

         INTEGER(C_INT) FUNCTION fteik_h5io_writeVelocityModel16iF(h5fl, velName, &
                                                                   nz, nx, ny,    &
                                                                   vel)           &
         BIND(C, NAME='fteik_h5io_writeVelocityModel16iF')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT64_T), INTENT(IN), VALUE :: h5fl
         CHARACTER(C_CHAR), INTENT(IN) :: velName(*)
         INTEGER(C_INT), INTENT(IN), VALUE :: nx, ny, nz
         INTEGER(C_INT16_T), INTENT(IN) :: vel((nx-1)*(ny-1)*(nz-1))
         END FUNCTION

    END INTERFACE


    CONTAINS
!----------------------------------------------------------------------------------------!
!                                 Begin the Code                                         !
!----------------------------------------------------------------------------------------!
!>    @brief Initializes the HDF5 archive.
!>    @param[in] fileName   Name of HDF5 archive.
!>    @result 0 indicates success.
!>    @ingroup h5io
      INTEGER(C_INT) FUNCTION fteik_h5io_initializeF(fileName) &
      RESULT(ierr) BIND(C, NAME='fteik_h5io_initializeF')
      USE FTEIK_MODEL64F, ONLY : nz, nx, ny, z0, x0, y0, dx, dy, dz, lhaveModel
      USE ISO_C_BINDING
      IMPLICIT NONE
      CHARACTER(C_CHAR), INTENT(IN) :: fileName(*)
      INTEGER(C_SIZE_T) lenos
      IF (.NOT.lhaveModel) THEN
         WRITE(*,*) 'fteik_h5io_initializeF: File not initialized'
         ierr =-1
         RETURN
      ENDIF
      fname(:) = C_NULL_CHAR
      linitH5FL = .FALSE.
      ierr = fteik_h5io_initialize(fileName,   &
                                   nz, nx, ny, &
                                   dz, dx, dy, &
                                   z0, x0, y0)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'Error initializing file'
         RETURN
      ENDIF
      lenos = strlenF(fileName)
      fname(1:lenos) = fileName(1:lenos)
      linitH5FL = .TRUE.
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Finalizes the HDF5 archive and writes the XDMF file.
!>    @result 0 indicates success.
!>    @ingroup h5io
      INTEGER(C_INT) FUNCTION fteik_h5io_finalizeF() &
      RESULT(ierr) BIND(C, NAME='fteik_h5io_finalizeF')
      ierr = fteik_h5io_finalize()
      !ierr = eik_h5io_closeFileF(h5fl) 
      linitH5FL = .FALSE.
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Writes the level of each node in the grid.
!>
!>    @param[in] nx        Number of x grid points in travel time field.
!>    @param[in] ny        Number of y grid points in travel time field.
!>    @param[in] ngrd      Number of grid points in travel time field.
!>    @param[in] nLevels   Number of levels in level scheduling method.
!>    @param[in] levelPtr  Maps to the first node of the level'th level.  This is
!>                         a vector of dimension [nLevels+1]. 
!>    @param[in] ijkv      Contains the (i, j, k, level) of each node in the solver.
!>                         This is a vector of dimension [4 x ngrd] with leading
!>                         dimension 4.
!>    @param[out] levels   The level of each grid point in the mesh.  This is a 
!>                         vector of dimension [nz x ny x nx = ngrd] with leading
!>                         dimension nz.
!>    @ingroup h5io
      PURE SUBROUTINE fteik_h5io_ijkv2levelsF(nx, ny, ngrd, nLevels, &
                                             levelPtr, ijkv, levels) &
      BIND(C, NAME='fteik_h5io_ijkv2levelsF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: ngrd, nLevels, nx, ny
      INTEGER(C_INT), INTENT(IN) :: levelPtr(nLevels+1), ijkv(4*ngrd)
      INTEGER(C_INT16_T), INTENT(OUT) :: levels(ngrd)
      INTEGER(C_INT) level, indx, ix, iy, iz, node
      levels(:) = 0
      DO 1 level=1,nLevels
         DO 2 node=levelPtr(level),levelPtr(level+1)-1
            iz = ijkv(4*(node-1)+1)
            ix = ijkv(4*(node-1)+2)
            iy = ijkv(4*(node-1)+3)
            indx = (iz - 1)*nx*ny + (iy - 1)*nx + ix
            levels(indx) = INT2(level)
    2    CONTINUE
    1 CONTINUE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Convenience function to write all the levels schedules in the solver
!>           for all 8 sweeps to the archive.
!>    @result 0 indicates success.
!>    @ingroup h5io
      INTEGER(C_INT) FUNCTION fteik_h5io_writeLevelSchedulesF( ) &
      RESULT(ierr) BIND(C, NAME='fteik_h5io_writeLevelSchedulesF')
      USE FTEIK_SOLVER64F, ONLY : ijkv1, ijkv2, ijkv3, ijkv4,    &
                                  ijkv5, ijkv6, ijkv7, ijkv8,    &
                                  levelPtr, nLevels
      USE FTEIK_MODEL64F, ONLY : nz, nx, ny, ngrd, lhaveModel, lis3dModel
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT16_T), ALLOCATABLE :: levels(:)
      INTEGER(C_INT64_T) h5fl
      INTEGER(C_INT) ierr2, sweep
      ierr = 0
      IF (.NOT.linitH5FL) THEN
         WRITE(*,*) 'fteik_h5io_writeLevelSchedulesF: File not initialized'
         ierr =-1 
         RETURN
      ENDIF
      IF (.NOT.lhaveModel) THEN
         WRITE(*,*) 'fteik_h5io_writeLevelSchedulesF: Grid not set'
         ierr =-1 
         RETURN
      ENDIF
      IF (.NOT.lis3dModel) THEN
         WRITE(*,*) 'fteik_h5io_writeLevelSchedulesF: Not yet done for 2D'
         ierr =-1
         RETURN
      ENDIF
      ! Open file 
      h5fl = fteik_h5io_openFileReadWriteF(fname)
      ALLOCATE(levels(ngrd))
      DO 1 sweep=1,8
         IF (sweep == 1) THEN 
            CALL fteik_h5io_ijkv2levelsF(nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv1, levels)
         ELSEIF (sweep == 2) THEN
            CALL fteik_h5io_ijkv2levelsF(nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv2, levels)
         ELSEIF (sweep == 3) THEN
            CALL fteik_h5io_ijkv2levelsF(nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv3, levels)
         ELSEIF (sweep == 4) THEN
            CALL fteik_h5io_ijkv2levelsF(nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv4, levels)
         ELSEIF (sweep == 5) THEN 
            CALL fteik_h5io_ijkv2levelsF(nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv5, levels)
         ELSEIF (sweep == 6) THEN
            CALL fteik_h5io_ijkv2levelsF(nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv6, levels)
         ELSEIF (sweep == 7) THEN
            CALL fteik_h5io_ijkv2levelsF(nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv7, levels)
         ELSEIF (sweep == 8) THEN
            CALL fteik_h5io_ijkv2levelsF(nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv8, levels)
         ENDIF
         ierr2 = fteik_h5io_writeLevelScheduleF(h5fl, sweep, nz, nx, ny, levels)
         IF (ierr2 /= 0) THEN
            WRITE(*,*) 'fteik_h5io_writeLevelSchedulesF: Error sweep schedule:', sweep
            ierr = ierr + 1
         ENDIF
    1 CONTINUE
      DEALLOCATE(levels)
      ierr2 = fteik_h5io_closeFileF(h5fl)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Writes the travel time field in the solver to the archive.
!>    @param[in] dataName   Name of travel time field to write.
!>    @result 0 indicates success.
!>    @ingroup h5io
      INTEGER(C_INT) FUNCTION fteik_h5io_writeTravelTimesF(dataName) &
      RESULT(ierr) BIND(C, NAME='fteik_h5io_writeTravelTimesF')
      !USE FTEIK_MODEL64F, ONLY : fteik_model_grid2indexF !fteik_model_grid2indexF
      USE FTEIK_MODEL64F, ONLY : nz, nx, ny, ngrd, lis3dModel
      USE FTEIK_SOLVER64F, ONLY : ttimes3d => ttimes
      USE FTEIK_SOLVER64F, ONlY : lhaveTimes3d => lhaveTimes
      USE FTEIK2D_SOLVER64F, ONLY : ttimes2d => ttimes
      USE FTEIK2D_SOLVER64F, ONLY : lhaveTimes2d => lhaveTimes
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTERFACE
         PURE INTEGER(C_INT)                                       &
         FUNCTION fteik_model_grid2indexF(i, j, k, nz, nzx)        &
         BIND(C, NAME='fteik_model_grid2indexF')                   &
         RESULT(grid2indexF)
         !!$OMP DECLARE SIMD(fteik_model_grid2indexF) UNIFORM(nz, nzx)
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nzx 
         END FUNCTION
      END INTERFACE
      CHARACTER(C_CHAR), INTENT(IN) :: dataName(*)
      INTEGER(C_INT64_T) h5fl
      REAL(C_FLOAT), ALLOCATABLE :: ttout(:)
      INTEGER(C_INT) ierr2, indx, ix, iy, iz, jndx, jz, nzx
      ierr = 0
      IF (.NOT.linitH5FL) THEN
         WRITE(*,*) 'fteik_h5io_writeTravelTimesF: File not initialized'
         ierr =-1 
         RETURN
      ENDIF
      IF (lis3dModel) THEN
         IF (.NOT.lhaveTimes3d) THEN
            WRITE(*,*) 'fteik_h5io_writeTravelTimesF: 3D travel times not yet computed'
            ierr =-1
            RETURN
         ENDIF
         nzx = nz*nx
         ALLOCATE(ttout(ngrd))
         ttout(:) = 0.0
         DO 1 iz=1,nz
            DO 2 iy=1,ny
               DO 3 ix=1,nx
                  !jz = nz + 1 - iz ! flip coordinate syste
                  indx = fteik_model_grid2indexF(iz, ix, iy, nz, nzx)
                  jndx = (iz - 1)*nx*ny + (iy - 1)*nx + ix
                  !print *, jndx
                  ttout(jndx) = REAL(ttimes3d(indx))
    3          CONTINUE
    2       CONTINUE
    1    CONTINUE
      ELSE
         IF (.NOT.lhaveTimes2d) THEN
            WRITE(*,*) 'fteik_h5io_writeTravelTimesF: 2D travel times not yet computed'
            ierr =-1 
            RETURN
         ENDIF
         nzx = nz*nx
         ALLOCATE(ttout(ngrd))
         ttout(:) = 0.0 
         DO 11 iz=1,nz
            DO 12 ix=1,nx
               jz = nz + 1 - iz ! flip coordinate system
               indx = fteik_model_grid2indexF(jz, ix, 1, nz, nzx)
               jndx = (iz - 1)*nx + ix
               !print *, jndx
               ttout(jndx) = REAL(ttimes2d(indx))
   12       CONTINUE
   11    CONTINUE
      ENDIF
      ! Open and write the file
      h5fl = fteik_h5io_openFileReadWriteF(fname)
      ierr = fteik_h5io_writeTravelTimes32fF(h5fl, dataName, nz, nx, ny, ttout)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_h5io_writeTravelTimesF: Failed to travel times'
      ENDIF
      DEALLOCATE(ttout)
      ! Close it
      ierr2 = fteik_h5io_closeFileF(h5fl)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Writes the velocity model in the solver to the archive.
!>    @param[in] dataName   Name of velocity file to write.
!>    @ingroup h5io
      INTEGER(C_INT) FUNCTION fteik_h5io_writeVelocityModelF(dataName) &
                     RESULT(ierr) BIND(C, NAME='fteik_h5io_writeVelocityModelF')
      USE FTEIK_MODEL64F, ONLY : slow, ncell, nx, ny, nz, nzm1, &
                                 nzm1_nxm1, lhaveModel, lis3dModel
      !USE FTEIK_MODEL64F, ONLY : fteik_model_velGrid2indexF
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTERFACE
         PURE INTEGER(C_INT)                                             &
         FUNCTION fteik_model_velGrid2indexF(i, j, k, nzm1, nzm1_nxm1)   &
         BIND(C, NAME='fteik_model_velGrid2indexF')                      &
         RESULT(velGrid2IndexF)
         !!$OMP DECLARE SIMD(fteik_model_velGrid2indexF) UNIFORM(nzm1, nzm1_nxm1)
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nzm1, nzm1_nxm1
         END FUNCTION
      END INTERFACE
      CHARACTER(C_CHAR), INTENT(IN) :: dataName(*)
      INTEGER(C_INT64_T) h5fl
      INTEGER(C_INT16_T), ALLOCATABLE :: vout(:)
      INTEGER ierr2, indx, ix, iy, iz, jndx, jz
      ierr = 0
      IF (.NOT.linitH5FL) THEN
         WRITE(*,*) 'fteik_h5io_writeVelocityModelF: File not initialized'
         ierr =-1
         RETURN
      ENDIF
      IF (.NOT.lhaveModel) THEN
         WRITE(*,*) 'fteik_h5io_writeVelocityModelF: Velocity model not set'
         ierr =-1
         RETURN
      ENDIF
      ! First dimension in topology is slowest changing (z)
      ALLOCATE(vout(ncell))
      vout(:) = 0
      IF (lis3dModel) THEN
         DO 1 iz=1,nz-1
            DO 2 iy=1,ny-1
               DO 3 ix=1,nx-1
                  !jz = iz !nz-iz ! flip coordinate syste
                  indx = fteik_model_velGrid2indexF(iz, ix, iy, nzm1, nzm1_nxm1) 
                  jndx = (iz-1)*(nx-1)*(ny-1) + (iy-1)*(nx-1) + ix
                  vout(jndx) = INT2(INT(1.d0/slow(indx)  + 0.5d0))
                  !print *, vout(jndx), 1.d0/slow(indx)
    3          CONTINUE
    2       CONTINUE
    1    CONTINUE
      ELSE
         iy = 1
         DO 11 iz=1,nz-1
            DO 12 ix=1,nx-1
               jz = nz-iz ! flip coordinate syste
               indx = fteik_model_velGrid2indexF(jz, ix, iy, nzm1, nzm1_nxm1) 
               jndx = (iz-1)*(nx-1) + ix
               vout(jndx) = INT2(INT(1.d0/slow(indx)  + 0.5d0))
               !print *, vout(jndx), 1.d0/slow(indx)
   12       CONTINUE
   11    CONTINUE
      ENDIF
      ! Open and write the file
      h5fl = fteik_h5io_openFileReadWriteF(fname) 
      ierr = fteik_h5io_writeVelocityModel16iF(h5fl, dataName, nz, nx, ny, vout)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_h5io_writeVelocityModelF: Failed to write model'
      ENDIF
      DEALLOCATE(vout)
      ! Close it
      ierr2 = fteik_h5io_closeFileF(h5fl)
      RETURN
      END

END MODULE !FTEIK_H5IO64F
