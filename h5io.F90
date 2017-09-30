
      INTEGER(C_INT) FUNCTION fteik_h5io_initializeF(fileName) &
                     RESULT(ierr) BIND(C, NAME='fteik_h5io_initializeF')
      USE FTEIK_UTILS64F, ONLY : lhaveGrid, lhaveGridSpacing
      USE FTEIK_H5IO, ONLY : fteik_h5io_initialize, strlenF, fname, linitH5FL
      USE ISO_C_BINDING
      CHARACTER(C_CHAR), INTENT(IN) :: fileName(*)
      INTEGER(C_SIZE_T) lenos
      IF (.NOT.lhaveGrid .OR. .NOT.lhaveGridSpacing) THEN
         WRITE(*,*) 'fteik_h5io_initializeF: File not initialized'
         ierr =-1
         RETURN
      ENDIF
      fname(:) = C_NULL_CHAR
      linitH5FL = .FALSE.
      ierr = fteik_h5io_initialize(fileName)
      IF (ierr /= 0) WRITE(*,*) 'Error initializing file'
      lenos = strlenF(fileName)
      fname(1:lenos) = fileName(1:lenos)
      linitH5FL = .TRUE.
      RETURN
      END

      SUBROUTINE fteik_h5io_ijkv2levelsF(nz, nx, ny, ngrd, nLevels, &
                                         levelPtr, ijkv, levels) &
                 BIND(C, NAME='fteik_h5io_ijkv2levelsF')
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: ngrd, nLevels , nz, nx, ny
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
            levels(indx) = level
    2    CONTINUE
    1 CONTINUE
      RETURN
      END

      INTEGER(C_INT) FUNCTION fteik_h5io_writeLevelSchedulesF( ) &
                     RESULT(ierr) BIND(C, NAME='fteik_h5io_writeLevelSchedulesF')
      USE FTEIK_UTILS64F, ONLY : ijkv1, ijkv2, ijkv3, ijkv4,    &
                                 ijkv5, ijkv6, ijkv7, ijkv8,    &
                                 levelPtr, nLevels, nz, nx, ny, &
                                 ngrd, lhaveGrid
      USE FTEIK_H5IO, ONLY : fteik_h5io_ijkv2levelsF, fteik_h5io_writeLevelScheduleF, &
                             fteik_h5io_openFileReadWriteF, fteik_h5io_closeFileF
      USE FTEIK_H5IO, ONLY : fname, linitH5FL
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
      IF (.NOT.lhaveGrid) THEN
         WRITE(*,*) 'fteik_h5io_writeLevelSchedulesF: Grid not set'
         ierr =-1 
         RETURN
      ENDIF
      ! Open file 
      h5fl = fteik_h5io_openFileReadWriteF(fname)
      ALLOCATE(levels(ngrd))
      DO 1 sweep=1,8
         IF (sweep == 1) THEN 
            CALL fteik_h5io_ijkv2levelsF(nz, nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv1, levels)
         ELSEIF (sweep == 2) THEN
            CALL fteik_h5io_ijkv2levelsF(nz, nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv2, levels)
         ELSEIF (sweep == 3) THEN
            CALL fteik_h5io_ijkv2levelsF(nz, nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv3, levels)
         ELSEIF (sweep == 4) THEN
            CALL fteik_h5io_ijkv2levelsF(nz, nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv4, levels)
         ELSEIF (sweep == 5) THEN 
            CALL fteik_h5io_ijkv2levelsF(nz, nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv5, levels)
         ELSEIF (sweep == 6) THEN
            CALL fteik_h5io_ijkv2levelsF(nz, nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv6, levels)
         ELSEIF (sweep == 7) THEN
            CALL fteik_h5io_ijkv2levelsF(nz, nx, ny, ngrd, nLevels,  &
                                         levelPtr, ijkv7, levels)
         ELSEIF (sweep == 8) THEN
            CALL fteik_h5io_ijkv2levelsF(nz, nx, ny, ngrd, nLevels,  &
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

      INTEGER(C_INT) FUNCTION fteik_h5io_writeTravelTimesF(dataName) &
                     RESULT(ierr) BIND(C, NAME='fteik_h5io_writeTravelTimesF')
      USE FTEIK_UTILS64F, ONLY : ttimes, nzx, nz, nx, ny, ngrd, &
                                 lhaveGrid, lhaveTravelTimes
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      USE FTEIK_H5IO, ONLY : fteik_h5io_writeTravelTimes32fF, &
                             fteik_h5io_openFileReadWriteF, fteik_h5io_closeFileF
      USE FTEIK_H5IO, ONLY : fname, linitH5FL
      USE ISO_C_BINDING
      IMPLICIT NONE
      CHARACTER(C_CHAR), INTENT(IN) :: dataName(*)
      INTEGER(C_INT64_T) h5fl
      REAL(C_FLOAT), ALLOCATABLE :: ttout(:)
      INTEGER(C_INT) ierr2, indx, ix, iy, iz, jndx
      ierr = 0
      IF (.NOT.linitH5FL) THEN
         WRITE(*,*) 'fteik_h5io_writeTravelTimesF: File not initialized'
         ierr =-1 
         RETURN
      ENDIF
      IF (.NOT.lhaveTravelTimes) THEN
         WRITE(*,*) 'fteik_h5io_writeTravelTimesF: Travel times not yet computed'
         ierr =-1
         RETURN
      ENDIF
      IF (.NOT.lhaveGrid) THEN
         WRITE(*,*) 'fteik_h5io_writeTravelTimesF: Grid not set'
         ierr =-1
         RETURN
      ENDIF
      ALLOCATE(ttout(ngrd))
      ttout(:) = 0.0
      DO 1 iz=1,nz
         DO 2 iy=1,ny
            DO 3 ix=1,nx
               !jz = nz + 1 - iz ! flip coordinate syste
               indx = grid2indexF(iz, ix, iy, nz, nzx)
               jndx = (iz - 1)*nx*ny + (iy - 1)*nx + ix
               !print *, jndx
               ttout(jndx) = (ttimes(indx))
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
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

      INTEGER(C_INT) FUNCTION fteik_h5io_writeVelocityModelF(dataName) &
                     RESULT(ierr) BIND(C, NAME='fteik_h5io_writeVelocityModelF')
      USE FTEIK_UTILS64F, ONLY : slow, ncell, nx, ny, nz, nzm1, nzm1_nxm1, &
                                 lhaveSlownessModel, lhaveGrid !, dx, dy, dz
      USE FTEIK_UTILS64F, ONLY : velGrid2indexF
      USE FTEIK_H5IO, ONLY : fteik_h5io_writeVelocityModel16iF, &
                             fteik_h5io_openFileReadWriteF, fteik_h5io_closeFileF
      USE FTEIK_H5IO, ONLY : fname, linitH5FL
      USE ISO_C_BINDING
      IMPLICIT NONE
      CHARACTER(C_CHAR), INTENT(IN) :: dataName(*)
      INTEGER(C_INT64_T) h5fl
      INTEGER(C_INT16_T), ALLOCATABLE :: vout(:)
      INTEGER ierr2, indx, ix, iy, iz, jndx!, jz
      ierr = 0
      IF (.NOT.linitH5FL) THEN
         WRITE(*,*) 'fteik_h5io_writeVelocityModelF: File not initialized'
         ierr =-1
         RETURN
      ENDIF
      IF (.NOT.lhaveSlownessModel) THEN
         WRITE(*,*) 'fteik_h5io_writeVelocityModelF: Velocity model not set'
         ierr =-1
         RETURN
      ENDIF
      IF (.NOT.lhaveGrid) THEN
         WRITE(*,*) 'fteik_h5io_writeVelocityModelF: Grid not set'
         ierr =-1
         RETURN
      ENDIF
      ! First dimension in topology is slowest changing (z)
      ALLOCATE(vout(ncell))
      vout(:) = 0
      DO 1 iz=1,nz-1
         DO 2 iy=1,ny-1
            DO 3 ix=1,nx-1
               !jz = iz !nz-iz ! flip coordinate syste
               indx = velGrid2indexF(iz, ix, iy, nzm1, nzm1_nxm1) 
               jndx = (iz-1)*(nx-1)*(ny-1) + (iy-1)*(nx-1) + ix
               vout(jndx) = INT(1.d0/slow(indx)  + 0.5d0)
               !print *, vout(jndx), 1.d0/slow(indx)
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
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
