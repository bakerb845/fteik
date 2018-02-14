MODULE FTEIK_RECEIVER64F
   USE FTEIK_MODEL64F, ONLY : fteik_model_grid2indexF
   USE ISO_C_BINDING
   !> Receiver index in z.  This is a vector of dimension [nrec]. 
   INTEGER(C_INT), PROTECTED, ALLOCATABLE, SAVE :: zri(:)
   !> Receiver index in x.  This is a vector of dimension [nrec].
   INTEGER(C_INT), PROTECTED, ALLOCATABLE, SAVE :: xri(:)
   !> Receiver index in y.  This is a vector of dimension [nrec]. 
   INTEGER(C_INT), PROTECTED, ALLOCATABLE, SAVE :: yri(:)
   !> z distance from DBLE(zri).  This is for the linear interpolation.
   REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, SAVE :: zdr(:)
   !> x distance from DBLE(xri).  This is for the linear interpolation.
   REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, SAVE :: xdr(:)
   !> y distance from DBLE(yri) . This is for the linear interpolation.
   REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, SAVE :: ydr(:)
   !> Flag indicating whether or not the module has been initialized.
   LOGICAL(C_BOOL), PRIVATE, SAVE :: linit = .FALSE. 
   !> Number of receivers.
   INTEGER(C_INT), PROTECTED, SAVE :: nrec = 0
   !> Controls verbosity.
   INTEGER(C_INT), PROTECTED, SAVE :: verbose = 0 
   INTEGER(C_INT), PARAMETER :: INTERP_NEAREST = 0      !< Nearest neighbor interpolation.
   INTEGER(C_INT), PARAMETER :: INTERP_LINEAR = 1       !< Linear interpolation.
   INTEGER(C_INT), PARAMETER :: INTERP_HIGHACCURACY = 2 !< Not programmed.
   ! Subroutine availability
   PUBLIC :: fteik_receiver_initialize64f
   PUBLIC :: fteik_receiver_getNumberOfReceivers
   PUBLIC :: fteik_receiver_getTravelTimes64f
   PUBLIC :: fteik_receiver_setVerobosity
   PUBLIC :: fteik_receiver_free
   PRIVATE :: fteik_model_grid2indexF
   CONTAINS
!----------------------------------------------------------------------------------------!
!                                      Begin the Code                                    !
!----------------------------------------------------------------------------------------!
!>
!>    @brief Sets the receiver indices in the grid.  This will make for a quick 
!>           linear interpolation.
!>           https://github.com/bottero/IMCMCrun/blob/master/src/functions.cpp
!>
!>    @param[in] nrecIn     Number of receivers.  If this is -1 then this function will
!>                          simply deallocate/reset the receiver information.
!>    @param[in] z          Receiver positions in z (meters).
!>    @param[in] x          Receiver positions in x (meters).
!>    @param[in] y          Receiver positions in y (meters).
!>    @param[in] verboseIn  Controls the verbosity.  < 1 is quiet.
!>
!>    @param[out] ierr      0 indicate success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_receiver_initialize64f(nrecIn, z, x, y, verboseIn, ierr) &
                 BIND(C, NAME='fteik_receiver_initialize64f')
      USE FTEIK_MODEL64F, ONLY : lis3dModel, dx, dy, dz, nx, ny, nz, x0, y0, z0
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: nrecIn, verboseIn
      REAL(C_DOUBLE), INTENT(IN) :: z(nrecIn), x(nrecIn), y(nrecIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) xr, yr, zr
      INTEGER(C_INT) irec, ix, iy, iz
      ! Initialize 
      ierr = 0
      CALL fteik_receiver_free()
      CALL fteik_receiver_setVerobosity(verboseIn)
      ! Check that the model has been initialized
      IF (nx < 1 .OR. ny < 1 .OR. nz < 1) THEN
         WRITE(*,*) 'fteik_receiver_initialize64f: Grid not yet set'
         ierr = 1
         RETURN
      ENDIF
      ! Verify the inputs
      IF (nrecIn < 1) THEN
         WRITE(*,*) 'fteik_receiver_initialize64f: No receivers'
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
            WRITE(*,*) 'fteik_receiver_initialize64f: z position is out of bounds', &
                       z(irec), z0, z0 + dz*DBLE(nz - 1)
            ierr = ierr + 1
         ENDIF
         IF (x(irec) < x0 .OR. x(irec) > x0 + dx*DBLE(nx - 1)) THEN
            WRITE(*,*) 'fteik_receiver_initialize64f: x position is out of bounds', &
                       x(irec), x0, x0 + dx*DBLE(nx - 1)
            ierr = ierr + 1
         ENDIF
         IF (y(irec) < y0 .OR. y(irec) > y0 + dy*DBLE(ny - 1)) THEN
            WRITE(*,*) 'fteik_receiver_initialize64f: y position is out of bounds', &
                       y(irec), y0, y0 + dy*DBLE(ny - 1)
            ierr = ierr + 1
         ENDIF
         ! Get the grid point position
         zr = (z(irec) - z0)/dz
         xr = (x(irec) - x0)/dx
         yr = (y(irec) - y0)/dy
         iz = INT(zr + 1.d0) ! convert to Fortran grid point
         ix = INT(xr + 1.d0)
         iy = INT(yr + 1.d0)
         IF (.NOT.lis3dModel) iy = 1
         ! Interpolation goes from (iz,ix,iy) to (iz+1,ix+1,iy+1).  This is an issue when
         ! the receiver is at the edge.  In this case; shift it to model bounds - 1.
         ! Additionally, being extra generous, I'll force the receiver into the model
         ! if the receiver is out of bounds.
         zri(irec) = MAX(1, MIN(nz-1, iz))
         xri(irec) = MAX(1, MIN(nx-1, ix))
         yri(irec) = MAX(1, MIN(ny-1, iy))
         zdr(irec) = (z(irec) - (z0 + dz*DBLE(zri(irec)-1)))!/dz
         xdr(irec) = (x(irec) - (x0 + dx*DBLE(xri(irec)-1)))!/dx
         ydr(irec) = (y(irec) - (y0 + dy*DBLE(yri(irec)-1)))!/dy
         IF (.NOT.lis3dModel) ydr(irec) = 0.d0
         IF (verbose > 0) THEN
            IF (lis3dModel) THEN
               WRITE(*,900) z(irec), x(irec), y(irec)
               WRITE(*,901) iz, ix, iy               
            ELSE
               WRITE(*,910) z(irec), x(irec)
               WRITE(*,911) iz, ix, iy
            ENDIF
         ENDIF
    1 CONTINUE
      linit = .TRUE.
  900 FORMAT(' fteik_receiver_initialize64f: Original receiver coordinates (z,x,y)=', &
             3F12.2, ' (m)')
  901 FORMAT(' fteik_receiver_initialize64f: Interpolation grid point (iz0,ix0,iy0)=', &
             3I3)
  910 FORMAT(' fteik_receiver_initialize64f: Original receiver coordinates (z,x)=', &
             2F12.2, ' (m)')
  911 FORMAT(' fteik_receiver_initialize64f: Interpolation grid point (iz0,ix0)=', &
             2I2)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the verbosity on the module.
!>
!>    @param[in] verboseIn   Verbosity level to set.  Less than 1 is quiet.
!>
      SUBROUTINE fteik_receiver_setVerobosity(verboseIn) &
      BIND(C, NAME='fteik_receiver_setVerbosity')
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: verboseIn
      verbose = verboseIn
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility function for determining the number of receivers set in the model.
!>
!>    @param[out] nrecOut  Number of receivers set in module.
!>    @param[out] ierr     0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_receiver_getNumberOfReceivers(nrecOut, ierr)  &
      BIND(C, NAME='fteik_receiver_getNumberOfReceivers')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: nrecOut, ierr
      ierr = 0
      nrecOut = 0
      IF (.NOT.linit) THEN
         WRITE(*,*) 'fteik_receiver_getNumberOfReceivers: Module not initialized'
         ierr = 1
         RETURN
      ENDIF
      nrecOut = nrec
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gets the travel times at each receiver via linear interpolation
!>           (e.g., https://en.wikipedia.org/wiki/Trilinear_interpolation)
!>
!>    @param[in] nrecIn    Workspace of ttr.
!>    @param[in] ngrd      Number of grid points in travel-time field.
!>    @param[in] ttimes    Travel-time field (seconds).  This has dimension [ngrd].
!>
!>    @param[out] ttr      Travel-times (seconds) at each receiver.  This is a vector of
!>                         dimension [nrecIn] however only the first nrec values are 
!>                         interpolated.  Excess values will be set to FTEIK_HUGE.
!>    @param[out] ierr     0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_receiver_getTravelTimes64f(nrecIn, ngrd,       &
                                                  ttimes, ttr, ierr)  &
      BIND(C, NAME='fteik_receiver_getTravelTimes64f')
      USE FTEIK_MODEL64F, ONLY : lis3dModel, dz, dx, dy, nz, nzx, dx, dy, dz
      USE FTEIK_CONSTANTS64F, ONLY : one, FTEIK_HUGE
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nrecIn, ngrd
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(ngrd)
      REAL(C_DOUBLE), INTENT(OUT) :: ttr(nrecIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) c0, c1, c00, c01, c10, c11, one_m_zd, one_m_xd, &
                     tt000, tt100, tt010, tt110, tt001, tt101, tt011, tt111, &
                     tt00, tt10, tt01, tt11, xd, yd, zd 
      INTEGER(C_INT) i, i0, i1, i2, i3, i4, i5, i6, i7
      ierr = 0
      IF (.NOT.linit) THEN
         WRITE(*,*) 'fteik_receiver_getTravelTimes64fF: Module not initialized'
         ierr = 1
         RETURN
      ENDIF
      IF (nrecIn < 1 .AND. verbose > 0) THEN
         WRITE(*,*) 'fteik_receiver_getTravelTimes64fF: No receivers'
         RETURN
      ENDIF
      IF (nrecIn /= nrec .AND. verbose > 0) THEN
         WRITE(*,*) 'fteik_receiver_getTravelTimes64fF: Warning nrecIn /= nrec', &
                    nrecIn, nrec
         ttr(1:nrecIn) = FTEIK_HUGE
      ENDIF
      IF (lis3dModel) THEN
         DO 1 i=1,MIN(nrec, nrecIn)
            ! extract travel-times in cache-friendly way
            i0 = fteik_model_grid2indexF(zri(i),   xri(i),   yri(i),   nz, nzx)
            i1 = fteik_model_grid2indexF(zri(i)+1, xri(i),   yri(i),   nz, nzx)
            i2 = fteik_model_grid2indexF(zri(i),   xri(i)+1, yri(i),   nz, nzx)
            i3 = fteik_model_grid2indexF(zri(i)+1, xri(i)+1, yri(i),   nz, nzx)
            i4 = fteik_model_grid2indexF(zri(i),   xri(i),   yri(i)+1, nz, nzx)
            i5 = fteik_model_grid2indexF(zri(i)+1, xri(i),   yri(i)+1, nz, nzx)
            i6 = fteik_model_grid2indexF(zri(i),   xri(i)+1, yri(i)+1, nz, nzx)
            i7 = fteik_model_grid2indexF(zri(i)+1, xri(i)+1, yri(i)+1, nz, nzx)
            tt000 = ttimes(i0)
            tt100 = ttimes(i1)
            tt010 = ttimes(i2)
            tt110 = ttimes(i3)
            tt001 = ttimes(i4)
            tt101 = ttimes(i5)
            tt011 = ttimes(i6)
            tt111 = ttimes(i7)
            ! perform trilinear interpolation
            zd = zdr(i)/dz ! zd = (z - z0)/(z1 - z0) means this is in [0,1]
            xd = xdr(i)/dx ! xd = (x - x0)/(x1 - x0) means this is in [0,1]
            yd = ydr(i)/dy ! yd = (y - y0)/(y1 - y0) means this is in [0,1]
            ! step 1 - interpolate in first direction (z)
            one_m_zd = one - zd
            c00 = tt000*one_m_zd + tt100*zd
            c01 = tt001*one_m_zd + tt101*zd
            c10 = tt010*one_m_zd + tt110*zd
            c11 = tt011*one_m_zd + tt111*zd
            ! step 2 - interpolate in second direction (x)
            one_m_xd = one - xd
            c0 = c00*one_m_xd + c10*xd
            c1 = c01*one_m_xd + c11*xd
            ! step 3 - interpolate in third direction (y)
            ttr(i) = c0*(one - yd) + c1*yd 
!if (ttr(i) < min(tt000, tt100, tt010, tt110, tt001, tt101, tt011, tt111) .or. &
!    ttr(i) > max(tt000, tt100, tt010, tt110, tt001, tt101, tt011, tt111)) then
!print *, 'weird'
!endif
    1     CONTINUE
      ELSE
         DO 2 i=1,MIN(nrec, nrecIn)
            ! extract travel-times in cache-friendly way
            i0 = fteik_model_grid2indexF(zri(i),   xri(i),   yri(i), nz, 0)
            i1 = fteik_model_grid2indexF(zri(i)+1, xri(i),   yri(i), nz, 0)
            i2 = fteik_model_grid2indexF(zri(i),   xri(i)+1, yri(i), nz, 0)
            i3 = fteik_model_grid2indexF(zri(i)+1, xri(i)+1, yri(i), nz, 0)
            tt00 = ttimes(i0)
            tt10 = ttimes(i1)
            tt01 = ttimes(i2)
            tt11 = ttimes(i3)
            zd = zdr(i)/dz      ! zd = (z - z0)/(z1 - z0) means this is in [0,1]
            xd = xdr(i)/dx      ! xd = (x - x0)/(x1 - x0) means this is in [0,1]
            one_m_zd = one - zd
            one_m_xd = one - xd
            ! perform bilinear interpolation
            ttr(i) = tt00*one_m_zd*one_m_xd &
                   + tt10*zd      *one_m_xd &
                   + tt01*one_m_zd*xd       &
                   + tt11*zd      *xd
    2    CONTINUE
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases memory in the receiver module and sets all variables to 0.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_receiver_free( ) &
                 BIND(C, NAME='fteik_receiver_free')
      USE ISO_C_BINDING
      IMPLICIT NONE
      nrec = 0 
      verbose = 0
      IF (ALLOCATED(zri)) DEALLOCATE(zri)
      IF (ALLOCATED(xri)) DEALLOCATE(xri)
      IF (ALLOCATED(yri)) DEALLOCATE(yri)
      IF (ALLOCATED(zdr)) DEALLOCATE(zdr)
      IF (ALLOCATED(xdr)) DEALLOCATE(xdr)
      IF (ALLOCATED(ydr)) DEALLOCATE(ydr)
      linit = .FALSE.
      RETURN
      END
!----------------------------------------------------------------------------------------!
!                                    Begin the MPI                                       !
!----------------------------------------------------------------------------------------!
#if defined(FTEIK_FORTRAN_USE_MPI) 
      SUBROUTINE fteik_receiver_broadcastF(root, comm, mpierr) &
                 BIND(C, NAME='fteik_receiver_broadcastF')
      USE MPI
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: comm, root 
      INTEGER(C_INT), INTENT(OUT) :: mpierr
      INTEGER(C_INT) myid
      CALL MPI_Comm_rank(comm, myid, mpierr)
      CALL MPI_Bcast(nrec,    1, MPI_INT, root, comm, mpierr) 
      IF (nrec < 1) RETURN 
      CALL MPI_Bcast(verbose, 1, MPI_INT, root, comm, mpierr)
      IF (myid /= root) THEN 
         CALL fteik_receiver_free()
         ALLOCATE(zri(nrec))
         ALLOCATE(xri(nrec))
         ALLOCATE(yri(nrec))
         ALLOCATE(zdr(nrec))
         ALLOCATE(xdr(nrec))
         ALLOCATE(ydr(nrec))
      ENDIF
      CALL MPI_Bcast(zri, nrec, MPI_INT, root, comm, mpierr)
      CALL MPI_Bcast(xri, nrec, MPI_INT, root, comm, mpierr)
      CALL MPI_Bcast(yri, nrec, MPI_INT, root, comm, mpierr)
      CALL MPI_Bcast(zdr, nrec, MPI_DOUBLE, root, comm, mpierr)
      CALL MPI_Bcast(xdr, nrec, MPI_DOUBLE, root, comm, mpierr)
      CALL MPI_Bcast(ydr, nrec, MPI_DOUBLE, root, comm, mpierr) 
      END SUBROUTINE
#endif
!----------------------------------------------------------------------------------------!
!                                     End the Code                                       !
!----------------------------------------------------------------------------------------!
END MODULE
