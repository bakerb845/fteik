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
      SUBROUTINE fteik_receiver_initialize64fF(nrecIn, z, x, y, ierr) &
                 BIND(C, NAME='fteik_receiver_initialize64fF')
      USE FTEIK_RECEIVER64F, ONLY : fteik_receiver_finalizeF
      USE FTEIK_RECEIVER64F, ONLY : nrec, xri, yri, zri, xdr, ydr, zdr
      USE FTEIK_MODEL64F, ONLY : dx, dy, dz, nx, ny, nz, x0, y0, z0
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: nrecIn
      REAL(C_DOUBLE), INTENT(IN) :: z(nrecIn), x(nrecIn), y(nrecIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) xr, yr, zr
      INTEGER(C_INT) irec, ix, iy, iz
      ! Initialize 
      ierr = 0
      CALL fteik_receiver_finalizeF() 
      !nrec = 0
      !IF (ALLOCATED(zri)) DEALLOCATE(zri)
      !IF (ALLOCATED(xri)) DEALLOCATE(xri)
      !IF (ALLOCATED(yri)) DEALLOCATE(yri)
      !IF (ALLOCATED(zdr)) DEALLOCATE(zdr)
      !IF (ALLOCATED(xdr)) DEALLOCATE(xdr)
      !IF (ALLOCATED(ydr)) DEALLOCATE(ydr)
      ! Check that the model has been initialized
      IF (nx < 1 .OR. ny < 1 .OR. nz < 1) THEN
         WRITE(*,*) 'fteik_receiver_initialize64fF: Grid not yet set'
         ierr = 1
         RETURN
      ENDIF
      ! Verify the inputs
      IF (nrecIn < 1) THEN
         WRITE(*,*) 'fteik_receiver_initialize64fF: No receivers'
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
            WRITE(*,*) 'fteik_receiver_initialize64fF: z position is out of bounds', &
                       z(irec)
            ierr = ierr + 1
         ENDIF
         IF (x(irec) < x0 .OR. x(irec) > x0 + dx*DBLE(nx - 1)) THEN
            WRITE(*,*) 'fteik_receiver_initialize64fF: x position is out of bounds', &
                       x(irec)
            ierr = ierr + 1
         ENDIF
         IF (y(irec) < y0 .OR. y(irec) > y0 + dy*DBLE(ny - 1)) THEN
            WRITE(*,*) 'fteik_receiver_initialize64fF: y position is out of bounds', &
                       y(irec)
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
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility function for determining the number of receivers set in the model.
!>
!>    @result Number of receivers set in module.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      INTEGER(C_INT) FUNCTION fteik_receiver_getNumberOfReceiversF( )     &
      RESULT(nrecOut) BIND(C, NAME='fteik_receiver_getNumberOfReceiversF')
      USE FTEIK_RECEIVER64F, ONLY : nrec
      USE ISO_C_BINDING
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
!>
!>    @param[out] ttr      Travel-times (seconds) at each receiver.  This is a vector of
!>                         dimension [nrecIn] however only the first nrec values are 
!>                         interpolated.  Excess values will be set to FTEIK_HUGE.
!>    @param[out] ierr     0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_receiver_getTravelTimes64fF(nrecIn, ttr, ierr) &
      BIND(C, NAME='fteik_receiver_getTravelTimes64fF')
      USE FTEIK_UTILS64F, ONLY : grid2indexF
      USE FTEIK_UTILS64F, ONLY : ttimes, lhaveTravelTimes
      USE FTEIK_UTILS64F, ONLY : one, FTEIK_HUGE
      USE FTEIK_RECEIVER64F, ONLY : nrec, xdr, ydr, zdr, xri, yri, zri 
      USE FTEIK_MODEL64F, ONLY : dz, dx, dy, nz, nzx, dx, dy, dz
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nrecIn
      REAL(C_DOUBLE), INTENT(OUT) :: ttr(nrecIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) c0, c1, c00, c01, c10, c11, one_m_zd, one_m_xd, &
                     tt000, tt100, tt010, tt110, tt001, tt101, tt011, tt111, xd, yd, zd 
      INTEGER(C_INT) i, i0, i1, i2, i3, i4, i5, i6, i7
      ierr = 0
      IF (nrec < 1) THEN
         WRITE(*,*) 'fteik_receiver_getTravelTimes64fF: No receivers'
         ierr = 1
         RETURN
      ENDIF
      IF (nrecIn < nrec) THEN
         WRITE(*,*) 'fteik_receiver_getTravelTimes64fF: Error nrecIn < nrec'
         ierr = 1
         RETURN
      ENDIF
      IF (nrecIn > nrec) ttr(nrec+1:nrecIn) = FTEIK_HUGE
      IF (.NOT.lhaveTravelTimes) THEN
         WRITE(*,*) 'fteik_receiver_getTravelTimes64fF: Travel-times not yet computed'
         ierr = 1
         RETURN
      ENDIF
      DO 1 i=1,nrec
         ! extract travel-times in cache-friendly way
         i0 = grid2indexF(zri(i),   xri(i),   yri(i),   nz, nzx)
         i1 = grid2indexF(zri(i)+1, xri(i),   yri(i),   nz, nzx)
         i2 = grid2indexF(zri(i),   xri(i+1), yri(i),   nz, nzx)
         i3 = grid2indexF(zri(i)+1, xri(i+1), yri(i),   nz, nzx)
         i4 = grid2indexF(zri(i),   xri(i),   yri(i)+1, nz, nzx)
         i5 = grid2indexF(zri(i)+1, xri(i),   yri(i)+1, nz, nzx)
         i6 = grid2indexF(zri(i),   xri(i+1), yri(i)+1, nz, nzx)
         i7 = grid2indexF(zri(i)+1, xri(i+1), yri(i)+1, nz, nzx)
         tt000 = ttimes(i0)
         tt100 = ttimes(i1)
         tt010 = ttimes(i2)
         tt110 = ttimes(i3)
         tt001 = ttimes(i4)
         tt101 = ttimes(i5)
         tt011 = ttimes(i6)
         tt111 = ttimes(i7)
         ! perform trilinear interpolation
         zd = zdr(i)/dz
         xd = xdr(i)/dx
         yd = ydr(i)/dy
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
    1 CONTINUE
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
      SUBROUTINE fteik_receiver_finalizeF( ) &
                 BIND(C, NAME='fteik_receiver_finalizeF')
      USE FTEIK_RECEIVER64F, ONLY : nrec, xri, yri, zri, xdr, ydr, zdr
      USE ISO_C_BINDING
      IMPLICIT NONE
      nrec = 0 
      IF (ALLOCATED(zri)) DEALLOCATE(zri)
      IF (ALLOCATED(xri)) DEALLOCATE(xri)
      IF (ALLOCATED(yri)) DEALLOCATE(yri)
      IF (ALLOCATED(zdr)) DEALLOCATE(zdr)
      IF (ALLOCATED(xdr)) DEALLOCATE(xdr)
      IF (ALLOCATED(ydr)) DEALLOCATE(ydr)
      RETURN
      END
