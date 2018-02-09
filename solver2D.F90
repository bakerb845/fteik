MODULE FTEIK2D_SOLVER64F
  USE FTEIK_CONSTANTS64F, ONLY : zero
  USE ISO_C_BINDING
  IMPLICIT NONE
  !> Holds the travel-times (seconds).  This has dimension [ngrd].
  REAL(C_DOUBLE), PROTECTED, ALLOCATABLE, SAVE :: ttimes(:)
  !!!!DIR$ ATTRIBUTES ALIGN: 64 :: ttimes
  !> Defines the number of levels.  For 2D that is nz + nx
  INTEGER(C_INT), PROTECTED, SAVE :: nLevels = 0
  !> Defines the max level size.
  INTEGER(C_INT), PROTECTED, SAVE :: maxLevelSize = 0
  !> Controls the verbosity
  INTEGER(C_INT), PROTECTED, SAVE :: verbose = 0
  !> Defines the transition from the spherical to the Cartesian solver solver
  !> during the initialization phase.  This has units of grid points.
  REAL(C_DOUBLE), PROTECTED, SAVE :: epsS2C = zero
  !> Defines the number of Gauss-Seidel iterations.
  INTEGER(C_INT), PROTECTED, SAVE :: nsweep = 0
  !> Flag indicating whether or not the travel times were computed.
  LOGICAL(C_BOOL), PROTECTED, SAVE :: lhaveTimes = .FALSE.
  !> 
  LOGICAL(C_BOOL), PARAMETER :: lis3d = .FALSE.
  !> Private variables for the local solver
  INTEGER(C_INT), PARAMETER :: alignment = 64
  REAL(C_DOUBLE), SAVE, PRIVATE :: zsa, xsa, szero, szero2, eps
  INTEGER(C_INT), SAVE, PRIVATE :: zsi, xsi
  REAL(C_DOUBLE), SAVE, PRIVATE :: dx, dxi, dz, dzi, dz2i, dx2i, dz2i_dx2i, &
                                   dz2i_p_dx2i, dz2i_p_dx2i_inv
  REAL(C_DOUBLE), SAVE, PRIVATE, ALLOCATABLE :: slocWork(:), ttvecWork(:), tupdWork(:)
  LOGICAL(C_BOOL), SAVE, PRIVATE, ALLOCATABLE :: lupdWork(:)
  ! Label the subroutines/functions
  PUBLIC :: fteik_solver2d_setVelocityModel64f
  PUBLIC :: fteik_solver2d_initialize64f
  PUBLIC :: fteik_solver2d_free
  PUBLIC :: fteik_solver2d_solveSourceLSM
  PUBLIC :: fteik_solver2d_solveSourceFSM
  PUBLIC :: fteik_solver2d_getTravelTimeField64f
  PUBLIC :: fteik_solver2d_getTravelTimes64f
  PUBLIC :: fteik_solver2d_setNumberOfSweeps
  PUBLIC :: fteik_solver2d_setReceivers64f
  PUBLIC :: fteik_solver2d_setSources64f
  PUBLIC :: fteik_solver2d_setSphereToCartEpsilon
  PUBLIC :: fteik_solver2d_getNumberOfReceivers
  PRIVATE :: fteik_solver2d_isUpdateNodeF
  PRIVATE :: fteik_solver2d_prefetchTravelTimes64fF
  PRIVATE :: fteik_solver2d_prefetchSlowness64fF
  PRIVATE :: fteik_solver2d_updateTravelTimes64fF
  PRIVATE :: fteik_localSolver2d_noInit64fF
  PRIVATE :: fteik_localSolver2d_init64fF
  PRIVATE :: fteik_localSolver2d_tAna64fF
  PRIVATE :: fteik_localSolver2d_tAnaD64fF
  PRIVATE :: fteik_localSolver2d_initialize64f
  CONTAINS
  !--------------------------------------------------------------------------------------!
  !                                      Begin the Code                                  !
  !--------------------------------------------------------------------------------------!
      SUBROUTINE fteik_solver2d_initialize64f(nzIn, nxIn,       &
                                              z0In, x0In,       &
                                              dzIn, dxIn,       &   
                                              nsweepIn, epsIn,  &
                                              verboseIn, ierr)  &
                 BIND(C, NAME='fteik_solver2d_initialize64f')
      USE ISO_C_BINDING
      USE FTEIK_MEMORY, ONLY : padLength64F
      USE FTEIK_MODEL64F, ONLY : fteik_model_intializeGeometry, ngrd
      USE FTEIK_CONSTANTS64F, ONLY : one, zero, FALSE
      IMPLICIT NONE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x0In, z0In
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dxIn, dzIn
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: epsIn
      INTEGER(C_INT), VALUE, INTENT(IN) :: nxIn, nzIn
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsweepIn, verboseIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      IF (verboseIn > 0) WRITE(*,*) 'fteik_solver_initialize64f: Initializing...'
      CALL fteik_solver2d_free()
      CALL fteik_solver2d_setVerobosity(verboseIn)
      CALL fteik_model_intializeGeometry(lis3d,            &
                                         nzIn, nxIn, 1,    &
                                         dzIn, dxIn, one,  &
                                         z0In, x0In, zero, &
                                         ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver2d_initialize64f: Error setting grid geometry'
         RETURN
      ENDIF
      CALL fteik_solver2d_setNumberOfSweeps(nsweepIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver2d_initialize64f: Failed to set max number of sweeps'
         RETURN
      ENDIF
      CALL fteik_solver2d_setSphereToCartEpsilon(epsIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver2d_initialize64f: Failed to set epsilon'
         RETURN
      ENDIF
      ! Compute the number of levels
      nLevels = nxIn + nzIn
      maxLevelSize = MIN(nxIn, nzIn) + padLength64F(alignment, MIN(nxIn, nzIn))
      ! Set space for the travel-times
      IF (.NOT.ALLOCATED(ttimes)) ALLOCATE(ttimes(ngrd))
      ALLOCATE(lupdWork(maxLevelSize))
      ALLOCATE(slocWork(4*maxLevelSize))
      ALLOCATE(ttvecWork(4*maxLevelSize))
      ALLOCATE(tupdWork(maxLevelSize))
      lupdWork(:) = FALSE
      slocWork(:) = zero
      ttvecWork(:) = zero
      tupdWork(:) = zero
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the verbosity on the module.
!>
!>    @param[in] verboseIn   Verbosity level to set.  Less than 1 is quiet.
!>
      SUBROUTINE fteik_solver2d_setVerobosity(verboseIn) &
      BIND(C, NAME='fteik_solver2d_setVerbosity')
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: verboseIn
      verbose = verboseIn
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases all memory associated with the solver and resets all variables.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solver2d_free()   &
      BIND(C, NAME='fteik_solver2d_free')
      USE FTEIK_RECEIVER64F, ONLY : fteik_receiver_free
      USE FTEIK_SOURCE64F, ONLY : fteik_source_free
      USE FTEIK_MODEL64F, ONLY : fteik_model_free
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      lhaveTimes = .FALSE.
      epsS2C = zero
      nsweep = 0
      nLevels = 0
      verbose = 0
      IF (ALLOCATED(ttimes))     DEALLOCATE(ttimes)
      IF (ALLOCATED(lupdWork))   DEALLOCATE(lupdWork)
      IF (ALLOCATED(slocWork))   DEALLOCATE(slocWork)
      IF (ALLOCATED(ttvecWork))  DEALLOCATE(ttvecWork)
      IF (ALLOCATED(tupdWork))   DEALLOCATE(tupdWork)
      CALL fteik_receiver_free()
      CALL fteik_source_free()
      CALL fteik_model_free()
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility routine to return the number of receivers.
!>
!>    @param[out] nrec   Number of receivers.
!>    @param[out] ierr   0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver2d_getNumberOfReceivers(nrec, ierr) &
      BIND(C, NAME='fteik_solver2d_getNumberOfReceivers')
      USE FTEIK_RECEIVER64F, ONLY : fteik_receiver_getNumberOfReceivers
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(OUT) :: nrec, ierr
      ierr = 0 
      CALL fteik_receiver_getNumberOfReceivers(nrec, ierr)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Determines if this is an update node in the level-set fast sweeping method.
!>
!>    @param[in] linit      If true then this is an initialization step and we can
!>                          use a subset of the grid based on the source position. \n
!>                          Otherwise, this is general sweeping during the Gauss-Seidel
!>                          iteration. 
!>    @param[in] level      Current level in level-set.
!>    @param[in] sweep      Sweep number.  This is in the range [1,4].
!>    @param[in] nz         Number of z grid points in travel time field.
!>    @param[in] nx         Number of x grid points in travel time field.
!>    @param[in] i1         Start index of level set.
!>    @param[in] i2         End indx of level set.
!>    @param[in] zsi        Source index in z.  This is only used if linit is true. 
!>    @param[in] xsi        Source index in x.  This is only used if linit is true.
!>
!>    @param[out] lupd      If true then the k'th node in the level is to be updated.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      PURE SUBROUTINE fteik_solver2d_isUpdateNodeF(linit, level, sweep, nz, nx, &
                                                   i1, i2, zsi, xsi, lupd)
      USE FTEIK_CONSTANTS64F, ONLY : TRUE, FALSE
      USE ISO_C_BINDING
      IMPLICIT NONE
      LOGICAL(C_BOOL), VALUE, INTENT(IN) :: linit
      INTEGER(C_INT), VALUE, INTENT(IN) :: i1, i2, level, nx, nz, sweep, xsi, zsi
      LOGICAL(C_BOOL), DIMENSION(:), INTENT(OUT) :: lupd
      INTEGER(C_INT) i, ix, iz
      IF (.NOT.linit) THEN
         IF (sweep == 1) THEN
            !$OMP SIMD
            DO i=i1,i2
               ix = i
               iz = level - i
               lupd(i+1-i1) = FALSE
               IF (ix > 1 .AND. iz > 1) lupd(i+1-i1) = TRUE
            ENDDO
         ELSEIF (sweep == 2) THEN
            !$OMP SIMD
            DO i=i1,i2
               ix = nx + 1 - i  !flip ix
               iz = level - i
               lupd(i+1-i1) = FALSE
               IF (ix < nx .AND. iz > 1) lupd(i+1-i1) = TRUE
            ENDDO
         ELSEIF (sweep == 3) THEN
            !$OMP SIMD
            DO i=i1,i2
               ix = i
               iz = nz + 1 - (level - i) ! flip iz
               lupd(i+1-i1) = FALSE
               IF (ix > 1 .AND. iz < nz) lupd(i+1-i1) = TRUE
            ENDDO
         ELSE
            !$OMP SIMD
            DO i=i1,i2
               ix = nx + 1 - i           ! flip ix
               iz = nz + 1 - (level - i) ! flip iz
               lupd(i+1-i1) = FALSE
               IF (ix < nx .AND. iz < nz) lupd(i+1-i1) = TRUE
            ENDDO
         ENDIF
      ELSE
         IF (sweep == 1) THEN
            !$OMP SIMD
            DO i=i1,i2
               ix = i
               iz = level - i
               lupd(i+1-i1) = FALSE
               IF (ix > MAX(1, xsi-1) .AND. iz > MAX(1, zsi-1)) lupd(i+1-i1) = TRUE
            ENDDO
         ELSEIF (sweep == 2) THEN 
            !$OMP SIMD
            DO i=i1,i2
               ix = nx + 1 - i  !flip ix
               iz = level - i
               lupd(i+1-i1) = FALSE
               IF (ix < MIN(nx, xsi+2) .AND. iz > MAX(1, zsi-1)) lupd(i+1-i1) = TRUE
            ENDDO
         ELSEIF (sweep == 3) THEN 
            !$OMP SIMD
            DO i=i1,i2
               ix = i
               iz = nz + 1 - (level - i) ! flip iz
               lupd(i+1-i1) = FALSE
               IF (ix > MAX(1, xsi-1) .AND. iz < MIN(nz, zsi+2)) lupd(i+1-i1) = TRUE 
            ENDDO
         ELSE
            !$OMP SIMD
            DO i=i1,i2
               ix = nx + 1 - i           ! flip ix
               iz = nz + 1 - (level - i) ! flip iz
               lupd(i+1-i1) = FALSE
               IF (ix < MIN(nx, xsi+2) .AND. iz < MIN(nz, zsi+2)) lupd(i+1-i1) = TRUE 
            ENDDO
         ENDIF
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Prefetches the slownesses for the local solver.
!>
!>    @param[in] level   Current level in level-set.
!>    @param[in] sweep   Sweep number.  This is in the range [1,4].
!>    @param[in] nz      Number of z grid points in travel time field.
!>    @param[in] nx      Number of x grid points in travel time field.
!>    @param[in] i1      Start index of level set.
!>    @param[in] i2      End indx of level set.
!>    @param[in] lupd    If true then the k'th node in the level is to be updated.
!>    @param[in] slow    Slowness field (s/m) defined at cells.  This is a vector
!>                       of dimension [nz-1 x nx-1] with leading dimenesion [nz-1].
!>
!>    @param[out] sloc   Slownesses (s/m) for finite difference stencil in the V, WE, and
!>                       home position respetively.  This is a vector of dimension
!>                       [4 x maxLevelSize] with leading dimension 4. 
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      PURE SUBROUTINE fteik_solver2d_prefetchSlowness64fF(level, sweep, nz, nx, &
                                                          i1, i2,               &
                                                          lupd, slow, sloc)     &
      BIND(C, NAME='fteik_solver2d_prefetchSlowness64fF')
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: i1, i2, level, nx, nz, sweep 
      REAL(C_DOUBLE), DIMENSION(:), INTENT(IN) :: slow
      LOGICAL(C_BOOL), DIMENSION(:), INTENT(IN) :: lupd
      REAL(C_DOUBLE), DIMENSION(:), INTENT(OUT) :: sloc
      INTEGER(C_INT) i, ix, ixcell, indx1, indx2, indx3, indx4, indx5, iz, izcell
      ! Extract slownesses
      IF (sweep == 1) THEN
         !$OMP SIMD
         DO i=i1,i2
            ix = i 
            iz = level - i 
            izcell = iz - 1
            ixcell = ix - 1
            indx1 = (MAX(ix-1, 1) - 1)*(nz - 1) + izcell ! V Plane
            indx2 = (MIN(ix,nx-1) - 1)*(nz - 1) + izcell ! V Plane
            indx3 = (ixcell - 1)*(nz - 1) + MAX(iz-1,1)  ! WE Plane
            indx4 = (ixcell - 1)*(nz - 1) + MIN(iz,nz-1) ! WE Plane
            indx5 = (ixcell - 1)*(nz - 1) + izcell
            !sloc(i+1-i1) = zero
            sloc(4*(i-i1)+1:4*(i-i1)+4) = zero
            IF (lupd(i+1-i1)) THEN
               !sloc(i+1-i1) = slow((ix - 2)*(nz - 1) + iz - 1)
               sloc(4*(i-i1)+1) = MIN(slow(indx1), slow(indx2)) 
               sloc(4*(i-i1)+2) = MIN(slow(indx3), slow(indx4))
               sloc(4*(i-i1)+3) = slow(indx5)
            ENDIF
         ENDDO
      ELSEIF (sweep == 2) THEN
         !$OMP SIMD
         DO i=i1,i2
            ix = nx + 1 - i  !flip ix
            iz = level - i
            izcell = iz - 1
            ixcell = ix
            indx1 = (MAX(ix-1, 1) - 1)*(nz - 1) + izcell ! V Plane
            indx2 = (MIN(ix,nx-1) - 1)*(nz - 1) + izcell ! V Plane
            indx3 = (ixcell - 1)*(nz - 1) + MAX(iz-1,1)  ! WE Plane
            indx4 = (ixcell - 1)*(nz - 1) + MIN(iz,nz-1) ! WE Plane
            indx5 = (ixcell - 1)*(nz - 1) + izcell
            !sloc(i+1-i1) = zero
            sloc(4*(i-i1)+1:4*(i-i1)+4) = zero
            IF (lupd(i+1-i1)) THEN
               !sloc(i+1-i1) = slow((ix - 1)*(nz - 1) + iz - 1)
               sloc(4*(i-i1)+1) = MIN(slow(indx1), slow(indx2)) 
               sloc(4*(i-i1)+2) = MIN(slow(indx3), slow(indx4))
               sloc(4*(i-i1)+3) = slow(indx5)
            ENDIF
         ENDDO
      ELSEIF (sweep == 3) THEN
         !$OMP SIMD
         DO i=i1,i2
            ix = i 
            iz = nz + 1 - (level - i) ! flip iz
            izcell = iz
            ixcell = ix - 1
            indx1 = (MAX(ix-1, 1) - 1)*(nz - 1) + izcell ! V Plane
            indx2 = (MIN(ix,nx-1) - 1)*(nz - 1) + izcell ! V Plane
            indx3 = (ixcell - 1)*(nz - 1) + MAX(iz-1,1)  ! WE Plane
            indx4 = (ixcell - 1)*(nz - 1) + MIN(iz,nz-1) ! WE Plane
            indx5 = (ixcell - 1)*(nz - 1) + izcell
            !sloc(i+1-i1) = zero
            sloc(4*(i-i1)+1:4*(i-i1)+4) = zero
            IF (lupd(i+1-i1)) THEN
               !sloc(i+1-i1) = slow((ix - 2)*(nz - 1) + iz)
               sloc(4*(i-i1)+1) = MIN(slow(indx1), slow(indx2)) 
               sloc(4*(i-i1)+2) = MIN(slow(indx3), slow(indx4))
               sloc(4*(i-i1)+3) = slow(indx5)
            ENDIF
         ENDDO 
      ELSE
         !$OMP SIMD
         DO i=i1,i2
            ix = nx + 1 - i           ! flip ix
            iz = nz + 1 - (level - i) ! flip iz
            izcell = iz
            ixcell = ix
            indx1 = (MAX(ix-1, 1) - 1)*(nz - 1) + izcell ! V Plane
            indx2 = (MIN(ix,nx-1) - 1)*(nz - 1) + izcell ! V Plane
            indx3 = (ixcell - 1)*(nz - 1) + MAX(iz-1,1)  ! WE Plane
            indx4 = (ixcell - 1)*(nz - 1) + MIN(iz,nz-1) ! WE Plane
            indx5 = (ixcell - 1)*(nz - 1) + izcell
            !sloc(i+1-i1) = zero
            sloc(4*(i-i1)+1:4*(i-i1)+4) = zero
            IF (lupd(i+1-i1)) THEN 
               !sloc(i+1-i1) = slow((ix - 1)*(nz - 1) + iz)
               sloc(4*(i-i1)+1) = MIN(slow(indx1), slow(indx2))
               sloc(4*(i-i1)+2) = MIN(slow(indx3), slow(indx4))
               sloc(4*(i-i1)+3) = slow(indx5)
            ENDIF
         ENDDO 
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Prefetches the travel times for the local solver.
!>
!>    @param[in] level   Current level in level-set.
!>    @param[in] sweep   Sweep number.  This is in the range [1,4].
!>    @param[in] nz      Number of z grid points in travel time field.
!>    @param[in] nx      Number of x grid points in travel time field.
!>    @param[in] i1      Start index of level set.
!>    @param[in] i2      End indx of level set.
!>    @param[in] lupd    If true then the k'th node in the level is to be updated.
!>    @param[in] tt      Travel time field (s) defined at nodes.  This is a vector
!>                       of dimension [nz x nx] with leading dimension [nz].
!>
!>    @param[out] ttvec  Travel times for finite differencing and updating.  
!>                       This is a vector of dimension [4, nnodes].  Each
!>                       four-tuple is packed tv, te, tev, tt.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      PURE SUBROUTINE fteik_solver2d_prefetchTravelTimes64fF(level, sweep, nz, nx, &
                                                             i1, i2,               &
                                                             lupd, tt, ttvec)      &
      BIND(C, NAME='fteik_solver2d_prefetchTravelTimes64fF')
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: i1, i2, level, nx, nz, sweep
      REAL(C_DOUBLE), DIMENSION(:), INTENT(IN) :: tt
      LOGICAL(C_BOOL), DIMENSION(:), INTENT(IN) :: lupd
      REAL(C_DOUBLE), DIMENSION(:), INTENT(OUT) :: ttvec
      REAL(C_DOUBLE) t, te, tev, tv
      INTEGER(C_INT) i, ix, iz
      IF (sweep == 1) THEN
         !$OMP SIMD
         DO i=i1,i2
            ix = i
            iz = level - i
            tev = FTEIK_HUGE
            te  = FTEIK_HUGE
            tv  = FTEIK_HUGE
            t   = FTEIK_HUGE
            IF (lupd(i+1-i1)) THEN
               tev = tt((ix - 2)*nz + iz - 1)
               te  = tt((ix - 2)*nz + iz)
               tv  = tt((ix - 1)*nz + iz - 1)
               t   = tt((ix - 1)*nz + iz)
            ENDIF
            ttvec(4*(i-i1)+1) = tv
            ttvec(4*(i-i1)+2) = te
            ttvec(4*(i-i1)+3) = tev
            ttvec(4*(i-i1)+4) = t
         ENDDO
      ELSEIF (sweep == 2) THEN
         !$OMP SIMD
         DO i=i1,i2
            ix = nx + 1 - i  !flip ix
            iz = level - i
            tev = FTEIK_HUGE
            te  = FTEIK_HUGE
            tv  = FTEIK_HUGE
            t   = FTEIK_HUGE
            IF (lupd(i+1-i1)) THEN
               tv  = tt((ix - 1)*nz + iz - 1)
               t   = tt((ix - 1)*nz + iz) 
               tev = tt(ix*nz + iz - 1)
               te  = tt(ix*nz + iz)
            ENDIF
            ttvec(4*(i-i1)+1) = tv
            ttvec(4*(i-i1)+2) = te
            ttvec(4*(i-i1)+3) = tev
            ttvec(4*(i-i1)+4) = t
         ENDDO
      ELSEIF (sweep == 3) THEN
         !$OMP SIMD
         DO i=i1,i2
            ix = i 
            iz = nz + 1 - (level - i) ! flip iz
            tev = FTEIK_HUGE
            te  = FTEIK_HUGE
            tv  = FTEIK_HUGE
            t   = FTEIK_HUGE
            IF (lupd(i+1-i1)) THEN
               te  = tt((ix - 2)*nz + iz)
               tev = tt((ix - 2)*nz + iz + 1)
               t   = tt((ix - 1)*nz + iz)
               tv  = tt((ix - 1)*nz + iz + 1)
            ENDIF
            ttvec(4*(i-i1)+1) = tv
            ttvec(4*(i-i1)+2) = te
            ttvec(4*(i-i1)+3) = tev
            ttvec(4*(i-i1)+4) = t
         ENDDO
      ELSE
         !$OMP SIMD
         DO i=i1,i2
            ix = nx + 1 - i           ! flip ix
            iz = nz + 1 - (level - i) ! flip iz
            tev = FTEIK_HUGE
            te  = FTEIK_HUGE
            tv  = FTEIK_HUGE
            t   = FTEIK_HUGE
            IF (lupd(i+1-i1)) THEN
               t   = tt((ix - 1)*nz + iz)
               tv  = tt((ix - 1)*nz + iz + 1)
               te  = tt(ix*nz + iz)
               tev = tt(ix*nz + iz + 1)
            ENDIF
            ttvec(4*(i-i1)+1) = tv
            ttvec(4*(i-i1)+2) = te
            ttvec(4*(i-i1)+3) = tev
            ttvec(4*(i-i1)+4) = t
         ENDDO
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Updates the travel time at each grid point in the sweep.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      PURE SUBROUTINE fteik_solver2d_updateTravelTimes64fF(level, sweep, nz, nx, &
                                                           i1, i2, lupd,         &
                                                           tupd, tt)
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: i1, i2, level, nx, nz, sweep
      REAL(C_DOUBLE), DIMENSION(:), INTENT(IN) :: tupd
      LOGICAL(C_BOOL), DIMENSION(:), INTENT(IN) :: lupd
      REAL(C_DOUBLE), DIMENSION(:), INTENT(INOUT) :: tt

      INTEGER i, indx, ix, iz
      IF (sweep == 1) THEN
         !$OMP SIMD
         DO i=i1,i2
            ix = i
            iz = level - i 
            indx = (ix - 1)*nz + iz
            !print *, level, ix, iz, indx, sngl(tt(indx)), lupd(i+1-i1)
            IF (lupd(i+1-i1)) THEN
               tt(indx) = MIN(tt(indx), tupd(i+1-i1))
            ENDIF
         ENDDO
      ELSEIF (sweep == 2) THEN
         !$OMP SIMD
         DO i=i1,i2
            ix = nx + 1 - i  !flip ix
            iz = level - i
            indx = (ix - 1)*nz + iz
            IF (lupd(i+1-i1)) THEN
               tt(indx) = MIN(tt(indx), tupd(i+1-i1))
            ENDIF
         ENDDO
      ELSEIF (sweep == 3) THEN
         !$OMP SIMD
         DO i=i1,i2
            ix = i 
            iz = nz + 1 - (level - i) ! flip iz
            indx = (ix - 1)*nz + iz
            IF (lupd(i+1-i1)) THEN
               tt(indx) = MIN(tt(indx), tupd(i+1-i1))
            ENDIF
         ENDDO
      ELSE
         !$OMP SIMD
         DO i=i1,i2
            ix = nx + 1 - i           ! flip ix
            iz = nz + 1 - (level - i) ! flip iz
            indx = (ix - 1)*nz + iz
            IF (lupd(i+1-i1)) THEN
               tt(indx) = MIN(tt(indx), tupd(i+1-i1)) 
            ENDIF
         ENDDO
      ENDIF
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Solves the eikonal equation with the level-set fast sweeping method.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver2d_solveSourceLSM(isrc, ierr)   &
                 BIND(C, NAME='fteik_solver2d_solveSourceLSM')
      USE FTEIK_SOURCE64F, ONLY : nsrc
      USE FTEIK_MODEL64F, ONLY : lhaveModel, nx, nz, slow
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, TRUE, FALSE
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) ts4(4), t0, t1
      INTEGER(C_INT) dest(4), i1, i2, kiter, level, nnodes, sweep
      INTEGER(C_INT) iz, ix
      ierr = 0
      lhaveTimes = FALSE
      IF (.NOT.lhaveModel) THEN
         WRITE(*,*) 'fteik_solver2d_solveSourceLSMF: Model not yet set'
         ierr = 1
         GOTO 500
      ENDIF
      IF (isrc < 1 .OR. isrc > nsrc) THEN
         WRITE(*,*) 'fteik_solver2d_solveSourceLSMF: Invalid source number:',isrc, 1, nsrc
         ierr = 1
         GOTO 500
      ENDIF
      CALL fteik_localSolver2d_initialize64f(isrc, dest, ts4, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver2d_solveSourceLSMF: Error setting mesh constants'
         GOTO 500
      ENDIF
      t0 = 0.d0
      CALL CPU_TIME(t0)
      ttimes(:) = FTEIK_HUGE
      ttimes(dest(1:4)) = ts4(1:4)
!print *, dest
!print *, ts4
!return
      ! Initialization
      ! First sweeping: Top->Bottom ; West->East
      DO ix=MAX(2,xsi),nx
         DO iz=MAX(2,zsi),nz
            CALL fteik_localSolver2d_init64fF(1, nz, nx, iz, ix, &
                                              slow, ttimes, tupdWork(1))
            ttimes((ix-1)*nz + iz) = MIN(ttimes((ix-1)*nz + iz), tupdWork(1))
         ENDDO
      ENDDO
      ! Second sweeping: Top->Bottom ; East->West
      DO ix=xsi+1,1,-1
         DO iz=MAX(2,zsi),nz
            CALL fteik_localSolver2d_init64fF(2, nz, nx, iz, ix, &
                                              slow, ttimes, tupdWork(1))
            ttimes((ix-1)*nz + iz) = MIN(ttimes((ix-1)*nz + iz), tupdWork(1))
         ENDDO
      ENDDO
      ! Third sweep: Bottom->Top ; West->East
      DO ix=MAX(2,xsi),nx
         DO iz=zsi+1,1,-1
            CALL fteik_localSolver2d_init64fF(3, nz, nx, iz, ix, &
                                              slow, ttimes, tupdWork(1))
            ttimes((ix-1)*nz + iz) = MIN(ttimes((ix-1)*nz + iz), tupdWork(1))
         ENDDO
      ENDDO
      ! Fourth sweeping: Bottom->Top ; East->West
      DO ix=xsi+1,1,-1
         DO iz=zsi+1,1,-1
            CALL fteik_localSolver2d_init64fF(4, nz, nx, iz, ix, &
                                              slow, ttimes, tupdWork(1))
            ttimes((ix-1)*nz + iz) = MIN(ttimes((ix-1)*nz + iz), tupdWork(1))
         ENDDO
      ENDDO
!     ! Initialize
!     DO sweep=1,4
!        ! Loop on the levels
!        DO level=2,nlevels
!           i1 = MAX(1, level - nz)
!           i2 = MIN(nx, level - 1) 
!           nnodes = i2 - i1 + 1
!           ! Get the update nodes
!           CALL fteik_solver2d_isUpdateNodeF(TRUE, level, sweep, nz, nx, &
!                                             i1, i2, zsi, xsi, lupdWork)
!           ! Prefetch the slowness
!           CALL fteik_solver2d_prefetchSlowness64fF(level, sweep, nz, nx,   &
!                                                    i1, i2,                 &
!                                                    lupdWork, slow, slocWork)
!           ! Prefetch travel times
!           CALL fteik_solver2d_prefetchTravelTimes64fF(level, sweep, nz, nx,      &
!                                                       i1, i2,                    &
!                                                       lupdWork, ttimes, ttvecWork)
!           ! Compute the candidate travel times at each node in level set
!           !CALL fteik_localSolver2d_noInit64fF(nnodes, ttvecWork, slocWork, tupdWork)
!           ! And update
!           !CALL fteik_solver2d_updateTravelTimes64fF(level, sweep, nz, nx, &
!           !                                          i1, i2, lupdWork,     &
!           !                                          tupdWork, ttimes)
!        ENDDO
!     ENDDO
      ! Number of Gauss-Seidel iterations
      DO kiter=1,nsweep
         ! Loop on sweeping directions
         DO sweep=1,4
            ! Loop on the levels
            DO level=2,nlevels
               i1 = MAX(1, level - nz)
               i2 = MIN(nx, level - 1)
               nnodes = i2 - i1 + 1
               ! Get the update nodes
               CALL fteik_solver2d_isUpdateNodeF(FALSE, level, sweep, nz, nx, &
                                                 i1, i2, zsi, xsi, lupdWork)
               ! Prefetch the slowness
               CALL fteik_solver2d_prefetchSlowness64fF(level, sweep, nz, nx,   &
                                                        i1, i2,                 &
                                                        lupdWork, slow, slocWork)
               ! Prefetch travel times
               CALL fteik_solver2d_prefetchTravelTimes64fF(level, sweep, nz, nx,      &
                                                           i1, i2,            &
                                                           lupdWork, ttimes, ttvecWork)
               ! Compute the candidate travel times at each node in level set
               CALL fteik_localSolver2d_noInit64fF(nnodes, ttvecWork, slocWork, tupdWork)
               ! And update
               CALL fteik_solver2d_updateTravelTimes64fF(level, sweep, nz, nx, &
                                                         i1, i2, lupdWork, &
                                                         tupdWork, ttimes)
            ENDDO
            !do i2=1,nx*nz
            !   print *, i2, ttimes(i2)
            !enddo
            !print *, minval(ttimes), maxval(ttimes)
            IF (verbose > 3) WRITE(*,901) kiter, sweep, MINVAL(ttimes), MAXVAL(ttimes)
         ENDDO
         IF (verbose > 2) WRITE(*,900) kiter, MINVAL(ttimes), MAXVAL(ttimes)
      ENDDO
      IF (verbose > 2) THEN
         CALL CPU_TIME(t1)
         WRITE(*,902) t1 - t0 
      ENDIF
      lhaveTimes = TRUE
  500 CONTINUE
      IF (ierr /= 0) THEN
         ttimes(:) = FTEIK_HUGE
         lhaveTimes = FALSE
      ENDIF
  900 FORMAT(' fteik_solver2d_solveSourceLSM: (iteration,ttmin,ttmax)=', &
             I4, 2F14.8) 
  901 FORMAT(' fteik_solver2d_solveSourceLSM: (iteration,sweep,ttmin,ttmax)=', &
             2I4, 2F14.8)
  902 FORMAT(' fteik_solver2d_solveSourceLSM: Solver time in seconds=', F14.8)
!print *, lhaveTimes
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE fteik_solver2d_solveSourceFSM(isrc, ierr) &
      BIND(C, NAME='fteik_solver2d_solveSourceFSM')
      USE FTEIK_SOURCE64F, ONLY : nsrc
      USE FTEIK_MODEL64F, ONLY : lhaveModel, nx, nz, slow
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, TRUE, FALSE
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) ts4(4), t0
      INTEGER(C_INT) dest(4), i, i1, ix, iz, j, j1, kiter, level, sgntx, sgntz, sgnvx, sgnvz 
      INTEGER(C_INT) indx1, indx2, indx3, indx4, indx5
      ierr = 0 
      lhaveTimes = FALSE
      IF (.NOT.lhaveModel) THEN
         WRITE(*,*) 'fteik_solver2d_solveSourceFSM: Model not yet set'
         ierr = 1
         GOTO 500
      ENDIF
      IF (isrc < 1 .OR. isrc > nsrc) THEN
         WRITE(*,*) 'fteik_solver2d_solveSourceFSM: Invalid source number:',isrc, 1, nsrc
         ierr = 1
         GOTO 500
      ENDIF
      CALL fteik_localSolver2d_initialize64f(isrc, dest, ts4, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver2d_solveSourceFSM: Error setting mesh constants'
         GOTO 500
      ENDIF
      CALL CPU_TIME(t0)
      ttimes(:) = FTEIK_HUGE
      ttimes(dest(1:4)) = ts4(1:4)
      ! Initialization
      ! First sweeping: Top->Bottom ; West->East
      sgntz = 1
      sgntx = 1
      sgnvz = 1
      sgnvx = 1
      DO ix=MAX(2,xsi),nx
         DO iz=MAX(2,zsi),nz
            CALL fteik_localSolver2d_init64fF(1, nz, nx, iz, ix, &
                                              slow, ttimes, tupdWork(1))
            ttimes((ix-1)*nz + iz) = MIN(ttimes((ix-1)*nz + iz), tupdWork(1))
         ENDDO
      ENDDO
      ! Second sweeping: Top->Bottom ; East->West
      sgntz = 1
      sgntx =-1
      sgnvz = 1
      sgnvx = 0
      DO ix=MIN(nx-1,xsi+1),1,-1
         DO iz=MAX(2,zsi),nz
            CALL fteik_localSolver2d_init64fF(2, nz, nx, iz, ix, &
                                              slow, ttimes, tupdWork(1))
            ttimes((ix-1)*nz + iz) = MIN(ttimes((ix-1)*nz + iz), tupdWork(1))
         ENDDO
      ENDDO
      ! Third sweep: Bottom->Top ; West->East
      sgntz =-1
      sgntx = 1
      sgnvz = 0
      sgnvx = 1
      DO ix=MAX(2,xsi),nx
         DO iz=MIN(nz-1,zsi+1),1,-1
            CALL fteik_localSolver2d_init64fF(3, nz, nx, iz, ix, &
                                              slow, ttimes, tupdWork(1))
            ttimes((ix-1)*nz + iz) = MIN(ttimes((ix-1)*nz + iz), tupdWork(1))
         ENDDO
      ENDDO
      ! Fourth sweeping: Bottom->Top ; East->West
      sgntz =-1
      sgntx =-1
      sgnvz = 0
      sgnvx = 0
      DO ix=MIN(nx-1,xsi+1),1,-1
         DO iz=MIN(nz-1,zsi+1),1,-1
            CALL fteik_localSolver2d_init64fF(4, nz, nx, iz, ix, &
                                              slow, ttimes, tupdWork(1))
            ttimes((ix-1)*nz + iz) = MIN(ttimes((ix-1)*nz + iz), tupdWork(1))
         ENDDO
      ENDDO 
      ! Gauss-Seidel iterations 
      DO kiter=1,nsweep
         ! First sweeping: Top->Bottom ; West->East
         sgntz = 1 
         sgntx = 1 
         sgnvz = 1 
         sgnvx = 1 
         DO j = 2, nx
            DO i = 2, nz
               i1 = i - sgnvz 
               j1 = j - sgnvx

               ix = i
               iz = i
               level = iz + ix
               ttvecWork(1) = ttimes((j - 1)*nz         + i - sgntz)
               ttvecWork(2) = ttimes((j - sgntx - 1)*nz + i)
               ttvecWork(3) = ttimes((j - sgntx - 1)*nz + i - sgntz)
               ttvecWork(4) = ttimes((j - 1)*nz         + i) 
               indx1 = (MAX(j-1,1)  - 1)*(nz - 1) + i - sgnvz ! V Plane
               indx2 = (MIN(j,nx-1) - 1)*(nz - 1) + i - sgnvz ! V Plane
               indx3 = (j - sgnvx - 1)*(nz - 1) + MAX(i-1,1)  ! WE Plane
               indx4 = (j - sgnvx - 1)*(nz - 1) + MIN(i,nz-1) ! WE Plane
               indx5 = (j - sgnvx - 1)*(nz - 1) + i - sgnvz
               slocWork(1) = MIN(slow(indx1), slow(indx2))
               slocWork(2) = MIN(slow(indx3), slow(indx4))
               slocWork(3) = slow(indx5)
               !slocWork(1) = slow((j1 - 1)*(nz - 1) + i1)
               CALL fteik_localSolver2d_noInit64fF(1, ttvecWork, slocWork, tupdWork)
               ttimes((j-1)*nz + i) = MIN(ttimes((j-1)*nz + i), tupdWork(1))
            ENDDO
         ENDDO
!do i=1,nx*nz
!print *, 'ref', i, ttimes(i)
!enddo
print *, 'p1:', minval(ttimes), maxval(ttimes)
         ! Second sweeping: Top->Bottom ; East->West
         sgntz = 1 
         sgntx =-1
         sgnvz = 1 
         sgnvx = 0 
         DO j = nx-1, 1, -1
            DO i = 2, nz
               i1 = i - sgnvz 
               j1 = j - sgnvx

               ix = i
               iz = i
               level = iz + ix
!              CALL fteik_solver2d_prefetchTravelTimes64fF(level, 2, nz, nx,  &
!                                                          1, ix, ix,         &
!                                                          lupd, ttimes, ttvec)
!              CALL fteik_solver2d_prefetchSlowness64fF(level, 2, nz, nx, &
!                                                       1, ix, ix,        &
!                                                       lupd, slow, sloc)
!print *, sngl(ttvec(1:4))
               ttvecWork(1) = ttimes((j - 1)*nz         + i - sgntz)
               ttvecWork(2) = ttimes((j - sgntx - 1)*nz + i)
               ttvecWork(3) = ttimes((j - sgntx - 1)*nz + i - sgntz)
               ttvecWork(4) = ttimes((j - 1)*nz         + i)
!print *, sngl(ttvec(1:4))
!pause
               indx1 = (MAX(j-1,1)  - 1)*(nz - 1) + i - sgnvz ! V Plane
               indx2 = (MIN(j,nx-1) - 1)*(nz - 1) + i - sgnvz ! V Plane
               indx3 = (j - sgnvx - 1)*(nz - 1) + MAX(i-1,1)  ! WE Plane
               indx4 = (j - sgnvx - 1)*(nz - 1) + MIN(i,nz-1) ! WE Plane
               indx5 = (j - sgnvx - 1)*(nz - 1) + i - sgnvz
               slocWork(1) = MIN(slow(indx1), slow(indx2))
               slocWork(2) = MIN(slow(indx3), slow(indx4))
               slocWork(3) = slow(indx5)
               !slocWork(1) = slow((j1 - 1)*(nz - 1) + i1) 
               CALL fteik_localSolver2d_noInit64fF(1, ttvecWork, slocWork, tupdWork)
               ttimes((j-1)*nz + i) = MIN(ttimes((j-1)*nz + i), tupdWork(1))
            ENDDO
         ENDDO
print *, 'p2:', minval(ttimes), maxval(ttimes)
         ! Third sweep: Bottom->Top ; West->East
         sgntz =-1
         sgntx = 1
         sgnvz = 0
         sgnvx = 1
         DO j = 2, nx
            DO i = nz-1, 1, -1
               i1 = i - sgnvz 
               j1 = j - sgnvx
               ttvecWork(1) = ttimes((j - 1)*nz         + i - sgntz)
               ttvecWork(2) = ttimes((j - sgntx - 1)*nz + i) 
               ttvecWork(3) = ttimes((j - sgntx - 1)*nz + i - sgntz)
               ttvecWork(4) = ttimes((j - 1)*nz         + i)
               indx1 = (MAX(j-1,1)  - 1)*(nz - 1) + i - sgnvz ! V Plane
               indx2 = (MIN(j,nx-1) - 1)*(nz - 1) + i - sgnvz ! V Plane
               indx3 = (j - sgnvx - 1)*(nz - 1) + MAX(i-1,1)  ! WE Plane
               indx4 = (j - sgnvx - 1)*(nz - 1) + MIN(i,nz-1) ! WE Plane
               indx5 = (j - sgnvx - 1)*(nz - 1) + i - sgnvz
               slocWork(1) = MIN(slow(indx1), slow(indx2))
               slocWork(2) = MIN(slow(indx3), slow(indx4))
               slocWork(3) = slow(indx5)
               !slocWork(1) = slow((j1 - 1)*(nz - 1) + i1)
               CALL fteik_localSolver2d_noInit64fF(1, ttvecWork, slocWork, tupdWork)
               ttimes((j-1)*nz + i) = MIN(ttimes((j-1)*nz + i), tupdWork(1))
            ENDDO
         ENDDO
print *, 'p3:', minval(ttimes), maxval(ttimes)
         ! Fourth sweeping: Bottom->Top ; East->West
         sgntz =-1
         sgntx =-1
         sgnvz = 0
         sgnvx = 0
         DO j = nx-1, 1, -1
            DO i = nz-1, 1, -1
               i1 = i - sgnvz 
               j1 = j - sgnvx
               ttvecWork(1) = ttimes((j - 1)*nz         + i - sgntz)
               ttvecWork(2) = ttimes((j - sgntx - 1)*nz + i) 
               ttvecWork(3) = ttimes((j - sgntx - 1)*nz + i - sgntz)
               ttvecWork(4) = ttimes((j - 1)*nz         + i)
               indx1 = (MAX(j-1,1)  - 1)*(nz - 1) + i - sgnvz ! V Plane
               indx2 = (MIN(j,nx-1) - 1)*(nz - 1) + i - sgnvz ! V Plane
               indx3 = (j - sgnvx - 1)*(nz - 1) + MAX(i-1,1)  ! WE Plane
               indx4 = (j - sgnvx - 1)*(nz - 1) + MIN(i,nz-1) ! WE Plane
               indx5 = (j - sgnvx - 1)*(nz - 1) + i - sgnvz
               slocWork(1) = MIN(slow(indx1), slow(indx2))
               slocWork(2) = MIN(slow(indx3), slow(indx4))
               slocWork(3) = slow(indx5)
               !slocWork(1) = slow((j1 - 1)*(nz - 1) + i1)
               CALL fteik_localSolver2d_noInit64fF(1, ttvecWork, slocWork, tupdWork)
               ttimes((j-1)*nz + i) = MIN(ttimes((j-1)*nz + i), tupdWork(1))
            ENDDO
         ENDDO
print *, 'p4:', minval(ttimes), maxval(ttimes)
      ENDDO
      lhaveTimes = TRUE
  500 CONTINUE
      IF (ierr /= 0) THEN
         ttimes(:) = FTEIK_HUGE
         lhaveTimes = FALSE
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Returns the travel time field.
!>
!>    @param[in] ngin    Number of input grid points.  This must equal [nz x nx].
!>
!>    @param[out] ttout  Travel time field at the grid points.  This is a [nz x nx]
!>                       vector with leading dimension nz.
!>    @param[out] ierr   0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver2d_getTravelTimeField64f(ngin, ttout, ierr) &
      BIND(C, NAME='fteik_solver2d_getTravelTimeField64f')
      USE FTEIK_MODEL64F, ONLY : ngrd
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: ngin
      REAL(C_DOUBLE), INTENT(OUT) :: ttout(ngin)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      IF (.NOT.lhaveTimes) THEN 
         WRITE(*,*) 'fteik_solver2d_getTravelTimeField64f: Travel times not yet computed'
         ierr = 1
         RETURN
      ENDIF
      IF (ngin /= ngrd) THEN
         WRITE(*,*) 'fteik_solver2d_getTravelTimeField64f: No receivers'
         RETURN
      ENDIF
      ttout(:) = ttimes(:)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Extracts the travel times at the receivers.
!>
!>    @param[in] nrec   Number of receivers.
!>
!>    @param[out] ttr   Travel times (seconds) at the receivers.  This is a vector of
!>                      dimension [nrec].
!>    @param[out] ierr  0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver2d_getTravelTimes64f(nrec, ttr, ierr) &
      BIND(C, NAME='fteik_solver2d_getTravelTimes64f')
      USE FTEIK_RECEIVER64F, ONLY : fteik_receiver_getTravelTimes64f
      USE FTEIK_MODEL64F, ONLY : ngrd
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nrec
      REAL(C_DOUBLE), INTENT(OUT) :: ttr(nrec)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      IF (nrec <= 0) THEN
         WRITE(*,*) 'fteik_solver2d_getTravelTimes64f: No receivers'
         RETURN
      ENDIF
      IF (.NOT.lhaveTimes) THEN
         WRITE(*,*) 'fteik_solver2d_getTravelTimes64f: Travel times not yet computed'
         ierr = 1
         RETURN
      ENDIF
      CALL fteik_receiver_getTravelTimes64f(nrec, ngrd, ttimes, ttr, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver2d_getTravelTimes64f: Error getting travel times'
         ierr = 1
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief 2D local solver for Cartesian frame only.
!>
!>    @param[in] n       Number of nodes to update. 
!>    @param[in] ttvec   The local travel times (seconds) surrounding the point.  This is
!>                       a vector of dimension [4 x n] with leading dimemsion 4.  Each
!>                       four-tuple is packed [tv, te, tev, tt].
!>    @param[in] sloc    Slowness (s/m) at the finite difference points.  This is a vector
!>                       of dimension [n].
!>    @param[out] tupd   Updated travel times (seconds) at each node (seconds).  This is
!>                       a vector of dimension [n].
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      PURE SUBROUTINE fteik_localSolver2d_noInit64fF(n, ttvec, sloc, tupd)
      USE ISO_C_BINDING
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, four, chunkSize
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: n
      REAL(C_DOUBLE), INTENT(IN) :: ttvec(4*n), sloc(4*n)
      REAL(C_DOUBLE), INTENT(OUT) :: tupd(n)
      REAL(C_DOUBLE) four_sref2, s1, s2, s3, sref2, t1_2d, t12min, ta, tb, tab, &
                     tab2, temtv, te, tev, tv, tt
      INTEGER(C_INT) i
      ! Loop on nodes in update
      DO i=1,n
         s1 = sloc(4*(i-1)+1)
         s2 = sloc(4*(i-1)+2)
         s3 = sloc(4*(i-1)+3)
         tv  = ttvec(4*(i-1)+1)
         te  = ttvec(4*(i-1)+2)
         tev = ttvec(4*(i-1)+3)
         tt  = ttvec(4*(i-1)+4)
         temtv = te - tv
         ! 1D operators (refracted times)
         t12min = MIN(tv + dz*s1, te + dx*s2)
         ! 2D operator
         t1_2d = FTEIK_HUGE
         IF (temtv < dz*s3 .AND. -temtv < dx*s3) THEN
            sref2 = s3*s3
            ta = tev + temtv
            tb = tev - temtv
            tab = ta - tb
            tab2 = tab*tab
            four_sref2 = four*sref2
            t1_2d = ( (tb*dz2i + ta*dx2i) &
                     + SQRT(four_sref2*dz2i_p_dx2i - dz2i_dx2i*tab2) )*dz2i_p_dx2i_inv
         ENDIF
         tupd(i) = MIN(t12min, t1_2d)
      ENDDO
      RETURN
      END

      SUBROUTINE fteik_localSolver2d_init64fF(sweep, nz, nx, iz, ix, &
                                              slow, tt, tupd)
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, four
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: sweep, nz, nx, iz, ix
      REAL(C_DOUBLE), INTENT(OUT) :: tupd
      REAL(C_DOUBLE), DIMENSION(:), INTENT(IN) :: tt, slow
      REAL(C_DOUBLE) apoly, bpoly, cpoly, dpoly, &
                     ta, tb, taue, tauev, tauv, tdiag, tv, te, tev, t, &
                     t1, t2, t3, t1d, t2d, &
                     sref1, sref2, sref3, sgnrx, sgnrz, &
                     sgnrx_txc, sgnrz_tzc, sgnrx_txc_dxi, sgnrz_tzc_dzi, t0c, txc, tzc
      INTEGER(C_INT) i1, indx, indxv, indxe, indxev, j1
      INTEGER(C_INT) indx1, indx2, indx3, indx4, indx5
      INTEGER(C_INT) sgntx, sgntz, sgnvx, sgnvz
      !LOGICAL, PARAMETER :: linitk = .TRUE.
      IF (sweep == 1) THEN
         sgntz = 1
         sgntx = 1
         sgnvz = 1
         sgnvx = 1
      ELSE IF (sweep == 2) THEN
         sgntz = 1
         sgntx =-1
         sgnvz = 1
         sgnvx = 0
      ELSE IF (sweep == 3) THEN
         sgntz =-1
         sgntx = 1
         sgnvz = 0
         sgnvx = 1
      ELSE
         sgntz =-1
         sgntx =-1
         sgnvz = 0
         sgnvx = 0
      ENDIF
      sgnrz = DBLE(sgntz)
      sgnrx = DBLE(sgntx)
      indxv  = (ix - 1)*nz         + iz - sgntz
      indxe  = (ix - sgntx - 1)*nz + iz
      indxev = (ix - sgntx - 1)*nz + iz - sgntz 
      indx   = (ix - 1)*nz         + iz
      tv  = tt(indxv)  !dble( tt(i-sgntz,j,k) )
      te  = tt(indxe)  !dble( tt(i,j-sgntx,k) )
      tev = tt(indxev) !dble( tt(i-sgntz,j-sgntx,k) )
      t   = tt(indx)
      ! Extract slownesses
      i1 = iz - sgnvz
      j1 = ix - sgnvx
      indx1 = (MAX(ix-1,1)  - 1)*(nz - 1) + i1 ! V Plane
      indx2 = (MIN(ix,nx-1) - 1)*(nz - 1) + i1 ! V Plane
      indx3 = (j1 - 1)*(nz - 1) + MAX(iz-1,1)  ! WE Plane
      indx4 = (j1 - 1)*(nz - 1) + MIN(iz,nz-1) ! WE Plane
      indx5 = (j1 - 1)*(nz - 1) + i1
      sref1 = MIN(slow(indx1), slow(indx2))
      sref2 = MIN(slow(indx3), slow(indx4))
      sref3 = slow(indx5)
      t1d = MIN(tv + dz*sref1, te + dx*sref2) ! min(V, WE) planes
      ! 2D and diagonal operators
      CALL fteik_localSolver2d_tAnaD64fF(t0c,              &
                                         tzc, txc, iz, ix, &
                                         dz, dx, zsa, xsa, &
                                         szero)
      tauv = tv   - fteik_localSolver2d_tAna64fF(iz-sgntz,   ix,        &
                                                 dz, dx, zsa, xsa, szero)
      taue = te   - fteik_localSolver2d_tAna64fF(iz, ix-sgntx,          &
                                                 dz, dx, zsa, xsa, szero)
      tauev = tev - fteik_localSolver2d_tAna64fF(iz-sgntz,   ix-sgntx,  &
                                                 dz, dx, zsa, xsa, szero)
      sgnrz_tzc = sgnrz*tzc
      sgnrx_txc = sgnrx*txc
      sgnrz_tzc_dzi = sgnrz_tzc*dzi
      sgnrx_txc_dxi = sgnrx_txc*dxi
      ! Diagonal operator
      tdiag = tev + sref3*SQRT(dx*dx + dz*dz)
      ! Choose spherical or plane wave
      t1 = FTEIK_HUGE
      t2 = FTEIK_HUGE
      t3 = FTEIK_HUGE
      ! Choose spherical or plane wave; first test for Plane wave
      IF ( ( ABS(iz - zsi) > epsS2C .OR. ABS(ix - xsi) > epsS2C ) ) THEN
         ! 4 Point operator, if possible otherwise do three points
         IF (tv <= te + dx*sref3 .AND. te <= tv + dz*sref3 .AND. &
             te >= tev .AND. tv >= tev ) THEN
            ta = tev + te - tv
            tb = tev - te + tv
            t1 = ( ( tb * dz2i + ta * dx2i ) + SQRT( 4.d0 * sref3*sref3 * ( dz2i + dx2i ) &
                 - dz2i * dx2i * ( ta - tb ) * ( ta - tb ) ) ) / ( dz2i + dx2i )
         ! Two 3 point operators
         ELSEIF (( te - tev ) <= dz*dz*sref3/SQRT(dx*dx + dz*dz) .AND. te > tev) THEN
            t2 = te + dx*SQRT(sref3*sref3 - ((te - tev)/dz)**2)
         ELSEIF (( tv - tev ) <= dx*dx*sref3/SQRT(dx*dx + dz*dz) .AND. tv > tev) THEN
            t3 = tv + dz*SQRT(sref3*sref3 - ((tv - tev)/dx)**2)
         ENDIF
      ELSE
         ! Do spherical operator if conditions ok
         IF ( tv < te + dx*sref3 .AND. te < tv + dz*sref3 .AND. &
             te >= tev .AND. tv >= tev) THEN
            ta = tauev + taue - tauv   ! X
            tb = tauev - taue + tauv   ! Z
            apoly = dz2i + dx2i
            bpoly = 4.d0 * ( sgnrx * txc * dxi + sgnrz * tzc * dzi ) &
                  - 2.d0 * ( ta * dx2i + tb * dz2i )
            cpoly = ( ta*ta * dx2i ) + ( tb*tb * dz2i ) &
                  - 4.d0 * ( sgnrx * txc * dxi * ta + sgnrz * tzc * dzi * tb ) &
                  + 4.d0 * ( szero*szero - sref3*sref3)
            dpoly = bpoly*bpoly - 4.d0 * apoly * cpoly
            IF (dpoly >= 0.d0) t1 = 0.5d0*(SQRT(dpoly) - bpoly )/apoly + t0c
            IF (t1 < tv .OR. t1 < te) t1 = FTEIK_HUGE
         ENDIF
      ENDIF
      t2d  = MIN(t1, t2, t3)
      tupd = MIN(t, t1d, t2d)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Initializes parameters for the local solver and computes initial travel 
!>           times around the source.
!>
!>    @param[in] isrc   Source number.
!>
!>    @param[out] dest  Indices in travel time field where ts4 should be inserted.
!>    @param[out] ts4   Analytic travel times computed at grid points around source.
!>    @param[out] ierr  0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!> 
      SUBROUTINE fteik_localSolver2d_initialize64f(isrc, dest, ts4, ierr)
      USE FTEIK_SOURCE64F, ONLY : fteik_source_getSolverInfo64fF
      USE FTEIK_MODEL64F, ONLY : fteik_model_grid2indexF
      USE FTEIK_MODEL64F, ONLY : nz, nzx
      USE FTEIK_MODEL64F, ONLY : dzIn => dz
      USE FTEIK_MODEL64F, ONLY : dxIn => dx
      USE FTEIK_CONSTANTS64F, ONLY : one
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      REAL(C_DOUBLE), INTENT(OUT) :: ts4(4)
      INTEGER(C_INT), INTENT(OUT) :: dest(4), ierr
      REAL(C_DOUBLE) dz2, dx2
      REAL(C_DOUBLE) ysa
      INTEGER(C_INT) ysi 
      ! Things to copy from model
      ierr = 0
      dz = dzIn
      dx = dxIn
      ! Things to grab from the source
      CALL fteik_source_getSolverInfo64fF(isrc,          &
                                          zsi, xsi, ysi, &
                                          zsa, xsa, ysa, &
                                          szero, szero2, &
                                          ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_localSolver2d_initialize64f: Problem with source'
         RETURN
      ENDIF
      eps = epsS2C*szero*MIN(dz, dx)
      ! Some constants
      dz2 = dzIn*dzIn
      dx2 = dxIn*dxIn
      dzi = one/dzIn
      dxi = one/dxIn
      dz2i = one/dz2
      dx2i = one/dx2
      ! Some more constants 
      dz2i = one/dz2
      dx2i = one/dx2
      dz2i_dx2i = dz2i*dx2i
      dz2i_p_dx2i = dz2i + dx2i
      dz2i_p_dx2i_inv = one/(dz2i + dx2i)
      ! Initialize points around source
      ts4(1) = fteik_localSolver2d_tAna64fF(zsi,   xsi,   dz, dx, zsa, xsa, szero)
      ts4(2) = fteik_localSolver2d_tAna64fF(zsi+1, xsi,   dz, dx, zsa, xsa, szero)
      ts4(3) = fteik_localSolver2d_tAna64fF(zsi,   xsi+1, dz, dx, zsa, xsa, szero)
      ts4(4) = fteik_localSolver2d_tAna64fF(zsi+1, xsi+1, dz, dx, zsa, xsa, szero)
      dest(1) = fteik_model_grid2indexF(zsi,   xsi  , 1, nz, nzx)
      dest(2) = fteik_model_grid2indexF(zsi+1, xsi  , 1, nz, nzx)
      dest(3) = fteik_model_grid2indexF(zsi  , xsi+1, 1, nz, nzx)
      dest(4) = fteik_model_grid2indexF(zsi+1, xsi+1, 1, nz, nzx)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the number of Gauss-Seidel iterations  in the fast sweeping method.
!>
!>    @param[in] nsweepIn    Number of sweeps.  This cannot be negative and 1 is
!>                           usually sufficient.
!>
!>    @param[out] ierr       0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver2d_setNumberOfSweeps(nsweepIn, ierr) &
      BIND(C, NAME='fteik_solver2d_setNumberOfSweeps')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsweepIn 
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0 
      nsweep = 0 
      IF (nsweepIn < 0) THEN
         WRITE(*,*) 'fteik_solver2d_setNumberOfSweeps: nsweep must be positive', nsweep
         ierr = 1 
         RETURN
      ENDIF
      nsweep = nsweepIn
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the receivers in the model.
!>
!>    @param[in] nrec    Number of receivers to set.  
!>    @param[in] zrec    z locations (meters) of receivers.  This is a vector of
!>                       of dimension [nrec].
!>    @param[in] xrec    x locations (meters) of receivers.  This is a vector of
!>                       of dimension [nrec].
!>
!>    @param[out] ierr   0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver2d_setReceivers64f(nrec, zrec, xrec, ierr) &
      BIND(C, NAME='fteik_solver2d_setReceivers64f')
      USE FTEIK_MODEL64F, ONLY : y0
      USE FTEIK_RECEIVER64F, ONLY : fteik_receiver_initialize64f
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nrec
      REAL(C_DOUBLE), INTENT(IN) :: zrec(nrec), xrec(nrec)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE), ALLOCATABLE :: yrec(:)
      ! It's actually safe to have no receivers
      ierr = 0
      IF (nrec < 1) THEN
         WRITE(*,*) 'fteik_solver_setReceivers64fF: No receivers to set'
         RETURN
      ENDIF
      ALLOCATE(yrec(MAX(nrec, 1)))
      yrec(:) = y0
      CALL fteik_receiver_initialize64f(nrec, zrec, xrec, yrec, verbose, ierr)
      IF (ierr /= 0) WRITE(*,*) 'fteik_solver_setReceivers64f: Failed to set receivers'
      IF (ALLOCATED(yrec)) DEALLOCATE(yrec)
      RETURN
      END
!                                                                                        !!========================================================================================!
!                                                                                        !!>    @brief Initializes the source(s) on the solver.
!> 
!>    @param[in] nsrc     Number of sources.
!>    @param[in] zsrc     z locations (meters) of source.  This is a vector of dimension
!>                        [nsrc].
!>    @param[in] xsrc     x locations (meters) of source.  This is a vector of dimension
!>                        [nsrc].
!>
!>    @param[out] ierr    0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver2d_setSources64f(nsrc, zsrc, xsrc, ierr) &
      BIND(C, NAME='fteik_solver2d_setSources64f')
      USE FTEIK_MODEL64F, ONLY : y0
      USE FTEIK_SOURCE64F, ONLY : fteik_source_initialize64f
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsrc
      REAL(C_DOUBLE), INTENT(IN) :: zsrc(nsrc), xsrc(nsrc)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE), ALLOCATABLE :: ysrc(:)
      ALLOCATE(ysrc(MAX(1, nsrc)))
      ysrc(:) = y0
      CALL fteik_source_initialize64f(nsrc, zsrc, xsrc, ysrc, verbose, ierr)
      IF (ierr /= 0) WRITE(*,*) 'fteik_solver_setSources64f: Failed to set source'
      IF (ALLOCATED(ysrc)) DEALLOCATE(ysrc)
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
      SUBROUTINE fteik_solver2d_setSphereToCartEpsilon(epsIn, ierr) &
      BIND(C, NAME='fteik_solver2d_setSphereToCartEpsilon')
      USE ISO_C_BINDING
      USE FTEIK_MODEL64F, ONLY : nz, nx
      USE FTEIK_CONSTANTS64F, ONLY : zero
      IMPLICIT NONE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: epsIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      epsS2C = zero
      IF (nx < 1 .OR. nz < 1) THEN
         WRITE(*,*) 'fteik_solver_setSphereToCartEpsilon: Grid not yet set'
         ierr = 1
         RETURN
      ENDIF
      IF (INT(epsIn) > nz .OR. INT(epsIn) > nx) THEN
         IF (INT(epsIn) > nz) THEN
            WRITE(*,*) 'fteik_solver2d_setSphereToCartEpsilon: eps bigger than nz', &
                       INT(epsIn), nz
         ENDIF
         IF (INT(epsIn) > nx) THEN
            WRITE(*,*) 'fteik_solver2d_setSphereToCartEpsilon: eps bigger than nx', &
                       INT(epsIn), nx
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
!>    @brief Sets the velocity model on the solver.
!>
!>    @param[in] ncell    Number of cells in velocity model.  This should be 
!>                        (nz-1)*(nx-1).
!>    @param[in] vel      Velocity model (meters/second) in model cells.  This is
!>                        a vector of dimension [ncell] whose fastest direction is
!>                        z and whose slowest direction is y.
!>
!>    @param[out] ierr    0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver2d_setVelocityModel64f(ncell, vel, ierr) &
      BIND(C, NAME='fteik_solver2d_setVelocityModel64f')
      USE FTEIK_MODEL64F, ONLY : fteik_model_setVelocityModel64f
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: ncell
      REAL(C_DOUBLE), INTENT(IN) :: vel(ncell)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      CALL fteik_model_setVelocityModel64f(ncell, vel, ierr)
      IF (ierr /= 0) &
      WRITE(*,*) 'fteik_solver2d_setVelocityModel64f: Error setting velocity model'
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Compute analytic travel time at the point (i, j) in a homogeneous model.
!>
!>    @param[in] i      iz'th grid point (Fortran numbered).
!>    @param[in] j      ix'th grid point (Fortran numbered).
!>    @param[in] dz     Grid spacing (meters) in z.
!>    @param[in] dx     Grid spacing (meters) in x.
!>    @param[in] zsa    Source offset (meters) in z.
!>    @param[in] xsa    Source offset (meters) in x.
!>    @param[in] szero  Slowness at source (s/m).  
!>
!>    @result The travel time from the source at (zsa,xsa) to the (i,j)'th grid
!>            point given a constant slowness around the source.
!>
!>    @author Keurfon Luu, Mark Noble, Alexandrine Gesret, and Ben Baker.
!>
!>    @version 2
!>
!>    @date January 2018
!>
!>    @copyright MIT
!>
      PURE REAL(C_DOUBLE)                                    &
      FUNCTION fteik_localSolver2d_tAna64fF(i, j, dz, dx,    &
                                            zsa, xsa, szero)
      !$OMP DECLARE SIMD(fteik_localSolver2d_tAna64fF) &
      !$OMP UNIFORM(dz, dx, zsa, xsa, szero) 
      USE ISO_C_BINDING
      IMPLICIT NONE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dz, dx, zsa, xsa, szero
      INTEGER(C_INT), VALUE, INTENT(IN) :: i, j 
      REAL(C_DOUBLE) diffx, diffz
      diffz = (DBLE(i) - zsa)*dz
      diffx = (DBLE(j) - xsa)*dx
      fteik_localSolver2d_tAna64fF = szero*HYPOT(diffx, diffz)
      !t_ana = vzero * ( ( ( dfloat(i) - zsa ) * dz )**2.d0 &
      !                + ( ( dfloat(j) - xsa ) * dx )**2.d0 )**0.5d0
      RETURN
      END FUNCTION
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Compute derivative of analytic travel time and derivative of times
!>           at point (i, j, k) in a homogeneous model.
!>
!>    @param[in] i      iz'th grid point (Fortran numbered).
!>    @param[in] j      ix'th grid point (Fortran numbered).
!>    @param[in] dz     Grid spacing (meters) in z.
!>    @param[in] dx     Grid spacing (meters) in x.
!>    @param[in] zsa    Source offset (meters) in z.
!>    @param[in] xsa    Source offset (meters) in x.
!>    @param[in] szero  Slowness at source (s/m).  
!>
!>    @param[out] t_anad   Derivative of analytic travel time at point (i,j).
!>    @param[out] tzc      Derivative of analytic travel time in z.
!>    @param[out] txc      Derivative of analytic travel time in x.
!>
!>    @author Keurfon Luu, Mark Noble, Alexandrine, Gesret, and Ben Baker.
!>
!>    @version 2
!>
!>    @date January 2018
!>
!>    @copyright MIT
!>
      PURE SUBROUTINE fteik_localSolver2d_tAnaD64fF(t_anad,           &
                                                    tzc, txc, i, j,   &
                                                    dz, dx, zsa, xsa, &
                                                    szero)
      !$OMP DECLARE SIMD(fteik_localSolver2d_tAnaD64fF) &
      !$OMP UNIFORM(dz, dx, zsa, xsa, szero) 
      USE ISO_C_BINDING
      IMPLICIT NONE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dz, dx, zsa, xsa, szero
      INTEGER(C_INT), VALUE, INTENT(IN) :: i, j
      REAL(C_DOUBLE), INTENT(OUT) :: t_anad, tzc, txc
      REAL(C_DOUBLE) d0, d0i_szero, diffx, diffx2, diffz, diffz2, sqrtd0
      REAL(C_DOUBLE), PARAMETER :: zero = 0.d0
      diffz = (DBLE(i) - zsa)*dz
      diffx = (DBLE(j) - xsa)*dx
      diffz2 = diffz*diffz
      diffx2 = diffx*diffx
      d0 = diffz2 + diffx2

      d0i_szero = zero
      t_anad = zero
      IF (d0 > zero) THEN
         sqrtd0 = SQRT(d0)
         t_anad = szero*sqrtd0
         d0i_szero = szero/sqrtd0
      ENDIF
      tzc = d0i_szero*diffz
      txc = d0i_szero*diffx
!     d0 = ( ( dfloat(i) - zsa ) *dz )**2.d0 &
!          + ( ( dfloat(j) - xsa ) *dx )**2.d0
!     t_anad = vzero * (d0**0.5d0)
!     if ( d0 .gt. 0.d0 ) then
!       tzc = ( d0**(-0.5d0) ) * ( dfloat(i) - zsa ) * dz * vzero
!       txc = ( d0**(-0.5d0) ) * ( dfloat(j) - xsa ) * dx * vzero
!     else
!       tzc = 0.d0
!       txc = 0.d0
!     end if
      RETURN
      END SUBROUTINE 



END MODULE
