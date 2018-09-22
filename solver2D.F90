!> @defgroup solver2d 2D Eikonal Solver
!> @brief 2D eikonal equation solver.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE FTEIK2D_SOLVER64F
  USE FTEIK_CONSTANTS64F, ONLY : zero
  USE FTEIK_CONSTANTS64F, ONLY : FTEIK_NATURAL_ORDERING, & 
                                 FTEIK_ZX_ORDERING, &
                                 FTEIK_XZ_ORDERING
  USE FTEIK_MEMORY, ONLY : padLength64F
  USE FTEIK_MODEL64F, ONLY : fteik_model_initializeGeometry, &
                             fteik_model_setVelocityModel64f, &
                             fteik_model_free
  USE FTEIK_RECEIVER64F, ONLY : fteik_receiver_getNumberOfReceivers, &
                                fteik_receiver_getTravelTimes64f,    &
                                fteik_receiver_initialize64f,        &
                                fteik_receiver_free
  USE FTEIK_SOURCE64F, ONLY : fteik_source_initialize64f, &
                              fteik_source_free
  USE FTEIK_RAYS64F, ONLY : fteik_rays_initialize, &
                            fteik_rays_free
  USE ISO_FORTRAN_ENV
  USE ISO_C_BINDING
  IMPLICIT NONE
  !> Holds the travel times (seoncds) at receivers.  This has dimension [nrec x nsrc].
  DOUBLE PRECISION, PROTECTED, ALLOCATABLE, SAVE :: ttimesRec(:)
  !> Holds the travel times (seconds).  This has dimension [ngrd].
  DOUBLE PRECISION, PROTECTED, ALLOCATABLE, SAVE :: ttimes(:)
  !!!!DIR$ ATTRIBUTES ALIGN: 64 :: ttimes
  !> The gradient of the travel times (s/m) in x.  This has dimension [ngrd].
  DOUBLE PRECISION, PROTECTED, ALLOCATABLE, SAVE :: tgradx(:)
  !> The gradient of the travel times (s/m) in z.  This has dimension [ngrd].
  DOUBLE PRECISION, PROTECTED, ALLOCATABLE, SAVE :: tgradz(:)
  
  !> Gauss-Seidel method converges when update is less than tol
  DOUBLE PRECISION, PROTECTED, SAVE :: convTol = zero 
  !> Defines the number of levels.  For 2D that is nz + nx
  INTEGER, PROTECTED, SAVE :: nLevels = 0
  !> Defines the max level size.
  INTEGER, PROTECTED, SAVE :: maxLevelSize = 0
  !> Controls the verbosity
  INTEGER, PROTECTED, SAVE :: verbose = 0
  !> Defines the transition from the spherical to the Cartesian solver solver
  !> during the initialization phase.  This has units of grid points.
  DOUBLE PRECISION, PROTECTED, SAVE :: epsS2C = zero
  !> Defines the number of Gauss-Seidel iterations.
  INTEGER, PROTECTED, SAVE :: nsweep = 1
  !> Flag indicating whether or not the travel times were computed.
  LOGICAL, PROTECTED, SAVE :: lhaveTimes = .FALSE.
  !> Flag indicating whether or not the gradients were computed.
  LOGICAL, PROTECTED, SAVE :: lhaveGradient = .FALSE.
  !> Enforces that this is for 2D only.
  LOGICAL(C_BOOL), PRIVATE, PARAMETER :: lis3d = .FALSE.
  !> Private variables for the local solver
  INTEGER(C_SIZE_T), PARAMETER :: alignment = 64
  !DOUBLE PRECISION, SAVE, PRIVATE :: zsa, xsa, szero
  !INTEGER, SAVE, PRIVATE :: zsi, xsi
  DOUBLE PRECISION, SAVE, PRIVATE, ALLOCATABLE :: slocWork(:), ttvecWork(:), &
                                                  tupdWork(:), ttold(:),     &
                                                  gxWork(:), gzWork(:)
  LOGICAL(C_BOOL), SAVE, PRIVATE, ALLOCATABLE :: lupdWork(:)
  INTEGER, PRIVATE, SAVE :: srcNumber = -1
  ! Label the subroutines/functions
  PUBLIC :: fteik_solver2d_setVelocityModel64f
  PUBLIC :: fteik_solver2d_setNodalVelocityModel64f
  PUBLIC :: fteik_solver2d_setCellVelocityModel64f 
  PUBLIC :: fteik_solver2d_initialize64f
  PUBLIC :: fteik_solver2d_free
  PUBLIC :: fteik_solver2d_solveSourceLSM
  PUBLIC :: fteik_solver2d_solveSourceFSM
  PUBLIC :: fteik_solver2d_getTravelTimeField64f
  PUBLIC :: fteik_solver2d_getTravelTimeFieldCell64f
  PUBLIC :: fteik_solver2d_getTravelTimes64f
  PUBLIC :: fteik_solver2d_setNumberOfSweeps
  PUBLIC :: fteik_solver2d_setReceivers64f
  PUBLIC :: fteik_solver2d_setSources64f
  PUBLIC :: fteik_solver2d_setSphereToCartEpsilon
  PUBLIC :: fteik_solver2d_getNumberOfReceivers
  PUBLIC :: fteik_solver2d_setConvergenceTolerance
  PRIVATE :: fteik_solver2d_isUpdateNode
  PRIVATE :: fteik_solver2d_prefetchTravelTimes64f
  PRIVATE :: fteik_solver2d_prefetchSlowness64f
  PRIVATE :: fteik_solver2d_updateTravelTimes64f
  !PRIVATE :: fteik_localSolver2d_noInit64f
  !PRIVATE :: fteik_localSolver2d_init64f

  ! Fortran only routines
  !PUBLIC :: fteik_localSolver2d_initialize64f

  PRIVATE :: padLength64F
  PRIVATE :: fteik_model_initializeGeometry
  PRIVATE :: fteik_model_setVelocityModel64f
  PRIVATE :: fteik_model_free
  PRIVATE :: fteik_receiver_initialize64f
  PRIVATE :: fteik_receiver_getNumberOfReceivers
  PRIVATE :: fteik_receiver_getTravelTimes64f
  PRIVATE :: fteik_receiver_free
  PRIVATE :: fteik_source_initialize64f
  PRIVATE :: fteik_source_free
  PRIVATE :: solver2d_extractTravelTimes
  CONTAINS
!----------------------------------------------------------------------------------------!
!                                      Begin the Code                                    !
!----------------------------------------------------------------------------------------!
!>    @brief Initializes the 2D eikonal solver.
!>    @param[in] nzIn       Number of grid points in z.  This must be at least 3.
!>    @param[in] nxIn       Number of grid points in x.  This must be at least 3.
!>    @param[in] z0In       z origin (meters).
!>    @param[in] x0In       x origin (meters).
!>    @param[in] nsweepIn   This is the number of sweeps or iterations of the 
!>                          Gauss-Seidel method to be performed.  One is generally
!>                          sufficient.
!>    @param[in] epsIn      This is the radius, in number of grid points, around source
!>                          where the spherical approximation finite difference stencils
!>                          will be used.  This cannot be negative.
!>    @param[in] convTolIn  The Gauss-Seidel sweeping will stop if the updates
!>                          to the travel time field are less than convTol (seconds).
!>    @param[in] verboseIn  Controls the verbosity.  Less than 1 is quiet. 
!>    @param[out] ierr      0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_initialize64f(nzIn, nxIn,           &
                                              z0In, x0In,           &
                                              dzIn, dxIn,           &
                                              nsweepIn, epsIn,      &
                                              convTolIn, verboseIn, &
                                              ierr)                 &
      BIND(C, NAME='fteik_solver2d_initialize64f')
      USE FTEIK_MODEL64F, ONLY : ngrd
      USE FTEIK_CONSTANTS64F, ONLY : one, zero, FALSE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x0In, z0In
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dxIn, dzIn
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: epsIn, convtolIn
      INTEGER(C_INT), VALUE, INTENT(IN) :: nxIn, nzIn
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsweepIn, verboseIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER, PARAMETER :: nyIn = 1
      DOUBLE PRECISION, PARAMETER :: y0In = zero
      DOUBLE PRECISION, PARAMETER :: dyIn = one
      IF (verboseIn > 0) WRITE(OUTPUT_UNIT,800)
      CALL fteik_solver2d_free()
      CALL fteik_solver2d_setVerobosity(verboseIn)
      CALL fteik_solver2d_setConvergenceTolerance(convTolIn)
      CALL fteik_model_initializeGeometry(lis3d,            &
                                          nzIn, nxIn, nyIn, &
                                          z0In, x0In, y0In, &
                                          dzIn, dxIn, dyIn,  &
                                          ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,901)
         RETURN
      ENDIF
      CALL fteik_solver2d_setNumberOfSweeps(nsweepIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,902)
         RETURN
      ENDIF
      CALL fteik_rays_initialize(ierr)
      IF (ierr /= 0) THEN 
         WRITE(ERROR_UNIT,903)
         RETURN
      ENDIF 
      CALL fteik_solver2d_setSphereToCartEpsilon(epsIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,904)
         RETURN
      ENDIF
      ! Compute the number of levels
      nLevels = nxIn + nzIn
      maxLevelSize = MIN(nxIn, nzIn) + padLength64F(alignment, MIN(nxIn, nzIn))
      ! Set space for the travel-times
      IF (.NOT.ALLOCATED(ttimes)) ALLOCATE(ttimes(ngrd))
      IF (.NOT.ALLOCATED(tgradx)) ALLOCATE(tgradx(ngrd))
      IF (.NOT.ALLOCATED(tgradz)) ALLOCATE(tgradz(ngrd))
      ALLOCATE(lupdWork(maxLevelSize))
      ALLOCATE(slocWork(4*maxLevelSize))
      ALLOCATE(ttvecWork(4*maxLevelSize))
      ALLOCATE(tupdWork(maxLevelSize))
      ALLOCATE(gxWork(maxLevelSize))
      ALLOCATE(gzWork(maxLevelSize))
      ALLOCATE(ttold(ngrd))
      lupdWork(:) = FALSE
      slocWork(:) = zero
      ttvecWork(:) = zero
      tupdWork(:) = zero
      gxWork(:) = zero
      gzWork(:) = zero
  800 FORMAT('fteik_solver2d_initialize64f: Initializing...')
  901 FORMAT('fteik_solver2d_initialize64f: Error setting grid geometry')
  902 FORMAT('fteik_solver2d_initialize64f: Failed to set max number of sweeps')
  903 FORMAT('fteik_solver2d_initialize64f: Failed to initialize ray tracer')
  904 FORMAT('fteik_solver2d_initialize64f: Failed to set epsilon')
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the verbosity on the module.
!>    @param[in] verboseIn   Verbosity level to set.  Less than 1 is quiet.
!>    @ingroup solver2d
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
!>    @brief Sets the convergence tolerance for the Gauss-Seidel sweeping.
!>    @param[in] convTolIn   The Gauss-Seidel sweeping will terminate when the updated
!>                           travel time perturbation is less than convTolIn (seconds).
!>                           Note that if convTolin <= 0 then it will be disabled.
!>
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_setConvergenceTolerance(convTolIn) &
      BIND(C, NAME='fteik_solver2d_setConvergenceTolerance')
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: convTolIn
      convTol = convTolIn
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases all memory associated with the solver and resets all variables.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_free()   &
      BIND(C, NAME='fteik_solver2d_free')
      !USE FTEIK_CONSTANTS64F, ONLY : zero
      lhaveTimes = .FALSE.
      lhaveGradient = .FALSE.
      epsS2C = zero
      convTol = zero 
      nsweep = 0
      nLevels = 0
      verbose = 0
      srcNumber =-1
      IF (ALLOCATED(ttimesRec))  DEALLOCATE(ttimesRec)
      IF (ALLOCATED(ttimes))     DEALLOCATE(ttimes)
      IF (ALLOCATED(tgradx))     DEALLOCATE(tgradx)
      IF (ALLOCATED(tgradz))     DEALLOCATE(tgradz)
      IF (ALLOCATED(ttold))      DEALLOCATE(ttold)
      IF (ALLOCATED(lupdWork))   DEALLOCATE(lupdWork)
      IF (ALLOCATED(slocWork))   DEALLOCATE(slocWork)
      IF (ALLOCATED(ttvecWork))  DEALLOCATE(ttvecWork)
      IF (ALLOCATED(tupdWork))   DEALLOCATE(tupdWork)
      IF (ALLOCATED(gxWork))     DEALLOCATE(gxWork)
      IF (ALLOCATED(gzWork))     DEALLOCATE(gzWork)
      CALL fteik_receiver_free()
      CALL fteik_source_free()
      CALL fteik_model_free()
      CALL fteik_rays_free()
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility routine to return the number of receivers.
!>    @param[out] nrec   Number of receivers.
!>    @param[out] ierr   0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_getNumberOfReceivers(nrec, ierr) &
      BIND(C, NAME='fteik_solver2d_getNumberOfReceivers')
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(OUT) :: nrec, ierr
      ierr = 0 
      CALL fteik_receiver_getNumberOfReceivers(nrec, ierr)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the ray paths to receivers from the current source in the
!>           corresponding travel time field.
!>    @param[out] ierr   0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_computeRaysToReceivers(ierr) &
      BIND(C, NAME='fteik_solver2d_computeRaysToReceivers')
      USE FTEIK_RECEIVER64F, ONLY : nrec, xrec, zrec
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      IF (nrec < 1) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      CALL fteik_solver2d_computeRaysToPoints(nrec, xrec, zrec, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_uNIT,901)
         ierr = 1
         RETURN
      ENDIF
  900 FORMAT('fteik_solver2d_computeRaysToReceivers: No receivers set')
  901 FORMAT('fteik_solver2d_computeRaysToReceivers: Failed to compute rays')
      RETURN
      END
!>    @brief Computes the ray paths from the current source to the given (xr,zr)
!>           points in the travel time field.
!>    @param[in] np     Number of points.
!>    @param[in] xp     x destinations of rays.  This has dimension [nr].
!>    @param[in] zp     z destinations of rays.  This has dimension [nr].
!>    @param[out] ierr  0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_computeRaysToPoints(np, xp, zp, ierr) &
      BIND(C, NAME='fteik_solver2d_computeRaysToPoints')
      USE FTEIK_MODEL64F, ONLY : ngrd
      USE FTEIK_SOURCE64F, ONLY : xstrue, zstrue, xsrc, zsrc
      USE FTEIK_RAYS64F, ONLY : fteik_rays_setRayDestinations2D, &
                                fteik_rays_setSource, fteik_rays_trace
      USE FTEIK_LOCALSOLVER2D64F, ONLY : xsi, zsi, zsa, xsa, szero
      INTEGER(C_INT), VALUE, INTENT(IN) :: np
      REAL(C_DOUBLE), DIMENSION(np), INTENT(IN) :: xp, zp
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      IF (.NOT.lhaveTimes) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      ! Set the source information
      CALL fteik_rays_setSource(zstrue(srcNumber), xstrue(srcNumber), 0.d0, &
                                zsrc(srcNumber),   xsrc(srcNumber),   0.d0, & 
                                zsa, xsa, 0.d0,                             &
                                zsi, xsi, 0,                                &
                                szero)
      ! Set the ray destinations 
      CALL fteik_rays_setRayDestinations2D(np, xp, zp, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,901)
         ierr = 1
         RETURN
      ENDIF
      ! Compute the raypaths
      CALL fteik_rays_trace(ngrd, ttimes, tgradx, tgradz, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,902)
         ierr = 1
         RETURN
      ENDIF
  900 FORMAT('fteik_solver2d_computeRaysToPoints: Travel time field not yet computed')
  901 FORMAT('fteik_solver2d_computeRaysToPoints: Failed to set ray origins')
  902 FORMAT('fteik_solver2d_computeRaysToPoints: Failed to compute rays')
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
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_isUpdateNode(linit, level, sweep, nz, nx, &
                                             i1, i2, zsi, xsi, lupd)
      USE FTEIK_CONSTANTS64F, ONLY : TRUE, FALSE
      LOGICAL, INTENT(IN) :: linit
      INTEGER, INTENT(IN) :: i1, i2, level, nx, nz, sweep, xsi, zsi
      LOGICAL(C_BOOL), DIMENSION(:), INTENT(OUT) :: lupd
      INTEGER i, ix, iz
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
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_prefetchSlowness64f(level, sweep, nz, nx, &
                                                    i1, i2,               &
                                                    lupd, slow, sloc)
      USE FTEIK_CONSTANTS64F, ONLY : zero
      INTEGER, INTENT(IN) :: i1, i2, level, nx, nz, sweep 
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: slow
      LOGICAL(C_BOOL), DIMENSION(:), INTENT(IN) :: lupd
      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: sloc
      INTEGER i, ix, ixcell, indx1, indx2, indx3, indx4, indx5, iz, izcell
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
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_prefetchTravelTimes64f(level, sweep, nz, nx, &
                                                       i1, i2,               &
                                                       lupd, tt, ttvec)
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE
      INTEGER, INTENT(IN) :: i1, i2, level, nx, nz, sweep
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tt
      LOGICAL(C_BOOL), DIMENSION(:), INTENT(IN) :: lupd
      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: ttvec
      DOUBLE PRECISION t, te, tev, tv
      INTEGER i, ix, iz
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
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_updateTravelTimes64f(level, sweep, nz, nx, &
                                                     i1, i2, lupd,         &
                                                     tupd, tt)
      INTEGER, INTENT(IN) :: i1, i2, level, nx, nz, sweep
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tupd
      LOGICAL(C_BOOL), DIMENSION(:), INTENT(IN) :: lupd
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: tt

      INTEGER i, indx, ix, iz
      IF (sweep == 1) THEN
         !$OMP SIMD
         DO i=i1,i2
            ix = i
            iz = level - i 
            indx = (ix - 1)*nz + iz
            !print *, level, ix, iz, indx, sngl(tt(indx)), lupd(i+1-i1)
            IF (lupd(i+1-i1)) tt(indx) = MIN(tt(indx), tupd(i+1-i1))
         ENDDO
      ELSEIF (sweep == 2) THEN
         !$OMP SIMD
         DO i=i1,i2
            ix = nx + 1 - i  !flip ix
            iz = level - i
            indx = (ix - 1)*nz + iz
            IF (lupd(i+1-i1)) tt(indx) = MIN(tt(indx), tupd(i+1-i1))
         ENDDO
      ELSEIF (sweep == 3) THEN
         !$OMP SIMD
         DO i=i1,i2
            ix = i 
            iz = nz + 1 - (level - i) ! flip iz
            indx = (ix - 1)*nz + iz
            IF (lupd(i+1-i1)) tt(indx) = MIN(tt(indx), tupd(i+1-i1))
         ENDDO
      ELSE
         !$OMP SIMD
         DO i=i1,i2
            ix = nx + 1 - i           ! flip ix
            iz = nz + 1 - (level - i) ! flip iz
            indx = (ix - 1)*nz + iz
            IF (lupd(i+1-i1)) tt(indx) = MIN(tt(indx), tupd(i+1-i1)) 
         ENDDO
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Convenience function to solve all travel time fields for all sources with
!>           the level-set fast sweeping method.
!>    @param[out] ierr  0 indicates success.
!>    @ingroup solve2d
      SUBROUTINE fteik_solver2d_solveLSM(ierr)    &
      BIND(C, NAME='fteik_solver2d_solveLSM')
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE
      USE FTEIK_MODEL64F, ONLY : lhaveModel
      USE FTEIK_SOURCE64F, ONLY : nsrc
      USE FTEIK_RECEIVER64F, ONLY : nrec
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) isrc
      ierr = 0
      IF (.NOT.lhaveModel) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      IF (nsrc < 1) THEN
         IF (verbose > 0) WRITE(ERROR_UNIT,901)
         RETURN
      ENDIF
      IF (nsrc > 0 .AND. nrec > 0) THEN
         IF (.NOT.ALLOCATED(ttimesRec)) ALLOCATE(ttimesRec(nsrc*nrec))
         ttimesRec(:) = FTEIK_HUGE
      ENDIF
      DO isrc=1,nsrc
         CALL fteik_solver2d_solveSourceLSM(isrc, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,902) isrc
            RETURN 
         ENDIF
      ENDDO
  900 FORMAT('fteik_solver2d_solveLSM: Model not yet set')
  901 FORMAT('fteik_solver2d_solveLSM: No sources')
  902 FORMAT('fteik_solver2d_solveLSM: Error solving source ', I0)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Solves the eikonal equation with the level-set fast sweeping method.
!>    @param[in] isrc   The source number for which the eikonal equation will be solved.
!>                      This must be in the range [1,nsrc].
!>    @param[out] ierr  0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_solveSourceLSM(isrc, ierr)   &
                 BIND(C, NAME='fteik_solver2d_solveSourceLSM')
      USE FTEIK_SOURCE64F, ONLY : nsrc
      USE FTEIK_MODEL64F, ONLY : lhaveModel, ngrd, nx, nz, slow
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, TRUE, FALSE
      USE FTEIK_LOCALSOLVER2D64F, ONLY : xsi, zsi
      USE FTEIK_LOCALSOLVER2D64F, ONLY : fteik_localSolver2d_initialize64f, &
                                         fteik_localSolver2d_init64f,       &
                                         fteik_localSolver2d_noInit64f,     &
                                         fteik_localSolver2d_getSweepSigns
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION ts4(4), t0, t1, tdiff
      INTEGER dest(4), i, i1, i2, indx, kiter, level, nnodes, sweep
      INTEGER iz, ix, sgntz, sgntx, sgnvz, sgnvx
      ierr = 0
      lhaveTimes = .FALSE.
      lhaveGradient = .FALSE.
      srcNumber =-1
      IF (.NOT.lhaveModel) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         GOTO 500
      ENDIF
      IF (isrc < 1 .OR. isrc > nsrc) THEN
         WRITE(ERROR_UNIT,901) isrc, nsrc
         ierr = 1
         GOTO 500
      ENDIF
      CALL fteik_localSolver2d_initialize64f(isrc, epsS2C, dest, ts4, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,902)
         GOTO 500
      ENDIF
      t0 = 0.d0
      CALL CPU_TIME(t0)
      ttimes(:) = FTEIK_HUGE
      tgradx(:) = 0.d0
      tgradz(:) = 0.d0 
      ttimes(dest(1:4)) = ts4(1:4)
!print *, dest
!print *, ts4
!return
      ! Initialization
      ! First sweeping: Top->Bottom ; West->East
      CALL fteik_localSolver2d_getSweepSigns(1, sgntz, sgntx, sgnvz, sgnvx, ierr)
      DO ix=MAX(2,xsi),nx
         DO iz=MAX(2,zsi),nz
            CALL fteik_localSolver2d_init64f(sgntz, sgntx,   &
                                             sgnvz, sgnvx,   &
                                             nz, nx, iz, ix, &
                                             slow, ttimes, tupdWork(1))
            ttimes((ix-1)*nz + iz) = MIN(ttimes((ix-1)*nz + iz), tupdWork(1))
         ENDDO
      ENDDO
      ! Second sweeping: Top->Bottom ; East->West
      CALL fteik_localSolver2d_getSweepSigns(2, sgntz, sgntx, sgnvz, sgnvx, ierr)
      DO ix=xsi+1,1,-1
         DO iz=MAX(2,zsi),nz
            CALL fteik_localSolver2d_init64f(sgntz, sgntx,   &
                                             sgnvz, sgnvx,   &
                                             nz, nx, iz, ix, &
                                             slow, ttimes, tupdWork(1))
            ttimes((ix-1)*nz + iz) = MIN(ttimes((ix-1)*nz + iz), tupdWork(1))
         ENDDO
      ENDDO
      ! Third sweep: Bottom->Top ; West->East
      CALL fteik_localSolver2d_getSweepSigns(3, sgntz, sgntx, sgnvz, sgnvx, ierr)
      DO ix=MAX(2,xsi),nx
         DO iz=zsi+1,1,-1
            CALL fteik_localSolver2d_init64f(sgntz, sgntx,   &
                                             sgnvz, sgnvx,   &
                                             nz, nx, iz, ix, &
                                             slow, ttimes, tupdWork(1))
            ttimes((ix-1)*nz + iz) = MIN(ttimes((ix-1)*nz + iz), tupdWork(1))
         ENDDO
      ENDDO
      ! Fourth sweeping: Bottom->Top ; East->West
      CALL fteik_localSolver2d_getSweepSigns(4, sgntz, sgntx, sgnvz, sgnvx, ierr)
      DO ix=xsi+1,1,-1
         DO iz=zsi+1,1,-1
            CALL fteik_localSolver2d_init64f(sgntz, sgntx,   &
                                             sgnvz, sgnvx,   &
                                             nz, nx, iz, ix, &
                                             slow, ttimes, tupdWork(1))
            ttimes((ix-1)*nz + iz) = MIN(ttimes((ix-1)*nz + iz), tupdWork(1))
         ENDDO
      ENDDO
!print *, 'hey', minval(ttimes), maxval(ttimes)
!     ! Initialize
!     DO sweep=1,4
!        ! Loop on the levels
!        DO level=2,nlevels
!           i1 = MAX(1, level - nz)
!           i2 = MIN(nx, level - 1) 
!           nnodes = i2 - i1 + 1
!           ! Get the update nodes
!           CALL fteik_solver2d_isUpdateNode(.TRUE., level, sweep, nz, nx, &
!                                            i1, i2, zsi, xsi, lupdWork)
!           ! Prefetch the slowness
!           CALL fteik_solver2d_prefetchSlowness64f(level, sweep, nz, nx,   &
!                                                   i1, i2,                 &
!                                                   lupdWork, slow, slocWork)
!           ! Prefetch travel times
!           CALL fteik_solver2d_prefetchTravelTimes64f(level, sweep, nz, nx,      &
!                                                      i1, i2,                    &
!                                                      lupdWork, ttimes, ttvecWork)
!           ! Compute the candidate travel times at each node in level set
!           !CALL fteik_localSolver2d_noInit64f(nnodes, ttvecWork, slocWork, tupdWork)
!           ! And update
!           !CALL fteik_solver2d_updateTravelTimes64f(level, sweep, nz, nx, &
!           !                                         i1, i2, lupdWork,     &
!           !                                         tupdWork, ttimes)
!        ENDDO
!     ENDDO
      IF (convTol > zero) ttold(1:ngrd) = ttimes(1:ngrd)
      ! Number of Gauss-Seidel iterations
      DO kiter=1,nsweep
         ! Loop on sweeping directions
         DO sweep=1,4
            CALL fteik_localSolver2d_getSweepSigns(sweep, sgntz, sgntx, sgnvz, sgnvx, &
                                                   ierr)
            ! Loop on the levels
            DO level=2,nlevels
               i1 = MAX(1, level - nz)
               i2 = MIN(nx, level - 1)
               nnodes = i2 - i1 + 1
               ! Get the update nodes
               CALL fteik_solver2d_isUpdateNode(.FALSE., level, sweep, nz, nx, &
                                                i1, i2, zsi, xsi, lupdWork)
               ! Prefetch the slowness
               CALL fteik_solver2d_prefetchSlowness64f(level, sweep, nz, nx,   &
                                                       i1, i2,                 &
                                                       lupdWork, slow, slocWork)
               ! Prefetch travel times
               CALL fteik_solver2d_prefetchTravelTimes64f(level, sweep, nz, nx,      &
                                                          i1, i2,                    &
                                                          lupdWork, ttimes, ttvecWork)
               ! Compute the candidate travel times at each node in level set
               CALL fteik_localSolver2d_noInit64f(nnodes, ttvecWork,  &
                                                  slocWork, tupdWork)
               ! And update
               CALL fteik_solver2d_updateTravelTimes64f(level, sweep, nz, nx,    &
                                                        i1, i2, lupdWork,        &
                                                        tupdWork, ttimes)
            ENDDO
            !do i2=1,nx*nz
            !   print *, i2, ttimes(i2)
            !enddo
            !print *, minval(ttimes), maxval(ttimes)
            IF (verbose > 3) &
            WRITE(OUTPUT_UNIT,801) kiter, sweep, MINVAL(ttimes), MAXVAL(ttimes)
         ENDDO
         IF (verbose > 2) WRITE(OUTPUT_UNIT,800) kiter, MINVAL(ttimes), MAXVAL(ttimes)
         ! Check convergence
         tdiff = FTEIK_HUGE
         IF (nsweep > 1 .AND. convTol > zero) THEN
            !$OMP SIMD REDUCTION(MAX:tdiff)
            DO i=1,ngrd
               tdiff = MAX(tdiff, ABS(ttimes(i) - ttold(i)))
               ttold(i) = ttimes(i)
            ENDDO
            IF (verbose > 3) WRITE(OUTPUT_UNIT,803) tdiff 
         ENDIF
         IF (tdiff < convTol) THEN
            EXIT
         ENDIF
      ENDDO
      lhaveTimes = .TRUE.
      srcNumber = isrc
      ! Differenate the travel times
      CALL fteik_solver2d_finiteDifferenceGradient()
      lhaveGradient = .TRUE.
      ! Extract the travel time solutions
      CALL solver2d_extractTravelTimes(isrc, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,903)
      ENDIF 
      IF (verbose > 2) THEN
         CALL CPU_TIME(t1)
         WRITE(OUTPUT_UNIT,802) t1 - t0 
      ENDIF
  500 CONTINUE
      IF (ierr /= 0) THEN
         ttimes(:) = FTEIK_HUGE
         lhaveTimes = .FALSE.
         lhaveGradient = .FALSE.
         srcNumber =-1
      ENDIF
  800 FORMAT('fteik_solver2d_solveSourceLSM: (iteration,ttmin,ttmax)=', &
             I4, 2F14.8) 
  801 FORMAT('fteik_solver2d_solveSourceLSM: (iteration,sweep,ttmin,ttmax)=', &
             2I4, 2F14.8)
  802 FORMAT('fteik_solver2d_solveSourceLSM: Solver time in seconds=', F14.8)
  803 FORMAT('fteik_solver2d_solveSourceLSM: Max perturbation=', F16.10)
  900 FORMAT('fteik_solver2d_solveSourceLSM: Model not yet set')
  901 FORMAT('fteik_solver2d_solveSourceLSM: Source number =',I0, &
             ' must be in range [1,',I0,']')
  902 FORMAT('fteik_solver2d_solveSourceLSM: Error setting mesh constants')
  903 FORMAT('fteik_solver2d_solveSourceLSM: Failed to get travel times')
!print *, lhaveTimes
      RETURN
      END

      SUBROUTINE fteik_solver2d_finiteDifferenceGradient()
      USE FTEIK_LOCALSOLVER2D64F, ONLY : xsi, zsi
      USE FTEIK_MODEL64F, ONLY : dx, dz, nx, nz
      INTEGER ix, iz, m1, m2, m3
      ! Do the heavy lifting and forward difference the field 
      DO ix=1,nx
         DO iz=2,nz-1
            m1 = (ix-1)*nz + iz - 1
            m2 = (ix-1)*nz + iz
            m3 = (ix-1)*nz + iz + 1
            !tgradz(m2) = (ttimes(m3) - ttimes(m1))/(2.d0*dz)
            tgradz(m2) = (ttimes(m2) - ttimes(m1))/dz
         ENDDO
      ENDDO
      DO ix=2,nx-1
         DO iz=1,nz
            m1 = (ix-2)*nz + iz
            m2 = (ix-1)*nz + iz
            m3 = (ix-0)*nz + iz
            !tgradx(m2) = (ttimes(m3) - ttimes(m1))/(2.d0*dx)
            tgradx(m2) = (ttimes(m2) - ttimes(m1))/dx
         ENDDO
      ENDDO
      ! Top/bottom boundaries
      DO ix=1,nx
         ! Forward difference top boundary
         m1 = (ix-1)*nz + 1     ! home
         m2 = (ix-1)*nz + 2 ! ahead
         tgradz(m1) = (ttimes(m2) - ttimes(m1))/dz
         ! Backwards difference bottom boundary
         m1 = (ix-1)*nz + nz - 1 ! behind
         m2 = (ix-1)*nz + nz     ! home 
         tgradz(m2) = (ttimes(m2) - ttimes(m1))/dz
      ENDDO
      ! Left/right boundaries
      DO iz=1,nz
         ! Forward difference left
         m1 = (1-1)*nz + iz ! home
         m2 = (2-1)*nz + iz ! ahead
         tgradx(m1) = (ttimes(m2) - ttimes(m1))/dx
         ! Backwards difference the righ boundary
         m1 = (nx-2)*nz + iz ! behind
         m2 = (nx-1)*nz + iz ! home 
         tgradx(m2) = (ttimes(m2) - ttimes(m1))/dx
      ENDDO
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Solves the eikonal equation with the fast sweeping method.
!>    @param[in] isrc   The source number for which the eikonal equation will be solved.
!>                      This must be in the range [1,nsrc].
!>    @param[out] ierr  0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_solveSourceFSM(isrc, ierr) &
      BIND(C, NAME='fteik_solver2d_solveSourceFSM')
      USE FTEIK_SOURCE64F, ONLY : nsrc
      USE FTEIK_MODEL64F, ONLY : lhaveModel, nx, nz, slow
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, TRUE, FALSE
      USE FTEIK_LOCALSOLVER2D64F, ONLY : xsi, zsi
      USE FTEIK_LOCALSOLVER2D64F, ONLY : fteik_localSolver2d_initialize64f, &
                                         fteik_localSolver2d_init64f,       &
                                         fteik_localSolver2d_noInit64f
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION ts4(4), t0
      INTEGER dest(4), i, i1, ix, iz, j, j1, kiter, level, sgntx, sgntz, sgnvx, sgnvz 
      INTEGER indx1, indx2, indx3, indx4, indx5
      ierr = 0 
      lhaveTimes = .FALSE.
      lhaveGradient = .FALSE.
      srcNumber =-1
      IF (.NOT.lhaveModel) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         GOTO 500
      ENDIF
      IF (isrc < 1 .OR. isrc > nsrc) THEN
         WRITE(ERROR_UNIT,901) isrc, nsrc
         ierr = 1
         GOTO 500
      ENDIF
      CALL fteik_localSolver2d_initialize64f(isrc, epsS2C, dest, ts4, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,902)
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
            CALL fteik_localSolver2d_init64f(sgntz, sgntx, &
                                             sgnvz, sgnvx,  &
                                             nz, nx, iz, ix, &
                                             slow, ttimes, tupdWork(1))
!print *, iz, ix, ttimes((ix-1)*nz + iz), tupdWork(1)
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
            CALL fteik_localSolver2d_init64f(sgntz, sgntx,   &
                                             sgnvz, sgnvx,   &
                                             nz, nx, iz, ix, &
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
            CALL fteik_localSolver2d_init64f(sgntz, sgntx,   &
                                             sgnvz, sgnvx,   &
                                             nz, nx, iz, ix, &
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
            CALL fteik_localSolver2d_init64f(sgntz, sgntx,   &
                                             sgnvz, sgnvx,   &
                                             nz, nx, iz, ix, &
                                             slow, ttimes, tupdWork(1))
            ttimes((ix-1)*nz + iz) = MIN(ttimes((ix-1)*nz + iz), tupdWork(1))
         ENDDO
      ENDDO 
!print *, 'begin', minval(ttimes), maxval(ttimes)
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
               CALL fteik_localSolver2d_noInit64f(1, ttvecWork, slocWork, tupdWork)
               ttimes((j-1)*nz + i) = MIN(ttimes((j-1)*nz + i), tupdWork(1))
            ENDDO
         ENDDO
!do i=1,nx*nz
!print *, 'ref', i, ttimes(i)
!enddo
!print *, 'p1:', minval(ttimes), maxval(ttimes)
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
!              CALL fteik_solver2d_prefetchTravelTimes64f(level, 2, nz, nx,  &
!                                                         1, ix, ix,         &
!                                                         lupd, ttimes, ttvec)
!              CALL fteik_solver2d_prefetchSlowness64f(level, 2, nz, nx, &
!                                                      1, ix, ix,        &
!                                                      lupd, slow, sloc)
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
               CALL fteik_localSolver2d_noInit64f(1, ttvecWork, slocWork, tupdWork)
               ttimes((j-1)*nz + i) = MIN(ttimes((j-1)*nz + i), tupdWork(1))
            ENDDO
         ENDDO
!print *, 'p2:', minval(ttimes), maxval(ttimes)
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
               CALL fteik_localSolver2d_noInit64f(1, ttvecWork, slocWork, tupdWork)
               ttimes((j-1)*nz + i) = MIN(ttimes((j-1)*nz + i), tupdWork(1))
            ENDDO
         ENDDO
!print *, 'p3:', minval(ttimes), maxval(ttimes)
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
               CALL fteik_localSolver2d_noInit64f(1, ttvecWork, slocWork, tupdWork)
               ttimes((j-1)*nz + i) = MIN(ttimes((j-1)*nz + i), tupdWork(1))
            ENDDO
         ENDDO
!print *, 'p4:', minval(ttimes), maxval(ttimes)
      ENDDO
      lhaveTimes = .TRUE.
      srcNumber =-1
  500 CONTINUE
      IF (ierr /= 0) THEN
         ttimes(:) = FTEIK_HUGE
         lhaveTimes = .FALSE.
      ENDIF
  900 FORMAT('fteik_solver2d_solveSourceFSM: Model not yet set')
  901 FORMAT('fteik_solver2d_solveSourceFSM: Source number =',I0, &
             ' must be in range [1,',I0,']')
  902 FORMAT('fteik_solver2d_solveSourceFSM: Error setting mesh constants')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Returns the nodal travel time field.
!>    @param[in] ngin    Number of input grid points.  This must equal [nz x nx].
!>    @param[in] order   Desired ordering of output.
!>    @param[in] order   If order == FTEIK_NATURAL_ORDERING or FTEIK_ZX_ORDERING then
!>                       ttimes will be [nz x nx] with leading dimension nz.
!>    @param[in] order   If order == FTEIK_XZ_ORDERING then ttimes will be [nx x nz]
!>                       with leading dimension nx.
!>    @param[out] ttout  Travel time field at the grid points.  This is a [nz x nx]
!>                       vector with leading dimension nz.
!>    @param[out] ierr   0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_getTravelTimeField64f(ngin, order, ttout, ierr) &
      BIND(C, NAME='fteik_solver2d_getTravelTimeField64f')
      USE FTEIK_MODEL64F, ONLY : ngrd, nx, nz
      INTEGER(C_INT), VALUE, INTENT(IN) :: ngin, order
      REAL(C_DOUBLE), INTENT(OUT) :: ttout(ngin)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER indx, ix, iz, jndx
      ierr = 0
      IF (.NOT.lhaveTimes) THEN 
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      IF (ngin /= ngrd) THEN
         WRITE(ERROR_UNIT,901) ngin, ngrd
         ierr = 1
         RETURN
      ENDIF
      IF (order == FTEIK_XZ_ORDERING) THEN
         DO iz=1,nz
            !$OMP SIMD
            DO ix=1,nx
               indx = (ix - 1)*nz + iz
               jndx = (iz - 1)*nx + ix
               ttout(jndx) = ttimes(indx)
            ENDDO
         ENDDO
      ELSE
         IF (order /= FTEIK_NATURAL_ORDERING .AND. order /= FTEIK_ZX_ORDERING) THEN
            WRITE(OUTPUT_UNIT,800)
         ENDIF
         ttout(:) = ttimes(:)
      ENDIF
  800 FORMAT('fteik_solver2d_getTravelTimeField64f: Defaulting to natural order')
  900 FORMAT('fteik_solver2d_getTravelTimeField64f: Travel times not yet computed')
  901 FORMAT('fteik_solver2d_getTravelTimeField64f: ngin=', I0, ' /= ngrd =', I0)
      RETURN
      END
!>    @brief Returns the cell-based travel time field.  The values are computed by
!>           averaging the nodal values of the travel time field.
!>    @param[in] ngin    Number of input grid points.  This must equal [nz-1 x nx-1].
!>    @param[in] order   Desired ordering of output.
!>    @param[in] order   If order == FTEIK_NATURAL_ORDERING or FTEIK_ZX_ORDERING then
!>                       ttimes will be [nz-1 x nx-1] with leading dimension nz-1.
!>    @param[in] order   If order == FTEIK_XZ_ORDERING then ttimes will be [nx-1 x nz-1]
!>                       with leading dimension nx-1.
!>    @param[out] ttout  Travel time field at the grid points.  This is a [nz-1 x nx-1]
!>                       vector with leading dimension nz-1.
!>    @param[out] ierr   0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_getTravelTimeFieldCell64f(ncin, order, ttout, ierr) &
      BIND(C, NAME='fteik_solver2d_getTravelTimeFieldCell64f')
      USE FTEIK_MODEL64F, ONLY : ncell, nx, nz
      INTEGER(C_INT), VALUE, INTENT(IN) :: ncin, order
      REAL(C_DOUBLE), DIMENSION(ncin), INTENT(OUT) :: ttout
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER indx(4), ix, iz, jndx
      ierr = 0
      IF (.NOT.lhaveTimes) THEN 
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      IF (ncin /= ncell) THEN 
         WRITE(ERROR_UNIT,901) ncin, ncell
         ierr = 1
         RETURN
      ENDIF
      IF (order == FTEIK_XZ_ORDERING) THEN 
         DO iz=1,nz-1
            !$OMP SIMD
            DO ix=1,nx-1
               indx(1) = (ix - 1)*nz + iz 
               indx(2) = (ix - 1)*nz + iz + 1
               indx(3) = (ix - 0)*nz + iz
               indx(4) = (ix - 0)*nz + iz + 1
               jndx = (iz - 1)*(nx - 1) + ix 
               ttout(jndx) = 0.25d0*SUM(ttimes(indx(1:4))) ! Compute average
            ENDDO
         ENDDO
      ELSE 
         IF (order /= FTEIK_NATURAL_ORDERING .AND. order /= FTEIK_ZX_ORDERING) THEN 
            WRITE(OUTPUT_UNIT,800)
         ENDIF
         ttout(:) = ttimes(:)
      ENDIF
  800 FORMAT('fteik_solver2d_getTravelTimeFieldCell64f: Defaulting to natural order')
  900 FORMAT('fteik_solver2d_getTravelTimeFieldCell64f: Travel times not yet computed')
  901 FORMAT('fteik_solver2d_getTravelTimeFieldCell64f: ncin=', I0, ' /= ncell =', I0)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Extracts the travel times at the receivers.
!>
!>    @param[in] ldr    Leading dimension of ttr.  This must >= nrec.
!>
!>    @param[out] ttr   Travel times (seconds) at the receivers for all sources.
!>                      This is a vector of dimension [ldr x nsrc] with leading
!>                      dimension ldr.
!>    @param[out] ierr  0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_getTravelTimes64f(ldr, ttr, ierr) &
      BIND(C, NAME='fteik_solver2d_getTravelTimes64f')
      USE FTEIK_SOURCE64F, ONLY : nsrc
      USE FTEIK_RECEIVER64F, ONLY : nrec 
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldr
      REAL(C_DOUBLE), INTENT(OUT) :: ttr(ldr*nsrc)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i1, i2, isrc, j1, j2
      ierr = 0
      IF (nrec < 1) THEN
         WRITE(OUTPUT_UNIT,800)
         RETURN
      ENDIF
      IF (nsrc < 1) THEN
         WRITE(ERROR_UNIT,901)
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.lhaveTimes) THEN
         WRITE(ERROR_UNIT,902)
         ierr = 1
         RETURN
      ENDIF
      IF (ldr < nrec) THEN
         WRITE(ERROR_UNIT,903) ldr, nrec
         ierr = 1
         RETURN
      ENDIF
      IF (ldr > nrec) ttr(1:ldr*nsrc) = 0.d0
      DO isrc=1,nsrc
         i1 = (isrc - 1)*ldr + 1
         i2 = i1 + nrec - 1
         j1 = (isrc - 1)*nrec + 1
         j2 = j1 + nrec - 1
         ttr(i1:i2) = ttimesRec(j1:j2)
         !print *, ttr(i1:i2)
      ENDDO
!     CALL fteik_receiver_getTravelTimes64f(nrec, ngrd, ttimes, ttr, ierr)
!     IF (ierr /= 0) THEN
!        WRITE(*,*) 'fteik_solver2d_getTravelTimes64f: Error getting travel times'
!        ierr = 1
!     ENDIF
  800 FORMAT('fteik_solver2d_getTravelTimes64f: No receivers')
  901 FORMAT('fteik_solver2d_getTravelTimes64f: No sources')
  902 FORMAT('fteik_solver2d_getTravelTimes64f: Travel times not yet computed')
  903 FORMAT('fteik_solver2d_getTravelTimes64f: ldr = ', I0, ' < nrec = ', I0)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE solver2d_extractTravelTimes(isrc, ierr)
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE
      USE FTEIK_MODEL64F, ONLY : ngrd
      USE FTEIK_RECEIVER64F, ONLY : nrec
      USE FTEIK_SOURCE64F, ONLY : nsrc
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i1, i2
      ierr = 0
      IF (nrec < 1) RETURN
      IF (nsrc < 1) RETURN
      IF (.NOT.lhaveTimes) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      i1 = (isrc - 1)*nrec + 1
      i2 = i1 + nrec - 1
      CALL fteik_receiver_getTravelTimes64f(nrec, ngrd, ttimes, ttimesRec(i1:i2), ierr) 
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,901)
         ttimesRec(i1:i2) = FTEIK_HUGE
      ENDIF
  900 FORMAT('solver2d_extractTravelTimes: Travel times not yet computed')
  901 FORMAT('solver2d_extractTravelTimes: Failed to get travel times')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the number of Gauss-Seidel iterations  in the fast sweeping method.
!>
!>    @param[in] nsweepIn    Number of sweeps.  This cannot be negative and 1 is
!>                           usually sufficient.
!>    @param[out] ierr       0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_setNumberOfSweeps(nsweepIn, ierr) &
      BIND(C, NAME='fteik_solver2d_setNumberOfSweeps')
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsweepIn 
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      nsweep = 1
      IF (nsweepIn < 0) THEN
         WRITE(ERROR_UNIT,900) nsweep
         ierr = 1
         RETURN
      ENDIF
      nsweep = nsweepIn
  900 FORMAT('fteik_solver2d_setNumberOfSweeps: nsweep=',I0,' must be positive')
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the receivers in the model.
!>
!>    @param[in] nrecIn  Number of receivers to set.  
!>    @param[in] zrec    z locations (meters) of receivers.  This is a vector of
!>                       of dimension [nrecIn].
!>    @param[in] xrec    x locations (meters) of receivers.  This is a vector of
!>                       of dimension [nrecIn].
!>    @param[out] ierr   0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_setReceivers64f(nrecIn, zrec, xrec, ierr) &
      BIND(C, NAME='fteik_solver2d_setReceivers64f')
      USE FTEIK_RECEIVER64F, ONLY : nrec
      USE FTEIK_SOURCE64F, ONLY : nsrc
      USE FTEIK_MODEL64F, ONLY : y0
      INTEGER(C_INT), VALUE, INTENT(IN) :: nrecIn
      REAL(C_DOUBLE), INTENT(IN) :: zrec(nrecIn), xrec(nrecIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE :: yrec(:)
      ! It's actually safe to have no receivers
      ierr = 0
      IF (nrecIn < 1) THEN
         WRITE(OUTPUT_UNIT,800)
         RETURN
      ENDIF
      ALLOCATE(yrec(MAX(nrecIn, 1)))
      yrec(:) = y0
      CALL fteik_receiver_initialize64f(nrecIn, zrec, xrec, yrec, verbose, ierr)
      IF (ierr /= 0) WRITE(ERROR_UNIT,901)
      IF (ALLOCATED(yrec)) DEALLOCATE(yrec)
      ! May be ready to allocate space for travel times 
      IF (ALLOCATED(ttimesRec)) DEALLOCATE(ttimesRec)
      IF (nsrc > 0) THEN
         IF (ALLOCATED(ttimesRec)) DEALLOCATE(ttimesRec)
         ALLOCATE(ttimesRec(nrec*nsrc))
      ENDIF
  800 FORMAT('fteik_solver2d_setReceivers64f: No receivers to set')
  901 FORMAT('fteik_solver2d_setReceivers64f: Failed to set receivers')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the source(s) on the solver.
!> 
!>    @param[in] nsrc     Number of sources.
!>    @param[in] zsrc     z locations (meters) of source.  This is a vector of dimension
!>                        [nsrc].
!>    @param[in] xsrc     x locations (meters) of source.  This is a vector of dimension
!>                        [nsrc].
!>    @param[out] ierr    0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_setSources64f(nsrc, zsrc, xsrc, ierr) &
      BIND(C, NAME='fteik_solver2d_setSources64f')
      USE FTEIK_MODEL64F, ONLY : y0
      USE FTEIK_RECEIVER64F, ONLY : nrec
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsrc
      REAL(C_DOUBLE), INTENT(IN) :: zsrc(nsrc), xsrc(nsrc)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE :: ysrc(:)
      ALLOCATE(ysrc(MAX(1, nsrc)))
      ysrc(:) = y0
      CALL fteik_source_initialize64f(nsrc, zsrc, xsrc, ysrc, verbose, ierr)
      IF (ierr /= 0) WRITE(ERROR_UNIT,900)
      IF (ALLOCATED(ysrc)) DEALLOCATE(ysrc)
      ! May be ready to allocate space for travel times 
      IF (ALLOCATED(ttimesRec)) DEALLOCATE(ttimesRec)
      IF (nrec > 0) THEN
         IF (ALLOCATED(ttimesRec)) DEALLOCATE(ttimesRec)
         ALLOCATE(ttimesRec(nrec*nsrc))
      ENDIF
  900 FORMAT('fteik_solver2d_setSources64f: Failed to set source')
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
!>    @date July 2017
!>    @copyright CeCILL-3
!>    @ingroup solver2d
!>    This is now a subroutine and incorporated with the fteik Fortran modules.
      SUBROUTINE fteik_solver2d_setSphereToCartEpsilon(epsIn, ierr) &
      BIND(C, NAME='fteik_solver2d_setSphereToCartEpsilon')
      USE FTEIK_MODEL64F, ONLY : nz, nx
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE FTEIK_RAYS64F, ONLY : fteik_rays_setEpsS2C
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: epsIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      epsS2C = zero
      IF (nx < 1 .OR. nz < 1) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      IF (INT(epsIn) > nz .OR. INT(epsIn) > nx) THEN
         IF (INT(epsIn) > nz) WRITE(ERROR_UNIT,901) INT(epsIn), nz
         IF (INT(epsIn) > nx) WRITE(ERROR_UNIT,902) INT(epsIn), nx
         ierr = 1
         RETURN
      ENDIF
      epsS2C = epsIn
      CALL fteik_rays_setEpsS2C(epsS2C)
  900 FORMAT('fteik_solver2d_setSphereToCartEpsilon: Grid not yet set')
  901 FORMAT('fteik_solver2d_setSphereToCartEpsilon: eps=', I0, ' bigger than nz =',I0)
  902 FORMAT('fteik_solver2d_setSphereToCartEpsilon: eps=', I0, ' bigger than nx =',I0)
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the nodal velocity model.
!>
!>    @param[in] ng     Number of grid points in velocity.  This should be equal to
!>                      [nz x nx].
!>    @param[in] order  Defines the column major order of vel.
!>    @param[in] order  If order == FTEIK_XZ_ORDERING then vel has dimension [nx x nz]
!>                      where nx is the leading dimension and the model is 2D.
!>    @param[in] order  If order == FTEIK_ZX_ORDERING or FTEIK_NATURAL_ORDERING
!>                      then vel has dimension [nz x nx] where nz is the leading
!>                      dimension and the model is 2D.
!>    @param[in] vel    Nodal velocity model whose leading dimension are given
!>                      by order.
!>    @param[out] ierr  0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_setNodalVelocityModel64f(ng, order, vel, ierr)  &
      BIND(C, NAME='fteik_solver2d_setNodalVelocityModel64f')
      USE FTEIK_MODEL64F, ONLY : fteik_model_setNodalVelocityModel64f
      INTEGER(C_INT), VALUE, INTENT(IN) :: ng, order
      REAL(C_DOUBLE), INTENT(IN) :: vel(ng)
      INTEGER(C_INT), INTENT(OUT) :: ierr 
      CALL fteik_model_setNodalVelocityModel64f(ng, order, vel, ierr)
      IF (ierr /= 0) WRITE(ERROR_UNIT,900)
  900 FORMAT('fteik_solver2d_setNodalVelocityModel64f: Error setting velocity model')
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the cell-based velocity model.
!>
!>    @param[in] ng     Number of grid points in velocity.  This should be equal to
!>                      [nz-1 x nx-1].
!>    @param[in] order  Defines the column major order of vel. \n
!>                      If order == FTEIK_XZ_ORDERING then vel has dimension [nx-1 x nz-1]
!>                      where nx-1 is the leading dimension. \n
!>                      If order == FTEIK_ZX_ORDERING or FTEIK_NATURAL_ORDERING
!>                      then vel has dimension [nz-1 x nx-1] where nz-1 is the leading
!>                      dimension.
!>    @param[in] vel    Cell-based velocity model whose leading dimension are given
!>                      by order.
!>    @param[out] ierr  0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_setCellVelocityModel64f(nc, order, vel, ierr)  &
      BIND(C, NAME='fteik_solver2d_setCellVelocityModel64f')
      USE FTEIK_MODEL64F, ONLY : fteik_model_setCellVelocityModel64f 
      INTEGER(C_INT), VALUE, INTENT(IN) :: nc, order
      REAL(C_DOUBLE), INTENT(IN) :: vel(nc)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      CALL fteik_model_setCellVelocityModel64f(nc, order, vel, ierr)
      IF (ierr /= 0) WRITE(ERROR_UNIT,900)
  900 FORMAT('fteik_solver2d_setCellVelocityModel64f: Error setting velocity model')
      RETURN
      END 

!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the velocity model on the solver.
!>    @param[in] ncell    Number of cells in velocity model.  This should be 
!>                        (nz-1)*(nx-1).
!>    @param[in] vel      Velocity model (meters/second) in model cells.  This is
!>                        a vector of dimension [ncell] whose fastest direction is
!>                        z and whose slowest direction is y.
!>    @param[out] ierr    0 indicates success.
!>    @ingroup solver2d
      SUBROUTINE fteik_solver2d_setVelocityModel64f(ncell, vel, ierr) &
      BIND(C, NAME='fteik_solver2d_setVelocityModel64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ncell
      REAL(C_DOUBLE), INTENT(IN) :: vel(ncell)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      CALL fteik_model_setVelocityModel64f(ncell, vel, ierr)
      IF (ierr /= 0) WRITE(ERROR_UNIT,900)
  900 FORMAT('fteik_solver2d_setVelocityModel64f: Error setting velocity model')
      RETURN
      END
END MODULE
