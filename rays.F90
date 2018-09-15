!> @defgroup rays Raypaths
!> @ingroup solver2d
!> @ingroup solver3d
!> @brief Utilities for computing raypaths in the travel time field.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE FTEIK_RAYS64F
      USE ISO_FORTRAN_ENV
      USE ISO_C_BINDING
      TYPE rayPath
         !> Defines the (x,y,z) along a ray segment.
         DOUBLE PRECISION, ALLOCATABLE :: xyzs(:)
         !> The cell indices that the ray visits.
         INTEGER, ALLOCATABLE :: cells(:)
         !> Number of segments.
         INTEGER :: nseg
      END TYPE
      !> Defines the ray trajectories.
      TYPE(rayPath), ALLOCATABLE, PROTECTED, SAVE :: rayPaths(:)
      !> The x origin of the rays.  This has dimension [nrays].
      DOUBLE PRECISION, ALLOCATABLE, PROTECTED, SAVE :: xr0(:)
      !> The y origin of the rays.  This has dimension [nrays].
      DOUBLE PRECISION, ALLOCATABLE, PROTECTED, SAVE :: yr0(:)
      !> The z origin of the rays.  This has dimension [nrays].
      DOUBLE PRECISION, ALLOCATABLE, PROTECTED, SAVE :: zr0(:)
      !> The number of ray paths.
      INTEGER, PROTECTED, SAVE :: nrays = 0
      !> Defines whether the module is initialized.
      LOGICAL, PRIVATE, SAVE :: linit = .FALSE.

      PUBLIC :: fteik_rays_initialize64f
      PUBLIC :: fteik_rays_free

      PRIVATE :: fteik_rays_trace2d
      PRIVATE :: fteik_rays_trace3d
      PRIVATE :: fteik_rays_resetRayPath
      PRIVATE :: fteik_rays_setRayOrigins
      PRIVATE :: computeGradient2D 
      CONTAINS
!========================================================================================!
!                                       Begin the code                                   !
!========================================================================================!
      SUBROUTINE fteik_rays_initialize(ierr) &
      BIND(C, NAME='fteik_rays_initialize')
      USE FTEIk_MODEL64F, ONLY : ngrd
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      CALL fteik_rays_free()
      ! Get the grid information
      IF (ngrd < 1) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      linit = .TRUE.
  900 FORMAT('fteik_ray_initialize64f: Model not yet set')
      RETURN
      END SUBROUTINE

      SUBROUTINE fteik_rays_free() &
      BIND(C, NAME='fteik_rays_free')
      CALL fteik_rays_resetRayPath()
      IF (ALLOCATED(xr0)) DEALLOCATE(xr0)
      IF (ALLOCATED(yr0)) DEALLOCATE(yr0)
      IF (ALLOCATED(zr0)) DEALLOCATE(zr0)
      nrays = 0
      linit = .FALSE.
      RETURN
      END

!>    @brief Defines the ray end points.
!>    @param[in] nr     Number of ray-path origins.
!>    @param[in] xr     Defines the x origin (meters) of the ir'th ray-path.
!>                      This has dimension [nr].
!>    @param[in] yr     Defines the y origin (meters) of the ir-th ray-path.
!>                      This has diemnsion [nr].  This is optional.
!>    @param[in] zr     Defines the z origin (meters) of the ir-th ray-path.
!>                      This has dimension [nr].
!>    @param[out] ierr  0 indicates success.
!>    @ingroup rays
      SUBROUTINE fteik_rays_setRayDestinations(nr, xr, zr, ierr, yr)
      USE FTEIK_MODEL64F, ONLY : xmin => x0
      USE FTEIK_MODEL64F, ONLY : ymin => y0
      USE FTEIK_MODEL64F, ONLY : zmin => z0
      USE FTEIK_MODEL64F, ONLY : dx, dy, dz, lis3dModel, nx, ny, nz
      INTEGER, INTENT(IN) :: nr
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: xr, zr 
      DOUBLE PRECISION, DIMENSION(:), OPTIONAL, INTENT(IN) :: yr
      INTEGER, INTENT(OUT) :: ierr
      DOUBLE PRECISION :: xmax, ymax, zmax
      INTEGER ir
      ! Ensure the
      ierr = 0
      nrays = 0
      CALL fteik_rays_resetRayPath()
      IF (ALLOCATED(xr0)) DEALLOCATE(xr0)
      IF (ALLOCATED(yr0)) DEALLOCATE(yr0)
      IF (ALLOCATED(zr0)) DEALLOCATE(zr0)
      IF (.NOT.linit) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      IF (nr < 1) THEN
         WRITE(ERROR_UNIT,901)
         ierr = 1
         RETURN
      ENDIF
      IF (lis3dModel .AND. .NOT.PRESENT(yr)) THEN
         WRITE(ERROR_UNIT,902)
         ierr = 1
         RETURN
      ENDIF
      zmax = zmin + DBLE(nz - 1)*dz
      xmax = xmin + DBLE(nx - 1)*dx
      ymax = ymin
      IF (lis3dModel) ymax = ymin + DBLE(ny - 1)*dy
      ! Ensure the origins are in the model
      DO ir=1,nr
         IF (zr(ir) < zmin .OR. zr(ir) > zmax) THEN
            WRITE(ERROR_UNIT,905) ir, xr(ir), xmin, xmax
            ierr = 1
         ENDIF
         IF (xr(ir) < xmin .OR. xr(ir) > xmax) THEN
            WRITE(ERROR_UNIT,906) ir, xr(ir), xmin, xmax
            ierr = 1
         ENDIF
         IF (lis3dModel .AND. (yr(ir) < ymin .OR. yr(ir) > ymax)) THEN
            WRITE(ERROR_UNIT,907) ir, yr(ir), ymin, ymax
            ierr = 1
         ENDIF
      ENDDO
      IF (ierr /= 0) RETURN
      IF (ALLOCATED(xr0)) DEALLOCATE(xr0)
      IF (ALLOCATED(yr0)) DEALLOCATE(yr0)
      IF (ALLOCATED(zr0)) DEALLOCATE(zr0)
      nrays = nr
      ALLOCATE(xr0(nrays))
      ALLOCATE(yr0(nrays))
      ALLOCATE(zr0(nrays))
      zr0(1:nr) = zr(1:nr)
      xr0(1:nr) = xr(1:nr)
      IF (lis3dModel) THEN
         yr0(1:nr) = yr(1:nr)
      ELSE
         yr0(1:nr) = 0.d0
      ENDIF
  900 FORMAT('fteik_rays_setRayDestinations: Module not initialized')
  901 FORMAT('fteik_rays_setRayDestinations: No origins!')
  902 FORMAT('fteik_rays_setRayDestinations: y ray origins must be given')
  905 FORMAT('fteik_rays_setRayDestinations: zr(', I0, ')=', E12.5, &
             ' must be in range [', E12.5, E12.5, ']')
  906 FORMAT('fteik_rays_setRayDestinations: xr(', I0, ')=', E12.5, &
             ' must be in range [', E12.5, E12.5, ']')
  907 FORMAT('fteik_rays_setRayDestinations: yr(', I0, ')=', E12.5, &
             ' must be in range [', E12.5, E12.5, ']')
      RETURN
      END SUBROUTINE

      SUBROUTINE fteik_rays_resetRayPath()
      INTEGER i
      IF (ALLOCATED(rayPaths)) THEN
         DO i=1,SIZE(rayPaths)
            IF (ALLOCATED(rayPaths(i)%xyzs))  DEALLOCATE(rayPaths(i)%xyzs)
            IF (ALLOCATED(rayPaths(i)%cells)) DEALLOCATE(rayPaths(i)%cells)
            rayPaths(i)%nseg = 0 
         ENDDO
         DEALLOCATE(rayPaths)
      ENDIF
      RETURN
      END

!>    @brief Sets the ray-path destinations.
      SUBROUTINE fteik_rays_setRayDestinations3D(nr, xr, yr, zr, ierr) &
      BIND(C, NAME='fteik_rays_setRayDestinations3D')
      INTEGER(C_INT), VALUE, INTENT(IN) :: nr
      REAL(C_DOUBLE), DIMENSION(nr), INTENT(IN) :: xr, yr, zr
      INTEGER(C_INT), INTENT(OUT) :: ierr
      CALL fteik_rays_setRayDestinations(nr, xr, zr, ierr, yr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1 
      ENDIF
  900 FORMAT('fteik_rays_setRayDestinations2D: Failed to set origins')
      RETURN
      END SUBROUTINE

!>    @brief Sets the ray-path origins.
      SUBROUTINE fteik_rays_setRayDestinations2D(nr, xr, zr, ierr) &
      BIND(C, NAME='fteik_rays_setDestinations2D')
      INTEGER(C_INT), VALUE, INTENT(IN) :: nr
      REAL(C_DOUBLE), DIMENSION(nr), INTENT(IN) :: xr, zr
      INTEGER(C_INT), INTENT(OUT) :: ierr
      CALL fteik_rays_setRayDestinations(nr, xr, zr, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
      ENDIF
  900 FORMAT('fteik_rays_setRayOrigins2D: Failed to set origins')
      RETURN
      END SUBROUTINE

      SUBROUTINE fteik_rays_trace(ngrdIn, ttimes, xs, ys, zs, ierr) &
      BIND(C, NAME='fteik_rays_trace')
      USE FTEIK_MODEL64F, ONLY : lis3dModel, ngrd
      INTEGER(C_INT), VALUE, INTENT(IN) :: ngrdIn
      REAL(C_DOUBLE), INTENT(IN) :: ttimes(ngrdIn), xs, ys, zs
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      IF (nrays < 1) THEN
         WRITE(ERROR_UNIT,901)
         ierr = 1
         RETURN
      ENDIF
      IF (ngrdIn /= ngrd) THEN
         WRITE(ERROR_UNIT,902) ngrdIn, ngrd
         ierr = 1
         RETURN
      ENDIF
      IF (lis3dModel) THEN
         CALL fteik_rays_trace3d(xs, ys, zs, ttimes, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,905)
            ierr = 1
         ENDIF
      ELSE
         CALL fteik_rays_trace2d(xs, zs, ttimes, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,906)
            ierr = 1
         ENDIF
      ENDIF
  901 FORMAT('fteik_rays_trace: No ray origins are set')
  902 FORMAT('fteik_rays_trace: ngrdIn = ', I0, ' expecting ngrd = ', I0)
  905 FORMAT('fteik_rays_trace: Failed the 3D ray tracing')
  906 FORMAT('fteik_rays_trace: Failed the 2D ray tracing')
      RETURN
      END

      SUBROUTINE fteik_rays_trace3d(xs, ys, zs, ttimes, ierr)
      DOUBLE PRECISION, INTENT(IN) :: xs, ys, zs
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ttimes
      INTEGER, INTENT(OUT) :: ierr
      ierr = 1
      RETURN
      END

      SUBROUTINE fteik_rays_trace2d(xs, zs, ttimes, ierr)
      USE FTEIK_MODEL64F, ONLY : dx, dz, x0, z0, nx, nz
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: xs, zs
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ttimes
      INTEGER, INTENT(OUT) :: ierr
      DOUBLE PRECISION :: dtmag, dtdx, dtdz, t1, t2, t3, t4, x, z 
      DOUBLE PRECISION :: dxs02, dzs02, dxs12, dzs12, &
                          dx1ds, dx3ds, eta, extra, l, x0r, xi, z0r, x1r, z1r
      INTEGER icell, igrd, igrd1, igrd2, igrd3, igrd4, ir, k, ix, ixs, iz, izs, maxk
      ! Want to move just enough to pop into the next element
      DOUBLE PRECISION, PARAMETER :: pert = 1.d-5
      LOGICAL lexit
      DOUBLE PRECISION, PARAMETER :: rad02 = 1.01d0*1.01d0 ! elements squared
      ierr = 0
      maxk = nx*nz
      extra = MAX(dx, dz)*1.d-5
      !deta_dx = 0
      !dxi_dz = 0 
      ixs = MIN(MAX(1, INT((xs - x0)/dx + 0.5d0)), nx - 1)
      izs = MIN(MAX(1, INT((zs - z0)/dz + 0.5d0)), nz - 1)
      DO ir=1,nrays
         ! Begin tracing the ray at the ray origin 
         x0r = xr0(ir) - x0
         z0r = zr0(ir) - z0
         DO k=1,maxk
            ! We're definitely in the grid
            ix = MIN(MAX(1, INT(x0r/dx + 0.5d0)), nx - 1)
            iz = MIN(MAX(1, INT(z0r/dz + 0.5d0)), nz - 1)
!           print *, dble(ix-1)*dx, x0r, dble(ix)*dx
!           print *, dble(iz-1)*dz, z0r, dble(iz)*dz

            ! Compute the local coordinates in reference element.
            xi  = (2.d0*x0r - DBLE(ix - 1)*dx - DBLE(ix)*dx)/dx
            eta = (2.d0*z0r - DBLE(iz - 1)*dz - DBLE(iz)*dz)/dz
!print *, xi, eta
            ! Compute the gradient of the travel time at this point.
            igrd  = (ix - 1)*nz + iz

            igrd1 = (ix - 1)*nz + iz 
            igrd2 = (ix - 0)*nz + iz
            igrd3 = (ix - 0)*nz + iz + 1
            igrd4 = (ix - 1)*nz + iz + 1

            icell = (ix - 1)*(nz - 1) + iz
            ! Get finite element ordering
            t1 = ttimes(igrd1)
            t4 = ttimes(igrd4)
            t2 = ttimes(igrd2)
            t3 = ttimes(igrd3)
            CALL computeGradient2D(dx, dz,         &
                                   xi, eta,        &
                                   t1, t2, t3, t4, &
                                   dtdx, dtdz)
!           ! Compute the forward difference
!           dtdx = (t2 - t1)/dx ! Note +x increases +west
!           dtdz = (t4 - t1)/dz ! Note +z increases +down
!print *, dtdx, dtdz
!print *, (t2 - t1)/dx, (t1 - t4)/dz
!print *, 'ttimes:'
!print *, t4, t3
!print *, t1, t2
            ! Compute the slowness vector p = \{ \frac{dT}{dx}, \frac{dT}{dz} \}.
            ! Compute dT/dx = dT/dxi dxi/dx + dT/deta deta/dx = dT/dxi dxi/dx
            ! Compute dT/dz = dT/dxi dxi/dz + dT/deta deta/dz = dT/deta deta/dz
            ! where both equations have been simplified b/c the elements are rectangles.

            !(eta - 1)*(t4 - t1) + (eta + 1)*(t3 - t2)/(2.d0/(2.d0*dx)
            !(1 - xi)*(t1 - t2)*xi + (t3 - t4)*xi - t1 - t2 + t3 + t4)/(2.d0*dz)

            ! Normalize.  To obtain a unit normal in the slowness direction.  Technically,
            ! the slowness. p, satisfies the eikonal equation so that
            ! |p| = |\nabla T| = \frac{1}{c} e.g., Chapman pg. 140 Eqn. 5.1.16.
            ! So in effect, these are the direction cosines of the slowness.
            dtmag = HYPOT(dtdx, dtdz) !DSQRT(dtdx**2 + dtdz**2)
            dx1ds = dtdx/dtmag
            dx3ds = dtdz/dtmag
            ! The gradient points in the direction of max increase of the travel time
            ! field.  However, we are tracing from a point back to the source so we 
            ! need to go in the opposite direction.
            dx1ds =-dx1ds
            dx3ds =-dx3ds
!print *, 'direction:', dx1ds, dx3ds
            ! The element is piecewise linear so the gradient is constant.  Let's step
            ! just far enough so that we enter the next element.  This means we solve
            ! for the length ds by considering a line beginning at x0, finishing at x1,
            ! and whose slope is dx1/dx.  Solving x1 = x0 + dx1/ds*l for the distance
            ! in l  yields l = (x1 - x0)/(dx1/ds).
            l = 0.d0
            IF (DABS(dx1ds) > DABS(dx3ds)) THEN
               IF (dx1ds < 0.d0) THEN
                  l = x0r + DBLE(ix - 1)*dx - pert*dx 
               ELSE ! Move far enough to the east.
                  l = x0r + DBLE(ix)*dx + pert*dx
               ENDIF
            ELSE
               ! Move far enough up so that we move into the next element.
               IF (dx3ds < 0.d0) THEN
                  l = z0r + DBLE(iz - 1)*dz - pert*dz
               ELSE ! Move far enough down.
                  l = z0r + DBLE(iz)*dz + pert*dz 
               ENDIF
            ENDIF
            ! Move this distance along the segment
            x1r = x0r + dx1ds*dx!/DSQRT(2.d0)/28.d0 !l
            z1r = z0r + dx3ds*dz!/DSQRT(2.d0)/28.d0 !l
            ! If moving from (x0r,z0r) to the source is a smaller step then move to the
            ! source.  Otherwise, move to the next point and continue.
            dxs02 = (x0r - xs)**2
            dzs02 = (z0r - zs)**2
            dxs12 = (x1r - xs)**2
            dzs12 = (z1r - zs)**2
            lexit = .FALSE.
            IF (dxs02 + dzs02 < rad02*(dx*dx + dz*dz) + dxs12 + dzs12) THEN
               x1r = xs
               z1r = zs
               lexit = .TRUE.
            ENDIF 
            ! Now update
            !nseg = nseg + 1
            print *, k, x1r, z1r
            IF (lexit) EXIT
!pause
            x0r = x1r
            z0r = z1r
!           IF (dtdx < 0.d0) THEN
!              xpert = x - (x0 + DBLE(ix - 1)*dx) - tol
!           ELSE
!            ! xpert = 
!           ENDIF
!           IF (dtdz < 0.d0) THEN

!           ELSE

!           ENDIF
            
  
              
!           dtdx = 0.25d0*( (1.d0 - xi)* 
!           dtdz = 0.25d0

            !0.25d0* 
         ENDDO
!pause
      ENDDO
      RETURN
      END SUBROUTINE

      PURE SUBROUTINE computeGradient2D(dx, dz,          &
                                        xi, eta,         &
                                        t1, t2, t3, t4,  &
                                        dtdx, dtdz)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: dx, dz, xi, eta, t1, t2, t3, t4
      DOUBLE PRECISION, INTENT(OUT) :: dtdx, dtdz
      dtdx = 0.5d0*(eta*t1 - eta*t2 + eta*t3 - eta*t4 - t1 + t2 + t3 - t4)/dx
      dtdz = 0.5d0*(t1*xi  - t2*xi  + t3*xi  - t4*xi  - t1 - t2 + t3 + t4)/dz
      RETURN
      END

END MODULE
