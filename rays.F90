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
         !> Defines the length of each segment.  
         DOUBLE PRECISION, ALLOCATABLE :: l(:)
         !> The travel time (seconds) of the ray.
         DOUBLE PRECISION :: time
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

      !> Number of grid points.
      INTEGER, PRIVATE, SAVE :: ngrd = 0
      !> Number of grid points in z.
      INTEGER, PRIVATE, SAVE :: nz = 0
      !> Number of grid points in x.
      INTEGER, PRIVATE, SAVE :: nx = 0
      !> Number of grid points in y.
      INTEGER, PRIVATE, SAVE :: ny = 0

      DOUBLE PRECISION, PRIVATE, SAVE :: z0 = 0.d0
      DOUBLE PRECISION, PRIVATE, SAVE :: x0 = 0.d0
      DOUBLE PRECISION, PRIVATE, SAVE :: y0 = 0.d0 

      DOUBLE PRECISION, PRIVATE, SAVE :: dz = 0.d0
      DOUBLE PRECISION, PRIVATE, SAVE :: dx = 0.d0
      DOUBLE PRECISION, PRIVATE, SAVE :: dy = 0.d0


      DOUBLE PRECISION, PRIVATE, SAVE :: epsS2C = 0.d0
      !> Flag indicating that we are to be doing 3D ray tracing.
      LOGICAL, PRIVATE, SAVE :: ltrace3D = .TRUE.
      !> True source offset (meters) in z.
      DOUBLE PRECISION, PRIVATE, SAVE :: zsTrue = 0.d0
      !> True source offset (meters) in x.
      DOUBLE PRECISION, PRIVATE, SAVE :: xsTrue = 0.d0
      !> True source offset (meters) in y.
      DOUBLE PRECISION, PRIVATE, SAVE :: ysTrue = 0.d0
      !> The source grid point (meters) in z.
      DOUBLE PRECISION, PRIVATE, SAVE :: zsrc = 0.d0
      !> The source grid point (meters) in z.
      DOUBLE PRECISION, PRIVATE, SAVE :: xsrc = 0.d0
      !> The source grid point (meters) in y.
      DOUBLE PRECISION, PRIVATE, SAVE :: ysrc = 0.d0
      !> Slowness at source.
      DOUBLE PRECISION, PRIVATE, SAVE :: szero
      !> Source grid point in z.
      DOUBLE PRECISION, PRIVATE, SAVE :: zsa = 0.d0
      !> Source grid point in x.
      DOUBLE PRECISION, PRIVATE, SAVE :: xsa = 0.d0
      !> Source grid point in y.
      DOUBLE PRECISION, PRIVATE, SAVE :: ysa = 0.d0
      !> The source index in z.
      INTEGER, PRIVATE, SAVE :: zsi = 0
      !> The source index in x.
      INTEGER, PRIVATE, SAVE :: xsi = 0
      !> The source index in y.
      INTEGER, PRIVATE, SAVE :: ysi = 0

      !> Solve the ray path equation
      !>   \f$ \frac{\textbf{x}}{ds} = \frac{\textbf{s}}{u} \f$
      !> with the (forward) Euler method.
      INTEGER, PARAMETER :: FTEIK_RAYS_ODE_EULER = 1
      !> Solve the ray path equation
      !>   \f$ \frac{\textbf{x}}{ds} = \frac{\textbf{s}}{u} \f$
      !> with the (forward) midpoint rule. 
      INTEGER, PARAMETER :: FTEIK_RAYS_ODE_MIDPOINT = 2
      !> Solve the ray path equation
      !>   \f$ \frac{\textbf{x}}{ds} = \frac{\textbf{s}}{u} \f$
      !> with the Runge-Kutta 4.
      INTEGER, PARAMETER :: FTEIK_RAYS_ODE_RK4 = 3

      !> The integrator used to solve the raypath equation.
      INTEGER, PRIVATE, SAVE :: odeIntegrator = FTEIK_RAYS_ODE_EULER

      ! Accessible from Fortran
      PUBLIC :: fteik_rays_initialize64f
      PUBLIC :: fteik_rays_free
      PUBLIC :: fteik_rays_setSource
      PUBLIC :: fteik_rays_setEpsS2C
      PUBLIC :: fteik_rays_setODEIntegrator
      PUBLIC :: fteik_rays_setRayDestinations3D
      PUBLIC :: fteik_rays_setRayDestinations2D

      ! Private routines
      PRIVATE :: fteik_rays_trace2d
      PRIVATE :: fteik_rays_trace3d
      PRIVATE :: fteik_rays_resetRayPath
      PRIVATE :: fteik_rays_setRayOrigins
      PRIVATE :: updateMidpoint2D
      PRIVATE :: updateRK42D
      PRIVATE :: updateEuler2D
      PRIVATE :: intersect2d
      PRIVATE :: evalf2d
      CONTAINS
!========================================================================================!
!                                       Begin the code                                   !
!========================================================================================!
!>    @brief Initializes the ray tracer.
!>    @param[out] ierr   0 indicates success.
!>    @ingroup rays
      SUBROUTINE fteik_rays_initialize(ierr)
      USE FTEIK_MODEL64F, ONLY : ngrdModel => ngrd
      USE FTEIK_MODEL64F, ONLY : nzModel => nz 
      USE FTEIK_MODEL64F, ONLY : nxModel => nx
      USE FTEIK_MODEL64F, ONLY : nyModel => ny
      USE FTEIK_MODEL64F, ONLY : dzModel => dz
      USE FTEIK_MODEL64F, ONLY : dxModel => dx
      USE FTEIK_MODEL64F, ONLY : dyModel => dy
      USE FTEIK_MODEL64F, ONLY : z0Model => z0
      USE FTEIK_MODEL64F, ONLY : x0Model => x0
      USE FTEIK_MODEL64F, ONLY : y0Model => y0
      USE FTEIK_MODEL64F, ONLY : lis3dModel
      INTEGER, INTENT(OUT) :: ierr
      ierr = 0
      CALL fteik_rays_free()
      ! Determine if the grid has been initialized.
      IF (ngrdModel < 1) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      IF (dxModel < EPSILON(1.d0) .OR. dzModel < EPSILON(1.d0)) THEN
         WRITE(ERROR_UNIT,905)
         ierr = 1
         RETURN
      ENDIF
      ! Copy the grid information.
      ngrd = ngrdModel
      nz = nzModel
      nx = nxModel
      ny = MAX(1, nyModel)
      dz = dzModel
      dx = dxModel
      dy = dyModel
      z0 = z0Model
      x0 = x0Model
      y0 = y0Model
      ltrace3D = lis3dModel
      linit = .TRUE.
  900 FORMAT('fteik_ray_initialize64f: Model not yet set')
  905 FORMAT('fteik_ray_initialize64f: Only dx > 0 and dz > 0 considered')
      RETURN
      END SUBROUTINE

      SUBROUTINE fteik_rays_free()
      CALL fteik_rays_resetRayPath()
      IF (ALLOCATED(xr0)) DEALLOCATE(xr0)
      IF (ALLOCATED(yr0)) DEALLOCATE(yr0)
      IF (ALLOCATED(zr0)) DEALLOCATE(zr0)
      odeIntegrator = FTEIK_RAYS_ODE_EULER
      nrays = 0
      ngrd = 0
      nz = 0
      nx = 0
      ny = 0
      dz = 0.d0
      dx = 0.d0
      dy = 0.d0 
      z0 = 0.d0
      x0 = 0.d0
      y0 = 0.d0
      epsS2C = 0.d0
      ltrace3d = .TRUE.
      linit = .FALSE.
      RETURN
      END

!>    @brief The defines the integrator used to solve the ray-path equation 
!>           \f$
!>              \frac{d \textbf{x}}{ds} = \textbf{s}{u}
!>           \f$
!>           where the slowness vector, $\textbf{s} = \nabla T(x,y,z) \f$,
!>           the slowness, \f$ u \f$, is the reciprocal of the velocity,
!>           \textbf{x} describes the position of the ray,
!>           and \f$ ds \f$ is the length of the segment.
!>    @param[in] integrator   FTEIK_RAYS_ODE_EULER will use Euler's method (fast).
!>    @param[in] integrator   FTEIK_RAYS_ODE_MIDPOINT will use the midpoint method
!>                            which is likely more accurate than Euler but slower.
!>    @param[in] integrator   FTEIK_RAYS_ODE_RK4 will use Runge-Kutta 4 which is
!>                            the slowest but likely the most accurate.
!>    @param[out] ierr  0 indicates success.
!>    @ingroup rays
      SUBROUTINE fteik_rays_setODEIntegrator(integrator, ierr) &
      BIND(C, NAME='fteik_rays_setODEIntegrator')
      INTEGER(C_INT), VALUE, INTENT(IN) :: integrator
      INTEGER(C_INT), INTENT(OUT) :: ierr
      odeIntegrator = FTEIK_RAYS_ODE_EULER  
      IF (ode < FTEIK_RAYS_ODE_EULER .OR. ode > FTEIK_RAYS_ODE_RK4) THEN
         WRITE(ERROR_UNIT,900) ode
         ierr = 1
         RETURN
      ENDIF
      odeIntegrator = integrator
  900 FORMAT('fteik_rays_setODEIntegrator: Invalid integrator = ', I0, ' using Euler')
      RETURN
      END

      SUBROUTINE fteik_rays_setEpsS2C(epsS2Cin)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: epsS2Cin
      epsS2C = epsS2Cin
      RETURN
      END

      SUBROUTINE fteik_rays_setSource(zsTrueIn,  xsTrueIn,  ysTrueIn, &
                                      zsrcIn,    xsrcIn,    ysrcIn,   &
                                      zsaIn,     xsaIn,     ysaIn,    &
                                      zsiIn,     xsiIn,     ysiIn,    &
                                      szeroIn)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: zsTrueIn, xsTrueIn, ysTrueIn, &
                                      zsrcIn, xsrcIn, ysrcIn,       &
                                      zsaIn, xsaIn, ysaIn
      INTEGER, INTENT(IN) :: zsiIn, xsiIn, ysiIn
      DOUBLE PRECISION, INTENT(IN) :: szeroIn
      zsTrue = zsTrueIn - z0
      xsTrue = xsTrueIn - x0
      ysTrue = ysTrueIn - y0
      zsrc = zsrcIn
      xsrc = xsrcIn
      ysrc = ysrcIn
      zsa = zsaIn
      xsa = xsaIn
      ysa = ysaIn
      zsi = zsiIn
      xsi = xsiIn
      ysi = ysiIn
   
      szero = szeroIn
      RETURN
      END SUBROUTINE
!========================================================================================!
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
      IMPLICIT NONE
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
      IF (ltrace3D .AND. .NOT.PRESENT(yr)) THEN
         WRITE(ERROR_UNIT,902)
         ierr = 1
         RETURN
      ENDIF
      zmax = z0 + DBLE(nz - 1)*dz
      xmax = x0 + DBLE(nx - 1)*dx
      ymax = y0
      IF (ltrace3D) ymax = y0 + DBLE(ny - 1)*dy
      ! Ensure the origins are in the model
      DO ir=1,nr
         IF (zr(ir) < z0 .OR. zr(ir) > zmax) THEN
            WRITE(ERROR_UNIT,905) ir, xr(ir), z0, zmax
            ierr = 1
         ENDIF
         IF (xr(ir) < x0 .OR. xr(ir) > xmax) THEN
            WRITE(ERROR_UNIT,906) ir, xr(ir), x0, xmax
            ierr = 1
         ENDIF
         IF (ltrace3D) THEN
            IF (yr(ir) < y0 .OR. yr(ir) > ymax) THEN
               WRITE(ERROR_UNIT,907) ir, yr(ir), y0, ymax
               ierr = 1
            ENDIF
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
      IF (ltrace3D) THEN
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
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Resets the ray paths.
!>    @ingroup rays
      SUBROUTINE fteik_rays_resetRayPath()
      INTEGER i
      IF (ALLOCATED(rayPaths)) THEN
         DO i=1,SIZE(rayPaths)
            IF (ALLOCATED(rayPaths(i)%xyzs))  DEALLOCATE(rayPaths(i)%xyzs)
            IF (ALLOCATED(rayPaths(i)%l))     DEALLOCATE(rayPaths(i)%l)
            IF (ALLOCATED(rayPaths(i)%cells)) DEALLOCATE(rayPaths(i)%cells)
            rayPaths(i)%time = 0.d0
            rayPaths(i)%nseg = 0 
         ENDDO
         DEALLOCATE(rayPaths)
      ENDIF
      RETURN
      END
!>    @brief Sets the ray-path destinations.
!>    @param[in] nr     Number of rays to compute.
!>    @param[in] xr     The ray origins (meters) in x.
!>    @param[in] yr     The ray origins (meters) in y.
!>    @param[in] zr     The ray origins (meters) in z.
!>    @param[out] ierr  0 indicates success.
!>    @ingroup rays
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
  900 FORMAT('fteik_rays_setRayDestinations3D: Failed to set origins')
      RETURN
      END SUBROUTINE
!>    @brief Sets the ray-path origins.
!>    @param[in] nr     Number of rays to compute.
!>    @param[in] zr     The ray origins (meters) in z.
!>    @param[out] ierr  0 indicates success.
!>    @ingroup rays
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

      SUBROUTINE fteik_rays_trace(ngrdIn, ttimes, tgradx, tgradz, ierr, tgrady)
      INTEGER, INTENT(IN) :: ngrdIn
      DOUBLE PRECISION, DIMENSION(ngrdIn), INTENT(IN) :: ttimes, tgradx, tgradz
      DOUBLE PRECISION, OPTIONAL, DIMENSION(ngrdIn), INTENT(IN) :: tgrady
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
      IF (ltrace3D) THEN
         CALL fteik_rays_trace3d(ttimes, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,905)
            ierr = 1
         ENDIF
      ELSE
         CALL fteik_rays_trace2d(ttimes, tgradx, tgradz, ierr)
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

      SUBROUTINE fteik_rays_trace3d(ttimes, ierr)
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ttimes
      INTEGER, INTENT(OUT) :: ierr
      ierr = 1
      RETURN
      END

      SUBROUTINE fteik_rays_trace2d(ttimes, tgradx, tgradz, ierr)
!     IMPLICIT NONE
      USE FTEIK_MODEL64F, ONLY : slow
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ttimes, tgradx, tgradz
      INTEGER, INTENT(OUT) :: ierr
      DOUBLE PRECISION :: x0r, z0r
      INTEGER icell, ir, ix, ixs, iz, izs, jbeg, &
              k, nk, maxk, nseg
double precision ds, trs
      DOUBLE PRECISION, ALLOCATABLE :: xpts(:), zpts(:), xray(:), zray(:), rayl(:)
      LOGICAL lterminate
      ! Want to move just enough to pop into the next element
      DOUBLE PRECISION, PARAMETER :: pert = 1.d-8
      DOUBLE PRECISION, PARAMETER :: tol = 1.d-8  ! 0 tolerance for vertica line
      DOUBLE PRECISION, PARAMETER :: rad02 = 1.01d0*1.01d0 ! elements squared
      INTEGER ix0, iz0, l, sgnx, sgnz
      INTEGER, ALLOCATABLE :: cells(:)
double precision t(2,2), dx1ds, dx3ds, dsr, tn, xn, zn, xn1, zn1, x1, z1, x2, z2
      ierr = 0
      maxk = nx*nz
      extra = MAX(dx, dz)*1.d-5
      CALL fteik_rays_resetRayPath()
      ALLOCATE(rayPaths(nrays))
      ALLOCATE(xpts(maxk))
      ALLOCATE(zpts(maxk))
      ALLOCATE(xray(maxk))
      ALLOCATE(zray(maxk))
      ALLOCATE(rayl(maxk))
      ALLOCATE(cells(maxk))
      !deta_dx = 0
      !dxi_dz = 0 
      ixs = MIN(MAX(1, INT((xsrc - x0)/dx + 1.d0)), nx - 1)
      izs = MIN(MAX(1, INT((zsrc - z0)/dz + 1.d0)), nz - 1)
      DO ir=1,nrays
         ! Reduce the ray position so that the origin is (0,0)
         x0r = xr0(ir) - x0
         z0r = zr0(ir) - z0
         xpts(1) = x0r
         zpts(1) = z0r
         nk = 0 ! Counter for the number of segments
ds = MIN(dx, dz)/2.d0!/3.d0
         ! Set the counters to zero
         tn = 0.d0
         tne = 0.d0
       
         ! Initialize the position at the receiver 
         xn = x0r
         zn = z0r
         lterminate = .FALSE.
         DO k=1,maxk
            ix = MIN(MAX(1, INT(xn/dx + 1.d0)), nx - 1) 
            iz = MIN(MAX(1, INT(zn/dz + 1.d0)), nz - 1)
            icell = (ix - 1)*(nz - 1) + iz 
            x1 = DBLE(ix-1)*dx; x2 = DBLE(ix)*dx
            z1 = DBLE(iz-1)*dz; z2 = DBLE(iz)*dz
            ! Find the position of the current point
            igrd11 = (ix - 1)*nz + iz
            igrd12 = (ix - 1)*nz + iz + 1
            igrd21 = (ix - 0)*nz + iz
            igrd22 = (ix - 0)*nz + iz + 1
            t(1,1) = ttimes(igrd11)
            t(1,2) = ttimes(igrd12)
            t(2,1) = ttimes(igrd21)
            t(2,2) = ttimes(igrd22)
            ! Update the soln to dx/ds =-\nabla s/u at this point (Shearer pg 85 Eq 4.61).
            ! Note, the negative is b/c we are tracing from receiver to source.
            !  x_{n+1} = x_n + ds*(-u(x_n,z_n)*dT(x_n)/dx)
            !  z_{n_1} = z_n + ds*(-u(x_n,z_n)*dT(z_n)/dz)
            IF (odeIntegrator == FTEIK_RAYS_ODE_MIDPOINT) THEN
               CALL updateMidpoint2D(nx, nz,               &   
                                     ds, dx, dz, xn, zn,   &
                                     slow, tgradx, tgradz, &
                                     xn1, zn1)   
            ELSEIF (odeIntegrator == FTEIK_RAYS_ODE_RK4) THEN
               CALL updateRK42D(nx, nz,               &
                                ds, dx, dz, xn, zn,   &
                                slow, tgradx, tgradz, &
                                xn1, zn1)
            ELSE
               CALL updateEuler2D(nx, nz,               &
                                  ds, dx, dz, xn, zn,   &
                                  slow, tgradx, tgradz, &
                                  xn1, zn1)
            ENDIF
            ! Is the ray in the ballpark of the source?
            IF ( (xn1 - (xsrc - x0))**2 + (zn1 - (zsrc - z0))**2 < (dx**2 + dz**2)/2.d0 ) THEN
               lterminate = .TRUE.
               xn1 = xsrc - x0
               zn1 = zsrc - z0
               nk = k + 1
            ENDIF
            ! Update the travel time along the ray
            ix = MIN(MAX(1, INT(xn/dx + 1.d0)), nx - 1)
            iz = MIN(MAX(1, INT(zn/dz + 1.d0)), nz - 1) 
            icell = (ix - 1)*(nz - 1) + iz
            ! Tabulate the travel time with a Huygen's principle concept where the next
            ! point in the field is the result of a point source at an earlier time.
            tne = tne - (bilinear(xn1, zn1, x1, z1, x2, z2, t) - &
                         bilinear(xn, zn,   x1, z1, x2, z2, t))
            ! Tabulate the travel time with a distance*velocity calculation
            tn  = tn  + DSQRT( (xn1 - xn)**2 + (zn1 - zn)**2 )*slow(icell)
!if (k == 1) print *, xn, zn, 0.d0
!print *, xn1, zn1, 1.d0/slow(icell) !tn
            xn = xn1
            zn = zn1
            xpts(k+1) = xn
            zpts(k+1) = zn
            IF (lterminate) EXIT
         ENDDO
         IF (.NOT.lterminate) THEN
            WRITE(ERROR_UNIT,905) ir 
            ierr = 1
            CYCLE
         ENDIF
         ! Cheap trick to force ray to terminate
         xpts(nk+1) =-2.d0*dx
         zpts(nk+1) =-2.d0*dz
         ! Now tabulate the rays
         iseg = 1
         ! All rays start at the receiver
         xray(1) = xpts(1)
         zray(1) = zpts(1)
         ix0 = INT(xpts(1)/dx + 1.d0)
         iz0 = INT(zpts(1)/dz + 1.d0)
         cells(iseg) = (ix0 - 1)*(nz - 1) + iz0 
         dsr = 0.d0
         jbeg = 1
         DO k=1,nk
            lterminate = .FALSE.
            ! Loop along the ray-path until the end of a cell is reached
            DO j=jbeg,nk
               ix1 = INT(xpts(j+1)/dx + 1.d0)
               iz1 = INT(zpts(j+1)/dz + 1.d0)
               !! Reached the source - finish
               !IF (ix1 < 0 .OR. iz1 < 0) THEN
               !   iseg = iseg + 1
               !   xray(iseg) = xpts(nk)
               !   zray(iseg) = zpts(nk)
               !   cells(iseg) = (ix0 - 1)*(nz - 1) + iz0
               !   lterminate = .TRUE.
               !   EXIT
               !ENDIF
               IF (ix1 /= ix0 .OR. iz1 /= iz0 .OR. j == nk - 1) THEN
                  sgnx = ISIGN(1, ix1 - ix0)
                  sngz = ISIGN(1, iz1 - iz0)
                  ! Crossing multiple cells in one shot isn't done
                  IF (.NOT.lterminate .AND.                            &
                      (ABS(ix1 - ix0) > 1 .OR. ABS(iz1 - iz0) > 1)) THEN
                     WRITE(ERROR_UNIT,910)
                  ELSE
                     ! Find where the ray intersects the grid
                     iseg = iseg + 1
                     dx1ds = xpts(j+1) - xpts(j)
                     dx3ds = zpts(j+1) - zpts(j)
                     cells(iseg) = (ix0 - 1)*(nz - 1) + iz0
                     CALL intersect2d(xpts(j), zpts(j),       &
                                      dx1ds, dx3ds,           &
                                      dx, dz,                 &
                                      xray(iseg), zray(iseg), &
                                      tol, pert)
                     ! Just finish it off
                     IF (j + 1 == nk) THEN
                        xray(iseg) = xpts(nk)
                        zray(iseg) = zpts(nk)
                        lterminate = .TRUE.
                     ENDIF
                  ENDIF
                  ! Tabulate the length
                  rayl(iseg) = DSQRT( (xpts(jbeg) - xray(iseg-1))**2 &
                                    + (zpts(jbeg) - zray(iseg-1))**2 )
                  DO l=jbeg+1,j
                     rayl(iseg) = rayl(iseg) + DSQRT( (xpts(l) - xpts(l-1))**2 &
                                                    + (zpts(l) - zpts(l-1))**2)
                  ENDDO
                  rayl(iseg) = rayl(iseg) + DSQRT( (xray(iseg) - xpts(j))**2 &
                                                 + (zray(iseg) - zpts(j))**2 )

                  ix0 = ix1
                  iz0 = iz1
                  jbeg = j + 1
                  EXIT
               ENDIF
            ENDDO
            IF (lterminate) EXIT
         ENDDO
         IF (.NOT.lterminate) THEN
            WRITE(ERROR_UNIT,915)
            CYCLE
         ENDIF
         nseg = iseg
         ALLOCATE(rayPaths(ir)%xyzs(3*nseg))
         ALLOCATE(rayPaths(ir)%l(nseg))
         ALLOCATE(rayPaths(ir)%cells(nseg))
         rayPaths(ir)%time = 0.d0
         DO iseg=1,nseg
            rayPaths(ir)%xyzs(3*(iseg-1)+1) = xray(iseg)
            rayPaths(ir)%xyzs(3*(iseg-1)+1) = 0.d0 
            rayPaths(ir)%xyzs(3*(iseg-1)+3) = zray(iseg)
            rayPaths(ir)%l(iseg) = rayl(iseg)
            rayPaths(ir)%cells(iseg) = cells(iseg)
            rayPaths(ir)%time = rayPaths(ir)%time + rayl(iseg)*slow(cells(iseg))
            print *, xray(iseg), zray(iseg), rayl(iseg)!, dsqrt( (xray(iseg+1) - xray(iseg))**2 + (zray(iseg+1)-zray(iseg))**2 )
         ENDDO 
print *, tn, tne, iseg, rayPaths(ir)%time
pause
return

!pause
      ENDDO
      IF (ALLOCATED(xpts))  DEALLOCATE(xpts)
      IF (ALLOCATED(zpts))  DEALLOCATE(zpts)
      IF (ALLOCATED(xray))  DEALLOCATE(xray)
      IF (ALLOCATED(zray))  DEALLOCATE(zray)
      IF (ALLOCATED(rayl))  DEALLOCATE(rayl)
      IF (ALLOCATED(cells)) DEALLOCATE(cells)
  900 FORMAT('fteik_rays_trace2d: Ray path failed to terminate')
  905 FORMAT('fteik_rays_trace2d: Ray path ', I0, ' failed to terminate')
  910 FORMAT('fteik_rays_trace2d: Skipping cell is not yet considered')
  915 FORMAT('fteiK_rays_trace2d: Internal error discretizing ray')
      RETURN
      END SUBROUTINE

      SUBROUTINE updateRK42D(nx, nz,               &
                             ds, dx, dz, xn, zn,   &
                             slow, tgradx, tgradz, &
                             xn1, zn1)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nx, nz
      DOUBLE PRECISION, DIMENSION(nx*nz), INTENT(IN) :: tgradx, tgradz
      DOUBLE PRECISION, DIMENSION((nz-1)*(nx-1)), INTENT(IN) :: slow
      DOUBLE PRECISION, INTENT(IN) :: ds, dx, dz, xn, zn
      DOUBLE PRECISION, INTENT(OUT) :: xn1, zn1
      DOUBLE PRECISION, DIMENSION(4,2) :: kmat
      DOUBLE PRECISION, DIMENSION(2) :: dxrk
      DOUBLE PRECISION xk, zk
      DOUBLE PRECISION, PARAMETER :: dsscal(3) = [0.5d0, 0.5d0, 1.d0]
      INTEGER i
      ! Initialize by evaluating k1 = f(\textbf{x}_n) = \nabla T(\textbf{x}_n)
      CALL evalf2d(dx, dz, xn, zn,       &
                   slow, tgradx, tgradz, &
                   kmat(1,1), kmat(1,2))
      kmat(1,1:2) = ds*kmat(1,1:2)
      CALL updatePosition2D(nx, nz, dx, dz,                     &
                            xn, zn, 1.d0, kmat(1,1), kmat(1,2), &
                            xk, zk)
      ! Update:
      ! k1 = f(x_n)
      ! k2 = f(x_n + ds/2*k1)
      ! k3 = f(x_n + ds/2*k2)
      ! k4 = f(x_n + ds*k3)
      DO i=2,4
         ! Evaluate the next gradient at the updated position
         CALL evalf2d(dx, dz, xk, zk,       &
                      slow, tgradx, tgradz, &
                      kmat(i,1), kmat(i,2))
         kmat(i,1:2) = ds*dsscal(i-1)*kmat(i,1:2)
         CALL updatePosition2D(nx, nz, dx, dz,                     &
                               xn, zn, 1.d0, kmat(i,1), kmat(i,2), &
                               xk, zk)
      ENDDO
      ! x_{n+1} = x_n + ds/6*(k1 + 2*k2 + 2*k3 + k4)
      dxrk(1) = 1.d0/6.d0*(kmat(1,1) + 2.d0*(kmat(2,1) + kmat(3,1)) + kmat(4,1))
      dxrk(2) = 1.d0/6.d0*(kmat(1,2) + 2.d0*(kmat(2,2) + kmat(3,2)) + kmat(4,2))
      CALL updatePosition2D(nx, nz, dx, dz,                 &
                            xn, zn, 1.d0, dxrk(1), dxrk(2), &
                            xn1, zn1)
      RETURN
      END
!>    @brief Finds the next point of the ray by solving the ray equation
!>           \f$
!>               \frac{d \textbf{x}}{ds}
!>             =-\frac{1}{u(\textbf{x}_n)} \nabla T(\textbf{x})
!>           \$
!>           using the Euler method, i.e., 
!>           \f$
!>               \textbf{x}_{n+1}
!>             = \textbf{x}_n
!>             - ds \frac{1}{u(\textbf{x})_n} \nabla T(\textbf{x}_n)
!>           \f$.  For more, see https://en.wikipedia.org/wiki/Euler_method.
!>    @param[in] nx      Number of x grid points.
!>    @param[in] nz      Number of z grid points.
!>    @param[in] ds      The ray step length (meters).
!>    @param[in] dx      The grid spacing in x (meters).
!>    @param[in] dz      The grid spacing in z (meters).
!>    @param[in] xn      The starting position of the ray in x (meters). 
!>    @param[in] zn      The starting position of the ray in z (meters).
!>    @param[in] slow    The slowness in each cell (s/m).
!>    @param[in] tgradx  The travel time gradient in x (s/m).
!>    @param[in] tgradz  The travel time gradient in z (s/m).
!>    @param[out] xn1    The updated ray position in x (meters).
!>    @parma[out] zn1    The updated ray position in z (meters).
!>    @ingroup rays
      SUBROUTINE updateEuler2D(nx, nz,               &
                               ds, dx, dz, xn, zn,   &
                               slow, tgradx, tgradz, &
                               xn1, zn1)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nx, nz
      DOUBLE PRECISION, DIMENSION(nx*nz), INTENT(IN) :: tgradx, tgradz
      DOUBLE PRECISION, DIMENSION((nz-1)*(nx-1)), INTENT(IN) :: slow
      DOUBLE PRECISION, INTENT(IN) :: ds, dx, dz, xn, zn
      DOUBLE PRECISION, INTENT(OUT) :: xn1, zn1
      DOUBLE PRECISION dx1ds, dx3ds
      ! Evaluate the gradient at the current position \{x_n, z_n \}
      CALL evalf2d(dx, dz, xn, zn,       &
                   slow, tgradx, tgradz, &
                   dx1ds, dx3ds)
      CALL updatePosition2D(nx, nz, dx, dz,           &
                            xn, zn, ds, dx1ds, dx3ds, &
                            xn1, zn1)
      RETURN
      END
!========================================================================================!
!>    @brief Computes the update according to the midpoint rule
!>    \f$
!>        \textbf{x}_{n+1}
!>      = \textbf{x}_n
!>      + ds \nabla T \left ( \textbf{x}_n + \frac{1}{2} ds \nabla T(\textbf{x}_n) \right)
!>    \f$
!>    @param[in] nx      Number of grid points in x.
!>    @param[in] nz      Number of grid points in z.
!>    @param[in] ds      The proposed length of the ray in meters.
!>    @param[in] dx      The grid spacing (meters) in x.
!>    @param[in] dz      The grid spacing (meters) in z.
!>    @param[in] xn      The current position (meters) in x.
!>    @param[in] zn      The current position (meters) in z.
!>    @param[in] slow    The slowness (s/m) in each cell.
!>    @param[in] tgradx  The gradient of the velocity in x (m/s) at each node.
!>    @param[in] tgradz  The gradient of the velocity in z (m/s) at each node.
!>    @param[out] xn1    The updated position (meters) in x.
!>    @param[out] zn1    The updated position (meters) in z.
!>    @ingroup rays
      SUBROUTINE updateMidpoint2D(nx, nz,               &
                                  ds, dx, dz, xn, zn,   &
                                  slow, tgradx, tgradz, &
                                  xn1, zn1) 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nx, nz
      DOUBLE PRECISION, DIMENSION(nx*nz), INTENT(IN) :: tgradx, tgradz 
      DOUBLE PRECISION, DIMENSION((nz-1)*(nx-1)), INTENT(IN) :: slow
      DOUBLE PRECISION, INTENT(IN) :: ds, dx, dz, xn, zn
      DOUBLE PRECISION, INTENT(OUT) :: xn1, zn1 
      DOUBLE PRECISION ds2, dx1ds, dx3ds, xnew, znew
      ! Evaluate the gradient at the current position \{x_n, z_n \}
      CALL evalf2d(dx, dz, xn, zn,       &
                   slow, tgradx, tgradz, &
                   dx1ds, dx3ds) 
      ! Compute the midpoint \textbf{x}_m = \frac{1}{2} ds \nabla \nabla T(\textbf{x})
      ds2 = 0.5d0*ds
      CALL updatePosition2D(nx, nz, dx, dz,            &
                            xn, zn, ds2, dx1ds, dx3ds, &
                            xnew, znew)
      ! Compute the new position:
      ! \f$ \textbf{x}_{n+1} = \textbf{x}_n + ds \nabla T (\textbf{x}_m)
      CALL evalf2d(dx, dz, xnew, znew,   &
                   slow, tgradx, tgradz, &
                   dx1ds, dx3ds)
      CALL updatePosition2D(nx, nz, dx, dz,           &
                            xn, zn, ds, dx1ds, dx3ds, &
                            xn1, zn1)
      RETURN
      END

      PURE SUBROUTINE updatePosition2D(nx, nz, dx, dz,           &
                                       xn, zn, ds, dx1ds, dx3ds, &
                                       xn1, zn1)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nx, nz
      DOUBLE PRECISION, INTENT(IN) :: dx, dz, xn, zn, ds, dx1ds, dx3ds
      DOUBLE PRECISION, INTENT(OUT) :: xn1, zn1 
      DOUBLE PRECISION, PARAMETER :: eps = EPSILON(1.d0) 
      xn1 = xn + ds*dx1ds
      zn1 = zn + ds*dx3ds
      xn1 = MIN(MAX(eps, xn1), DBLE(nx-1)*dx - eps)
      zn1 = MIN(MAX(eps, zn1), DBLE(nz-1)*dz - eps)
      RETURN
      END
 
!>    @brief Evaluates the right hand side of the ray equations (Shearer pg. 85)
!>           \f$ 
!>               \frac{d \textbf{x}}{ds} =-\frac{\textbf{s}(\textbf{x})}{u(\textbf{x})}
!>           \f$
!>           where \f$ \textbf{x} \f$ is the unknown position, 
!>           \f$ \textbf{s} = \nabla T \f$ is the slowness vector, 
!>           and \f$ u \f$ is the slowness (inverse of seismic velocity). 
!>           Note, that this marches down the gradient because of the minus sign and
!>           is only appropriate for receiver-to-source ray tracing.
!>    @param[in] dx      Grid spacing in x (meters).
!>    @param[in] dz      Grid spacing in z (meters).
!>    @param[in] xn      The x coordinate of the current position (meters).
!>    @param[in] zn      The z coordinate of the current position (meters).
!>    @param[in] slow    The slowness (seconds/meter) evaluated in each cell.
!>    @param[in] tgradx  The travel time gradient in x (seconds/meter) evaluated at 
!>                       each grid point of the model.
!>    @param[in] tgradz  The travel time gradient in z (seconds/meter) evaluated at
!>                       each grid point of the model.
!>    @param[out] dx1ds  The change in x for a given change in ray length, ds.
!>    @param[out] dx3ds  The change in z for a given change in ray length, ds.
!>    @ingroup rays 
      SUBROUTINE evalf2d(dx, dz, xn, zn,       &
                         slow, tgradx, tgradz, & 
                         dx1ds, dx3ds)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: dx, dz, xn, zn
      DOUBLE PRECISION, INTENT(oUT) :: dx1ds, dx3ds
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: slow
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tgradx, tgradz 
      DOUBLE PRECISION, DIMENSION(2,2) :: tx, tz
      DOUBLE PRECISION dtdx, dtdz, x1, x2, z1, z2
      INTEGER icell, igrd11, igrd12, igrd21, igrd22, ix, iz
      ix = MIN(MAX(1, INT(xn/dx + 1.d0)), nx - 1)
      iz = MIN(MAX(1, INT(zn/dz + 1.d0)), nz - 1)
      icell = (ix - 1)*(nz - 1) + iz
      x1 = DBLE(ix-1)*dx; x2 = DBLE(ix)*dx
      z1 = DBLE(iz-1)*dz; z2 = DBLE(iz)*dz
      ! Find the position of the current point
      igrd11 = (ix - 1)*nz + iz
      igrd12 = (ix - 1)*nz + iz + 1
      igrd21 = (ix - 0)*nz + iz
      igrd22 = (ix - 0)*nz + iz + 1
      tx(1,1) = tgradx(igrd11)
      tx(1,2) = tgradx(igrd12)
      tx(2,1) = tgradx(igrd21)
      tx(2,2) = tgradx(igrd22)
      tz(1,1) = tgradz(igrd11)
      tz(1,2) = tgradz(igrd12)
      tz(2,1) = tgradz(igrd21)
      tz(2,2) = tgradz(igrd22)
      ! Interpolate the derivatives at the point
      dtdx = bilinear(xn, zn, x1, z1, x2, z2, tx)
      dtdz = bilinear(xn, zn, x1, z1, x2, z2, tz)
      ! Step down the gradient
      dx1ds =-dtdx/slow(icell)
      dx3ds =-dtdz/slow(icell) 
      RETURN
      END 

      DOUBLE PRECISION FUNCTION bilinear(x, z, x1, z1, x2, z2, t)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: t
      DOUBLE PRECISION, INTENT(IN) :: x, z, x1, z1, x2, z2
      DOUBLE PRECISION a0, a1, a2, a3, den
      DOUBLE PRECISION c11, c12, c21, c22
      den = (x1 - x2)*(z1 - z2)
      c11 = t(1,1)
      c12 = t(1,2)
      c21 = t(2,1)
      c22 = t(2,2)
      a0 = ( c11*x2*z2 - c12*x2*z1 - c21*x1*z2 + c22*x1*z1)/den
      a1 = (-c11*z2    + c12*z1    + c21*z2    - c22*z1   )/den
      a2 = (-c11*x2    + c12*x2    + c21*x1    - c22*x1   )/den
      a3 = ( c11       - c12       - c21       + c22      )/den
      bilinear = a0 + a1*x + a2*z + a3*x*z
      RETURN
      END

      SUBROUTINE intersect2d(x0, z0,         &
                             dx1ds, dx3ds,   &
                             dxgrid, dzgrid, &
                             x1, z1,         &
                             tol, pert)
      IMPLICIT NONE 
      DOUBLE PRECISION, INTENT(IN) :: x0, z0, dx1ds, dx3ds, dxgrid, dzgrid, tol, pert 
      DOUBLE PRECISION, INTENT(OUT) :: x1, z1
      DOUBLE PRECISION dxgridAbs, dzgridAbs, dxdz, dzdx, dxprop, dzprop, &
                       x1est, x2est, z1est, z2est, xdir, zdir 
      INTEGER ix, iz
      ! Compute cell index of point
      dxgridAbs = DABS(dxgrid)
      dzgridAbs = DABS(dzgrid)
      ix = INT(x0/dxgridAbs + 1.d0)
      iz = INT(z0/dzgridAbs + 1.d0)
      xdir = 1.d0 
      IF (dx1ds < 0.d0) xdir =-1.d0
      zdir = 1.d0 
      IF (dx3ds < 0.d0) zdir =-1.d0
      ! The proposed move is in increasing x.  pert/dxgrid normalizes for highly
      ! anisotropic grid spacings.
      IF (xdir > 0.d0) THEN
         dxprop = (DBLE(ix)*dxgridAbs - x0) + pert!/dxgridAbs ! will push into next cell
      ELSE ! Decreasing x
         dxprop = (x0 - DBLE(ix-1)*dxgridAbs) + pert!/dxgridAbs ! will push into next cell
      ENDIF
      ! The proposed move is in increasing z.  pert/dzgrid normalizes for highly
      ! anisotropic grid spacings.
      IF (zdir > 0.d0) THEN
         dzprop = (DBLE(iz)*dzgridAbs - z0) + pert!/dzgridAbs ! will push into next cell
      ELSE
         dzprop = (z0 - DBLE(iz-1)*dzgridAbs) + pert!/dzgridAbs ! will push into next cell
      ENDIF
      ! Edge cases - Vertical line 
      IF (DABS(1.d0 - DABS(dx3ds)) < tol) THEN
         ! Vertical line increasing positive down in solver coordinates
         x1 = x0
         z1 = z0 + zdir*dzprop
      ! Edge case - Horizontal line
      ELSE IF (DABS(1.d0 - DABS(dx1ds)) < tol) THEN
         ! Horizontal line increasing positive right in solver coordinates
         x1 = x0 + xdir*dxprop
         z1 = z0
      ELSE
         ! Make the move to the point that produces the smallest change in travel time.
         dzdx = dx3ds/dx1ds ! dz/dx = (ds/dx)/(ds/dz) = (ds*dz/ds*dx) = dz/dx
         x1est = x0 + xdir*dxprop
         z1est = z0 + zdir*DABS(dzdx*dxprop) ! dz/dx = Delta z/Delta x -> Delta z = Delta x*dz/dx

         dxdz = dx1ds/dx3ds
         x2est = x0 + xdir*DABS(dxdz*dzprop) ! dx/dz = Delta x/Delta z -> Delta x = Delta z*dx/dz
         z2est = z0 + zdir*dzprop
         IF ((x1est - x0)**2 + (z1est - z0)**2 < &
             (x2est - x0)**2 + (z2est - z0)**2) THEN
            x1 = x1est
            z1 = z1est
         ELSE
            x1 = x2est
            z1 = z2est
         ENDIF
      ENDIF
      RETURN
      END

!>    @brief Computes the gradient of the bilinear interpolant:
!>           \f$ f(x,y,z) = a_0 + a_1 x + a_2 y + a_3 z 
!>                        + a_4 x y  + a_5 x z + a_6 y z + a_7 x y z \f$.
!>           where the coefficients are given by:
!>           https://en.wikipedia.org/wiki/Trilinear_interpolation#Alternative_algorithm
!>    @param[in] x         x position (meters).
!>    @param[in] y         y position (meters).
!>    @param[in] z         z position (meters).
!>    @param[in] x0        North, west, upper grid point (meters).
!>    @param[in] y0        North, west, upper grid point (meters).
!>    @param[in] z0        North, west, upper grid point (meters).
!>    @param[in] x1        South, east, lower grid point (meters).
!>    @param[in] y1        South, east, lower grid point (meters).
!>    @param[in] z1        South, east, lower grid point (meters).
!>    @param[in] t         Travel-times (seconds) at grid points where the first index
!>                         corresponds to x0 and x1, respectively, and the second index
!>                         corresponds to y0 and y1, respectively, and the third index
!>                         corresponds to z0 and z1, respectively.
!>    @param[out] tint     Interpolated travel time (seconds).
!>    @param[out] dx1ds    Direction cosine in x1 direction.
!>    @param[out] dx2ds    Direction cosine in x2 direction.
!>    @param[out] dx3ds    Direction cosine in x3 direction.
!>    @param[in] lnegGrad  If true then compute the negative of the gradient to move
!>                         in the direction of max descent of the travel time field.
!>    @param[in] lnegGrad  If false then move in the direction of max increase of the
!>                         travel time field.
!>    @ingroup rays
      SUBROUTINE trilinearGrad(x, y, z,                   &
                               x0, y0, z0,                &
                               x1, y1, z1, t,             &
                               tint, dx1ds, dx2ds, dx3ds, &
                               lnegGrad)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: t
      DOUBLE PRECISION, INTENT(IN) :: x, y, z, x0, y0, z0, x1, y1, z1
      DOUBLE PRECISION, INTENT(OUT) :: tint, dx1ds, dx2ds, dx3ds
      LOGICAL, INTENT(IN) :: lnegGrad
      DOUBLE PRECISION c000, c100, c010, c110, c001, c101, c011, c111
      DOUBLE PRECISION a0, a1, a2, a3, a4, a5, a6, a7, den, dtdx, dtdy, dtdz, xnorm
      den = (x0 - x1)*(y0 - y1)*(z0 - z1)
      c000 = t(1,1,1)
      c100 = t(2,1,1)
      c010 = t(1,2,1)
      c110 = t(2,2,1)
      c001 = t(1,1,2)
      c101 = t(2,1,2)
      c011 = t(1,2,2)
      c111 = t(2,2,2)
      a0 =(-c000*x1*y1*z1 + c001*x1*y1*z0 + c010*x1*y0*z1 - c011*x1*y0*z0 &
          + c100*x0*y1*z1 - c101*x0*y1*z0 - c110*x0*y0*z1 + c111*x0*y0*z0)/den
      a1 =( c000*y1*z1    - c001*y1*z0    - c010*y0*z1    + c011*y0*z0    &
          - c100*y1*z1    + c101*y1*z0    + c110*y0*z1    - c111*y0*z0)/den
      a2 =( c000*x1*z1    - c001*x1*z0    - c010*x1*z1    + c011*x1*z0    &
          - c100*x0*z1    + c101*x0*z0    + c110*x0*z1    - c111*x0*z0)/den
      a3 =( c000*x1*y1    - c001*x1*y1    - c010*x1*y0    + c011*x1*y0    &
          - c100*x0*y1    + c101*x0*y1    + c110*x0*y0    - c111*x0*y0)/den
      a4 =(-c000*z1       + c001*z0       + c010*z1       - c011*z0       &
          + c100*z1       - c101*z0       - c110*z1       + c111*z0)/den
      a5 =(-c000*y1       + c001*y1       + c010*y0       - c011*y0       &
          + c100*y1       - c101*y1       - c110*y0       + c111*y0)/den
      a6 =(-c000*x1       + c001*x1       + c010*x1       - c011*x1       &
          + c100*x0       - c101*x0       - c110*x0       + c111*x0)/den
      a7 =( c000          - c001          - c010          + c011          &
          - c100          + c101          + c110          - c111)/den
      ! f(x,y,z) = a0 + a1*x + a2*y + a3*z + a4*x*y + a5*x*z + a6*y*z + a7*x*y*z
      tint = a0 + a1*x + a2*y + a3*z + a4*x*y + a5*x*z + a6*y*z + a7*x*y*z
      dtdx = a1 + a4*y + a5*z + a7*y*z
      dtdy = a2 + a4*x + a6*z + a7*x*z
      dtdz = a3 + a5*x + a6*y + a7*x*y
      ! Normalized vector requires hypotenuse is unity.
      xnorm = DSQRT(dtdx*dtdx + dtdy*dtdy + dtdz*dtdz)
      ! Compute direction cosines of Lay and Wallace 3.16 so that:
      ! (dx1/ds)^2 + (dx2/ds)^2 + (dx3/ds)^2 = 1^2
      IF (lnegGrad) THEN
         dx1ds =-dtdx/xnorm
         dx2ds =-dtdy/xnorm
         dx3ds =-dtdz/xnorm
      ELSE
         dx1ds = dtdx/xnorm
         dx2ds = dtdy/xnorm
         dx3ds = dtdz/xnorm
      ENDIF
      RETURN
      END

END MODULE
