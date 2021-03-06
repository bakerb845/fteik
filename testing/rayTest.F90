  INTEGER ierr
  CALL test_rays3d(ierr)
  IF (ierr /= 0) THEN
     ERROR STOP 'Error calling rays3d'
  ENDIF
stop
  CALL test_rays2d(ierr) 
  IF (ierr /= 0) THEN
     ERROR STOP 'Error calling rays2d'
  ENDIF
  STOP
  END

      SUBROUTINE test_rays3d(ierr)
      USE ISO_FORTRAN_ENV
      USE ISO_C_BINDING
      USE FTEIK3D_SOLVER64F, ONLY : fteik_solver3d_initialize, &
                                    fteik_solver3d_free,       &
                                    fteik_solver3d_setNodalVelocityModel64f, &
                                    fteik_solver3d_setSources64f, &
                                    fteik_solver3d_setReceivers64f, &
                                    fteik_solver3d_solveSourceLSM, &
                                    fteik_solver3d_traceRaysToPoints, &
                                    FTEIK_ZXY_ORDERING
      USE FTEIK_MODEL64F, ONLY : fteik_model_grid2indexF
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ierr
      INTEGER, PARAMETER :: nx = 51
      INTEGER, PARAMETER :: ny = 51
      INTEGER, PARAMETER :: nz = 51
      DOUBLE PRECISION, PARAMETER :: x0 = 0.d0
      DOUBLE PRECISION, PARAMETER :: y0 = 0.d0
      DOUBLE PRECISION, PARAMETER :: z0 = 0.d0
      DOUBLE PRECISION, PARAMETER :: x1 = 5.d3
      DOUBLE PRECISION, PARAMETER :: y1 = 4.d3
      DOUBLE PRECISION, PARAMETER :: z1 = 2.5d3
      DOUBLE PRECISION dx, dy, dz
      DOUBLE PRECISION, PARAMETER :: conv = 1.d-10
      DOUBLE PRECISION, PARAMETER :: eps = 3.d0
      DOUBLE PRECISION, ALLOCATABLE :: vel(:)
      DOUBLE PRECISION xr(1), yr(1), zr(1)
      DOUBLE PRECISION xsrc(1), ysrc(1), zsrc(1), v0, vmax, vgrad, z
      INTEGER igrd, isrc, ix, iy, iz, ngrd
      INTEGER, PARAMETER :: ngs_sweeps = 2
      INTEGER, PARAMETER :: verb = 3
      INTEGER, PARAMETER :: nsrc = 1
      INTEGER, PARAMETER :: nrec = 1
      ierr = 0
      dx = (x1 - x0)/DBLE(nx - 1)
      dy = (y1 - y0)/DBLE(ny - 1)
      dz = (z1 - z0)/DBLE(nz - 1)
      CALL fteik_solver3d_initialize(nz, nx, ny,      &
                                     z0, x0, y0,      &
                                     dz, dx, dy,      &
                                     ngs_sweeps, eps, &
                                     conv, verb,      &
                                     ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,*) 'test_rays3d: Failed to initialize solver'
         RETURN
      ENDIF
      ! Create a velocity model with a linear gradient.
      ngrd = nx*ny*nz
      ALLOCATE(vel(ngrd)); vel(:) = 0.d0
      v0 = 3.d3
      vmax = 5.d3
      vgrad = (vmax - v0)/(z1 - z0) !(m/s)/(m) = 1/s -> this is slotnick's `a'
      DO iy=1,ny
         DO ix=1,nx
            DO iz=1,nz
               z = DBLE(iz - 1)*dz
               igrd = fteik_model_grid2indexF(iz, ix, iy, nz, nz*nx)
               vel(igrd) = v0 + vgrad*(z - z0)
            ENDDO
         ENDDO
      ENDDO
      CALL fteik_solver3d_setNodalVelocityModel64f(ngrd, FTEIK_ZXY_ORDERING, &
                                                   vel, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,*) 'test_rays3d: Failed to set velocity model'
         RETURN
      ENDIF
      ! Create a source
      xsrc(1) = x0 + (x1 + x0)/2.d0 + 5.d0
      ysrc(1) = y0 + (y1 + y0)/2.d0
      zsrc(1) = z0 !+ (z1 + z0)/2.d0
      CALL fteik_solver3d_setSources64f(nsrc, zsrc, xsrc, ysrc, ierr)
      ! Set the receiver
      xr(1) = x0 + (x1 + x0)/8.d0
      yr(1) = y0 + (y1 + y0)/8.d0
      zr(1) = z0
      CALL fteik_solver3d_setReceivers64f(nrec, zr, xr, yr, ierr)
      ! Solve
      DO isrc=1,nsrc
         CALL fteik_solver3d_solveSourceLSM(isrc, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,*) 'test_rays3d: Failed to solve source', isrc
            RETURN
         ENDIF
         CALL fteik_solver3d_traceRaysToPoints(1, zr, xr, yr, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,*) 'test_rays3d: Failed to trace rays', isrc
            RETURN
         ENDIF
      ENDDO

      CALL fteik_solver3d_free()
      DEALLOCATE(vel)
      RETURN
      END

      SUBROUTINE test_rays2d(ierr)
      USE ISO_FORTRAN_ENV
      USE FTEIK2D_SOLVER64F, ONLY : fteik_solver2d_initialize,               &
                                    fteik_solver2d_setNodalVelocityModel64f, &
                                    fteik_solver2d_setSources64f,            &
                                    fteik_solver2d_solveSourceLSM,           &
                                    fteik_solver2d_traceRaysToPoints,        &
                                    fteik_solver2d_free, &
                                    fteik_solver2d_setReceivers64f
      USE FTEIK2D_SOLVER64F, ONLY : FTEIK_ZX_ORDERING, ttimes, ttimesRec
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ierr
      INTEGER, PARAMETER :: nx = 201
      INTEGER, PARAMETER :: nz = 101
      INTEGER, PARAMETER :: ncellx = nx - 1
      INTEGER, PARAMETER :: ncellz = nz - 1
      DOUBLE PRECISION, PARAMETER :: x0 = 0.d0
      DOUBLE PRECISION, PARAMETER :: z0 = 0.d0
      DOUBLE PRECISION, PARAMETER :: x1 = 5.0d3
      DOUBLE PRECISION, PARAMETER :: z1 = 2.5d3
      DOUBLE PRECISION dx, dz
      DOUBLE PRECISION, PARAMETER :: eps = 3.d0
      DOUBLE PRECISION, PARAMETER :: conv = 1.d-10
      DOUBLE PRECISION, ALLOCATABLE :: vel(:)
      INTEGER, PARAMETER :: ngs_sweeps = 4
      INTEGER, PARAMETER :: verb = 3
      INTEGER, PARAMETER :: nsrc = 1
      DOUBLE PRECISION z
      DOUBLE PRECISION xsrc(nsrc), zsrc(nsrc), xr(1), zr(1)
      DOUBLE PRECISION arg, arg1, arg2, xhalf, v0, t, vgrad, vmax, p, dmax, aoi, eta
      INTEGER isrc, iz, ix
      ierr = 0
      dx = (x1 - x0)/DBLE(nx - 1)
      dz = (z1 - z0)/DBLE(nz - 1)
      CALL fteik_solver2d_initialize(nz, nx,           &
                                     z0, x0,           &
                                     dz, dx,           &
                                     ngs_sweeps,  eps, &
                                     conv, verb,       &
                                     ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,*) 'test_rays2d: Failed to initialize solver'
         RETURN
      ENDIF
      ! Create a velocity model with a linear gradient.
      ALLOCATE(vel(nx*nz)); vel(:) = 0.d0
      v0 = 3.d3
      vmax = 5.d3
      vgrad = (vmax - v0)/(z1 - z0) !(m/s)/(m) = 1/s -> this is slotnick's `a'
      DO iz=1,nz
         DO ix=1,nx
            z = DBLE(iz - 1)*dz
            vel((ix-1)*nz + iz) = v0 + vgrad*(z - z0)
         ENDDO
      ENDDO
      CALL fteik_solver2d_setNodalVelocityModel64f(nx*nz, FTEIK_ZX_ORDERING, vel, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,*) 'test_rays2d: Failed to make vel'
         RETURN
      ENDIF
      ! TODO put source where analytic computation can be done
print *, dx, dz
      xsrc(1) = x0 + (x1 + x0)/2.d0 + 5.d0
      zsrc(1) = z0 !+ (z1 + z0)/2.d0
      CALL fteik_solver2d_setSources64f(nsrc, zsrc, xsrc, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,*) 'test_rays2d: Failed to set sources'
         RETURN
      ENDIF
      xr(1) = x0 + (x1 + x0)/8.d0
      zr(1) = z0
      CALL fteik_solver2d_setReceivers64f(1, zr, xr, ierr)
      DO isrc=1,nsrc
         CALL fteik_solver2d_solveSourceLSM(isrc, ierr) 
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,*) 'test_rays2d: Failed to solve source', isrc 
            RETURN 
         ENDIF
         CALL fteik_solver2d_traceRaysToPoints(1, xr, zr, ierr)
      ENDDO
! iteratively solve for ray parameter
arg = 1.d0 + 0.5d0/v0*(1.d0/v0)*(vgrad**2)*( (xr(1) - xsrc(1))**2 )
t = acosh(arg)/vgrad
print *, 'analytic, estimated', t, ttimesREc
do ix=1,20
  aoi = dble(ix)*3.14159267d0/180.d0/4.d0
  p = dsin(aoi)/v0
  eta = dsqrt(1.d0/(v0)**2 - p**2)
  t = ( dlog( (1.d0/v0 + eta)/p ) - eta/v0 )/vgrad + eta/(vgrad/v0)
  eta = dsqrt(1.d0/(vmax)**2 - p**2)
  t = t - ( dlog( (1.d0/vmax + eta)/p ) - eta/vmax )/vgrad - eta/(vgrad/vmax)
! arg2 = (1.d0/vmax + eta)/p - eta/vmax
! print *, t, arg, p, aoi
enddo
t = 1.d0/vgrad*dlog(arg)

print *, 'toa', datan2(634.7133d0 - 624.d0, 2.394d0)*180.d0/3.14159267d0
p = dsin(datan2(634.7133d0 - 624.d0, 2.394d0))/v0
xhalf = (xsrc(1) + xr(1))/2.d0
dmax = (1.d0 - p*v0)/(vgrad*p)
print *, dmax, v0/vgrad, 1.d0/(vgrad*p), 'hmax', 1.d0/(vgrad*p) - v0/vgrad, xhalf
print *, 'do this anyway'
do iz=1,nz
 do ix=1,nx
   write(44,900) x0 +(ix-1)*dx, z0 +(iz-1)*dz, ttimes((ix-1)*nz+iz)
 enddo
 write(44,*)
enddo
900 format(3e14.6)
      CALL fteik_solver2d_free()
      IF (ALLOCATED(vel)) DEALLOCATE(vel)
      RETURN
      END SUBROUTINE
