program test2d
USE ISO_FORTRAN_ENV
integer ierr
call compare2d(ierr)
IF (ierr /= 0) THEN
   WRITE(ERROR_UNIT,*) 'Failed compare2d test'
   ERROR STOP
ENDIF
stop
end

      SUBROUTINE compare2d(ierr)
      USE ISO_FORTRAN_ENV
      USE FTEIK2D
      USE FTEIK2D_SOLVER64F, ONLY : fteik_solver2d_initialize,               &
                                    fteik_solver2d_setCellVelocityModel64f,  &
                                    fteik_solver2d_setSources64f,            &
                                    fteik_solver2d_solveSourceFSM,           &
                                    fteik_solver2d_solveSourceLSM,           &
                                    fteik_solver2d_traceRaysToPoints,        &
                                    fteik_solver2d_free
      USE FTEIK2D_SOLVER64F, ONLY : FTEIK_ZX_ORDERING, ttimes, tgradx, tgradz
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
      DOUBLE PRECISION, PARAMETER :: xsrc(1) = (x1 - x0)/3.d0
      DOUBLE PRECISION, PARAMETER :: zsrc(1) = 2.d0*(z1 - z0)/3.d0
      DOUBLE PRECISION, PARAMETER :: dx = (x1 - x0)/DBLE(nx - 1)
      DOUBLE PRECISION, PARAMETER :: dz = (z1 - z0)/DBLE(nz - 1)
      DOUBLE PRECISION, ALLOCATABLE :: slow(:,:), tt(:,:), ttgrad(:,:,:), vel(:)
      DOUBLE PRECISION dmax, z
      INTEGER, PARAMETER :: ngs_sweep = 2
      INTEGER, PARAMETER :: verb = 3 
      DOUBLE PRECISION, PARAMETER :: eps = 5.d0
      DOUBLE PRECISION, PARAMETER :: conv = 1.d-10
      INTEGER iz, ix, indx
      ierr = 0 
      ALLOCATE(vel((nx-1)*(nz-1)))
      ALLOCATE(slow(nz-1,nx-1))
      ALLOCATE(tt(nz,nx))
      ALLOCATE(ttgrad(nz,nx,2))

      DO iz=1,nz-1
         DO ix=1,nx-1
            z = DBLE(iz - 1)*dz
            vel((ix-1)*(nz-1) + iz) = 3000.d0 + (5000.d0 - 3000.d0)/(z1 - z0)*(z - z0)
            slow(iz,ix) = 1.d0/vel((ix-1)*(nz-1) + iz)
         ENDDO
      ENDDO

      CALL solver2d(slow, tt, nz, nx, zsrc(1), xsrc(1), dz, dx, ngs_sweep, ttgrad)

      CALL fteik_solver2d_initialize(nz, nx,           &
                                     z0, x0,           &
                                     dz, dx,           &
                                     ngs_sweep,  eps,  &
                                     conv, verb,       &
                                     ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,*) 'compare2d: Failed to initialize solver'
         RETURN
      ENDIF
      CALL fteik_solver2d_setCellVelocityModel64f((nz-1)*(nx-1), FTEIK_ZX_ORDERING, vel, ierr)
      CALL fteik_solver2d_setSources64f(1, zsrc, xsrc, ierr)
      CALL fteik_solver2d_solveSourceLSM(1, ierr) 
      dmax = 0.d0
      DO iz=1,nz
         DO ix=1,nx
            indx = (ix-1)*nz + iz
            dmax = MAX(dmax, ABS(tt(iz,ix) - ttimes(indx)))
write(45,850)  dble(ix-1)*dx, dble(iz-1)*dz, ttimes(indx), tgradx(indx), tgradz(indx) !ttgrad(iz,ix,1), ttgrad(iz,ix,2) !tgradx(indx), tgradz(indx)
850 format(5(e15.7,1x))
         ENDDO
write(45,*) 
      ENDDO
close(45)
      IF (dmax > 1.d-10) THEN
         WRITE(ERROR_UNIT,*) 'Max residual exceeded', dmax
         ierr = 1
      ENDIF
      print *, dmax, minval(ttimes), maxval(ttimes)

      CALL fteik_solver2d_free()

      RETURN 
      END

