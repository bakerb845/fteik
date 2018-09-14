  INTEGER ierr
  CALL test_rays2d(ierr) 
  IF (ierr /= 0) THEN
     ERROR STOP 'Error calling rays2d'
  ENDIF
  STOP
  END

      SUBROUTINE test_rays2d(ierr)
      USE ISO_FORTRAN_ENV
      USE FTEIK2D_SOLVER64F, ONLY : fteik_solver2d_initialize64f,            &
                                    fteik_solver2d_setNodalVelocityModel64f, &
                                    fteik_solver2d_setSources64f,            &
                                    fteik_solver2d_solveSourceLSM,           &
                                    fteik_solver2d_free
      USE FTEIK2D_SOLVER64F, ONLY : FTEIK_ZX_ORDERING
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ierr
      INTEGER, PARAMETER :: nx = 101
      INTEGER, PARAMETER :: nz = 51
      INTEGER, PARAMETER :: ncellx = nx - 1
      INTEGER, PARAMETER :: ncellz = nz - 1
      DOUBLE PRECISION, PARAMETER :: x0 = 0.d0
      DOUBLE PRECISION, PARAMETER :: z0 = 0.d0
      DOUBLE PRECISION, PARAMETER :: x1 = 5.d3
      DOUBLE PRECISION, PARAMETER :: z1 = 2.5d3
      DOUBLE PRECISION dx, dz
      DOUBLE PRECISION, PARAMETER :: eps = 3.d0
      DOUBLE PRECISION, PARAMETER :: conv = 1.d-8
      DOUBLE PRECISION, ALLOCATABLE :: vel(:)
      INTEGER, PARAMETER :: ngs_sweeps = 3
      INTEGER, PARAMETER :: verb = 3
      INTEGER, PARAMETER :: nsrc = 1
      DOUBLE PRECISION z
      DOUBLE PRECISION xsrc(nsrc), zsrc(nsrc)
      INTEGER isrc, iz, ix
      ierr = 0
      dx = (x1 - x0)/DBLE(nx - 1)
      dz = (z1 - z0)/DBLE(nz - 1)
      CALL fteik_solver2d_initialize64f(nz, nx,           &
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
      DO iz=1,nz
         DO ix=1,nx
            z = DBLE(iz - 1)*dz
            vel((ix-1)*nz + iz) = 3000.d0 + (5000.d0 - 3000.d0)/(z1 - z0)*(z - z0)
         ENDDO
      ENDDO
      CALL fteik_solver2d_setNodalVelocityModel64f(nx*nz, FTEIK_ZX_ORDERING, vel, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,*) 'test_rays2d: Failed to make vel'
         RETURN
      ENDIF
      ! TODO put source where analytic computation can be done
      xsrc(1) = (x1 + x0)/2.d0
      zsrc(1) = (z1 + z0)/2.d0
      CALL fteik_solver2d_setSources64f(nsrc, zsrc, xsrc, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,*) 'test_rays2d: Failed to set sources'
         RETURN
      ENDIF
      DO isrc=1,nsrc
         CALL fteik_solver2d_solveSourceLSM(isrc, ierr) 
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,*) 'test_rays2d: Failed to solve source', isrc 
            RETURN 
         ENDIF
      ENDDO
 
      CALL fteik_solver2d_free()
      IF (ALLOCATED(vel)) DEALLOCATE(vel)
      RETURN
      END SUBROUTINE
