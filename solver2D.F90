MODULE FTEIK2D_SOLVER64F
  USE FTEIK_CONSTANTS64F, ONLY : zero
  USE ISO_C_BINDING
  IMPLICIT NONE
  !> Holds the travel-times (seconds).  This has dimension [ngrd].
  REAL(C_DOUBLE), PROTECTED, ALLOCATABLE, SAVE :: ttimes(:)
  !DIR$ ATTRIBUTES ALIGN: 64 :: ttimes
  LOGICAL(C_BOOL), PROTECTED, ALLOCATABLE, SAVE ::         &
                   lupd1(:), lupd2(:), lupd3(:), lupd4(:), &
                   lupd5(:), lupd6(:), lupd7(:), lupd8(:)
  !> Defines the number of levels.  For 2D that is nz + nx
  INTEGER(C_INT), PROTECTED, SAVE :: nLevels = 0
  !> Defines the max level size.
  INTEGER(C_INT), PROTECTED, SAVE :: maxLevelSize = 0
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
  REAL(C_DOUBLE), SAVE, PRIVATE :: dx, dz, dz2i, dx2i, dz2i_dx2i, dz2i_p_dx2i, dz2i_p_dx2i_inv
  CONTAINS

      SUBROUTINE fteik_solver2d_initialize64fF(nzIn, nxIn,            &
                                               z0In, x0In,            &
                                               dzIn, dxIn,            &   
                                               nsweepIn, epsIn, ierr) &
                 BIND(C, NAME='fteik_solver2d_initialize64fF')
      USE ISO_C_BINDING
      USE FTEIK_MEMORY, ONLY : padLength64F
      USE FTEIK_MODEL64F, ONLY : fteik_model_intializeGeometryF, ngrd
      USE FTEIK_CONSTANTS64F, ONLY : one, zero
      IMPLICIT NONE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x0In, z0In
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dxIn, dzIn
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: epsIn
      INTEGER(C_INT), VALUE, INTENT(IN) :: nxIn, nzIn
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsweepIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      WRITE(*,*) 'fteik_solver_initialize64fF: Initializing...'
      CALL fteik_solver2d_finalizeF()
      CALL fteik_model_intializeGeometryF(lis3d,            &
                                          nzIn, nxIn, 1,    &
                                          dzIn, dxIn, one,  &
                                          z0In, x0In, zero, &
                                          ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver2d_initialize64fF: Error setting grid geometry'
         RETURN
      ENDIF
      CALL fteik_solver2d_setNumberOfSweepsF(nsweepIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver2d_initialize64fF: Failed to set max number of sweeps'
         RETURN
      ENDIF
!     CALL fteik_solver_setSphereToCartEpsilonF(epsIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver2d_initialize64fF: Failed to set epsilon'
         RETURN
      ENDIF
!     CALL fteik_solver_computeGraphF(ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver2d_initialize64fF: Failed to compute graph'
         RETURN
      ENDIF
      ! Compute the number of levels
      nLevels = nxIn + nzIn
      maxLevelSize = padLength64F(alignment, MIN(nxIn, nzIn))
      ! Set space for the travel-times
      IF (.NOT.ALLOCATED(ttimes)) ALLOCATE(ttimes(ngrd))
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases all memory associated with the solver and resets all variables.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solver2d_finalizeF()             &
      BIND(C, NAME='fteik_solver2d_finalizeF')
      USE FTEIK_RECEIVER64F, ONLY : fteik_receiver_finalizeF
!     USE FTEIK_SOURCE64F, ONLY : fteik_source_finalizeF
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      lhaveTimes = .FALSE.
      epsS2C = zero
      nsweep = 0
      nLevels = 0
      IF (ALLOCATED(ttimes)) DEALLOCATE(ttimes)
      CALL fteik_receiver_finalizeF()
!     CALL fteik_source_finalizeF()
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
!>    @brief Solves the eikonal equation with the level-set fast sweeping method.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver2d_solveSourceLSMF(isrc, ierr)   &
                 BIND(C, NAME='fteik_solver2d_solveSourceLSMF')
      USE FTEIK_SOURCE64F, ONLY : nsrc
      USE FTEIK_MODEL64F, ONLY : lhaveModel, nx, nz, slow
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, TRUE, FALSE
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) ts4(4), t0
      INTEGER(C_INT) dest(4), i1, i2, kiter, level, sweep
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
      CALL fteik_localSolver2D_initialize64fF(isrc, dest, ts4, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver2d_solveSourceLSMF: Error setting mesh constants'
         GOTO 500
      ENDIF
      CALL CPU_TIME(t0)
      ttimes(:) = FTEIK_HUGE
      ttimes(dest(1:4)) = ts4(1:4)
      ! Number of Gauss-Seidel iterations
      DO kiter=1,nsweep
         ! Loop on sweeping directions
         DO sweep=1,4
            ! Loop on the levels
            DO level=1,nlevels
               i1 = MAX(1, level - nz) 
               i2 = MIN(nx, level - 1)
!              ! Get the update nodes
!              CALL fteik_solver2D_isUpdateNodeF(level, sweep, nx, nz, &
!                                                i1, i2, lupd)
!              ! Prefetch the slowness
!              CALL fteik_solver2D_prefetchSlowness64fF(level, sweep, nx, nz,    &
!                                                        i1, i2, lupd, slow, sloc)
!              ! Prefetch travel times
!              CALL fteik_solver2D_prefetchTravelTimes64fF(level, sweep, nx, nz, &
!                                                          i1, i2, lupd, ttimes, &
!                                                          ttvec)
            ENDDO
         ENDDO
      ENDDO
      lhaveTimes = .TRUE.
  500 CONTINUE
      IF (ierr /= 0) THEN
         ttimes(:) = FTEIK_HUGE
         lhaveTimes = .FALSE.
      ENDIF
      DO kiter=1,nsweep

      ENDDO
      RETURN
      END

      SUBROUTINE fteik_localSolver2D_prefetchSlowness( )

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
      PURE SUBROUTINE fteik_localSolver2D_noInit64fF(n, ttvec, sloc, tupd) &
      BIND(C, NAME='fteik_localSolver2D_noInit64fF')
      USE ISO_C_BINDING
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, four, chunkSize
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: n
      REAL(C_DOUBLE), INTENT(IN) :: ttvec(4*n), sloc(n)
      REAL(C_DOUBLE), INTENT(OUT) :: tupd(n)
      REAL(C_DOUBLE) four_sref2, sref2, t1_2d, t12min, ta, tb, tab, tab2, temtv, &
                     te, tev, tv, tt
      INTEGER(C_INT) i
      ! Loop on nodes in update
      DO i=1,n
         tv  = ttvec(4*(i-1)+1)
         te  = ttvec(4*(i-1)+2)
         tev = ttvec(4*(i-1)+3)
         tt  = ttvec(4*(i-1)+4)
         temtv = te - tv
         ! 1D operators (refracted times)
         t12min = MIN(tv + dz*sloc(i), te + dx*sloc(i))
         ! 2D operator
         t1_2d = FTEIK_HUGE
         IF (temtv < dz*sloc(i) .AND. -temtv < dx*sloc(i)) THEN
            sref2 = sloc(i)*sloc(i)
            ta = tev + temtv
            tb = tev - temtv
            tab = ta - tb
            tab2 = tab*tab
            four_sref2 = four*sref2
            t1_2d = ( (tb*dz2i + ta*dx2i) &
                     + SQRT(four_sref2*dz2i_p_dx2i - dz2i_dx2i*tab2) )*dz2i_p_dx2i_inv
         ENDIF
         t12min = MIN(t12min, t1_2d)
         tupd(i)   = MIN(tt, t12min)
      ENDDO
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Initializes parameters for the local solver and computes initial travel 
!>           times around the source.
!>
      SUBROUTINE fteik_localSolver2D_initialize64fF(isrc, dest, ts4, ierr) &
      BIND(C, NAME='fteik_localSolver2D_initialize64fF')
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
      REAL(C_DOUBLE) dz2, dx2, dzi, dxi, szero, szero2
      REAL(C_DOUBLE) zsa, xsa, ysa
      INTEGER(C_INT) zsi, xsi, ysi 
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
         WRITE(*,*) 'fteik_localSolver2D_initialize64fF: Problem with source'
         RETURN
      ENDIF
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
      ts4(1) = fteik_localSolver2D_tAna64fF(zsi,   xsi,   dz, dx, zsa, xsa, szero)
      ts4(2) = fteik_localSolver2D_tAna64fF(zsi+1, xsi,   dz, dx, zsa, xsa, szero)
      ts4(3) = fteik_localSolver2D_tAna64fF(zsi,   xsi+1, dz, dx, zsa, xsa, szero)
      ts4(4) = fteik_localSolver2D_tAna64fF(zsi+1, xsi+1, dz, dx, zsa, xsa, szero)
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
!>    @author Ben Baker
!>
!>    @copyright MIT 
!>
      SUBROUTINE fteik_solver2d_setNumberOfSweepsF(nsweepIn, ierr) &
      BIND(C, NAME='fteik_solver2d_setNumberOfSweepsF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsweepIn 
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0 
      nsweep = 0 
      IF (nsweepIn < 0) THEN
         WRITE(*,*) 'fteik_solver2d_setNumberOfSweepsF: nsweep must be positive', nsweep
         ierr = 1 
         RETURN
      ENDIF
      nsweep = nsweepIn
      RETURN
      END SUBROUTINE
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
      PURE REAL(C_DOUBLE)                                   &
      FUNCTION fteik_localSolver2D_tAna64fF(i, j, dz, dx,   &
                                           zsa, xsa, szero) &
      BIND(C, NAME='fteik_localSolver2D_tAna64fF')
      !$OMP DECLARE SIMD(fteik_localSolver2D_tAna64fF) &
      !$OMP UNIFORM(dz, dx, zsa, xsa, szero) 
      USE ISO_C_BINDING
      IMPLICIT NONE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dz, dx, zsa, xsa, szero
      INTEGER(C_INT), VALUE, INTENT(IN) :: i, j 
      REAL(C_DOUBLE) diffx, diffz
      diffz = (DBLE(i) - zsa)*dz
      diffx = (DBLE(j) - xsa)*dx
      fteik_localSolver2D_tAna64fF = szero*HYPOT(diffx, diffz)
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
      PURE SUBROUTINE fteik_localSolver2D_tAnaD64fF(t_anad,           &
                                                    tzc, txc, i, j,   &
                                                    dz, dx, zsa, xsa, &
                                                    szero)            &
      BIND(C, NAME='fteik_localSolver2D_tAnaD64fF')
      !$OMP DECLARE SIMD(fteik_localSolver2D_tAnaD64fF) &
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
