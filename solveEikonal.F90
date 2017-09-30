!>    @brief This is solves the eikonal equation with the fast sweeping method.
!>
!>    @param[out] ierr   0 indicates success.
!>
!>    @author Ben Baker
!>
      SUBROUTINE fteik_solveEikonalFSMF(ierr) &
                 BIND(C, NAME='fteik_solveEikonalFSMF')
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_initializeTravelTimesNearSourceF, &
                                 fteik_setSzeroF, fteik_setEpsSolverF, &
                                 fteik_getSweepSigns64fF, fteik_getTravelTimePermF, &
                                 fteik_getSlownessPermF, fteik_getSweepLimitsF, &
                                 fteik_localSolver64fF, &
                                 fteik_meshConstants64fF, grid2indexF
      USE FTEIK_UTILS64F, ONLY : fteik_getSlownessPermF, &
                                 fteik_getTravelTimePermF, &
                                 fteik_prefetchSlowness64fF, &
                                 fteik_prefetchTravelTimes64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_prefetchSweep1TravelTimes64fF, &
                                 fteik_prefetchSweep2TravelTimes64fF, &
                                 fteik_prefetchSweep3TravelTimes64fF, &
                                 fteik_prefetchSweep4TravelTimes64fF, &
                                 fteik_prefetchSweep5TravelTimes64fF, &
                                 fteik_prefetchSweep6TravelTimes64fF, &
                                 fteik_prefetchSweep7TravelTimes64fF, &
                                 fteik_prefetchSweep8TravelTimes64fF, &
                                 fteik_prefetchSweep1Slowness64fF, &
                                 fteik_prefetchSweep2Slowness64fF, &
                                 fteik_prefetchSweep3Slowness64fF, &
                                 fteik_prefetchSweep4Slowness64fF, &
                                 fteik_prefetchSweep5Slowness64fF, &
                                 fteik_prefetchSweep6Slowness64fF, &
                                 fteik_prefetchSweep7Slowness64fF, &
                                 fteik_prefetchSweep8Slowness64fF
      USE FTEIK_LOCALSOLVER64F, ONLY : fteik_localSolver_noInit64fF
      USE FTEIK_UTILS64F, ONLY : ttimes, slow, lhaveGrid, lhaveTravelTimes, nsweep, &
                                 nx, ny, nz, nzx, nzm1, nzm1_nxm1, xsi, ysi, zsi
      USE FTEIK_UTILS64F, ONLY : dxi, dyi, dzi
      USE FTEIK_UTILS64F, ONLY : FTEIK_HUGE
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) ttWork(8), slowWork(8), sgnrz, sgnrx, sgnry, &
                     sgnrz_dzi, sgnrx_dxi, sgnry_dyi, tupd
#ifdef DEBUG
      REAL(C_DOUBLE) slowWork1(8), ttWork1(8), tupd1
#endif
      INTEGER(C_INT) slowPerm(20), ttPerm(8)
      INTEGER(C_INT) i, j, k, kndx, &
                     kiter, sgntz, sgntx, sgnty, sgnvz, sgnvx, sgnvy, sweep
      LOGICAL(C_BOOL) linitk
real*8 t0, t1
      !DIR$ ATTRIBUTES ALIGN: 64:: slowWork, ttWork
      ierr = 0
      lhaveTravelTimes = .FALSE.
      ! some easy checks
      IF (.NOT. lhaveGrid) THEN
         WRITE(*,*) 'fteik_solveEikonalFSM64fF: Grid not yet set'
         ierr = 1
         RETURN
      ENDIF
      ! don't let mesh constants go out of scope
      CALL fteik_meshConstants64fF()
      ! initialize travel times to a big number
      ttimes(:) = FTEIK_HUGE
      ! set the slowness around the source for the local solver
      CALL fteik_setSzeroF(ierr) 
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solveEikonalFSM64fF: Error initializing szero'
         RETURN
      ENDIF
      ! set the spherical to cartesian transition in time 
      CALL fteik_setEpsSolverF(ierr) 
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solveEikonalFSM64fF: Error setting epsSolver'
         RETURN
      ENDIF
      ! Initialize the points around the source
      ttimes(:) = FTEIK_HUGE
      CALL fteik_initializeTravelTimesNearSourceF(ttimes, ierr)
      IF (ierr /= 0) THEN 
         WRITE(*,*) 'fteik_solveEikonalFSM64fF: Error initializing travel-time'
         RETURN
      ENDIF
call cpu_time(t0)
      ! iterative loop
      linitk = .TRUE.
      !-----------------------------------Sweep 1----------------------------------------!
      sweep = 1
      ! Get some stencil information for this sweep
      CALL fteik_getSweepSigns64fF(sweep,               &
                                   sgntz, sgntx, sgnty, &
                                   sgnvz, sgnvx, sgnvy, &
                                   sgnrz, sgnrx, sgnry)
      CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
      CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
      ! Some derivative items
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      DO k = MAX(2,ysi),ny
         DO j = MAX(2,xsi),nx
            DO i = MAX(2,zsi),nz
               ! prefetch the slownesses and traveltimes
               CALL fteik_prefetchSweep1Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                     slow, & !slowWork)
                                                     slowWork(1), slowWork(2), &
                                                     slowWork(3), slowWork(4), &
                                                     slowWork(5), slowWork(6), &
                                                     slowWork(7)) 
               CALL fteik_prefetchSweep1TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                        ttWork(1), ttWork(2), &
                                                        ttWork(3), ttWork(4), &
                                                        ttWork(5), ttWork(6), &
                                                        ttWorK(7), ttwork(8))
               ! apply the finite differencing
               tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &
                                            i, j, k,                       &
                                            1, 1, 1, & !sgntz, sgntx, sgnty,           &
                                            sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               ! update ttimes
               kndx = grid2indexF(i, j, k, nz, nzx)
               ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
            ENDDO
         ENDDO
      ENDDO
      !-----------------------------------Sweep 2----------------------------------------!
      sweep = 2
      ! Get some stencil information for this sweep
      CALL fteik_getSweepSigns64fF(sweep,               &
                                   sgntz, sgntx, sgnty, &
                                   sgnvz, sgnvx, sgnvy, &
                                   sgnrz, sgnrx, sgnry)
      CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
      CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
      ! Some derivative items
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      DO k = MAX(2,ysi),ny
         DO j = xsi+1,1,-1
            DO i = MAX(2,zsi),nz
               ! prefetch the slownesses and traveltimes
               CALL fteik_prefetchSweep2Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                     slow, & !slowWork)
                                                     slowWork(1), slowWork(2), &
                                                     slowWork(3), slowWork(4), &
                                                     slowWork(5), slowWork(6), &
                                                     slowWork(7))
               CALL fteik_prefetchSweep2TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                        ttWork(1), ttWork(2), &
                                                        ttWork(3), ttWork(4), &
                                                        ttWork(5), ttWork(6), &
                                                        ttWorK(7), ttwork(8))
               ! apply the finite differencing
               tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &   
                                            i, j, k,                       &
                                            1, -1, 1, & !sgntz, sgntx, sgnty,           &
                                            sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               ! update ttimes
               kndx = grid2indexF(i, j, k, nz, nzx)
               ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
            ENDDO
         ENDDO
      ENDDO
      !-----------------------------------Sweep 3----------------------------------------!
      sweep = 3
      ! Get some stencil information for this sweep
      CALL fteik_getSweepSigns64fF(sweep,               &
                                   sgntz, sgntx, sgnty, &
                                   sgnvz, sgnvx, sgnvy, &
                                   sgnrz, sgnrx, sgnry)
      CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
      CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
      ! Some derivative items
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      DO k = ysi+1,1,-1
         DO j = MAX(2,xsi),nx
            DO i = MAX(2,zsi),nz
               ! prefetch the slownesses and traveltimes
               CALL fteik_prefetchSweep3Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                     slow, & !slowWork)
                                                     slowWork(1), slowWork(2), &
                                                     slowWork(3), slowWork(4), &
                                                     slowWork(5), slowWork(6), &
                                                     slowWork(7))
               CALL fteik_prefetchSweep3TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                        ttWork(1), ttWork(2), &
                                                        ttWork(3), ttWork(4), &
                                                        ttWork(5), ttWork(6), &
                                                        ttWorK(7), ttwork(8))
               ! apply the finite differencing
               tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &
                                            i, j, k,                       &
                                            1, 1, -1, & !sgntz, sgntx, sgnty,           &
                                            sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               ! update ttimes
               kndx = grid2indexF(i, j, k, nz, nzx)
               ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
            ENDDO
         ENDDO
      ENDDO
      !-----------------------------------Sweep 4----------------------------------------!
      sweep = 4
      ! Get some stencil information for this sweep
      CALL fteik_getSweepSigns64fF(sweep,               &
                                   sgntz, sgntx, sgnty, &
                                   sgnvz, sgnvx, sgnvy, &
                                   sgnrz, sgnrx, sgnry)
      CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
      CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
      ! Some derivative items
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      DO k = ysi+1,1,-1
         DO j = xsi+1,1,-1
            DO i = MAX(2,zsi),nz
               ! prefetch the slownesses and traveltimes
               CALL fteik_prefetchSweep4Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                     slow, & !slowWork)
                                                     slowWork(1), slowWork(2), &
                                                     slowWork(3), slowWork(4), &
                                                     slowWork(5), slowWork(6), &
                                                     slowWork(7)) 
               CALL fteik_prefetchSweep4TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                        ttWork(1), ttWork(2), &
                                                        ttWork(3), ttWork(4), &
                                                        ttWork(5), ttWork(6), &
                                                        ttWorK(7), ttwork(8))
               ! apply the finite differencing
               tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &   
                                            i, j, k,                       &
                                            1, -1, -1, & !sgntz, sgntx, sgnty,           &
                                            sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               ! update ttimes
               kndx = grid2indexF(i, j, k, nz, nzx)
               ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
            ENDDO
         ENDDO
      ENDDO
      !-----------------------------------Sweep 5----------------------------------------!
      sweep = 5
      ! Get some stencil information for this sweep
      CALL fteik_getSweepSigns64fF(sweep,               &
                                   sgntz, sgntx, sgnty, &
                                   sgnvz, sgnvx, sgnvy, &
                                   sgnrz, sgnrx, sgnry)
      CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
      CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
      ! Some derivative items
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      DO k = MAX(2,ysi),ny
         DO j = MAX(2,xsi),nx
            DO i = zsi+1,1,-1
               ! prefetch the slownesses and traveltimes
               CALL fteik_prefetchSweep5Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                     slow, & !slowWork)
                                                     slowWork(1), slowWork(2), &
                                                     slowWork(3), slowWork(4), &
                                                     slowWork(5), slowWork(6), &
                                                     slowWork(7))
               CALL fteik_prefetchSweep5TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttwork)
                                                        ttWork(1), ttWork(2), &
                                                        ttWork(3), ttWork(4), &
                                                        ttWork(5), ttWork(6), &
                                                        ttWorK(7), ttwork(8))
               ! apply the finite differencing
               tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &   
                                            i, j, k,                       &
                                            -1, 1, 1, & !sgntz, sgntx, sgnty,           &
                                            sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               ! update ttimes
               kndx = grid2indexF(i, j, k, nz, nzx)
               ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
            ENDDO
         ENDDO
      ENDDO
      !-----------------------------------Sweep 6----------------------------------------!
      sweep = 6
      ! Get some stencil information for this sweep
      CALL fteik_getSweepSigns64fF(sweep,               &
                                   sgntz, sgntx, sgnty, &
                                   sgnvz, sgnvx, sgnvy, &
                                   sgnrz, sgnrx, sgnry)
      CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
      CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
      ! Some derivative items
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      DO k = MAX(2,ysi),ny
         DO j = xsi+1,1,-1
            DO i = zsi+1,1,-1
               ! prefetch the slownesses and traveltimes
               CALL fteik_prefetchSweep6Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                     slow, & !slowWork)
                                                     slowWork(1), slowWork(2), &
                                                     slowWork(3), slowWork(4), &
                                                     slowWork(5), slowWork(6), &
                                                     slowWork(7)) 
               CALL fteik_prefetchSweep6TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                        ttWork(1), ttWork(2), &
                                                        ttWork(3), ttWork(4), &
                                                        ttWork(5), ttWork(6), &
                                                        ttWorK(7), ttwork(8))
               ! apply the finite differencing
               tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &   
                                            i, j, k,                       &
                                            -1, -1, 1, & !sgntz, sgntx, sgnty,           &
                                            sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               ! update ttimes
               kndx = grid2indexF(i, j, k, nz, nzx)
               ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
            ENDDO
         ENDDO
      ENDDO
      !------------------------------------Sweep 7---------------------------------------!
      sweep = 7
      ! Get some stencil information for this sweep
      CALL fteik_getSweepSigns64fF(sweep,               &
                                   sgntz, sgntx, sgnty, &
                                   sgnvz, sgnvx, sgnvy, &
                                   sgnrz, sgnrx, sgnry)
      CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
      CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
      ! Some derivative items
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      DO k = ysi+1,1,-1
         DO j = MAX(2,xsi),nx
            DO i = zsi+1,1,-1
               ! prefetch the slownesses and traveltimes
               CALL fteik_prefetchSweep7Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                     slow, & !slowWork)
                                                     slowWork(1), slowWork(2), &
                                                     slowWork(3), slowWork(4), &
                                                     slowWork(5), slowWork(6), &
                                                     slowWork(7)) 
               CALL fteik_prefetchSweep7TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                        ttWork(1), ttWork(2), &
                                                        ttWork(3), ttWork(4), &
                                                        ttWork(5), ttWork(6), &
                                                        ttWorK(7), ttwork(8))
               ! apply the finite differencing
               tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &
                                            i, j, k,                       &
                                            -1, 1, -1, & !sgntz, sgntx, sgnty,           &
                                            sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               ! update ttimes
               kndx = grid2indexF(i, j, k, nz, nzx)
               ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
            ENDDO
         ENDDO
      ENDDO
      !-----------------------------------Sweep 8----------------------------------------!
      sweep = 8
      ! Get some stencil information for this sweep
      CALL fteik_getSweepSigns64fF(sweep,               &
                                   sgntz, sgntx, sgnty, &
                                   sgnvz, sgnvx, sgnvy, &
                                   sgnrz, sgnrx, sgnry)
      CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
      CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
      ! Some derivative items
      sgnrz_dzi = sgnrz*dzi
      sgnrx_dxi = sgnrx*dxi
      sgnry_dyi = sgnry*dyi
      DO k = ysi+1,1,-1
         DO j = xsi+1,1,-1
            DO i = zsi+1,1,-1
               ! prefetch the slownesses and traveltimes
               CALL fteik_prefetchSweep8Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                     slow, & !slowWork)
                                                     slowWork(1), slowWork(2), &
                                                     slowWork(3), slowWork(4), &
                                                     slowWork(5), slowWork(6), &
                                                     slowWork(7))

               CALL fteik_prefetchSweep8TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                        ttWork(1), ttWork(2), &
                                                        ttWork(3), ttWork(4), &
                                                        ttWork(5), ttWork(6), &
                                                        ttWorK(7), ttwork(8))
               ! apply the finite differencing
               tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &   
                                            i, j, k,                       &
                                            -1, -1, -1, & !sgntz, sgntx, sgnty,           &
                                            sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
               ! update ttimes
               kndx = grid2indexF(i, j, k, nz, nzx)
               ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
            ENDDO
         ENDDO
      ENDDO
      !----------------------------------------------------------------------------------!
      !                              Begin the iterative method                          !
      !----------------------------------------------------------------------------------!
      linitk = .FALSE.
      CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
      CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
      DO kiter=1,nsweep
         !-----------------------------------Sweep 1-------------------------------------!
         sweep = 1
         ! Get some stencil information for this sweep
         CALL fteik_getSweepSigns64fF(sweep,               &
                                      sgntz, sgntx, sgnty, &
                                      sgnvz, sgnvx, sgnvy, &
                                      sgnrz, sgnrx, sgnry)
         ! Some derivative items
         sgnrz_dzi = sgnrz*dzi
         sgnrx_dxi = sgnrx*dxi
         sgnry_dyi = sgnry*dyi
         DO k = 2,ny
            DO j = 2,nx
               DO i = 2,nz
                  ! prefetch the slownesses and traveltimes
                  CALL fteik_prefetchSweep1Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                        slow, & !slowWork)
                                                        slowWork(1), slowWork(2), &
                                                        slowWork(3), slowWork(4), &
                                                        slowWork(5), slowWork(6), &
                                                        slowWork(7)) 
                  CALL fteik_prefetchSweep1TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                           ttWork(1), ttWork(2), &
                                                           ttWork(3), ttWork(4), &
                                                           ttWork(5), ttWork(6), &
                                                           ttWorK(7), ttwork(8))
                  ! apply the finite differencing
                  tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &
                                               i, j, k,                       &
                                               1, 1, 1, & !sgntz, sgntx, sgnty,           &
                                               sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
#ifdef DEBUG
                  CALL fteik_prefetchSlowness64fF(i, j, k, nzm1, nzm1_nxm1,    &
                                                  sgnvz, sgnvx, sgnvy,         &
                                                  slowPerm, slow, slowWork1)
                  IF (MAXVAL(ABS(slowWork(1:7) - slowWork1(1:7))) > 0.d0) print *, 'error1'
                  CALL fteik_prefetchTravelTimes64fF(i, j, k, nz, nzx,     &
                                                     sgntz, sgntx, sgnty,  &    
                                                     ttPerm, ttimes, ttWork1)   
                  IF (MAXVAL(ABS(ttWork(1:8) - ttWork1(1:8))) > 0.d0) print *, 'error1t'
                  CALl fteik_localSolver_noInit64fF(1, &
                                                    ttWork1(1), ttWork1(2), &
                                                    ttWork1(3), ttWork1(4),  &
                                                    ttWork1(5), ttWork1(6),  &
                                                    ttwork1(7), &
                                                    slowWork1(1), slowWork1(2), & 
                                                    slowWork1(3), slowWork1(4), &
                                                    slowWork1(5), slowWork1(6), &
                                                    slowWork1(7),  ttWork1(8)) !tupd1)
                  tupd1 = ttWork1(8)
                  IF (ABS(tupd - tupd1) > 1.d-15) THEN
                     WRITE(*,*) 'failed1:', i, j, k, tupd, tupd1
                  ENDIF
#endif
                  ! update ttimes
                  kndx = grid2indexF(i, j, k, nz, nzx)
                  ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
               ENDDO
            ENDDO
         ENDDO
         !-----------------------------------Sweep 2-------------------------------------!
         sweep = 2
         ! Get some stencil information for this sweep
         CALL fteik_getSweepSigns64fF(sweep,               &
                                      sgntz, sgntx, sgnty, &
                                      sgnvz, sgnvx, sgnvy, &
                                      sgnrz, sgnrx, sgnry)
         CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
         CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
         ! Some derivative items
         sgnrz_dzi = sgnrz*dzi
         sgnrx_dxi = sgnrx*dxi
         sgnry_dyi = sgnry*dyi
         DO k = 2,ny
            DO j = nx-1,1,-1
               DO i = 2,nz
                  ! prefetch the slownesses and traveltimes
                  CALL fteik_prefetchSweep2Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                        slow, & !slowWork)
                                                        slowWork(1), slowWork(2), &
                                                        slowWork(3), slowWork(4), &
                                                        slowWork(5), slowWork(6), &
                                                        slowWork(7)) 
                  CALL fteik_prefetchSweep2TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                           ttWork(1), ttWork(2), &
                                                           ttWork(3), ttWork(4), &
                                                           ttWork(5), ttWork(6), &
                                                           ttWorK(7), ttwork(8))
                  ! apply the finite differencing
                  tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &   
                                               i, j, k,                       &
                                               1, -1, 1, & !sgntz, sgntx, sgnty,           &
                                               sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
#ifdef DEBUG
                  CALL fteik_prefetchSlowness64fF(i, j, k, nzm1, nzm1_nxm1,    &
                                                  sgnvz, sgnvx, sgnvy,         &
                                                  slowPerm, slow, slowWork1)
                  IF (MAXVAL(ABS(slowWork(1:7) - slowWork1(1:7))) > 0.d0) print *, 'error2'
                  CALL fteik_prefetchTravelTimes64fF(i, j, k, nz, nzx,     &
                                                     sgntz, sgntx, sgnty,  &    
                                                     ttPerm, ttimes, ttWork1)   
                  IF (MAXVAL(ABS(ttWork(1:8) - ttWork1(1:8))) > 0.d0) print *, 'error2t'
                  CALl fteik_localSolver_noInit64fF(1, &
                                                    ttWork1(1), ttWork1(2), &
                                                    ttWork1(3), ttWork1(4),  &
                                                    ttWork1(5), ttWork1(6),  &
                                                    ttwork1(7), &
                                                    slowWork1(1), slowWork1(2), & 
                                                    slowWork1(3), slowWork1(4), &
                                                    slowWork1(5), slowWork1(6), &
                                                    slowWork1(7),  ttWork1(8)) !tupd1)
                  tupd1 = ttWork1(8)
                  IF (ABS(tupd - tupd1) > 1.d-15) THEN
                     WRITE(*,*) 'failed8:', i, j, k, tupd, tupd1
                  ENDIF
#endif
                  ! update ttimes
                  kndx = grid2indexF(i, j, k, nz, nzx)
                  ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
               ENDDO
            ENDDO
         ENDDO
         !-----------------------------------Sweep 3-------------------------------------!
         sweep = 3
         ! Get some stencil information for this sweep
         CALL fteik_getSweepSigns64fF(sweep,               &
                                      sgntz, sgntx, sgnty, &
                                      sgnvz, sgnvx, sgnvy, &
                                      sgnrz, sgnrx, sgnry)
         CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
         CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
         ! Some derivative items
         sgnrz_dzi = sgnrz*dzi
         sgnrx_dxi = sgnrx*dxi
         sgnry_dyi = sgnry*dyi
         DO k = ny-1,1,-1
            DO j = 2,nx
               DO i = 2,nz
                  ! prefetch the slownesses and traveltimes
                  CALL fteik_prefetchSweep3Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                        slow, & !slowWork)
                                                        slowWork(1), slowWork(2), &
                                                        slowWork(3), slowWork(4), &
                                                        slowWork(5), slowWork(6), &
                                                        slowWork(7)) 
                  CALL fteik_prefetchSweep3TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                           ttWork(1), ttWork(2), &
                                                           ttWork(3), ttWork(4), &
                                                           ttWork(5), ttWork(6), &
                                                           ttWorK(7), ttwork(8))
                  ! apply the finite differencing
                  tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &
                                               i, j, k,                       &
                                               1, 1, -1, & !sgntz, sgntx, sgnty,           &
                                               sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
#ifdef DEBUG
                  CALL fteik_prefetchSlowness64fF(i, j, k, nzm1, nzm1_nxm1,    &
                                                  sgnvz, sgnvx, sgnvy,         &
                                                  slowPerm, slow, slowWork1)
                  IF (MAXVAL(ABS(slowWork(1:7) - slowWork1(1:7))) > 0.d0) print *, 'error3'
                  CALL fteik_prefetchTravelTimes64fF(i, j, k, nz, nzx,     &
                                                     sgntz, sgntx, sgnty,  &    
                                                     ttPerm, ttimes, ttWork1)   
                  IF (MAXVAL(ABS(ttWork(1:8) - ttWork1(1:8))) > 0.d0) print *, 'error3t'
                  CALl fteik_localSolver_noInit64fF(1, &
                                                    ttWork1(1), ttWork1(2), &
                                                    ttWork1(3), ttWork1(4),  &
                                                    ttWork1(5), ttWork1(6),  &
                                                    ttwork1(7), &
                                                    slowWork1(1), slowWork1(2), & 
                                                    slowWork1(3), slowWork1(4), &
                                                    slowWork1(5), slowWork1(6), &
                                                    slowWork1(7),  ttWork1(8)) !tupd1)
                  tupd1 = ttWork1(8)
                  IF (ABS(tupd - tupd1) > 1.d-15) THEN
                     WRITE(*,*) 'failed3:', i, j, k, tupd, tupd1
                  ENDIF
#endif
                  ! update ttimes
                  kndx = grid2indexF(i, j, k, nz, nzx)
                  ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
               ENDDO
            ENDDO
         ENDDO
         !-----------------------------------Sweep 4-------------------------------------!
         sweep = 4
         ! Get some stencil information for this sweep
         CALL fteik_getSweepSigns64fF(sweep,               &
                                      sgntz, sgntx, sgnty, &
                                      sgnvz, sgnvx, sgnvy, &
                                      sgnrz, sgnrx, sgnry)
         CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
         CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
         ! Some derivative items
         sgnrz_dzi = sgnrz*dzi
         sgnrx_dxi = sgnrx*dxi
         sgnry_dyi = sgnry*dyi
         DO k = ny-1,1,-1
            DO j = nx-1,1,-1
               DO i = 2,nz
                  ! prefetch the slownesses and traveltimes
                  CALL fteik_prefetchSweep4Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                        slow, & !slowWork)
                                                        slowWork(1), slowWork(2), &
                                                        slowWork(3), slowWork(4), &
                                                        slowWork(5), slowWork(6), &
                                                        slowWork(7)) 
                  CALL fteik_prefetchSweep4TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                           ttWork(1), ttWork(2), &
                                                           ttWork(3), ttWork(4), &
                                                           ttWork(5), ttWork(6), &
                                                           ttWorK(7), ttwork(8))
                  ! apply the finite differencing
                  tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &   
                                               i, j, k,                       &
                                               1, -1, -1, & !sgntz, sgntx, sgnty,           &
                                               sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
#ifdef DEBUG
                  CALL fteik_prefetchSlowness64fF(i, j, k, nzm1, nzm1_nxm1,    &
                                                  sgnvz, sgnvx, sgnvy,         &
                                                  slowPerm, slow, slowWork1)
                  IF (MAXVAL(ABS(slowWork(1:7) - slowWork1(1:7))) > 0.d0) print *, 'error4'
                  CALL fteik_prefetchTravelTimes64fF(i, j, k, nz, nzx,     &
                                                     sgntz, sgntx, sgnty,  &    
                                                     ttPerm, ttimes, ttWork1)   
                  IF (MAXVAL(ABS(ttWork(1:8) - ttWork1(1:8))) > 0.d0) print *, 'error4t'
                  CALl fteik_localSolver_noInit64fF(1, &
                                                    ttWork1(1), ttWork1(2), &
                                                    ttWork1(3), ttWork1(4),  &
                                                    ttWork1(5), ttWork1(6),  &
                                                    ttwork1(7), &
                                                    slowWork1(1), slowWork1(2), & 
                                                    slowWork1(3), slowWork1(4), &
                                                    slowWork1(5), slowWork1(6), &
                                                    slowWork1(7),  ttWork1(8)) !tupd1)
                  tupd1 = ttWork1(8)
                  IF (ABS(tupd - tupd1) > 1.d-15) THEN
                     WRITE(*,*) 'failed4:', i, j, k, tupd, tupd1
                  ENDIF
#endif
                  ! update ttimes
                  kndx = grid2indexF(i, j, k, nz, nzx)
                  ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
               ENDDO
            ENDDO
         ENDDO
         !-----------------------------------Sweep 5-------------------------------------!
         sweep = 5
         ! Get some stencil information for this sweep
         CALL fteik_getSweepSigns64fF(sweep,               &
                                      sgntz, sgntx, sgnty, &
                                      sgnvz, sgnvx, sgnvy, &
                                      sgnrz, sgnrx, sgnry)
         CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
         CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
         ! Some derivative items
         sgnrz_dzi = sgnrz*dzi
         sgnrx_dxi = sgnrx*dxi
         sgnry_dyi = sgnry*dyi
         DO k = 2,ny
            DO j = 2,nx
               DO i = nz-1,1,-1
                  ! prefetch the slownesses and traveltimes
                  CALL fteik_prefetchSweep5Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                        slow, & !slowWork)
                                                        slowWork(1), slowWork(2), &
                                                        slowWork(3), slowWork(4), &
                                                        slowWork(5), slowWork(6), &
                                                        slowWork(7)) 
                  CALL fteik_prefetchSweep5TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                           ttWork(1), ttWork(2), &
                                                           ttWork(3), ttWork(4), &
                                                           ttWork(5), ttWork(6), &
                                                           ttWorK(7), ttwork(8))
                  ! apply the finite differencing
                  tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &   
                                               i, j, k,                       &
                                               -1, 1, 1, & !sgntz, sgntx, sgnty,           &
                                               sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
#ifdef DEBUG
                  CALL fteik_prefetchSlowness64fF(i, j, k, nzm1, nzm1_nxm1,    &
                                                  sgnvz, sgnvx, sgnvy,         &
                                                  slowPerm, slow, slowWork1)
                  IF (MAXVAL(ABS(slowWork(1:7) - slowWork1(1:7))) > 0.d0) print *, 'error5'
                  CALL fteik_prefetchTravelTimes64fF(i, j, k, nz, nzx,     &
                                                     sgntz, sgntx, sgnty,  &    
                                                     ttPerm, ttimes, ttWork1)   
                  IF (MAXVAL(ABS(ttWork(1:8) - ttWork1(1:8))) > 0.d0) print *, 'error5t'
                  CALl fteik_localSolver_noInit64fF(1, &
                                                    ttWork1(1), ttWork1(2), &
                                                    ttWork1(3), ttWork1(4),  &
                                                    ttWork1(5), ttWork1(6),  &
                                                    ttwork1(7), &
                                                    slowWork1(1), slowWork1(2), & 
                                                    slowWork1(3), slowWork1(4), &
                                                    slowWork1(5), slowWork1(6), &
                                                    slowWork1(7),  ttWork1(8)) !tupd1)
                  tupd1 = ttWork1(8)
                  IF (ABS(tupd - tupd1) > 1.d-15) THEN
                     WRITE(*,*) 'failed5:', i, j, k, tupd, tupd1
                  ENDIF
#endif
                  ! update ttimes
                  kndx = grid2indexF(i, j, k, nz, nzx)
                  ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
               ENDDO
            ENDDO
         ENDDO
         !-----------------------------------Sweep 6-------------------------------------!
         sweep = 6
         ! Get some stencil information for this sweep
         CALL fteik_getSweepSigns64fF(sweep,               &
                                      sgntz, sgntx, sgnty, &
                                      sgnvz, sgnvx, sgnvy, &
                                      sgnrz, sgnrx, sgnry)
         CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
         CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
         ! Some derivative items
         sgnrz_dzi = sgnrz*dzi
         sgnrx_dxi = sgnrx*dxi
         sgnry_dyi = sgnry*dyi
         DO k = 2,ny
            DO j = nx-1,1,-1
               DO i = nz-1,1,-1
                  ! prefetch the slownesses and traveltimes
                  CALL fteik_prefetchSweep6Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                        slow, & !slowWork)
                                                        slowWork(1), slowWork(2), &
                                                        slowWork(3), slowWork(4), &
                                                        slowWork(5), slowWork(6), &
                                                        slowWork(7)) 
                  CALL fteik_prefetchSweep6TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                           ttWork(1), ttWork(2), &
                                                           ttWork(3), ttWork(4), &
                                                           ttWork(5), ttWork(6), &
                                                           ttWorK(7), ttwork(8))
                  ! apply the finite differencing
                  tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &   
                                               i, j, k,                       &
                                               -1, -1, 1, & !sgntz, sgntx, sgnty,           &
                                               sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
#ifdef DEBUG
                  CALL fteik_prefetchSlowness64fF(i, j, k, nzm1, nzm1_nxm1,    &
                                                  sgnvz, sgnvx, sgnvy,         &
                                                  slowPerm, slow, slowWork1)
                  IF (MAXVAL(ABS(slowWork(1:7) - slowWork1(1:7))) > 0.d0) print *, 'error6'
                  CALL fteik_prefetchTravelTimes64fF(i, j, k, nz, nzx,     &
                                                     sgntz, sgntx, sgnty,  &    
                                                     ttPerm, ttimes, ttWork1)   
                  IF (MAXVAL(ABS(ttWork(1:8) - ttWork1(1:8))) > 0.d0) print *, 'error6t'
                  CALl fteik_localSolver_noInit64fF(1, &
                                                    ttWork1(1), ttWork1(2), &
                                                    ttWork1(3), ttWork1(4),  &
                                                    ttWork1(5), ttWork1(6),  &
                                                    ttwork1(7), &
                                                    slowWork1(1), slowWork1(2), & 
                                                    slowWork1(3), slowWork1(4), &
                                                    slowWork1(5), slowWork1(6), &
                                                    slowWork1(7),  ttWork1(8)) !tupd1)
                  tupd1 = ttWork1(8)
                  IF (ABS(tupd - tupd1) > 1.d-15) THEN
                     WRITE(*,*) 'failed6:', i, j, k, tupd, tupd1
                  ENDIF
#endif
                  ! update ttimes
                  kndx = grid2indexF(i, j, k, nz, nzx)
                  ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
               ENDDO
            ENDDO
         ENDDO
         !------------------------------------Sweep 7------------------------------------!
         sweep = 7
         ! Get some stencil information for this sweep
         CALL fteik_getSweepSigns64fF(sweep,               &
                                      sgntz, sgntx, sgnty, &
                                      sgnvz, sgnvx, sgnvy, &
                                      sgnrz, sgnrx, sgnry)
         CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
         CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
         ! Some derivative items
         sgnrz_dzi = sgnrz*dzi
         sgnrx_dxi = sgnrx*dxi
         sgnry_dyi = sgnry*dyi
         DO k = ny-1,1,-1
            DO j = 2,nx
               DO i = nz-1,1,-1
                  ! prefetch the slownesses and traveltimes
                  CALL fteik_prefetchSweep7Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                        slow, & !slowWork)
                                                        slowWork(1), slowWork(2), &
                                                        slowWork(3), slowWork(4), &
                                                        slowWork(5), slowWork(6), &
                                                        slowWork(7)) 
                  CALL fteik_prefetchSweep7TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                           ttWork(1), ttWork(2), &
                                                           ttWork(3), ttWork(4), &
                                                           ttWork(5), ttWork(6), &
                                                           ttWorK(7), ttwork(8))
                  ! apply the finite differencing
                  tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &
                                               i, j, k,                       &
                                               -1, 1, -1, & !sgntz, sgntx, sgnty,           &
                                               sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
#ifdef DEBUG
                  CALL fteik_prefetchSlowness64fF(i, j, k, nzm1, nzm1_nxm1,    &
                                                  sgnvz, sgnvx, sgnvy,         &
                                                  slowPerm, slow, slowWork1)
                  IF (MAXVAL(ABS(slowWork(1:7) - slowWork1(1:7))) > 0.d0) print *, 'error7'
                  CALL fteik_prefetchTravelTimes64fF(i, j, k, nz, nzx,     &
                                                     sgntz, sgntx, sgnty,  &    
                                                     ttPerm, ttimes, ttWork1)   
                  IF (MAXVAL(ABS(ttWork(1:8) - ttWork1(1:8))) > 0.d0) print *, 'error7t'
                  CALl fteik_localSolver_noInit64fF(1, &
                                                    ttWork1(1), ttWork1(2), &
                                                    ttWork1(3), ttWork1(4),  &
                                                    ttWork1(5), ttWork1(6),  &
                                                    ttwork1(7), &
                                                    slowWork1(1), slowWork1(2), & 
                                                    slowWork1(3), slowWork1(4), &
                                                    slowWork1(5), slowWork1(6), &
                                                    slowWork1(7),  ttWork1(8)) !tupd1)
                  tupd1 = ttWork1(8)
                  IF (ABS(tupd - tupd1) > 1.d-15) THEN
                     WRITE(*,*) 'failed7:', i, j, k, tupd, tupd1
                  ENDIF
#endif
                  ! update ttimes
                  kndx = grid2indexF(i, j, k, nz, nzx)
                  ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
               ENDDO
            ENDDO
         ENDDO
         !-----------------------------------Sweep 8-------------------------------------!
         sweep = 8
         ! Get some stencil information for this sweep
         CALL fteik_getSweepSigns64fF(sweep,               &
                                      sgntz, sgntx, sgnty, &
                                      sgnvz, sgnvx, sgnvy, &
                                      sgnrz, sgnrx, sgnry)
         CALL fteik_getSlownessPermF(sweep, slowPerm, ierr)
         CALL fteik_getTravelTimePermF(sweep, ttPerm, ierr)
         ! Some derivative items
         sgnrz_dzi = sgnrz*dzi
         sgnrx_dxi = sgnrx*dxi
         sgnry_dyi = sgnry*dyi
         DO k = ny-1,1,-1
            DO j = nx-1,1,-1 
               DO i = nz-1,1,-1
                  ! prefetch the slownesses and traveltimes
                  CALL fteik_prefetchSweep8Slowness64fF(i, j, k, nz, nx, ny, nzm1, nzm1_nxm1, &
                                                        slow, & !slowWork)
                                                        slowWork(1), slowWork(2), &
                                                        slowWork(3), slowWork(4), &
                                                        slowWork(5), slowWork(6), &
                                                        slowWork(7)) 
                  CALL fteik_prefetchSweep8TravelTimes64fF(i, j, k, nz, nx, ny, nzx, ttimes, &!ttWork)
                                                           ttWork(1), ttWork(2), &
                                                           ttWork(3), ttWork(4), &
                                                           ttWork(5), ttWork(6), &
                                                           ttWorK(7), ttwork(8))
                  ! apply the finite differencing
                  tupd = fteik_localSolver64fF(ttWork, slowWork, linitk,      &   
                                               i, j, k,                       &
                                               -1, -1, -1, & !sgntz, sgntx, sgnty,           &
                                               sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
#ifdef DEBUG
                  CALL fteik_prefetchSlowness64fF(i, j, k, nzm1, nzm1_nxm1,    &
                                                  sgnvz, sgnvx, sgnvy,         &    
                                                  slowPerm, slow, slowWork1)
                  IF (MAXVAL(ABS(slowWork(1:7) - slowWork1(1:7))) > 0.d0) print *, 'error8'
                  CALL fteik_prefetchTravelTimes64fF(i, j, k, nz, nzx,     &
                                                     sgntz, sgntx, sgnty,  &    
                                                     ttPerm, ttimes, ttWork1)
                  IF (MAXVAL(ABS(ttWork(1:8) - ttWork1(1:8))) > 0.d0) print *, 'error8t'
                  CALl fteik_localSolver_noInit64fF(1, &
                                                    ttWork1(1), ttWork1(2), &
                                                    ttWork1(3), ttWork1(4),  &
                                                    ttWork1(5), ttWork1(6),  &
                                                    ttwork1(7), &
                                                    slowWork1(1), slowWork1(2), &
                                                    slowWork1(3), slowWork1(4), &
                                                    slowWork1(5), slowWork1(6), &
                                                    slowWork1(7),  ttWork1(8)) !tupd1)
                  tupd1 = ttWork1(8)
                  IF (ABS(tupd - tupd1) > 1.d-15) THEN
                     WRITE(*,*) 'failed8:', i, j, k, tupd, tupd1
                  ENDIF
#endif
                  ! update ttimes
                  kndx = grid2indexF(i, j, k, nz, nzx)
                  ttimes(kndx) = tupd !DMIN1(ttimes(kndx), tupd)
               ENDDO
            ENDDO
         ENDDO
      ENDDO ! end iterative loop on sweeps 
call cpu_time(t1)
print *, 'FSM time:', t1 - t0
print *, minval(ttimes), maxval(ttimes)
      lhaveTravelTimes = .TRUE.
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE fteik_solveEikonalLSMF(ierr) &
                 BIND(C, NAME='fteik_solveEikonalLSMF')
      USE ISO_C_BINDING
      USE FTEIK_CONSTANTS64F, ONLY : TRUE, FALSE
      USE FTEIK_UTILS64F, ONLY : ttimes, &
                                 nsweep, &
                                 lhaveGrid, lhaveTravelTimes 
      USE FTEIK_UTILS64F, ONLY : FTEIK_HUGE
      USE FTEIK_UTILS64F, ONLY : fteik_meshConstants64fF, fteik_setSzeroF, &
                                 fteik_setEpsSolverF, &
                                 fteik_initializeTravelTimesNearSourceF, &
                                 fteik_evaluateLSSweep164fF, &
                                 fteik_evaluateLSSweep264fF, &
                                 fteik_evaluateLSSweep64fF
      USE FTEIK_AUTOCODE, ONLY : fteik_evaluateSweep1LS64fF, &
                                 fteik_evaluateSweep2LS64fF, &
                                 fteik_evaluateSweep3LS64fF, &
                                 fteik_evaluateSweep4LS64fF, &
                                 fteik_evaluateSweep5LS64fF, &
                                 fteik_evaluateSweep6LS64fF, &
                                 fteik_evaluateSweep7LS64fF, &
                                 fteik_evaluateSweep8LS64fF
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) kiter, lsweep
      LOGICAL(C_BOOL) linitk
real*8 t0, t1

      ierr = 0 
      lhaveTravelTimes = .FALSE.
      ! some easy checks
      IF (.NOT. lhaveGrid) THEN
         WRITE(*,*) 'fteik_solveEikonalFSM64fF: Grid not yet set'
         ierr = 1 
         RETURN
      ENDIF
      ! don't let mesh constants go out of scope
      CALL fteik_meshConstants64fF()
      ! initialize travel times to a big number
      ttimes(:) = FTEIK_HUGE
      ! set the slowness around the source for the local solver
      CALL fteik_setSzeroF(ierr) 
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solveEikonalFSM64fF: Error initializing szero'
         RETURN
      ENDIF
      ! set the spherical to cartesian transition in time 
      CALL fteik_setEpsSolverF(ierr) 
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solveEikonalFSM64fF: Error setting epsSolver'
         RETURN
      ENDIF
      ! Initialize the points around the source
      ttimes(:) = FTEIK_HUGE
      CALL fteik_initializeTravelTimesNearSourceF(ttimes, ierr)
      IF (ierr /= 0) THEN 
         WRITE(*,*) 'fteik_solveEikonalFSM64fF: Error initializing travel-time'
         RETURN
      ENDIF
      linitk = .FALSE.

call cpu_time(t0)
      DO 101 lsweep=1,8
         If (lsweep == 1) THEN
            CALL fteik_evaluateSweep1LS64fF(TRUE, ttimes, ierr)
print*, minval(ttimes), maxval(ttimes)
         ELSEIF (lsweep == 2) THEN
            CALL fteik_evaluateSweep2LS64fF(TRUE, ttimes, ierr)
         ELSEIF (lsweep == 3) THEN
            CALL fteik_evaluateSweep3LS64fF(TRUE, ttimes, ierr)
         ELSEIF (lsweep == 4) THEN
            CALL fteik_evaluateSweep4LS64fF(TRUE, ttimes, ierr)
         ELSEIF (lsweep == 5) THEN
            CALL fteik_evaluateSweep5LS64fF(TRUE, ttimes, ierr)
         ELSEIF (lsweep == 6) THEN
            CALL fteik_evaluateSweep6LS64fF(TRUE, ttimes, ierr)
         ELSEIF (lsweep == 7) THEN
            CALL fteik_evaluateSweep7LS64fF(TRUE, ttimes, ierr)
         ELSE
            CALL fteik_evaluateSweep8LS64fF(TRUE, ttimes, ierr)
         ENDIF
  101 CONTINUE
      DO kiter=1,nsweep
         DO 102 lsweep=1,8
            If (lsweep == 1) THEN
               CALL fteik_evaluateSweep1LS64fF(FALSE, ttimes, ierr)
            ELSEIF (lsweep == 2) THEN
               CALL fteik_evaluateSweep2LS64fF(FALSE, ttimes, ierr)
            ELSEIF (lsweep == 3) THEN
               CALL fteik_evaluateSweep3LS64fF(FALSE, ttimes, ierr)
            ELSEIF (lsweep == 4) THEN
               CALL fteik_evaluateSweep4LS64fF(FALSE, ttimes, ierr)
            ELSEIF (lsweep == 5) THEN
               CALL fteik_evaluateSweep5LS64fF(FALSE, ttimes, ierr)
            ELSEIF (lsweep == 6) THEN
               CALL fteik_evaluateSweep6LS64fF(FALSE, ttimes, ierr)
            ELSEIF (lsweep == 7) THEN
               CALL fteik_evaluateSweep7LS64fF(FALSE, ttimes, ierr)
            ELSE
               CALL fteik_evaluateSweep8LS64fF(FALSE, ttimes, ierr)
            ENDIF
  102    CONTINUE
      ENDDO
call cpu_time(t1)
print *, 'execute:', t1 - t0
print *, minval(ttimes), maxval(ttimes)
print *, 'wait; did id o it?'
      lhaveTravelTimes = .TRUE.
      RETURN
      END

