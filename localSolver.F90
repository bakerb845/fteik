MODULE FTEIK_LOCALSOLVER64F
  USE ISO_C_BINDING
  USE FTEIK_CONSTANTS64F, ONLY : chunkSize
  INTEGER(C_INT), PRIVATE, SAVE :: zsi, xsi, ysi
  REAL(C_DOUBLE), PRIVATE, SAVE :: zsa, xsa, ysa  
  REAL(C_DOUBLE), PRIVATE, SAVE :: dz, dx, dy
  REAL(C_DOUBLE), PRIVATE, SAVE :: szero, szero2
  REAL(C_DOUBLE), PRIVATE, SAVE :: dxi, dyi, dzi
  REAL(C_DOUBLE), PRIVATE, SAVE :: dx2, dy2, dz2
  REAL(C_DOUBLE), PRIVATE, SAVE :: dx2i, dy2i, dz2i 
  REAL(C_DOUBLE), PRIVATE, SAVE :: dsum, dsum_nine, dsumi, epsSolver
  REAL(C_DOUBLE), PRIVATE, SAVE :: dz2dx2, dz2dy2, dx2dy2
  REAL(C_DOUBLE), PRIVATE, SAVE :: dz2i_dx2i, dz2i_dy2i, dx2i_dy2i
  REAL(C_DOUBLE), PRIVATE, SAVE :: dz2i_p_dx2i, dz2i_p_dy2i, dx2i_p_dy2i
  REAL(C_DOUBLE), PRIVATE, SAVE :: dz2i_p_dy2i_inv, dz2i_p_dx2i_inv, dx2i_p_dy2i_inv
  !REAL(C_DOUBLE), PRIVATE, SAVE :: t12min(chunkSize)
  !!DIR$ ATTRIBUTES ALIGN:64 :: t12min
  CONTAINS
!----------------------------------------------------------------------------------------!
!                                     Begin the Code                                     !
!----------------------------------------------------------------------------------------!
      PURE SUBROUTINE fteik_localSolver_noInit64fF(n, tv, te, tn, tev,         &
                                                   ten, tnv, tnve,             &
                                                   slow1, slow2, slow3, slow4, &
                                                   slow5, slow6, slow7,        &
                                                   tupd)                       &
      BIND(C, NAME='fteik_localSolver_noInit64fF')
      USE ISO_C_BINDING
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, chunkSize, half, four
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: n
      REAL(C_DOUBLE), INTENT(IN) :: tv(n), te(n), tn(n), tev(n),  &
                                    ten(n), tnv(n), tnve(n)!, tt0(n)
      REAL(C_DOUBLE), INTENT(IN) :: slow1(n), slow2(n), slow3(n), slow4(n), &
                                    slow5(n), slow6(n), slow7(n)
      REAL(C_DOUBLE), INTENT(INOUT) :: tupd(n)
      ! local variables
      REAL(C_DOUBLE) four_sref2, sref2,                                      &
                     t1, t2, t3, t3d,                      &
                     t1_2d, t2_2d, t3_2d,                                    &
                     ta, tb, tab, tab2, tac, tac2, tbc, tbc2, tc,            &
                     temtn, temtv, tvmtn, tmax
      INTEGER(C_INT) i
      REAL(C_DOUBLE) :: t12min(chunkSize)
      !DIR$ ATTRIBUTES ALIGN:64 :: t12min
      !DIR$ ASSUME_ALIGNED tv: 64
      !DIR$ ASSUME_ALIGNED te: 64
      !DIR$ ASSUME_ALIGNED tn: 64
      !DIR$ ASSUME_ALIGNED tev: 64
      !DIR$ ASSUME_ALIGNED ten: 64
      !DIR$ ASSUME_ALIGNED tnv: 64
      !DIR$ ASSUME_ALIGNED tnve: 64
      !DIR$ ASSUME_ALIGNED slow1: 64
      !DIR$ ASSUME_ALIGNED slow2: 64
      !DIR$ ASSUME_ALIGNED slow3: 64
      !DIR$ ASSUME_ALIGNED slow4: 64
      !DIR$ ASSUME_ALIGNED slow5: 64
      !DIR$ ASSUME_ALIGNED slow6: 64
      !DIR$ ASSUME_ALIGNED slow7: 64
      !DIR$ ASSUME_ALIGNED tupd: 64
      !!DIR$ ASSUME_ALIGNED tt0: 64

      !------------------1D operators, (refracted times),set times to BIG----------------!
      !DO i=1,n
      !   tupd(i) = tt0(i)
      !ENDDO
      ! V plane
      DO i=1,n
         !tupd(i) = MIN(tt0(i),  tv(i) + dz*slow1(i)) !t1 = tv + dz*slowLoc(1)
         !tupd(i) = MIN(tupd(i), tv(i) + dz*slow1(i))
         t12min(i) = tv(i) + dz*slow1(i)
      ENDDO
      ! WE plane
      DO i=1,n
         !tupd(i) = MIN(tupd(i), te(i) + dx*slow2(i)) !t2 = te + dx*slowLoc(2)
         t12min(i) = MIN(t12min(i), te(i) + dx*slow2(i))
      ENDDO
      ! NS plane
      DO i=1,n
         !tupd(i) = MIN(tupd(i), tn(i) + dy*slow3(i)) !t3 = tn + dy*slowLoc(3)
         t12min(i) = MIN(t12min(i), tn(i) + dy*slow3(i))
      ENDDO
      !--------------------------------------2D operators--------------------------------!
      ! ZX (VE) plane
      DO i=1,n
         t1_2d = FTEIK_HUGE
         !tvmte = tv(i) - te(i)
         temtv = te(i) - tv(i)
         !IF ( (tv(i) < te(i) + dx*slow4(i)) .AND. (te(i) < tv(i) + dz*slow4(i)) ) THEN
         IF (temtv < dz*slow4(i) .AND.-temtv < dx*slow4(i)) THEN
            sref2 = slow4(i)*slow4(i)
            !ta = tev(i) + te(i) - tv(i)
            !tb = tev(i) - te(i) + tv(i)
            ta = tev(i) + temtv !ta = tev(i) + te(i) - tv(i)
            tb = tev(i) - temtv !tb = tev(i) - te(i) + tv(i)
            tab = ta - tb
            tab2 = tab*tab
            four_sref2 = four*sref2
            t1_2d = ( (tb*dz2i + ta*dx2i) &
                     + SQRT(four_sref2*dz2i_p_dx2i - dz2i_dx2i*tab2) )*dz2i_p_dx2i_inv
         ENDIF
         !tupd(i) = MIN(tupd(i), t1_2d)
         t12min(i) = MIN(t12min(i), t1_2d)
      ENDDO 
      ! ZY (VN) plane
      DO i=1,n
         t2_2d = FTEIK_HUGE
         tvmtn = tv(i) - tn(i)
         !IF ( (tv(i) < tn(i) + dy*slow5(i)) .AND. (tn(i) < tv(i) + dz*slow5(i)) ) THEN
         IF (tvmtn < dy*slow5(i) .AND.-tvmtn < dz*slow5(i)) THEN
            sref2 = slow5(i)*slow5(i)
            !ta = tv(i) - tn(i) + tnv(i)
            !tb = tn(i) - tv(i) + tnv(i)
            ta = tvmtn + tnv(i) !ta = tv(i) - tn(i) + tnv(i)
            tb =-tvmtn + tnv(i) !tb = tn(i) - tv(i) + tnv(i)
            tab = ta - tb
            tab2 = tab*tab
            four_sref2 = four*sref2
            t2_2d = ( (ta*dz2i + tb*dy2i) &
                     + SQRT(four_sref2*dz2i_p_dy2i - dz2i_dy2i*tab2) )*dz2i_p_dy2i_inv
         ENDIF
         !tupd(i) = MIN(tupd(i), t2_2d)
         t12min(i) = MIN(t12min(i), t2_2d)
      ENDDO
      ! XY (EN) Plane
      DO i=1,n
         t3_2d = FTEIK_HUGE
         temtn = te(i) - tn(i)
         !IF ( (te(i) < tn(i) + dy*slow6(i)) .AND. (tn(i) < te(i) + dx*slow6(i)) ) THEN
         IF (temtn < dy*slow6(i) .AND.-temtn < dx*slow6(i)) THEN
            sref2 = slow6(i)*slow6(i)
            !ta = te(i) - tn(i) + ten(i)
            !tb = tn(i) - te(i) + ten(i)
            ta = temtn + ten(i) !ta = te(i) - tn(i) + ten(i)
            tb =-temtn + ten(i) !tb = tn(i) - te(i) + ten(i)
            tab = ta - tb
            tab2 = tab*tab
            four_sref2 = four*sref2
            t3_2d = ( (ta*dx2i + tb*dy2i) &
                     + SQRT(four_sref2*dx2i_p_dy2i - dx2i_dy2i*tab2) )*dx2i_p_dy2i_inv
         ENDIF
         !tupd(i) = MIN(tupd(i), t3_2d)
         t12min(i) = MIN(t12min(i), t3_2d)
      ENDDO
      !DO i=1,n
      !   t12min(i) = MIN(t1v(i), t2v(i))
      !ENDDO
      !-------------------------------------3D operators---------------------------------!
      DO i=1,n
         t3d = FTEIK_HUGE
         tmax = MAX(tv(i), te(i), tn(i))
         IF (t12min(i) > tmax) THEN !IF (tupd(i) > tmax) THEN
            sref2 = slow7(i)*slow7(i)
            ta = te(i) + half*(-tn(i) + ten(i) - tv(i) + tev(i)) - tnv(i) + tnve(i) ! X
            tb = tv(i) + half*(-tn(i) + tnv(i) - te(i) + tev(i)) - ten(i) + tnve(i) ! Z
            tc = tn(i) + half*(-te(i) + ten(i) - tv(i) + tnv(i)) - tev(i) + tnve(i) ! Y
            tab = ta - tb
            tbc = tb - tc
            tac = ta - tc
            t2 = sref2*dsum_nine
            tab2 = tab*tab
            tbc2 = tbc*tbc
            tac2 = tac*tac
            t3 = dz2dx2*tab2 + dz2dy2*tbc2 + dx2dy2*tac2
            IF (t2 >= t3) THEN
               t1 = tb*dz2i + ta*dx2i + tc*dy2i
               t3d = (t1 + SQRT(t2 - t3))*dsumi
            ENDIF
         ENDIF
         !tupd(i) = MIN(tupd(i), t3d)
         t3d = MIN(t3d, t12min(i))
         tupd(i) = MIN(tupd(i), t3d) !t12min(i), t3d)
      ENDDO 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      PURE REAL(C_DOUBLE)                                                  &
      FUNCTION fteik_localSolver_init64fF(tv, te, tn, tev,                 &
                                          ten, tnv, tnve, tt0,             &
                                          slow1, slow2, slow3, slow4,      &
                                          slow5, slow6, slow7,             &
                                          i, j, k,                         &
                                          sgntz, sgntx, sgnty,             &
                                          sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
      RESULT(tupd) BIND(C, NAME='fteik_localSolver_init64fF')
      USE ISO_C_BINDING
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, zero
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN), VALUE :: tv, te, tn, tev, ten, tnv, tnve, tt0
      REAL(C_DOUBLE), INTENT(IN), VALUE :: slow1, slow2, slow3, slow4, &
                                           slow5, slow6, slow7
      REAL(C_DOUBLE), INTENT(IN), VALUE :: sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, sgntx, sgnty, sgntz
      ! local variables
      REAL(C_DOUBLE) apoly, bpoly, cpoly, dpoly,                                &
                     four_sref2, sref2,                                         &
                     t0c, t1, t1d, t1d_t2d_min, t2, t3, t3d,                    &
                     t1_2d, t2_2d, t3_2d,                                       &
                     ta, tb, tab, tab2, tac, tac2, tbc, tbc2, tc,               &
                     taue, tauev, tauen, taun, taunv, taunve, tauv,             &
                     tmax, tmin, txc, tyc, tzc
      LOGICAL(C_BOOL) lcartSolver, ldo3D, ldoZX, ldoZXCart, ldoZXSphere,  &
                      ldoZY, ldoZYCart, ldoZYSphere, ldoXY, ldoXYCart, ldoXYSphere
         !------------------1D operators, (refracted times),set times to BIG-------------!
         ! V plane
         !t1d = MIN(tt0, tv + dz*slow1) !t1 = tv + dz*slowLoc(1)
         t1d = tv + dz*slow1 !t1 = tv + dz*slowLoc(1)
         ! WE plane
         t1d = MIN(t1d, te + dx*slow2) !t2 = te + dx*slowLoc(2)
         ! NS plane
         t1d = MIN(t1d, tn + dy*slow3) !t3 = tn + dy*slowLoc(3)
         ! Check time to see when to switch to plane approximation
         tmin = MIN(tv, te, tn)
         tmax = MAX(tv, te, tn)
         lcartSolver = .FALSE.
         IF (tmin > epsSolver) lcartSolver = .TRUE.
         ! Initialize the travel-time perturbation terms.  
         t0c = zero 
         tzc = zero
         txc = zero
         tyc = zero
         tauv = FTEIK_HUGE
         taue = FTEIK_HUGE
         taun = FTEIK_HUGE
         taun = FTEIK_HUGE
         tauev = FTEIK_HUGE
         tauen = FTEIK_HUGE
         taunv = FTEIK_HUGE
         taunve = FTEIK_HUGE
         ! I can't promise the travel-time perturbations will be used but they must
         ! at minimum be computed.  There may be a slight penalty incurred for extraneous
         ! computations but it is no worse than the original solver.
         IF (.NOT.lcartSolver) THEN
            ! This is the spherical solver.  It can only trigger if we are within the
            ! radius.  There may exist further times when the 2D operators do not apply.
            CALL fteik_localSolver_tAnaD64fF(t0c, tzc, txc, tyc, i, j, k,    &
                                             dz, dx, dy, zsa, xsa, ysa, szero)
            ! Convert times into pertubations
            tauv   = tv   - fteik_localSolver_tAna64fF(i-sgntz, j,       k,            &
                                                       dz, dx, dy, zsa, xsa, ysa, szero)
            taue   = te   - fteik_localSolver_tAna64fF(i,       j-sgntx, k,            &
                                                       dz, dx, dy, zsa, xsa, ysa, szero)
            taun   = tn   - fteik_localSolver_tAna64fF(i,       j,       k-sgnty,      &
                                                       dz, dx, dy, zsa, xsa, ysa, szero)
            tauev  = tev  - fteik_localSolver_tAna64fF(i-sgntz, j-sgntx, k,            &
                                                       dz, dx, dy, zsa, xsa, ysa, szero)
            tauen  = ten  - fteik_localSolver_tAna64fF(i,       j-sgntx, k-sgnty,      &
                                                       dz, dx, dy, zsa, xsa, ysa, szero)
            taunv  = tnv  - fteik_localSolver_tAna64fF(i-sgntz, j,       k-sgnty,      &
                                                       dz, dx, dy, zsa, xsa, ysa, szero)
            taunve = tnve - fteik_localSolver_tAna64fF(i-sgntz, j-sgntx, k-sgnty,      &
                                                       dz, dx, dy, zsa, xsa, ysa, szero)
         ENDIF
         !------------------------------------2D operators-------------------------------!
         ! These conditions will decide if the 2d operators are active
         ldoZX = .FALSE.
         ldoZY = .FALSE.
         ldoXY = .FALSE.
         IF ( (tv < te + dx*slow4) .AND. (te < tv + dz*slow4) ) ldoZX = .TRUE.
         IF ( (tv < tn + dy*slow5) .AND. (tn < tv + dz*slow5) ) ldoZY = .TRUE.
         IF ( (te < tn + dy*slow6) .AND. (tn < te + dx*slow6) ) ldoXY = .TRUE.
         ldoZXCart = .FALSE.
         ldoZYCart = .FALSE.
         ldoXYCart = .FALSE.
         IF (ldoZX .AND. (lcartSolver .OR. k /= ysi)) ldoZXCart = .TRUE.
         IF (ldoZY .AND. (lcartSolver .OR. j /= xsi)) ldoZYCart = .TRUE.
         IF (ldoXY .AND. (lcartSolver .OR. i /= zsi)) ldoXYCart = .TRUE.
         ldoZXSphere = .FALSE.
         ldoZYSphere = .FALSE.
         ldoXYSphere = .FALSE.
         IF (ldoZX .AND. .NOT.ldoZXCart) ldoZXSphere = .TRUE.
         IF (ldoZY .AND. .NOT.ldoZYCart) ldoZYSphere = .TRUE.
         IF (ldoXY .AND. .NOT.ldoXYCart) ldoXYSphere = .TRUE.
         ! ZX (VE) plane
         t1_2d = FTEIK_HUGE
         IF (ldoZXCart) THEN
            sref2 = slow4*slow4
            ta = tev + te - tv
            tb = tev - te + tv
            tab = ta - tb
            tab2 = tab*tab
            four_sref2 = 4.d0*sref2
            t1_2d = ( (tb*dz2i + ta*dx2i) &
                    + SQRT(four_sref2*dz2i_p_dx2i - dz2i_dx2i*tab2) )*dz2i_p_dx2i_inv
         ENDIF
         IF (ldoZXSphere) THEN
            sref2 = slow4*slow4
            ta = tauev + taue - tauv
            tb = tauev - taue + tauv
            apoly = dz2i_p_dx2i
            bpoly = 4.d0*(sgnrx_dxi*txc + sgnrz_dzi*tzc) - 2.d0*(ta*dx2i + tb*dz2i)
            cpoly = ((ta*ta)*dx2i) + ((tb*tb)*dz2i)            &
                  - 4.d0*(sgnrx_dxi*txc*ta + sgnrz_dzi*tzc*tb) &
                  + 4.d0*(szero2 - sref2 + tyc*tyc)
            dpoly = bpoly*bpoly - 4.d0*apoly*cpoly
            IF (dpoly >= 0.d0) t1_2d = (SQRT(dpoly) - bpoly)/(2.d0*apoly) + t0c
            IF (t1_2d < tv .OR. t1_2d < te) t1_2d = FTEIK_HUGE
         ENDIF
         t1d_t2d_min = MIN(t1d, t1_2d)
         ! ZY (VN) plane
         t2_2d = FTEIK_HUGE
         IF (ldoZYCart) THEN
            sref2 = slow5*slow5
            ta = tv - tn + tnv 
            tb = tn - tv + tnv 
            tab = ta - tb
            tab2 = tab*tab
            four_sref2 = 4.d0*sref2
            t2_2d = ( (ta*dz2i + tb*dy2i) &
                     + SQRT(four_sref2*dz2i_p_dy2i - dz2i_dy2i*tab2) )*dz2i_p_dy2i_inv
         ENDIF
         IF (ldoZYSphere) THEN
            sref2 = slow5*slow5
            ta = tauv - taun + taunv ! Z
            tb = taun - tauv + taunv ! Y
            apoly = dz2i_p_dy2i
            bpoly = 4.d0*(sgnry_dyi*tyc + sgnrz_dzi*tzc) - 2.d0*(ta*dz2i + tb*dy2i)
            cpoly = ((ta*ta)*dz2i) + ((tb*tb)*dy2i)            &
                  - 4.d0*(sgnrz_dzi*tzc*ta + sgnry_dyi*tyc*tb) &
                  + 4.d0*(szero2 - sref2 + txc*txc)
            dpoly = bpoly*bpoly - 4.d0*apoly*cpoly
            IF (dpoly >= 0.d0) t2_2d = (SQRT(dpoly) - bpoly)/(2.d0*apoly) + t0c
            IF (t2_2d < tv .OR. t2_2d < tn) t2_2d = FTEIK_HUGE
         ENDIF
         t1d_t2d_min = MIN(t1d_t2d_min, t2_2d)
         ! XY (EN) Plane
         t3_2d = FTEIK_HUGE
         IF (ldoXYCart) THEN
            sref2 = slow6*slow6
            ta = te - tn + ten 
            tb = tn - te + ten 
            tab = ta - tb
            tab2 = tab*tab
            four_sref2 = 4.d0*sref2
            t3_2d = ( (ta*dx2i + tb*dy2i) &
                     + SQRT(four_sref2*dx2i_p_dy2i - dx2i_dy2i*tab2) )*dx2i_p_dy2i_inv
         ENDIF
         IF (ldoXYSphere) THEN
            sref2 = slow6*slow6
            ta = taue - taun + tauen ! X
            tb = taun - taue + tauen ! Y
            apoly = dx2i_p_dy2i
            bpoly = 4.d0*(sgnry_dyi*tyc + sgnrx_dxi*txc) - 2.d0*(ta*dx2i + tb*dy2i)
            cpoly = ((ta*ta)*dx2i)+((tb*tb)*dy2i)              &
                  - 4.d0*(sgnrx_dxi*txc*ta + sgnry_dyi*tyc*tb) &
                  + 4.d0*(szero2 - sref2 + tzc*tzc)
            dpoly = bpoly*bpoly - 4.d0*apoly*cpoly
            IF (dpoly >= 0.d0) t3_2d = (SQRT(dpoly) - bpoly)/(2.d0*apoly) + t0c
            IF (t3_2d < te .OR. t3_2d < tn) t3_2d = FTEIK_HUGE
         ENDIF
         t1d_t2d_min = MIN(t1d_t2d_min, t3_2d)
         !------------------------------------3D operators-------------------------------!
         t3d = FTEIK_HUGE
         ldo3D = .FALSE.
         IF (t1d_t2d_min > tmax) ldo3D = .TRUE.
         IF (ldo3D .AND. lcartSolver) THEN 
            sref2 = slow7*slow7
            ta = te + 0.5d0*(-tn + ten - tv + tev) - tnv + tnve ! X
            tb = tv + 0.5d0*(-tn + tnv - te + tev) - ten + tnve ! Z
            tc = tn + 0.5d0*(-te + ten - tv + tnv) - tev + tnve ! Y
            tab = ta - tb
            tbc = tb - tc
            tac = ta - tc
            t2 = sref2*dsum_nine !dsum*9.d0
            tab2 = tab*tab
            tbc2 = tbc*tbc
            tac2 = tac*tac
            t3 = dz2dx2*tab2 + dz2dy2*tbc2 + dx2dy2*tac2
            IF (t2 >= t3) THEN
               t1 = tb*dz2i + ta*dx2i + tc*dy2i
               t3d = (t1 + SQRT(t2 - t3))*dsumi
            ENDIF
         ENDIF
         IF (ldo3D .AND. .NOT.lcartSolver) THEN
            sref2 = slow7*slow7
            ta = taue + 0.5d0*(-taun + tauen - tauv + tauev) - taunv + taunve ! X
            tb = tauv + 0.5d0*(-taun + taunv - taue + tauev) - tauen + taunve ! Z
            tc = taun + 0.5d0*(-taue + tauen - tauv + taunv) - tauev + taunve ! Y
            apoly = dx2i + dz2i + dy2i
            bpoly = 6.d0*(sgnrz_dzi*tzc + sgnry_dyi*tyc + sgnrx_dxi*txc) &
                  - 2.d0*(ta*dx2i + tb*dz2i + tc*dy2i)
            cpoly = ((ta*ta)*dx2i) + ((tb*tb)*dz2i) + ((tc*tc)*dy2i)               &
                  - 6.d0*( sgnrx_dxi*txc*ta + sgnrz_dzi*tzc*tb + sgnry_dyi*tyc*tc) &
                  + 9.d0*(szero2 - sref2)
            dpoly = bpoly*bpoly - 4.d0*apoly*cpoly
            IF (dpoly >= 0.d0) t3d = (SQRT(dpoly) - bpoly)/(2.d0*apoly) + t0c
            IF (t3d < te .OR. t3d < tn .OR. t3d < tv) t3d = FTEIK_HUGE
         ENDIF
         tupd = MIN(t1d_t2d_min, t3d)
         tupd = MIN(tt0, tupd)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Compute analytic travel time at the point (i, j, k) in a homogeneous model.
!>
!>    @param[in] i      iz'th grid point (Fortran numbered).
!>    @param[in] j      ix'th grid point (Fortran numbered).
!>    @param[in] k      iy'th grid point (Fortran numbered).
!>    @param[in] dz     Grid spacing (meters) in z.
!>    @param[in] dx     Grid spacing (meters) in x.
!>    @param[in] dy     Grid spacing (meters) in y.
!>    @param[in] zsa    Source offset (meters) in z.
!>    @param[in] xsa    Source offset (meters) in x.
!>    @param[in] ysa    Source offset (meters) in y.
!>    @param[in] szero  Slowness at source (s/m).  
!>
!>    @result The travel time from the source at (zsa,xsa,ysa) to the (i,j,k)'th grid
!>            point given a constant slowness around the source.
!>
!>    @author Mark Noble, Alexandrine, Gesret, and Ben Baker.
!>
!>    @version 2
!>
!>    @date July 2017
!>
!>    @copyright CeCILL-3
!>
!>
      PURE REAL(C_DOUBLE)                                        &
      FUNCTION fteik_localSolver_tAna64fF(i, j, k, dz, dx, dy,   &
                                          zsa, xsa, ysa, szero)  &
      BIND(C, NAME='fteik_localSolver_tAna64fF')
      !$OMP DECLARE SIMD(fteik_localSolver_tAna64fF) &
      !$OMP UNIFORM(dz, dx, dy, zsa, xsa, ysa, szero)
      USE ISO_C_BINDING
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN), VALUE :: dz, dx, dy, szero, zsa, xsa, ysa
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k
      REAL(C_DOUBLE) d0, diffz, diffz2, diffx, diffx2, diffy, diffy2
      diffz = (DBLE(i) - zsa)*dz
      diffx = (DBLE(j) - xsa)*dx
      diffy = (DBLE(k) - ysa)*dy
      diffz2 = diffz*diffz
      diffx2 = diffx*diffx
      diffy2 = diffy*diffy
      d0 = diffz2 + diffx2 + diffy2
      fteik_localSolver_tAna64fF = szero*SQRT(d0) !sqrt(diffz2 + diffx2 + diffy2);
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
!>    @param[in] k      iy'th grid point (Fortran numbered).
!>    @param[in] dz     Grid spacing (meters) in z.
!>    @param[in] dx     Grid spacing (meters) in x.
!>    @param[in] dy     Grid spacing (meters) in y.
!>    @param[in] zsa    Source offset (meters) in z.
!>    @param[in] xsa    Source offset (meters) in x.
!>    @param[in] ysa    Source offset (meters) in y.
!>    @param[in] szero  Slowness at source (s/m).  
!>
!>    @param[out] t_anad   Derivative of analytic travel time at point (i,j,k).
!>    @param[out] tzc      Derivative of analytic travel time in z.
!>    @param[out] txc      Derivative of analytic travel time in x.
!>    @param[out] tyc      Derivative of analytic travel time in y.
!>
!>    @author Mark Noble, Alexandrine, Gesret, and Ben Baker.
!>
!>    @version 2
!>
!>    @date July 2017
!>
!>    @copyright CeCILL-3
!>
      PURE SUBROUTINE fteik_localSolver_tAnaD64fF(t_anad,        &
                                                  tzc, txc, tyc, &
                                                  i, j, k,       &
                                                  dz, dx, dy,    &
                                                  zsa, xsa, ysa, &
                                                  szero)         &
                      BIND(C, NAME='fteik_localSolver_tAnaD64fF')
      !$OMP DECLARE SIMD(fteik_localSolver_tAnaD64fF) &
      !$OMP UNIFORM(dz, dx, dy, zsa, xsa, ysa, szero)
      USE ISO_C_BINDING
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN), VALUE :: dz, dx, dy, zsa, xsa, ysa, szero
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k
      REAL(C_DOUBLE), INTENT(OUT) :: t_anad, tzc, txc, tyc
      REAL(C_DOUBLE) d0, d0i_szero, diffz, diffz2, diffx, diffx2, diffy, diffy2, sqrtd0
      REAL(C_DOUBLE), PARAMETER :: zero = 0.d0
      diffz = (DBLE(i) - zsa)*dz
      diffx = (DBLE(j) - xsa)*dx
      diffy = (DBLE(k) - ysa)*dy
      diffz2 = diffz*diffz
      diffx2 = diffx*diffx
      diffy2 = diffy*diffy
      d0 = diffz2 + diffx2 + diffy2

      d0i_szero = zero
      t_anad = zero
      IF (d0 > zero) THEN
         sqrtd0 = SQRT(d0)
         t_anad = szero*sqrtd0
         d0i_szero = szero/sqrtd0
      ENDIF
      tzc = d0i_szero*diffz
      txc = d0i_szero*diffx
      tyc = d0i_szero*diffy
      RETURN
      END SUBROUTINE

!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility routine for computing some grid spacing constants use by the 
!>           finite difference stencils in the local solver.
!>
!>    @param[in] isrc     Source number.
!>    @param[in] epsS2C   Controls the transition from the spherical to the
!>                        Cartesian solver.
!>
!>    @param[out] dest    Maps from source 8 nodes near source to index in
!>                        ttimes.  This is a vector of dimension [8].
!>    @param[out] ts8     Travel times around source grid-point (seconds).
!>                        this is a vector of dimension [8].
!>    @param[out] ierr    0 indicates success.
!>
!>    @author Mark Noble, Alexandrine Gesret, and Ben Baker
!>
!>    @version 2
!>
!>    @date July 2017
!>
!>    Note, this now updates private variables in this module and computes a 
!>    few additional terms that are constantly used in the local solver.
!>
      SUBROUTINE fteik_localSolver_initialize64fF(isrc, epsS2C, dest, &
                                                  ts8, ierr)  &
      BIND(C, NAME='fteik_localSolver_initialize64fF')
      USE FTEIK_SOURCE64F, ONLY : fteik_source_getSolverInfo64fF
      USE FTEIK_MODEL64F, ONLY : fteik_model_grid2indexF
      USE FTEIK_MODEL64F, ONLY : nz, nzx
      USE FTEIK_MODEL64F, ONLY : dzIn => dz
      USE FTEIK_MODEL64F, ONLY : dxIn => dx
      USE FTEIK_MODEL64F, ONLY : dyIn => dy 
      USE FTEIK_CONSTANTS64F, ONLY : one, nine
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: epsS2C
      REAL(C_DOUBLE), INTENT(OUT) :: ts8(8)
      INTEGER(C_INT), INTENT(OUT) :: dest(8), ierr
      ! Things to copy from model
      dz = dzIn
      dx = dxIn
      dy = dyIn
      ! Things to grab from the source
      CALL fteik_source_getSolverInfo64fF(isrc,          &
                                          zsi, xsi, ysi, &
                                          zsa, xsa, ysa, &
                                          szero, szero2, &
                                          ierr) 
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_localSolver_meshConstants64fF: Problem with source'
         RETURN
      ENDIF
      ! Time after which solver switches from spherical to cartesian 
      epsSolver = epsS2C*MIN(dz, dx, dy)*szero
      ! Derivative products
      dz2 = dz*dz
      dx2 = dx*dx
      dy2 = dy*dy
      dzi = one/dz
      dxi = one/dx
      dyi = one/dy
      dz2i = one/dz2
      dx2i = one/dx2
      dy2i = one/dy2
      dsum = dz2i + dx2i + dy2i
      dsum_nine = dsum*nine
      dsumi = one/dsum
      dz2dx2 = one/(dz2*dx2)
      dz2dy2 = one/(dz2*dy2)
      dx2dy2 = one/(dx2*dy2)
      dz2i_dx2i = dz2i*dx2i
      dz2i_dy2i = dz2i*dy2i
      dx2i_dy2i = dx2i*dy2i
      dz2i_p_dy2i_inv = one/(dz2i + dy2i)
      dz2i_p_dx2i_inv = one/(dz2i + dx2i)
      dx2i_p_dy2i_inv = one/(dx2i + dy2i)
      dz2i_p_dx2i = dz2i + dx2i
      dz2i_p_dy2i = dz2i + dy2i
      dx2i_p_dy2i = dx2i + dy2i
      ! Compute local travel-times near source
      ts8(1) = fteik_localSolver_tAna64fF(zsi,   xsi,   ysi,              &
                                          dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(2) = fteik_localSolver_tAna64fF(zsi+1, xsi,   ysi,              &
                                          dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(3) = fteik_localSolver_tAna64fF(zsi,   xsi+1, ysi,              &
                                          dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(4) = fteik_localSolver_tAna64fF(zsi+1, xsi+1, ysi,              &
                                          dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(5) = fteik_localSolver_tAna64fF(zsi,   xsi,   ysi+1,            &
                                          dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(6) = fteik_localSolver_tAna64fF(zsi+1, xsi,   ysi+1,            &
                                          dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(7) = fteik_localSolver_tAna64fF(zsi,   xsi+1, ysi+1,            &
                                          dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(8) = fteik_localSolver_tAna64fF(zsi+1, xsi+1, ysi+1,            &
                                          dz, dx, dy, zsa, xsa, ysa, szero)
      dest(1) = fteik_model_grid2indexF(zsi,   xsi  , ysi, nz, nzx)
      dest(2) = fteik_model_grid2indexF(zsi+1, xsi  , ysi, nz, nzx)
      dest(3) = fteik_model_grid2indexF(zsi  , xsi+1, ysi, nz, nzx)
      dest(4) = fteik_model_grid2indexF(zsi+1, xsi+1, ysi, nz, nzx)
      dest(5) = fteik_model_grid2indexF(zsi,   xsi  , ysi+1, nz, nzx)
      dest(6) = fteik_model_grid2indexF(zsi+1, xsi  , ysi+1, nz, nzx)
      dest(7) = fteik_model_grid2indexF(zsi,   xsi+1, ysi+1, nz, nzx)
      dest(8) = fteik_model_grid2indexF(zsi+1, xsi+1, ysi+1, nz, nzx)
      RETURN
      END

      PURE INTEGER(C_INT)                                             &
      FUNCTION fteik_localSolver_grid2indexF(i, j, k, nz, nzx)        &
      BIND(C, NAME='fteik_localSolver_grid2indexF') RESULT(grid2indexF)
      !$OMP DECLARE SIMD(fteik_localSolver_grid2indexF) UNIFORM(nz, nzx)
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN), VALUE :: i, j, k, nz, nzx 
      grid2indexF = (k - 1)*nzx + (j - 1)*nz + i 
      RETURN
      END FUNCTION

!----------------------------------------------------------------------------------------!
!                                      End the Code                                      !
!----------------------------------------------------------------------------------------!
END MODULE
