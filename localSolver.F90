
      PURE REAL(C_DOUBLE)                                        &
      FUNCTION fteik_localSolverNoInit64fF(tv, te, tn, tev,      &  
                                     ten, tnv, tnve, tt0,        &
                                     slow1, slow2, slow3, slow4, &
                                     slow5, slow6, slow7)        &
      RESULT(tupd) BIND(C, NAME='fteik_localSolverNoInit64fF')
      !$OMP DECLARE SIMD(fteik_localSolverNoInit64fF)
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_tAnaD64fF, fteik_tAna64fF
      USE FTEIK_UTILS64F, ONLY : FTEIK_HUGE
      USE FTEIK_UTILS64F, ONLY : dx, dy, dz, &
                                 dx2i, dy2i, dz2i, &
                                 dsum, dsumi, &
                                 dz2dx2, dz2dy2, dx2dy2, &
                                 dz2i_dx2i, dz2i_dy2i, dx2i_dy2i, &
                                 dz2i_p_dx2i, dz2i_p_dy2i, dx2i_p_dy2i, &
                                 dz2i_p_dy2i_inv, dz2i_p_dx2i_inv, dx2i_p_dy2i_inv
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN), VALUE :: tv, te, tn, tev, ten, tnv, tnve, tt0 
      REAL(C_DOUBLE), INTENT(IN), VALUE :: slow1, slow2, slow3, slow4, &
                                           slow5, slow6, slow7
      ! local variables
      REAL(C_DOUBLE) four_sref2, sref2,                                      &
                     t1, t1d, t1d_t2d_min, t2, t3, t3d,                      &
                     t1_2d, t2_2d, t3_2d,                                    &
                     ta, tb, tab, tab2, tac, tac2, tbc, tbc2, tc,            &
                     tmax, tmin
      !--------------------1D operators, (refracted times),set times to BIG--------------!
      ! V plane
      t1d = MIN(tt0, tv + dz*slow1) !t1 = tv + dz*slowLoc(1)
      ! WE plane
      t1d = MIN(t1d, te + dx*slow2) !t2 = te + dx*slowLoc(2)
      ! NS plane
      t1d = MIN(t1d, tn + dy*slow3) !t3 = tn + dy*slowLoc(3)
      ! Check time to see when to switch to plane approximation
      tmin = MIN(tv, te, tn)
      tmax = MAX(tv, te, tn)
      !-------------------------------------2D operators---------------------------------!
      ! ZX (VE) plane
      t1_2d = FTEIK_HUGE
      IF ( (tv < te + dx*slow4) .AND. (te < tv + dz*slow4) ) THEN
         sref2 = slow4*slow4
         ta = tev + te - tv
         tb = tev - te + tv
         tab = ta - tb
         tab2 = tab*tab
         four_sref2 = 4.d0*sref2
         t1_2d = ( (tb*dz2i + ta*dx2i) &
                  + SQRT(four_sref2*dz2i_p_dx2i - dz2i_dx2i*tab2) )*dz2i_p_dx2i_inv
      ENDIF
      t1d_t2d_min = MIN(t1d, t1_2d)
      ! ZY (VN) plane
      t2_2d = FTEIK_HUGE
      IF ( (tv < tn + dy*slow5) .AND. (tn < tv + dz*slow5) ) THEN
         sref2 = slow5*slow5
         ta = tv - tn + tnv
         tb = tn - tv + tnv
         tab = ta - tb
         tab2 = tab*tab
         four_sref2 = 4.d0*sref2
         t2_2d = ( (ta*dz2i + tb*dy2i) &
                  + SQRT(four_sref2*dz2i_p_dy2i - dz2i_dy2i*tab2) )*dz2i_p_dy2i_inv
      ENDIF
      t1d_t2d_min = MIN(t1d_t2d_min, t2_2d)
      ! XY (EN) Plane
      t3_2d = FTEIK_HUGE
      IF ( (te < tn + dy*slow6) .AND. (tn < te + dx*slow6) ) THEN
         sref2 = slow6*slow6
         ta = te - tn + ten
         tb = tn - te + ten
         tab = ta - tb
         tab2 = tab*tab
         four_sref2 = 4.d0*sref2
         t3_2d = ( (ta*dx2i + tb*dy2i) &
                  + SQRT(four_sref2*dx2i_p_dy2i - dx2i_dy2i*tab2) )*dx2i_p_dy2i_inv
      ENDIF
      t1d_t2d_min = MIN(t1d_t2d_min, t3_2d)
      !-------------------------------------3D operators---------------------------------!
      t3d = FTEIK_HUGE
      IF (t1d_t2d_min > tmax) THEN
         sref2 = slow7*slow7
         ta = te + 0.5d0*(-tn + ten - tv + tev) - tnv + tnve ! X
         tb = tv + 0.5d0*(-tn + tnv - te + tev) - ten + tnve ! Z
         tc = tn + 0.5d0*(-te + ten - tv + tnv) - tev + tnve ! Y
         tab = ta - tb
         tbc = tb - tc
         tac = ta - tc
         t2 = sref2*dsum*9.d0
         tab2 = tab*tab
         tbc2 = tbc*tbc
         tac2 = tac*tac
         t3 = dz2dx2*tab2 + dz2dy2*tbc2 + dx2dy2*tac2
         IF (t2 >= t3) THEN
            t1 = tb*dz2i + ta*dx2i + tc*dy2i
            t3d = (t1 + SQRT(t2 - t3))*dsumi
         ENDIF
      ENDIF
      tupd = MIN(t1d_t2d_min, t3d)
      RETURN
      END
 

!                                                                                        !
!========================================================================================!
!                                                                                        !
      PURE REAL(C_DOUBLE)                                             &
      FUNCTION fteik_localSolverInit64fF(tv, te, tn, tev,             &
                                     ten, tnv, tnve, tt0,             &
                                     slow1, slow2, slow3, slow4,      &    
                                     slow5, slow6, slow7,             &
                                     i, j, k,                         &
                                     sgntz, sgntx, sgnty,             &
                                     sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
      RESULT(tupd) BIND(C, NAME='fteik_localSolverInit64fF')
      !$OMP DECLARE SIMD(fteik_localSolverInit64fF) &
      !$OMP UNIFORM(sgntz, sgntx, sgnty, sgnrz_dzi, sgnrx_dxi, sgnry_dyi)
      USE ISO_C_BINDING
      USE FTEIK_UTILS64F, ONLY : fteik_tAnaD64fF, fteik_tAna64fF
      USE FTEIK_UTILS64F, ONLY : FTEIK_HUGE
      USE FTEIK_UTILS64F, ONLY : dx, dy, dz, xsa, ysa, zsa, &
                                 dx2i, dy2i, dz2i, &
                                 dsum, dsumi, epsSolver, &
                                 dz2dx2, dz2dy2, dx2dy2, &
                                 dz2i_dx2i, dz2i_dy2i, dx2i_dy2i, &
                                 dz2i_p_dx2i, dz2i_p_dy2i, dx2i_p_dy2i, &
                                 dz2i_p_dy2i_inv, dz2i_p_dx2i_inv, dx2i_p_dy2i_inv, &
                                 szero, szero2
      USE FTEIK_UTILS64F, ONLY : xsi, ysi, zsi
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
      !--------------------1D operators, (refracted times),set times to BIG--------------!
      ! V plane
      t1d = MIN(tt0, tv + dz*slow1) !t1 = tv + dz*slowLoc(1)
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
      t0c = 0.d0
      tzc = 0.d0
      txc = 0.d0
      tyc = 0.d0
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
         CALL fteik_tAnaD64fF(t0c, tzc, txc, tyc, i, j, k,    &
                              dz, dx, dy, zsa, xsa, ysa, szero)
         ! Convert times into pertubations
         tauv   = tv   - fteik_tAna64fF(i-sgntz, j,       k,            &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
         taue   = te   - fteik_tAna64fF(i,       j-sgntx, k,            &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
         taun   = tn   - fteik_tAna64fF(i,       j,       k-sgnty,      &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
         tauev  = tev  - fteik_tAna64fF(i-sgntz, j-sgntx, k,            &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
         tauen  = ten  - fteik_tAna64fF(i,       j-sgntx, k-sgnty,      &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
         taunv  = tnv  - fteik_tAna64fF(i-sgntz, j,       k-sgnty,      &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
         taunve = tnve - fteik_tAna64fF(i-sgntz, j-sgntx, k-sgnty,      &
                                        dz, dx, dy, zsa, xsa, ysa, szero)
      ENDIF
      !-------------------------------------2D operators---------------------------------!
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
               + 4.d0*(szero2 - sref2 + tzc*tzc);
         dpoly = bpoly*bpoly - 4.d0*apoly*cpoly;
         IF (dpoly >= 0.d0) t3_2d = (SQRT(dpoly) - bpoly)/(2.d0*apoly) + t0c
         IF (t3_2d < te .OR. t3_2d < tn) t3_2d = FTEIK_HUGE
      ENDIF
      t1d_t2d_min = MIN(t1d_t2d_min, t3_2d)
      !-------------------------------------3D operators---------------------------------!
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
         t2 = sref2*dsum*9.d0
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
      RETURN
      END
