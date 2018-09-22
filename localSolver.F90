!> @defgroup localSolver2d 2D Local Solver
!> @brief The 2D local solver.
!> @ingroup solver2d
MODULE FTEIK_LOCALSOLVER2D64F
      USE ISO_FORTRAN_ENV
      IMPLICIT NONE

      DOUBLE PRECISION, PRIVATE, SAVE :: szero2
      DOUBLE PRECISION, SAVE, PRIVATE :: dx, dxi, dz, dzi, dz2i, dx2i, dz2i_dx2i, &
                                         dz2i_p_dx2i, dz2i_p_dx2i_inv, epsSolver
      DOUBLE PRECISION, PROTECTED, SAVE :: szero, zsa, xsa
      INTEGER, SAVE, PROTECTED :: zsi, xsi


      PUBLIC :: fteik_localSolver2d_tAnaD64f
      PUBLIC :: fteik_localSolver2d_tAna64f
      CONTAINS
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
!>    @version 2
!>    @date January 2018
!>    @copyright MIT
!>    @ingroup localSolver2d
      PURE DOUBLE PRECISION                                 &
      FUNCTION fteik_localSolver2d_tAna64f(i, j, dz, dx,    &    
                                           zsa, xsa, szero)
      !!$OMP DECLARE SIMD(fteik_localSolver2d_tAna64f) &
      !!$OMP UNIFORM(dz, dx, zsa, xsa, szero) 
      DOUBLE PRECISION, INTENT(IN) :: dz, dx, zsa, xsa, szero
      INTEGER, INTENT(IN) :: i, j 
      DOUBLE PRECISION diffx, diffz
      diffz = (DBLE(i) - zsa)*dz
      diffx = (DBLE(j) - xsa)*dx
      fteik_localSolver2d_tAna64f = szero*HYPOT(diffx, diffz)
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
!>    @version 2
!>    @date January 2018
!>    @copyright MIT
!>    @ingroup localSolver2d
      PURE SUBROUTINE fteik_localSolver2d_tAnaD64f(t_anad, tzc, txc, &
                                                   i, j,             &
                                                   dz, dx,           &
                                                   zsa, xsa,         &
                                                   szero)
      !!$OMP DECLARE SIMD(fteik_localSolver2d_tAnaD64f) &
      !!$OMP UNIFORM(dz, dx, zsa, xsa, szero) 
      DOUBLE PRECISION, INTENT(IN) :: dz, dx, zsa, xsa, szero
      INTEGER, INTENT(IN) :: i, j
      DOUBLE PRECISION, INTENT(OUT) :: t_anad, tzc, txc
      DOUBLE PRECISION d0, d0i_szero, diffx, diffx2, diffz, diffz2, sqrtd0
      DOUBLE PRECISION, PARAMETER :: zero = 0.d0
      diffz = (DBLE(i) - zsa)*dz
      diffx = (DBLE(j) - xsa)*dx
      diffz2 = diffz*diffz
      diffx2 = diffx*diffx
      d0 = diffz2 + diffx2

      d0i_szero = zero
      t_anad = zero
      IF (d0 > zero) THEN
         sqrtd0 = DSQRT(d0)
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
!>    @brief 2D local solver for Cartesian frame only.
!>    @param[in] n       Number of nodes to update. 
!>    @param[in] ttvec   The local travel times (seconds) surrounding the point.  This is
!>                       a vector of dimension [4 x n] with leading dimemsion 4.  Each
!>                       four-tuple is packed [tv, te, tev, tt].
!>    @param[in] sloc    Slowness (s/m) at the finite difference points.  This is a vector
!>                       of dimension [n].
!>    @param[out] tupd   Updated travel times (seconds) at each node (seconds).  This is
!>                       a vector of dimension [n].
!>    @ingroup localSolver2d
      PURE SUBROUTINE fteik_localSolver2d_noInit64f(n, ttvec, sloc, tupd)
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, four, chunkSize
      INTEGER, INTENT(IN) :: n 
      DOUBLE PRECISION, INTENT(IN) :: ttvec(4*n), sloc(4*n)
      DOUBLE PRECISION, INTENT(OUT) :: tupd(n)
      DOUBLE PRECISION four_sref2, s1, s2, s3, sref2, t1_2d, t12min, ta, tb, tab, &
                       tab2, temtv, te, tev, tv, tt
      INTEGER i
      ! Loop on nodes in update
      DO i=1,n
         s1 = sloc(4*(i-1)+1)
         s2 = sloc(4*(i-1)+2)
         s3 = sloc(4*(i-1)+3)
         tv  = ttvec(4*(i-1)+1)
         te  = ttvec(4*(i-1)+2)
         tev = ttvec(4*(i-1)+3)
         tt  = ttvec(4*(i-1)+4)
         temtv = te - tv 
         ! 1D operators (refracted times)
         t12min = MIN(tv + dz*s1, te + dx*s2)
         ! 2D operator
         t1_2d = FTEIK_HUGE
         IF (temtv < dz*s3 .AND. -temtv < dx*s3) THEN 
            sref2 = s3*s3
            ta = tev + temtv
            tb = tev - temtv
            tab = ta - tb 
            tab2 = tab*tab
            four_sref2 = four*sref2
            t1_2d = ( (tb*dz2i + ta*dx2i) &
                     + DSQRT(four_sref2*dz2i_p_dx2i - dz2i_dx2i*tab2) )*dz2i_p_dx2i_inv
         ENDIF
         tupd(i) = MIN(t12min, t1_2d)
      ENDDO
      RETURN
      END

      SUBROUTINE fteik_localSolver2d_getSweepSigns(sweep,        &
                                                   sgntz, sgntx, &
                                                   sgnvz, sgnvx, &
                                                   ierr)
      INTEGER, INTENT(IN) :: sweep
      INTEGER, INTENT(OUT) :: sgnvx, sgnvz, sgntx, sgntz, ierr
      ierr = 0
      IF (sweep == 1) THEN 
         sgntz = 1
         sgntx = 1
         sgnvz = 1
         sgnvx = 1
      ELSEIF (sweep == 2) THEN 
         sgntz = 1
         sgntx =-1
         sgnvz = 1
         sgnvx = 0
      ELSEIF (sweep == 3) THEN 
         sgntz =-1
         sgntx = 1
         sgnvz = 0
         sgnvx = 1
      ELSEIF (sweep == 4) THEN
         sgntz =-1
         sgntx =-1
         sgnvz = 0
         sgnvx = 0
      ELSE
         WRITE(ERROR_UNIT,900) sweep
         ierr = 1
      ENDIF
  900 FORMAT('fteik_localSolver2d_getSweepSigns: Invalid sweep = ', I0)
      RETURN
      END

      SUBROUTINE fteik_localSolver2d_init64f(sgntz, sgntx,   &
                                             sgnvz, sgnvx,   &
                                             nz, nx, iz, ix, &
                                             slow, tt, tupd)
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, four
      INTEGER, INTENT(IN) :: sgntz, sgntx, sgnvz, sgnvx, nz, nx, iz, ix
      DOUBLE PRECISION, INTENT(OUT) :: tupd
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tt, slow
      DOUBLE PRECISION apoly, bpoly, cpoly, dpoly, &
                       ta, tb, taue, tauev, tauv, tdiag, tv, te, tev, t, &
                       t1, t2, t3, t1d, t2d, &
                       sref1, sref2, sref3, sgnrx, sgnrz, &
                       sgnrx_txc, sgnrz_tzc, sgnrx_txc_dxi, sgnrz_tzc_dzi, t0c, txc, tzc
      INTEGER i1, indx, indxv, indxe, indxev, j1
      INTEGER indx1, indx2, indx3, indx4, indx5
      !LOGICAL, PARAMETER :: linitk = .TRUE.
      sgnrz = DBLE(sgntz)
      sgnrx = DBLE(sgntx)
      indxv  = (ix - 1)*nz         + iz - sgntz
      indxe  = (ix - sgntx - 1)*nz + iz
      indxev = (ix - sgntx - 1)*nz + iz - sgntz
      indx   = (ix - 1)*nz         + iz
      tv  = tt(indxv)  !dble( tt(i-sgntz,j,k) )
      te  = tt(indxe)  !dble( tt(i,j-sgntx,k) )
      tev = tt(indxev) !dble( tt(i-sgntz,j-sgntx,k) )
      t   = tt(indx)
      ! Extract slownesses
      i1 = iz - sgnvz
      j1 = ix - sgnvx
      indx1 = (MAX(ix-1,1)  - 1)*(nz - 1) + i1 ! V Plane
      indx2 = (MIN(ix,nx-1) - 1)*(nz - 1) + i1 ! V Plane
      indx3 = (j1 - 1)*(nz - 1) + MAX(iz-1,1)  ! WE Plane
      indx4 = (j1 - 1)*(nz - 1) + MIN(iz,nz-1) ! WE Plane
      indx5 = (j1 - 1)*(nz - 1) + i1
      sref1 = MIN(slow(indx1), slow(indx2))
      sref2 = MIN(slow(indx3), slow(indx4))
      sref3 = slow(indx5)
      t1d = MIN(tv + dz*sref1, te + dx*sref2) ! min(V, WE) planes
      ! 2D and diagonal operators
      CALL fteik_localSolver2d_tAnaD64f(t0c,              &
                                        tzc, txc, iz, ix, &
                                        dz, dx, zsa, xsa, &
                                        szero)
      tauv = tv   - fteik_localSolver2d_tAna64f(iz-sgntz,   ix,        &
                                                dz, dx, zsa, xsa, szero)
      taue = te   - fteik_localSolver2d_tAna64f(iz, ix-sgntx,          &
                                                dz, dx, zsa, xsa, szero)
      tauev = tev - fteik_localSolver2d_tAna64f(iz-sgntz,   ix-sgntx,  &
                                                dz, dx, zsa, xsa, szero)
      sgnrz_tzc = sgnrz*tzc
      sgnrx_txc = sgnrx*txc
      sgnrz_tzc_dzi = sgnrz_tzc*dzi
      sgnrx_txc_dxi = sgnrx_txc*dxi
      ! Diagonal operator
      tdiag = tev + sref3*DSQRT(dx*dx + dz*dz)
      ! Choose spherical or plane wave
      t1 = FTEIK_HUGE
      t2 = FTEIK_HUGE
      t3 = FTEIK_HUGE
      ! Choose spherical or plane wave; first test for Plane wave
      IF ( ( ABS(iz - zsi) > epsSolver .OR. ABS(ix - xsi) > epsSolver ) ) THEN
         ! 4 Point operator, if possible otherwise do three points
         IF (tv <= te + dx*sref3 .AND. te <= tv + dz*sref3 .AND. &
             te >= tev .AND. tv >= tev ) THEN
            ta = tev + te - tv
            tb = tev - te + tv
            t1 = ( ( tb * dz2i + ta * dx2i ) + DSQRT( 4.d0 * sref3*sref3*( dz2i + dx2i ) &
                 - dz2i * dx2i * ( ta - tb ) * ( ta - tb ) ) ) / ( dz2i + dx2i )
         ! Two 3 point operators
         ELSEIF (( te - tev ) <= dz*dz*sref3/DSQRT(dx*dx + dz*dz) .AND. te > tev) THEN
            t2 = te + dx*DSQRT(sref3*sref3 - ((te - tev)/dz)**2)
         ELSEIF (( tv - tev ) <= dx*dx*sref3/DSQRT(dx*dx + dz*dz) .AND. tv > tev) THEN
            t3 = tv + dz*DSQRT(sref3*sref3 - ((tv - tev)/dx)**2)
         ENDIF
      ELSE
         ! Do spherical operator if conditions ok
         IF ( tv < te + dx*sref3 .AND. te < tv + dz*sref3 .AND. &
             te >= tev .AND. tv >= tev) THEN
            ta = tauev + taue - tauv   ! X
            tb = tauev - taue + tauv   ! Z
            apoly = dz2i + dx2i
            bpoly = 4.d0 * ( sgnrx * txc * dxi + sgnrz * tzc * dzi ) &
                  - 2.d0 * ( ta * dx2i + tb * dz2i )
            cpoly = ( ta*ta * dx2i ) + ( tb*tb * dz2i ) &
                  - 4.d0 * ( sgnrx * txc * dxi * ta + sgnrz * tzc * dzi * tb ) &
                  + 4.d0 * ( szero*szero - sref3*sref3)
            dpoly = bpoly*bpoly - 4.d0 * apoly * cpoly
            IF (dpoly >= 0.d0) t1 = 0.5d0*(DSQRT(dpoly) - bpoly )/apoly + t0c
            IF (t1 < tv .OR. t1 < te) t1 = FTEIK_HUGE
         ENDIF
      ENDIF
      t2d  = MIN(t1, t2, t3)
      tupd = MIN(t, t1d, t2d)
      RETURN
      END
!>    @brief Initializes parameters for the local solver and computes initial travel 
!>           times around the source.
!>
!>    @param[in] isrc   Source number.
!>
!>    @param[out] dest  Indices in travel time field where ts4 should be inserted.
!>    @param[out] ts4   Analytic travel times computed at grid points around source.
!>    @param[out] ierr  0 indicates success.
!>    @ingroup localSolver2d
      SUBROUTINE fteik_localSolver2d_initialize64f(isrc, epsS2C, dest, ts4, ierr)
      USE FTEIK_SOURCE64F, ONLY : fteik_source_getSolverInfo64fF
      USE FTEIK_MODEL64F, ONLY : fteik_model_grid2indexF
      USE FTEIK_MODEL64F, ONLY : nz, nzx
      USE FTEIK_MODEL64F, ONLY : dzIn => dz
      USE FTEIK_MODEL64F, ONLY : dxIn => dx
      USE FTEIK_CONSTANTS64F, ONLY : one
      INTEGER, INTENT(IN) :: isrc
      DOUBLE PRECISION, INTENT(IN) :: epsS2C
      DOUBLE PRECISION, INTENT(OUT) :: ts4(4)
      INTEGER, INTENT(OUT) :: dest(4), ierr
      DOUBLE PRECISION dz2, dx2
      DOUBLE PRECISION ysa
      INTEGER ysi
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
         WRITE(ERROR_UNIT,900)
         RETURN
      ENDIF
      epsSolver = epsS2C !epsS2C*szero*MIN(dz, dx)
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
      ts4(1) = fteik_localSolver2d_tAna64f(zsi,   xsi,   dz, dx, zsa, xsa, szero)
      ts4(2) = fteik_localSolver2d_tAna64f(zsi+1, xsi,   dz, dx, zsa, xsa, szero)
      ts4(3) = fteik_localSolver2d_tAna64f(zsi,   xsi+1, dz, dx, zsa, xsa, szero)
      ts4(4) = fteik_localSolver2d_tAna64f(zsi+1, xsi+1, dz, dx, zsa, xsa, szero)
      dest(1) = fteik_model_grid2indexF(zsi,   xsi  , 1, nz, nzx)
      dest(2) = fteik_model_grid2indexF(zsi+1, xsi  , 1, nz, nzx)
      dest(3) = fteik_model_grid2indexF(zsi  , xsi+1, 1, nz, nzx)
      dest(4) = fteik_model_grid2indexF(zsi+1, xsi+1, 1, nz, nzx)
  900 FORMAT('fteik_localSolver2d_initialize64f: Problem with source')
      RETURN
      END

END MODULE
!========================================================================================!
!> @defgroup localSolver3d 3D Local Solver
!> @brief The 3D local solver.
!> @ingroup solver3d
MODULE FTEIK_LOCALSOLVER64F
  USE ISO_FORTRAN_ENV
  USE FTEIK_CONSTANTS64F, ONLY : chunkSize
  IMPLICIT NONE
  INTEGER, PRIVATE, SAVE :: zsi, xsi, ysi
  DOUBLE PRECISION, PRIVATE, SAVE :: zsa, xsa, ysa  
  DOUBLE PRECISION, PRIVATE, SAVE :: dz, dx, dy
  DOUBLE PRECISION, PRIVATE, SAVE :: szero, szero2
  DOUBLE PRECISION, PRIVATE, SAVE :: dxi, dyi, dzi
  DOUBLE PRECISION, PRIVATE, SAVE :: dx2, dy2, dz2
  DOUBLE PRECISION, PRIVATE, SAVE :: dx2i, dy2i, dz2i 
  DOUBLE PRECISION, PRIVATE, SAVE :: dsum, dsum_nine, dsumi, epsSolver
  DOUBLE PRECISION, PRIVATE, SAVE :: dz2dx2, dz2dy2, dx2dy2
  DOUBLE PRECISION, PRIVATE, SAVE :: dz2i_dx2i, dz2i_dy2i, dx2i_dy2i
  DOUBLE PRECISION, PRIVATE, SAVE :: dz2i_p_dx2i, dz2i_p_dy2i, dx2i_p_dy2i
  DOUBLE PRECISION, PRIVATE, SAVE :: dz2i_p_dy2i_inv, dz2i_p_dx2i_inv, dx2i_p_dy2i_inv

  PUBLIC :: fteik_localSolver3D_tAna64f
  PUBLIC :: fteik_localSolver3D_tAnaD64f
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
                                                   tupd)
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, chunkSize, half, four
      INTEGER, INTENT(IN) :: n
      DOUBLE PRECISION, INTENT(IN) :: tv(n), te(n), tn(n), tev(n),  &
                                      ten(n), tnv(n), tnve(n)!, tt0(n)
      DOUBLE PRECISION, INTENT(IN) :: slow1(n), slow2(n), slow3(n), slow4(n), &
                                      slow5(n), slow6(n), slow7(n)
      DOUBLE PRECISION, INTENT(INOUT) :: tupd(n)
      ! local variables
      DOUBLE PRECISION four_sref2, sref2,                                      &
                       t1, t2, t3, t3d,                      &
                       t1_2d, t2_2d, t3_2d,                                    &
                       ta, tb, tab, tab2, tac, tac2, tbc, tbc2, tc,            &
                       temtn, temtv, tvmtn, tmax
      INTEGER i
      DOUBLE PRECISION :: t12min(chunkSize)
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
                     + DSQRT(four_sref2*dz2i_p_dx2i - dz2i_dx2i*tab2) )*dz2i_p_dx2i_inv
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
                     + DSQRT(four_sref2*dz2i_p_dy2i - dz2i_dy2i*tab2) )*dz2i_p_dy2i_inv
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
                     + DSQRT(four_sref2*dx2i_p_dy2i - dx2i_dy2i*tab2) )*dx2i_p_dy2i_inv
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
               t3d = (t1 + DSQRT(t2 - t3))*dsumi
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
      PURE DOUBLE PRECISION                                                &
      FUNCTION fteik_localSolver_init64fF(tv, te, tn, tev,                 &
                                          ten, tnv, tnve, tt0,             &
                                          slow1, slow2, slow3, slow4,      &
                                          slow5, slow6, slow7,             &
                                          i, j, k,                         &
                                          sgntz, sgntx, sgnty,             &
                                          sgnrz_dzi, sgnrx_dxi, sgnry_dyi) &
      RESULT(tupd)
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, zero
      DOUBLE PRECISION, INTENT(IN) :: tv, te, tn, tev, ten, tnv, tnve, tt0
      DOUBLE PRECISION, INTENT(IN) :: slow1, slow2, slow3, slow4, &
                                           slow5, slow6, slow7
      DOUBLE PRECISION, INTENT(IN) :: sgnrz_dzi, sgnrx_dxi, sgnry_dyi
      INTEGER, INTENT(IN) :: i, j, k, sgntx, sgnty, sgntz
      ! local variables
      DOUBLE PRECISION apoly, bpoly, cpoly, dpoly,                                &
                       four_sref2, sref2,                                         &
                       t0c, t1, t1d, t1d_t2d_min, t2, t3, t3d,                    &
                       t1_2d, t2_2d, t3_2d,                                       &
                       ta, tb, tab, tab2, tac, tac2, tbc, tbc2, tc,               &
                       taue, tauev, tauen, taun, taunv, taunve, tauv,             &
                       tmax, tmin, txc, tyc, tzc
      LOGICAL lcartSolver, ldo3D, ldoZX, ldoZXCart, ldoZXSphere,  &
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
            CALL fteik_localSolver3D_tAnaD64f(t0c, tzc, txc, tyc, i, j, k,    &
                                              dz, dx, dy, zsa, xsa, ysa, szero)
            ! Convert times into pertubations
            tauv   = tv   - fteik_localSolver3D_tAna64f(i-sgntz, j,       k,            &
                                                        dz, dx, dy, zsa, xsa, ysa, szero)
            taue   = te   - fteik_localSolver3D_tAna64f(i,       j-sgntx, k,            &
                                                        dz, dx, dy, zsa, xsa, ysa, szero)
            taun   = tn   - fteik_localSolver3D_tAna64f(i,       j,       k-sgnty,      &
                                                        dz, dx, dy, zsa, xsa, ysa, szero)
            tauev  = tev  - fteik_localSolver3D_tAna64f(i-sgntz, j-sgntx, k,            &
                                                        dz, dx, dy, zsa, xsa, ysa, szero)
            tauen  = ten  - fteik_localSolver3D_tAna64f(i,       j-sgntx, k-sgnty,      &
                                                        dz, dx, dy, zsa, xsa, ysa, szero)
            taunv  = tnv  - fteik_localSolver3D_tAna64f(i-sgntz, j,       k-sgnty,      &
                                                        dz, dx, dy, zsa, xsa, ysa, szero)
            taunve = tnve - fteik_localSolver3D_tAna64f(i-sgntz, j-sgntx, k-sgnty,      &
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
                    + DSQRT(four_sref2*dz2i_p_dx2i - dz2i_dx2i*tab2) )*dz2i_p_dx2i_inv
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
            IF (dpoly >= 0.d0) t1_2d = (DSQRT(dpoly) - bpoly)/(2.d0*apoly) + t0c
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
                     + DSQRT(four_sref2*dz2i_p_dy2i - dz2i_dy2i*tab2) )*dz2i_p_dy2i_inv
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
            IF (dpoly >= 0.d0) t2_2d = (DSQRT(dpoly) - bpoly)/(2.d0*apoly) + t0c
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
                     + DSQRT(four_sref2*dx2i_p_dy2i - dx2i_dy2i*tab2) )*dx2i_p_dy2i_inv
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
            IF (dpoly >= 0.d0) t3_2d = (DSQRT(dpoly) - bpoly)/(2.d0*apoly) + t0c
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
               t3d = (t1 + DSQRT(t2 - t3))*dsumi
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
            IF (dpoly >= 0.d0) t3d = (DSQRT(dpoly) - bpoly)/(2.d0*apoly) + t0c
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
!>    @author Mark Noble, Alexandrine Gesret, and Ben Baker.
!>    @version 2
!>    @date July 2017
!>    @copyright CeCILL-3
      PURE DOUBLE PRECISION                                    &
      FUNCTION fteik_localSolver3D_tAna64f(i, j, k, dz, dx, dy,  &
                                           zsa, xsa, ysa, szero)
      !!$OMP DECLARE SIMD(fteik_localSolver_tAna64fF) &
      !!$OMP UNIFORM(dz, dx, dy, zsa, xsa, ysa, szero)
      DOUBLE PRECISION, INTENT(IN) :: dz, dx, dy, szero, zsa, xsa, ysa
      INTEGER, INTENT(IN) :: i, j, k
      DOUBLE PRECISION d0, diffz, diffz2, diffx, diffx2, diffy, diffy2
      diffz = (DBLE(i) - zsa)*dz
      diffx = (DBLE(j) - xsa)*dx
      diffy = (DBLE(k) - ysa)*dy
      diffz2 = diffz*diffz
      diffx2 = diffx*diffx
      diffy2 = diffy*diffy
      d0 = diffz2 + diffx2 + diffy2
      fteik_localSolver3D_tAna64f = szero*DSQRT(d0) !sqrt(diffz2 + diffx2 + diffy2);
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
!>    @version 2
!>    @date July 2017
!>    @copyright CeCILL-3
      PURE SUBROUTINE fteik_localSolver3D_tAnaD64f(t_anad,        &
                                                 tzc, txc, tyc, &
                                                 i, j, k,       &
                                                 dz, dx, dy,    &
                                                 zsa, xsa, ysa, &
                                                 szero)
      !!$OMP DECLARE SIMD(fteik_localSolver_tAnaD64fF) &
      !!$OMP UNIFORM(dz, dx, dy, zsa, xsa, ysa, szero)
      DOUBLE PRECISION, INTENT(IN) :: dz, dx, dy, zsa, xsa, ysa, szero
      INTEGER, INTENT(IN) :: i, j, k
      DOUBLE PRECISION, INTENT(OUT) :: t_anad, tzc, txc, tyc
      DOUBLE PRECISION d0, d0i_szero, diffz, diffz2, diffx, diffx2, diffy, diffy2, sqrtd0
      DOUBLE PRECISION, PARAMETER :: zero = 0.d0
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
         sqrtd0 = DSQRT(d0)
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
                                                  ts8, ierr)
      USE FTEIK_SOURCE64F, ONLY : fteik_source_getSolverInfo64fF
      USE FTEIK_MODEL64F, ONLY : fteik_model_grid2indexF
      USE FTEIK_MODEL64F, ONLY : nz, nzx
      USE FTEIK_MODEL64F, ONLY : dzIn => dz
      USE FTEIK_MODEL64F, ONLY : dxIn => dx
      USE FTEIK_MODEL64F, ONLY : dyIn => dy 
      USE FTEIK_CONSTANTS64F, ONLY : one, nine
      INTEGER, INTENT(IN) :: isrc
      DOUBLE PRECISION, INTENT(IN) :: epsS2C
      DOUBLE PRECISION, INTENT(OUT) :: ts8(8)
      INTEGER, INTENT(OUT) :: dest(8), ierr
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
         WRITE(ERROR_UNIT,900)
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
      ts8(1) = fteik_localSolver3D_tAna64f(zsi,   xsi,   ysi,              &
                                           dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(2) = fteik_localSolver3D_tAna64f(zsi+1, xsi,   ysi,              &
                                           dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(3) = fteik_localSolver3D_tAna64f(zsi,   xsi+1, ysi,              &
                                           dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(4) = fteik_localSolver3D_tAna64f(zsi+1, xsi+1, ysi,              &
                                           dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(5) = fteik_localSolver3D_tAna64f(zsi,   xsi,   ysi+1,            &
                                           dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(6) = fteik_localSolver3D_tAna64f(zsi+1, xsi,   ysi+1,            &
                                           dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(7) = fteik_localSolver3D_tAna64f(zsi,   xsi+1, ysi+1,            &
                                           dz, dx, dy, zsa, xsa, ysa, szero)
      ts8(8) = fteik_localSolver3D_tAna64f(zsi+1, xsi+1, ysi+1,            &
                                           dz, dx, dy, zsa, xsa, ysa, szero)
      dest(1) = fteik_model_grid2indexF(zsi,   xsi  , ysi, nz, nzx)
      dest(2) = fteik_model_grid2indexF(zsi+1, xsi  , ysi, nz, nzx)
      dest(3) = fteik_model_grid2indexF(zsi  , xsi+1, ysi, nz, nzx)
      dest(4) = fteik_model_grid2indexF(zsi+1, xsi+1, ysi, nz, nzx)
      dest(5) = fteik_model_grid2indexF(zsi,   xsi  , ysi+1, nz, nzx)
      dest(6) = fteik_model_grid2indexF(zsi+1, xsi  , ysi+1, nz, nzx)
      dest(7) = fteik_model_grid2indexF(zsi,   xsi+1, ysi+1, nz, nzx)
      dest(8) = fteik_model_grid2indexF(zsi+1, xsi+1, ysi+1, nz, nzx)
  900 FORMAT('fteik_localSolver3D_initialize64fF: Problem with source')
      RETURN
      END

      PURE INTEGER                                                    &
      FUNCTION fteik_localSolver_grid2indexF(i, j, k, nz, nzx)        &
      RESULT(grid2indexF)
      !!$OMP DECLARE SIMD(fteik_localSolver_grid2indexF) UNIFORM(nz, nzx)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i, j, k, nz, nzx 
      grid2indexF = (k - 1)*nzx + (j - 1)*nz + i 
      RETURN
      END FUNCTION

!----------------------------------------------------------------------------------------!
!                                      End the Code                                      !
!----------------------------------------------------------------------------------------!
END MODULE
