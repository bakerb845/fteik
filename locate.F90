#define ERRORMSG(msg) WRITE(0,'("[ERROR] at ",I4," in file ",/,A,/,"Error message: ",A)') __LINE__,__FILE__,msg

MODULE FTEIK_LOCATE
  USE ISO_C_BINDING
  IMPLICIT NONE
  !> Contains the largest epochal time for each gather.  This is a vector with 
  !> dimension [nEvents].
  REAL(C_DOUBLE), PROTECTED, ALLOCATABLE, SAVE :: tepoch(:)
  !> Observation weights (1/seconds).  This is a vector with dimension [nEvents x lntf].
  REAL(C_FLOAT), PROTECTED, POINTER, SAVE :: weights(:) 
  !> Observed pick times (seconds).  The largest epochal time has been removed to avoid
  !> numerical overflow.  This is a vector of dimension [nEvents x ldgrd].
  REAL(C_FLOAT), PROTECTED, POINTER, SAVE :: observations(:)
  !> Static corrections for each observation (that is each receiver and P or S
  !> velocity pair). This is a vector of dimension [ntf].
  REAL(C_FLOAT), PROTECTED, POINTER, SAVE :: staticCorrection(:)
  !> Holds the travel time fields (seconds).  This is a vector with 
  !> dimension [ntf x ldgrd].
  REAL(C_FLOAT), PROTECTED, POINTER, SAVE :: ttFields(:) 
  !> Flag indicating whether or not to skip this event/observation pair.  This is a vector
  !> of dimension [nEvents x lntf]
  LOGICAL(C_BOOL), PROTECTED, ALLOCATABLE, SAVE :: lskip(:)
  !> Number of events.
  INTEGER(C_INT), PROTECTED, SAVE :: nEvents = 0
  !> Number of travel time fields.
  INTEGER(C_INT), PROTECTED, SAVE :: ntf = 0
  !> Number of grid points.
  INTEGER(C_INT), PROTECTED, SAVE :: nGrd = 0
  !> Padded version of ntf which serves as a leading dimension.
  INTEGER(C_INT), PRIVATE, SAVE :: lntf = 0
  !> Padded version of nGrd which serves as a leading dimension.
  INTEGER(C_INT), PRIVATE, SAVE :: ldgrd = 0 
  !> Flag indicating whether or not locator has been initialized.
  LOGICAL(C_BOOL), PRIVATE, SAVE :: linit = .FALSE.
  !> Memory alignment.
  INTEGER(C_INT), PRIVATE, PARAMETER :: alignment = 64

  CONTAINS

      SUBROUTINE locate_locateF( ) BIND(C, NAME='locate_locateF')
      USE ISO_C_BINDING
     
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Initializes space on the locator.
!>
!>    @param[in] nEventsIn      Number of events to locate.
!>    @param[in] ntfIn          Number of travel-time fields that will be set.
!>    @param[in] ngrdIn         Number of grid-points in each travel-time field.
!>
!>    @param[out] ierr          0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_initializeF(nEventsIn, ntfIn, ngrdIn, ierr) &
      BIND(C, NAME='locate_initializeF')
      USE FTEIK_MEMORY
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: nEventsIn, ntfIn, ngrdIn
      INTEGER(C_INT), INTENT(OUT) :: ierr

      CALL locate_finalizeF()
      CALL locate_setNumberOfEventsF(nEventsIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,900)
 900     FORMAT('locate_initializeF: Failed to set number of events')
         RETURN
      ENDIF
      CALL locate_setNumberOfTravelTimeFieldsF(ntfIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,905)
 905     FORMAT('locate_initializeF: Failed to set number of travel time fields')
         RETURN
      ENDIF
      CALL locate_setNumberOfGridPointsF(ngrdIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,"('locate_initializeF: Failed to set number of grid points',A)")
         RETURN
      ENDIF
      ALLOCATE(tepoch(nEvents))
      lntf  = ntf  + padLength32F(alignment, ntf) 
      IF (MOD(4*lntf, alignment) /= 0) THEN
         WRITE(*,"('locate_initializeF: Failed to pad ntf',A)")
         ierr = 1
         RETURN
      ENDIF
      ALLOCATE(lskip(nEvents*lntf))
      ldgrd = ngrd + padLength32F(alignment, ngrd)
      IF (MOD(4*ldgrd, alignment) /= 0) THEN
         WRITE(*,"('locate_initializeF: Failed to pad ngrd',A)")
         ierr = 1
         RETURN
      ENDIF
      CALL allocate32f(weights,          alignment, nEvents*lntf)
      CALL allocate32f(observations,     alignment, nEvents*lntf)
      CALL allocate32f(ttFields,         alignment, ntf*ldgrd)
      CALL allocate32f(staticCorrection, alignment, ntf)
      tepoch(:) = 0.d0
      lskip(:) = .TRUE.
      weights(:) = 0.0
      observations(:) = 0.0
      staticCorrection(:) = 0.0
      ttFields(:) = 0.0
      linit = .TRUE.
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases memory on the locator and sets all variables to 0 or NULL.
!> 
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_finalizeF( ) &
      BIND(C, NAME='locate_finalizeF')
      USE FTEIK_MEMORY
      USE ISO_C_BINDING
      IF (ALLOCATED(tepoch)) DEALLOCATE(tepoch)
      IF (ALLOCATED(lskip))  DEALLOCATE(lskip)
      IF (linit) THEN
         CALL free32f(observations)
         CALL free32f(weights)
         CALL free32f(ttFields)
         CALL free32f(staticCorrection)
      ENDIF
      lntf = 0
      ldgrd = 0
      nGrd = 0
      nEvents = 0
      linit = .FALSE.
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the double precision travel time field on the solver.
!>
!>    @param[in] ngrdIn  Number of grid points in travel time field.  This should
!>                       equal ngrd.
!>    @param[in] itf     Travel time field number (Fortran indexed).
!>    @param[in] ttIn    Travel-time field (seconds) from receiver to all points
!>                       in the medium.  This is a vector of dimension [ngrdIn]. 
!>
!>    @param[out] ierr   0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_setTravelTimeField64fF(ngrdIn, itf, ttIn, ierr) &
      BIND(C, NAME='locate_setTravelTimeField64fF')
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: ngrdIn, itf
      REAL(C_DOUBLE), TARGET, INTENT(IN) :: ttIn(ngrdIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
#if defined(FTEIK_FORTRAN_USE_INTEL)
      INTERFACE
         INTEGER(C_INT) FUNCTION ippsConvert_64f32f(pSrc, pDst, len) &
         BIND(C, NAME='ippsConvert_64f32f')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT), VALUE, INTENT(IN) :: len
         REAL(C_DOUBLE), INTENT(IN) :: pSrc(len)
         REAL(C_FLOAT), INTENT(OUT) :: pDst(len)
         END FUNCTION
      END INTERFACE
#else
      REAL(C_FLOAT), POINTER :: ttptr(:)
      INTEGER(C_INT) igrd
#endif 
      INTEGER(C_INT) i1, i2
      ierr = 0
      IF (.NOT.linit) THEN
         WRITE(*,"('locate_setTravelTimeField64fF: Locator not initialized',A)")
         ierr = 1
         RETURN
      ENDIF
      IF (itf < 0 .OR. itf > ntf) THEN
         WRITE(*,900) itf
 900     FORMAT('locate_setTravelTimeField64fF: Field number=',I4,' must be in [1,',I3']')
         ierr = 1
         RETURN
      ENDIF
      IF (ngrdIn /= ngrd) THEN
         WRITE(*,905) ngrdIn, ngrd
 905     FORMAT('locate_setTravelTimeField64fF: ngrdIn=',I8,'not equal ngrd=', I8)
         ierr = 1
         RETURN
      ENDIF
      i1 = (itf - 1)*ldgrd + 1
      i2 = itf*ldgrd
#if defined(FTEIK_FORTRAN_USE_INTEL)
      ierr = ippsConvert_64f32f(ttIn, ttFields(i1:i2), ngrd)
      IF (ierr /= 0) THEN
         WRITE(*,910) itf
  910    FORMAT('locate_setTravelTimeField64fF: Failed to travel time field=', I4)
      ENDIF
#else
      ttptr => ttFields(i1:i2)
      !$OMP SIMD ALIGNED(ttptr: 64)
      DO 1 igrd=1,ngrd
         ttptr(igrd) = REAL(ttIn(igrd))
   1  CONTINUE
      NULLIFY(ttptr)
#endif
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the single precision travel time field on the solver.
!>
!>    @param[in] ngrdIn  Number of grid points in travel time field.  This should
!>                       equal ngrd.
!>    @param[in] itf     Travel time field number (Fortran indexed).
!>    @param[in] ttIn    Travel-time field (seconds) from receiver to all points
!>                       in the medium.  This is a vector of dimension [ngrdIn]. 
!>
!>    @param[out] ierr   0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_setTravelTimeField32fF(ngrdIn, itf, ttIn, ierr) &
      BIND(C, NAME='locate_setTravelTimeField32fF')
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: ngrdIn, itf
      REAL(C_FLOAT), TARGET, INTENT(IN) :: ttIn(ngrdIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_FLOAT), POINTER :: ttptr(:)
      INTEGER(C_INT) i1, i2, igrd
      ierr = 0 
      IF (.NOT.linit) THEN
         WRITE(*,"('locate_setTravelTimeField32fF: Locator not initialized',A)")
         ierr = 1 
         RETURN
      ENDIF
      IF (itf < 0 .OR. itf > ntf) THEN
         WRITE(*,900) itf 
 900     FORMAT('locate_setTravelTimeField32fF: Field number=',I4,' must be in [1,',I3']')
         ierr = 1 
         RETURN
      ENDIF
      IF (ngrdIn /= ngrd) THEN
         WRITE(*,905) ngrdIn, ngrd
 905     FORMAT('locate_setTravelTimeField32fF: ngrdIn=',I8,'not equal ngrd=', I8) 
         ierr = 1 
         RETURN
      ENDIF
      i1 = (itf - 1)*ldgrd + 1 
      i2 = itf*ldgrd
      ttptr => ttFields(i1:i2)
      !$OMP SIMD ALIGNED(ttptr: 64)
      DO 1 igrd=1,ngrd
         ttptr(igrd) = ttIn(igrd)
   1  CONTINUE
      NULLIFY(ttptr)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE locate_setObservations64fF(nEventsIn )
      INTEGER(C_INT), INTENT(IN) :: nEventsIn
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the number of events to locate.
!>
!>    @param[in] nEventsIn    Number of events to locate.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_setNumberOfEventsF(nEventsIn, ierr) &
      BIND(C, NAME='locate_setNumberOfEventsF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nEventsIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      nEvents = 0 
      ierr = 0
      IF (nEventsIn < 1) THEN
         WRITE(*, 900) nEventsIn
  900    FORMAT('locate_setNumberOfEventsf: nEvents=',I4, 'must be positive')
         ierr = 1 
         RETURN
      ENDIF
      nEvents = nEventsIn
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the number of travel time fields.
!>
!>    @param[in] nobs   Number of travel time fields. 
!>
!>    @param[out] ierr  0 indicate success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>                        
      SUBROUTINE locate_setNumberOfTravelTimeFieldsF(ntfIn, ierr) &
      BIND(C, NAME='locate_setNumberOfTravelTimeFieldsF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: ntfIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ntf = 0
      ierr = 0 
      IF (ntfIn < 1) THEN
         WRITE(*,900) ntfIn
  900    FORMAT('locate_setNumberOfTravelTimeFieldsF: ntf=', I4,' must be positive')
         ierr = 1 
         RETURN
      ENDIF
      ntf = ntfIn
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the number of grid points.
!>
!>    @param[in] ngrdIn     Number of grid points in travel-time field.
!>
!>    @param[out] ierr      0 indicates success.
!>
!>    @result 0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_setNumberOfGridPointsF(ngrdIn, ierr) &
      BIND(C, NAME='locate_setNumberOfGridPointsF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: ngrdIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ngrd = 0 
      ierr = 0 
      IF (ngrdIn < 1) THEN
         WRITE(*,900) ngrdIn
  900    FORMAT('locate_setNumberOfGridPointsF: Number of grid points=', &
                I8, ' must be positive')
         ierr = 1 
         RETURN
      ENDIF
      ngrd = ngrdIn
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the origin time for a least-squares inversion at each grid
!>           point in the model.  Following Moser and and Van Eck A5 the least-squares
!>           origin time for a diagonal weight matrix is computed by (A5):
!>           \f[
!>               t_0(\textbf{x}_0)
!>             = \sum_{i}^{n_{obs}}
!>               \frac{ \sum_i (\tau_i - (T(\textbf{x}_0; \textbf{x}_i) + t_i^s))}
!>                    { \sum_{i}^{n_{obs}} w_i }
!>           \f]
!>           where \f$ t_0(\textbf{x}_0) \f$ is the origin time,
!>           \f$ \tau_i \f$ is the travel time for the i'th observation, 
!>           \f$ T(\textbf{x}_0; \textbf{x}_i) \f$ are the traveltimes from the i'th
!>           observation to all points in the model, and \f$ t_i^2 \f$ is the static
!>           correction for the i'th observation.
!>           This minimizes the travel time L2 objective function:
!>           \f[
!>              \mathcal{C}(\textbf{x}_0)
!>             =\frac{1}{2} \sum_i^{n_{obs}}
!>               w_i \left (\tau_i - (T(\textbf{x}_0, \textbf{x}_i)
!>                         + t_i^s + t_0) \right )^2
!>           \f].
!>
!>    @param[in] evnmbr   Event number.
!>    @param[in] nobs     Number of observations.
!>    @param[in] tobs     Observations (seconds).  This is a vector of dimension [nobs].
!>    @param[in] tstatic  Static corrections (s) for the observations.  This is a
!>                        vector of dimension [nobs].
!>    @param[in] wts      Observation weight (1/second).  This is a vector of
!>                        dimension [nobs].
!>
!>    @param[out] t0      Origin time (seconds).  
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_computeL2OriginTime32fF(evnmbr, nobs, tstatic, &
                                                t0)             &
      BIND(C, NAME='locate_computeL2OriginTime32fF')
      USE FTEIK_CONSTANTS32F, ONLY : one, zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: evnmbr, nobs
      REAL(C_FLOAT), INTENT(IN) :: tstatic(nobs)
      REAL(C_FLOAT), INTENT(OUT) :: t0(ngrd)
      REAL(C_FLOAT) est, res, obs, toff, xnorm, wt
      REAL(C_FLOAT), POINTER :: ttptr(:)
      REAL(C_FLOAT), POINTER :: wtptr(:)
      REAL(C_FLOAT), POINTER :: obsptr(:)
      INTEGER(C_INT) i, igrd, k1, k2
      ! Initialize t0 
      !$OMP SIMD aligned(t0:64)
      DO 1 igrd=1,ngrd
         t0(igrd) = zero 
    1 CONTINUE
      k1 = (evnmbr - 1)*lntf + 1
      k2 = evnmbr*lntf
      obsptr = observations(k1:k2)
      wtptr = weights(k1:k2)
      ! Compute the sum of the weights - this is the normalization in Eqn A5
      xnorm = one/SUM(wtptr)
      ! Loop on the observations and stack in the residual
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP PRIVATE(est, i, k1, k2, obs, res, toff, ttptr, wt) &
      !$OMP SHARED(ldgrd, lskip, ngrd, nobs, obsptr, tstatic, ttFields, wtptr, xnorm) &
      !$OMP REDUCTION(+:t0)
      DO 2 i=1,nobs
         IF (lskip(i)) CYCLE
         wt = wtptr(i)
         toff = tstatic(i)
         obs = obsptr(i)
         k1 = (ldgrd - 1)*i + 1
         k2 = ldgrd*i
         ttptr => ttFields(k1:k2)
         !$OMP SIMD ALIGNED(t0, ttptr: 64)
         DO 3 igrd=1,ngrd
            est = ttptr(igrd) + toff
            res = xnorm*(obs - est)
            t0(igrd) = t0(igrd) + wt*res 
    3    CONTINUE
         NULLIFY(ttptr)
    2 CONTINUE 
      !$OMP END PARALLEL DO
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the objective function for a least-squres inversion at
!>           each grid point.  This minimizes the travel time L2 objective function:
!>           \f[
!>              \mathcal{C}(\textbf{x}_0)
!>             =\frac{1}{2} \sum_i^{n_{obs}}
!>               w_i \left (\tau_i - (T(\textbf{x}_0, \textbf{x}_i)
!>                         + t_i^s + t_0(\textbf{x}_0)) \right )^2
!>           \f]
!>           \f$ T(\textbf{x}_0; \textbf{x}_i) \f$ are the traveltimes from the i'th
!>           observation to all points in the model, \f$ t_0(\textbf{x}_0) \f$
!>           is the origin time at all points in the model, and \f$ t_i^2 \f$ is the 
!>           static correction for the i'th observation.
!>
!>    @param[in] ngrd     Number of grid points.
!>    @param[in] nobs     Number of observations.
!>    @param[in] tobs     Observations (seconds).  This is a vector of dimension [nobs].
!>    @param[in] tstatic  Static corrections (s) for the observations.
!>    @param[in] wts      Observation weight (1/second).  This is a vector of
!>                        dimension [nobs].
!>    @param[in] tt       Traveltimes from receiver to all grid-points (seconds).
!>                        This is a vector of dimension [ldg x nobs] where ldg >= ngrd
!>                        and mod(ldg, 64) must equal 0.
!>    @param[in] t0       Origin time (seconds).  This is a vector of dimension [nobs].
!>
!>    @param[out] tobj    Objective function at all grid points.  This is a vector of
!>                        dimension [ngrd].
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_computeObjectiveFunction32fF(ldg, ngrd, nobs, &
                                                     tobs, tstatic,   &
                                                     wts, tt, t0,     &
                                                     tobj)            &
      BIND(C, NAME='locate_computeObjectiveFunction32fF')
      USE FTEIK_CONSTANTS32F, ONLY : sqrtHalf, zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldg, ngrd, nobs
      REAL(C_FLOAT), INTENT(IN) :: tobs(nobs), tstatic(nobs), wts(nobs)
      REAL(C_FLOAT), INTENT(IN) :: t0(ngrd)
      REAL(C_FLOAT), TARGET, INTENT(IN) :: tt(ngrd)
      REAL(C_FLOAT), INTENT(OUT) :: tobj(ngrd)
      REAL(C_FLOAT) est, res, obs, toff, wt
      REAL(C_FLOAT), POINTER :: ttptr(:)
      INTEGER(C_INT) i, igrd, k1, k2
      !$OMP SIMD aligned(tobj: 64)
      DO 1 igrd=1,ngrd
         tobj(igrd) = zero
    1 CONTINUE 
      ! Now compute the travel times
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP PRIVATE(i, igrd, est, k1, k2, obs, res, toff, ttptr, wt) &
      !$OMP SHARED(ldg, ngrd, nobs, t0, tobs, tstatic, tt, wts) &
      !$OMP REDUCTION(+:tobj)
      DO 2 i=1,nobs
         wt = sqrtHalf*wts(i)
         toff = tstatic(i)
         obs = tobs(i)
         k1 = (ldg - 1)*i + 1 
         k2 = ldg*i
         ttptr => tt(k1:k2) 
         !$OMP SIMD ALIGNED(t0, ttptr, tobj: 64)
         DO 3 igrd=1,ngrd
            est = ttptr(igrd) + t0(igrd) + toff ! Eqn A1 + static correction
            res = wt*(obs - est)
            tobj(igrd) = tobj(igrd) + res*res
    3    CONTINUE
    2 CONTINUE 
      !$OMP END PARALLEL DO
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE locate_computeStatics32fF(ldo, nobs, nev, tobs, tt, &
                                           t0, wts, tstatic)    &
      BIND(C, NAME='locate_computeStatics32fF')
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: ldo, nev, nobs
      REAL(C_FLOAT), INTENT(IN) :: t0(nev)
      REAL(C_FLOAT), INTENT(IN), TARGET :: tobs(nev*ldo), wts(nev*ldo), &
                                           tt(nev*ldo)
      REAL(C_FLOAT), INTENT(OUT) :: tstatic(nobs)
      REAL(C_FLOAT), POINTER :: obsptr(:), ttptr(:), wtsptr(:)
      REAL(C_FLOAT) res, test, tori, ts
      INTEGER(C_INT) i, iv, k1, k2
      !$OMP SIMD aligned(tstatic: 64)
      DO 1 i=1,nobs
         tstatic(i) = zero
    1 CONTINUE 
      DO iv=1,nev
         tori = t0(iv)
         k1 = (iv - 1)*ldo + 1
         k2 = iv*ldo
         obsptr => tobs(k1:k2)
         ttptr  => tt(k1:k2)
         wtsptr => wts(k1:k2) 
         ts = zero
         !$OMP SIMD ALIGNED(tstatic, obsptr, ttptr, wtsptr: 64)
         DO i=1,nobs
            test = ttptr(i) + tori
            res = wtsptr(i)*(obsptr(i) - test)
            tstatic(i) = tstatic(i) + res
         ENDDO
      ENDDO
      RETURN
      END
END MODULE
