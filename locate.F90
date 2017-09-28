#define ERRORMSG(msg) WRITE(0,'("[ERROR] at ",I4," in file ",/,A,/,"Error message: ",A)') __LINE__,__FILE__,msg

MODULE FTEIK_LOCATE
  USE ISO_C_BINDING
  IMPLICIT NONE
  INTEGER(C_INT), PROTECTED, SAVE :: nEvents = 0

  CONTAINS

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
      SUBROUTINE fteik_locate_setNumberOfEventsF(nEventsIn, ierr) &
      BIND(C, NAME='fteik_locate_setNumberOfEventsF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nEventsIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      nEvents = 0 
      ierr = 0
      IF (nEvents < 1) THEN
         ERRORMSG("nEvents must be positive")
         ierr = 1 
         RETURN
      ENDIF
      nEvents = nEventsIn
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the origin time for a least-squares inversion at each grid
!>           point in the model.  Following Moser and and Van Eck A5 the least-squares
!>           origin time for a diagonal weight matrix is computed by (A5):
!>           \f[
!>               t_0(\textbf{x})
!>             =-\sum_{i}^{n_{obs}} \frac{ \sum_i (\tau_i - (T(\textbf{x}_0; \textbf{x}_i) + t_i^s))}
!>                                      { \sum_{i}^{n_{obs}} w_i }
!>           \f]
!>           This minimizes the travel time L2 objective function:
!>           \f[
!>              \mathcal{C}(\textbf{x}_0)
!>             =\frac{1}{2} \sum_i^{n_{obs}}
!>               w_i \left (\tau_i - (T(\textbf{x}_0, \textbf{x}_i)
!>                    + t_i^s + t_0) \right )^2
!>           \f].
!>
!>    @param[in] ngrd     Number of grid points.
!>    @param[in] nobs     Number of observations.
!>    @param[in] tobs     Observations (seconds).
!>    @param[in] tstatic  Static corrections (s) for the observations.
!>
      PURE SUBROUTINE locate_computeL2OriginTime32fF(ngrd, nobs,    &
                                                     tobs, tstatic, &
                                                     wts, tt,       &
                                                     t0, tpdf)      &
      BIND(C, NAME='locate_computeL2OriginTime32fF')
      USE FTEIK_CONSTANTS32F, ONLY : one, sqrtHalf, zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nobs, ngrd
      REAL(C_FLOAT), INTENT(IN) :: tobs(nobs), tstatic(nobs), wts(nobs)
      REAL(C_FLOAT), INTENT(IN) :: tt(ngrd)
      REAL(C_FLOAT), INTENT(OUT) :: t0(ngrd), tpdf(ngrd)
      REAL(C_FLOAT) est, res, obs, toff, xnorm, wt
      INTEGER(C_INT) i, igrd
      ! Compute the sum of the weights - this is the normalization in Eqn A5
      xnorm = one/SUM(wts)
      ! Initialize result
      t0(:) = zero 
      ! Loop on the observations and stack in the residual
      DO 1 i=1,nobs
         wt = wts(i)
         toff = tstatic(i)
         obs = tobs(i)
         !!$OMP SIMD
         DO 2 igrd=1,ngrd
            est = tt(igrd) + toff
            res = xnorm*(obs - est)
            t0(igrd) = t0(igrd) - wt*res 
    2    CONTINUE
    1 CONTINUE 
      ! Now compute the travel times
      DO 11 i=1,nobs
         wt = sqrtHalf*wts(i)
         toff = tstatic(i)
         obs = tobs(i)
         DO 12 igrd=1,ngrd
            est = tt(igrd) + t0(igrd) + toff ! Eqn A1 + static correction
            res = wt*(obs - est)
            tpdf(igrd) = tpdf(igrd) + res*res
   12    CONTINUE
   11 CONTINUE 
      RETURN
      END
END MODULE
