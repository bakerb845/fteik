MODULE FTEIK_LOCATE
  USE ISO_FORTRAN_ENV
  USE ISO_C_BINDING
#if defined(FTEIK_FORTRAN_USE_MPI)
  USE MPI_F08
#endif
  IMPLICIT NONE
  !> Contains the largest epochal time for each gather.  This is a vector with 
  !> dimension [nEvents].
  REAL(C_DOUBLE), PROTECTED, ALLOCATABLE, SAVE :: tepoch(:)
  !> Contains the reduced origin times if known a priori.  This a a vector of
  !> dimension [nEvents].
  REAL(C_FLOAT), PROTECTED, ALLOCATABLE, SAVE :: originTimes(:)
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
  !> Workspace array holding the objective function.  This is a vector of
  !> dimension [ngrd].
  REAL(C_FLOAT), PROTECTED, POINTER, SAVE :: objWork(:)
  !> Workspace array holding the origin time.  This is a vector of 
  !> dimension [MAX(blockSize, ngrd)].
  REAL(C_FLOAT), PROTECTED, POINTER, SAVE :: t0Work(:)
  !> Flag indicating whether or not to skip this event/observation pair.  This is a vector
  !> of dimension [nEvents x lntf]
  LOGICAL(C_BOOL), PROTECTED, ALLOCATABLE, TARGET, SAVE :: lskip(:)
  !> Flag indicating whether or not the origin time is known a priori.  This is a
  !> vector of dimension [nEvents].
  LOGICAL(C_BOOl), PROTECTED, ALLOCATABLE, SAVE :: lhaveOT(:)
  !> Number of events.
  INTEGER(C_INT), PROTECTED, SAVE :: nEvents = 0
  !> Number of travel time fields.
  INTEGER(C_INT), PROTECTED, SAVE :: ntf = 0
  !> Number of grid points.
  INTEGER(C_INT), PROTECTED, SAVE :: ngrd = 0
  !> Padded version of ntf which serves as a leading dimension.
  INTEGER(C_INT), PRIVATE, SAVE :: lntf = 0
  !> Padded version of ngrd which serves as a leading dimension.
  INTEGER(C_INT), PRIVATE, SAVE :: ldgrd = 0 
  !> Flag indicating whether or not locator has been initialized.
  LOGICAL(C_BOOL), PRIVATE, SAVE :: linit = .FALSE.
  !> Flag indicating whether or not to solve for the origin time.
  LOGICAL(C_BOOL), PRIVATE, SAVE :: lsolveOT = .TRUE.
  !> Flag indicating whether or not to save the origin time field.
  LOGICAL(C_BOOL), PRIVATE, SAVE :: lwantT0Field = .TRUE. 
  !> Memory alignment.
  INTEGER(C_SIZE_T), PRIVATE, PARAMETER :: alignment = 64
  !> Block size in grid search.
  INTEGER(C_INT), PRIVATE, PARAMETER :: blockSize = 512


  !> Total number of grid poitns in mesh
  INTEGER(C_INT), PRIVATE, SAVE :: ngrdAll
  !> Send counts in MPI scatterv.  This is a vector of dimension [nprocs].
  INTEGER(C_INT), PRIVATE, ALLOCATABLE, SAVE :: sendCounts(:)
  !> The displacements in MPI scatterv.  This is a vector of dimension [nprocs].
  INTEGER(C_INT), PRIVATE, ALLOCATABLE, SAVE :: displs(:)
  !> If true then the process is included in the locator. 
  LOGICAL(C_BOOL), PRIVATE, SAVE :: linLocator = .FALSE.
  !> If true then use the MPI locator.
  LOGICAL(C_BOOL), PRIVATE, SAVE :: luseMPI = .FALSE.
#if defined(FTEIK_FORTRAN_USE_MPI)
  !> Handle to the MPI communicator
  TYPE(MPI_Comm), PRIVATE, SAVE :: locatorComm !INTEGER, PRIVATE, SAVE :: locatorComm =-1
  !> Myid on solver.
  INTEGER, PRIVATE, SAVE :: mylocatorID = 0
  !> Number of processes.
  INTEGER, PRIVATE, SAVE :: nprocs = 1
  !> MPI error code.
  INTEGER, PRIVATE :: mpierr
#endif
  PUBLIC :: locate_initialize
  PUBLIC :: locate_finalize
  PUBLIC :: locate_setTravelTimeField64f
  PUBLIC :: locate_setTravelTimeField32f
  PUBLIC :: locate_setObservation64f
  PUBLIC :: locate_setNumberOfEvents
  PUBLIC :: locate_locateEvent

  PRIVATE :: locate_computeL2ObjectiveFunction32f
  CONTAINS

!     SUBROUTINE locate_locateF( ) BIND(C, NAME='locate_locateF')
!     USE ISO_C_BINDING
!     ! Loop on the events
!     DO 1 iev=1,nEevents
!   1 CONTINUE 
!     RETURN
!     END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Locates an event.
!>
!>    @param[in] evnmbr    Event number.
!>    @param[out] optIndx  Index in grid corresponding to the optimum.
!>    @param[out] t0Opt    The origin time (epochal seconds; UTC) at the optimal
!>                         grid point.
!>    @param[out] objOpt   Value of objective function at optimal grid point.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_locateEvent(evnmbr, optIndx, t0Opt, objOpt) &
      BIND(C, NAME='locate_locateEvent')
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: evnmbr
      REAL(C_DOUBLE), INTENT(OUT) :: t0Opt, objOpt
      INTEGER(C_INT), INTENT(OUT) :: optIndx
      optindx = 0
      ! Compute the objective function and potentially the origin time
      CALL locate_computeL2ObjectiveFunction32f(evnmbr)
      ! Find the optimal index (smallest value in \f$ L_2 \f$ objective function)
      optindx = MINLOC(objWork, 1)
      t0Opt  = DBLE(t0Work(optindx)) + tepoch(evnmbr)
      objOpt = DBLE(objWork(optindx))
      RETURN
      END
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
      SUBROUTINE locate_initialize(nEventsIn, ntfIn, ngrdIn, ierr) &
      BIND(C, NAME='locate_initialize')
      USE FTEIK_CONSTANTS64F, ONLY : FALSE, TRUE
      USE FTEIK_CONSTANTS32F, ONLY : zero
      USE FTEIK_MEMORY
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: nEventsIn, ntfIn, ngrdIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT64_T) nwork
      INTEGER(C_SIZE_T) nalloc 

      CALL locate_finalize()
      CALL locate_setNumberOfEvents(nEventsIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
 900     FORMAT('locate_initializeF: Failed to set number of events')
         RETURN
      ENDIF
      CALL locate_setNumberOfTravelTimeFields(ntfIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,905)
 905     FORMAT('locate_initializeF: Failed to set number of travel time fields')
         RETURN
      ENDIF
      CALL locate_setNumberOfGridPointsF(ngrdIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,"('locate_initializeF: Failed to set number of grid points',A)")
         RETURN
      ENDIF
      ALLOCATE(tepoch(nEvents))
      ALLOCATE(originTimes(nEvents))
      ALLOCATE(lhaveOT(nEvents))
      lntf  = ntf  + padLength32F(alignment, ntf) 
      IF (MOD(4*lntf, alignment) /= 0) THEN
         WRITE(ERROR_UNIT,"('locate_initializeF: Failed to pad ntf',A)")
         ierr = 1
         RETURN
      ENDIF
      ALLOCATE(lskip(nEvents*lntf))
      ldgrd = ngrd + padLength32F(alignment, ngrd)
      IF (MOD(4*ldgrd, alignment) /= 0) THEN
         WRITE(ERROR_UNIT,"('locate_initializeF: Failed to pad ngrd',A)")
         ierr = 1
         RETURN
      ENDIF
      nwork = ntf
      nwork = nwork*ldgrd
      IF (nwork > HUGE(1)) THEN
         WRITE(*,"('locate_initializeF: WARNING: Beware of integer overflow in MPI',A)")
      ENDIF
      nalloc = nEvents*lntf
      CALL allocate32f(weights,          alignment, nalloc)
      CALL allocate32f(observations,     alignment, nalloc)
      nalloc = ntf 
      nalloc = nalloc*ldgrd
      CALL allocate32f(ttFields,         alignment, nalloc)
      nalloc = ntf
      CALL allocate32f(staticCorrection, alignment, nalloc)
      nalloc = MAX(ngrd, blockSize) 
      CALL allocate32f(t0Work,           alignment, nalloc)
      ngrd = ngrd
      CALL allocate32f(objWork,          alignment, nalloc)
      tepoch(:) = 0.d0
      originTimes(:) = zero
      lhaveOT(:) = FALSE
      lskip(:) = TRUE
      weights(:) = zero 
      observations(:) = zero 
      staticCorrection(:) = zero 
      ttFields(:) = zero
      t0Work(:) = zero
      objWork(:) = zero
      linit = TRUE
      IF (.NOT.luseMPI) THEN
         ALLOCATE(sendCounts(1))
         ALLOCATE(displs(1))
         sendCounts(1) = ngrd 
         displs(1) = 0
         ngrdAll = ngrd
      ENDIF 
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
      SUBROUTINE locate_finalize( ) &
      BIND(C, NAME='locate_finalize')
      USE FTEIK_MEMORY
      USE ISO_C_BINDING
      IF (ALLOCATED(tepoch))      DEALLOCATE(tepoch)
      IF (ALLOCATED(originTimes)) DEALLOCATE(originTimes)
      IF (ALLOCATED(lskip))       DEALLOCATE(lskip)
      IF (ALLOCATED(lhaveOT))     DEALLOCATE(lhaveOT)
      IF (ALLOCATED(sendCounts))  DEALLOCATE(sendCounts)
      IF (ALLOCATED(displs))      DEALLOCATE(displs)
      IF (linit) THEN
         CALL free32f(observations)
         CALL free32f(weights)
         CALL free32f(ttFields)
         CALL free32f(staticCorrection)
         CALL free32f(t0Work)
         CALL free32f(objWork)
      ENDIF
      lntf = 0
      ldgrd = 0
      ngrd = 0
      ngrdAll = 0
      nEvents = 0
      linit = .FALSE.
#if defined(FTEIK_FORTRAN_USE_MPI)
      IF (linLocator) THEN
         CALL MPI_Comm_free(locatorComm, mpierr) 
         locatorComm = MPI_COMM_NULL !-1
         mylocatorID = 0
         nprocs = 0
      ENDIF
#endif
      luseMPI = .FALSE.
      linLocator = .FALSE.
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
      SUBROUTINE locate_setTravelTimeField64f(ngrdIn, itf, ttIn, ierr) &
      BIND(C, NAME='locate_setTravelTimeField64f')
#if defined(FTEIK_FORTRAN_USE_INTEL)
      USE IPPS_MODULE, ONLY : ippsConvert_64f32f
#endif
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: ngrdIn, itf
      REAL(C_DOUBLE), INTENT(IN) :: ttIn(ngrdIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_FLOAT), POINTER :: ttptr(:)
#if defined(FTEIK_FORTRAN_USE_INTEL)

#else
      INTEGER(C_INT) igrd
#endif 
      INTEGER(C_INT) i1, i2
      ierr = 0
      IF (.NOT.linit .AND. luseMPI) RETURN ! Quiet return
      IF (.NOT.linit) THEN
         WRITE(*,"('locate_setTravelTimeField64fF: Locator not initialized',A)")
         ierr = 1
         RETURN
      ENDIF
      IF (itf < 0 .OR. itf > ntf) THEN
         WRITE(*,900) itf
 900     FORMAT('locate_setTravelTimeField64f: Field number=',I4,' must be in [1,',I3']')
         ierr = 1
         RETURN
      ENDIF
      IF (ngrdIn /= ngrd) THEN
         WRITE(*,905) ngrdIn, ngrd
 905     FORMAT('locate_setTravelTimeField64f: ngrdIn=',I8,' not equal ngrd=', I8)
         ierr = 1
         RETURN
      ENDIF
      i1 = (itf - 1)*ldgrd + 1
      i2 = itf*ldgrd
      ttptr => ttFields(i1:i2)
      IF (.NOT.ASSOCIATED(ttptr)) THEN
         WRITE(*,"('locate_setTravelTimeField64fF: Failed to associate pointer',A)")
         ierr = 1
      ENDIF
#if defined(FTEIK_FORTRAN_USE_INTEL)
      ierr = ippsConvert_64f32f(ttIn, ttptr, ngrd) !ttFields(i1:i2), ngrd)
      IF (ierr /= 0) THEN
         WRITE(*,910) itf
  910    FORMAT('locate_setTravelTimeField64f: Failed to convert travel time field=', I4)
      ENDIF
#else
      !$OMP SIMD ALIGNED(ttptr: 64)
      DO 1 igrd=1,ngrd
         ttptr(igrd) = REAL(ttIn(igrd))
   1  CONTINUE
#endif
      NULLIFY(ttptr)
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
      SUBROUTINE locate_setTravelTimeField32f(ngrdIn, itf, ttIn, ierr) &
      BIND(C, NAME='locate_setTravelTimeField32f')
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: ngrdIn, itf
      REAL(C_FLOAT), INTENT(IN) :: ttIn(ngrdIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_FLOAT), POINTER :: ttptr(:)
      INTEGER(C_INT) i1, i2, igrd
      ierr = 0 
      IF (.NOT.linit) THEN
         WRITE(*,"('locate_setTravelTimeField32f: Locator not initialized',A)")
         ierr = 1 
         RETURN
      ENDIF
      IF (itf < 0 .OR. itf > ntf) THEN
         WRITE(*,900) itf 
 900     FORMAT('locate_setTravelTimeField32f: Field number=',I4,' must be in [1,',I3']')
         ierr = 1 
         RETURN
      ENDIF
      IF (ngrdIn /= ngrd) THEN
         WRITE(*,905) ngrdIn, ngrd
 905     FORMAT('locate_setTravelTimeField32f: ngrdIn=',I8,'not equal ngrd=', I8) 
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
!>    @brief This is a convenience routine for setting all observations in a catalog.
!>
!>
      SUBROUTINE locate_setObservations64f(nEventsIn )
      INTEGER(C_INT), INTENT(IN) :: nEventsIn
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the observations corresonding to event evnmbr.
!>
!>    @param[in] evnmbr     Event number to set.
!>    @param[in] nttimes    Number of travel times to set.
!>    @param[in] lhaveOTin  If true then the origin time is known.
!>    @param[in] t0In       If lhaveOTin is true then this is the origin time in
!>                          UTC epochal seconds.
!>    @param[in] obs2tf     This is a map from the observation number to the travel-time
!>                          field number.  This is a vector of dimension [nttimes].
!>    @param[in] pickTimes  Pick times (UTC epochal seconds).  This is a vector of
!>                          dimension [nttimes].
!>    @param[in] wts        These are the weights corresponding to the observations.
!>                          The observation will only be used if the weight is positive.
!>                          This is a vector of dimension [nttimes]. 
!> 
!>    @param[out] ierr      0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_setObservation64f(evnmbr, nttimes, lhaveOTin, &
                                          t0In, obs2tf,               &
                                          pickTimes, wts, ierr)       &
      BIND(C, NAME='locate_setObservation64f')
      USE FTEIK_CONSTANTS64F, ONLY : TRUE, FALSE
      USE FTEIK_CONSTANTS32F, ONLY : zero 
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: evnmbr, nttimes
      LOGICAL(C_BOOL), VALUE, INTENT(IN) :: lhaveOTin
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: t0In
      INTEGER(C_INT), INTENT(IN) :: obs2tf(nttimes)
      REAL(C_DOUBLE), INTENT(IN) :: pickTimes(nttimes), wts(nttimes)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) i, indx, k1, k2
      ierr = 0 
      IF (.NOT.linit) THEN
         WRITE(*,"('locate_setObservation64f: Locator not initialized',A)")
         ierr = 1 
         RETURN
      ENDIF
      IF (nttimes < 1 .OR. nttimes > ntf) THEN
         WRITE(*,900) nttimes, ntf
  900    FORMAT('locate_setObservation64f: nttimes=',I4,' must be in [1,',I3']')
         ierr = 1
         RETURN
      ENDIF
      IF (MINVAL(obs2tf) < 1) THEN
         WRITE(*,"('locate_setObservation64f: All values in obs2tf must be positive',A)")
         ierr = 1
         RETURN
      ENDIF
      IF (MAXVAL(obs2tf) > ntf) THEN
         WRITE(*,905) ntf
  905    FORMAT('locate_setObservation64f: All values in obs2tf must be <',I4)
         ierr = 1
         RETURN
      ENDIF
      IF (evnmbr < 1 .OR. evnmbr > nEvents) THEN
         WRITE(*,910) evnmbr, nEvents
  910    FORMAT('locate_setObservation64f: Event number=',I5,' must be in range [1,',I5)
         ierr = 1
         RETURN
      ENDIF
      ! Reduce the pick times to avoid overflow
      tepoch(evnmbr) = MAXVAL(pickTimes)
      ! Deal with the origin times
      originTimes(evnmbr) = zero
      IF (lhaveOTin) THEN
         lhaveOT(evnmbr) = TRUE
         originTimes(evnmbr) = REAL(t0In - tepoch(evnmbr))
      ENDIF
      ! Finally set the data
      k1 = (evnmbr - 1)*lntf + 1 
      k2 = evnmbr*lntf
      observations(k1:k2) = zero
      weights(k1:k2) = zero
      lskip(k1:k2) = TRUE
      DO 1 i=1,nttimes
         indx = k1 + obs2tf(i) - 1
         IF (wts(i) > 0.d0) THEN
            observations(indx)  = REAL(pickTimes(i) - tepoch(evnmbr)) ! Reduce time
            weights(indx) = REAL(wts(i))
            lskip(indx) = FALSE
         ENDIF
    1 CONTINUE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the number of events to locate.
!>
!>    @param[in] nEventsIn    Number of events to locate.
!>
!>    @param[out] ierr        0 indicates success. 
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_setNumberOfEvents(nEventsIn, ierr) &
      BIND(C, NAME='locate_setNumberOfEvents')
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
!>    @param[in] ntfIn  Number of travel time fields. 
!>
!>    @param[out] ierr  0 indicate success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>                        
      SUBROUTINE locate_setNumberOfTravelTimeFields(ntfIn, ierr) &
      BIND(C, NAME='locate_setNumberOfTravelTimeFields')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: ntfIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ntf = 0
      ierr = 0 
      IF (ntfIn < 1) THEN
         WRITE(*,900) ntfIn
  900    FORMAT('locate_setNumberOfTravelTimeFields: ntf=', I4,' must be positive')
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
!>    @brief Computes the origin time and objective function for a least-squres inversion
!>           at each grid point for the given event.  The origin time at each grid
!>           point, \f$ t_0(\textbf{x}_0) \f$, is found by computing
!>           \f[
!>               t_0(\textbf{x}_0)
!>             = \sum_{i}^{n_{obs}}
!>               \frac{ \sum_i (\tau_i - (T(\textbf{x}_0; \textbf{x}_i) + t_i^s))}
!>                    { \sum_{i}^{n_{obs}} w_i }
!>           \f]
!>           where \f$ \tau_i \f$ is the travel time for the i'th observation, 
!>           \f$ T(\textbf{x}_0; \textbf{x}_i) \f$ are the traveltimes from the i'th
!>           observation to all points in the model, and \f$ t_i^2 \f$ is the static
!>           correction for the i'th observation.
!>
!>           The origin time minimizes the travel time  L2 objective function:
!>           \f[
!>              \mathcal{C}(\textbf{x}_0)
!>             =\frac{1}{2} \sum_i^{n_{obs}}
!>               w_i \left (\tau_i - (T(\textbf{x}_0, \textbf{x}_i)
!>                         + t_i^s + t_0(\textbf{x}_0)) \right )^2
!>           \f]
!>
!>    @param[in] evnmbr   Event number.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE locate_computeL2ObjectiveFunction32f(evnmbr)
      USE FTEIK_CONSTANTS32F, ONLY : half, one, sqrtHalf, zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: evnmbr 
      REAL(C_FLOAT) est, res, obs, toff, wt, xnorm
      REAL(C_FLOAT), CONTIGUOUS, POINTER :: ttPtr(:)
      REAL(C_FLOAT), CONTIGUOUS, POINTER :: wtPtr(:)
      REAL(C_FLOAT), CONTIGUOUS, POINTER :: obsPtr(:)
      REAL(C_FLOAT), CONTIGUOUS, POINTER :: scPtr(:)
      REAL(C_FLOAT), CONTIGUOUS, POINTER :: objPtr(:)
      REAL(C_FLOAT), CONTIGUOUS, POINTER :: t0Ptr(:)
      LOGICAL(C_BOOL), CONTIGUOUS, POINTER :: skipPtr(:)
      INTEGER(C_INT) grd1, grd2, i, igrd, jgrd, k1, k2, ngrdLoc, obs1, obs2
      NULLIFY(ttPtr)
      NULLIFY(wtPtr)
      NULLIFY(obsPtr)
      NULLIFY(scPtr)
      NULLIFY(objPtr)
      NULLIFY(t0Ptr)
      NULLIFY(skipPtr)
      !$OMP SIMD aligned(objWork: 64)
      DO 1 igrd=1,ngrd
         objWork(igrd) = zero
    1 CONTINUE 
      ! Set my pointers
      obs1 = (evnmbr - 1)*lntf + 1 
      obs2 = evnmbr*lntf
      scPtr   => staticCorrection(1:ntf)
      obsPtr  => observations(obs1:obs2)
      scPtr   => staticCorrection(1:ntf)
      wtPtr   => weights(obs1:obs2)
      skipPtr => lskip(obs1:obs2)
      ! Compute the sum of the weights - this is the normalization in Eqn A5
      xnorm = one/SUM(wtPtr)
      ! Compute the origin times and travel times
      IF (lwantT0Field) THEN
         ! Loop on grid blocks
         !$OMP PARALLEL DO DEFAULT(NONE) &
         !$OMP PRIVATE(objPtr, t0Ptr, ttPtr, est, i, grd1, grd2, igrd, jgrd) &
         !$OMP PRIVATE(k1, k2, ngrdLoc, obs, res, toff, wt) &
         !$OMP FIRSTPRIVATE(xnorm) &
         !$OMP SHARED(ldgrd, lwantT0Field, ngrd, ntf, objWork, obsPtr, scPtr, skipPtr) &
         !$OMP SHARED(t0Work, ttFields, wtPtr)
         DO 11 igrd=1,ngrd,blockSize
            grd1 = igrd
            grd2 = MIN(igrd+blockSize-1, ngrd)
            ngrdLoc = grd2 - grd1 + 1
            IF (lwantT0Field) THEN
               t0Ptr => t0Work(grd1:grd2)
            ELSE
               t0Ptr => t0Work(1:blockSize)
            ENDIF
            objPtr => objWork(grd1:grd2)
            ! Loop on the observations and tabulate the origin time
            DO 12 i=1,ntf
               IF (skipPtr(i)) CYCLE
               wt = wtPtr(i)             ! Weight for this observation
               toff = scPtr(i)           ! Static correction for this observation
               obs  = obsPtr(i)          ! Travel-time field for this observation
               k1 = (i - 1)*ldgrd + grd1 ! Start chunk of observation's travel time field 
               k2 = (i - 1)*ldgrd + grd2 ! End chunk of observation's travel time field
               ttPtr  => ttFields(k1:k2) ! Chunk of travel-time field
               !$OMP SIMD ALIGNED(t0Ptr, ttptr: 64)
               DO 13 jgrd=1,ngrdLoc
                  est = ttPtr(jgrd) + toff           ! Estimate + static correction
                  res = xnorm*(obs - est)            ! Normalized residual
                  t0Ptr(jgrd) = t0Ptr(jgrd) + wt*res ! Stack the weighted residual
   13          CONTINUE
               NULLIFY(ttPtr)
   12       CONTINUE ! End computation of origin times
            ! Loop on the observations and tabulate the objective function
            DO 15 i=1,ntf
               IF (skipPtr(i)) CYCLE
               wt = sqrtHalf*wtPtr(i)    ! Weight for this observation
               toff = scPtr(i)           ! Static correction for this observation
               obs  = obsPtr(i)          ! Travel-time field for this observation
               k1 = (i - 1)*ldgrd + grd1 ! Start chunk of observation's travel time field 
               k2 = (i - 1)*ldgrd + grd2 ! End chunk of observation's travel time field
               ttPtr  => ttFields(k1:k2) ! Chunk of travel-time field
               !$OMP SIMD ALIGNED(ttPtr, t0Ptr, objPtr: 64) 
               DO 16 jgrd=1,ngrdLoc
                  est = ttPtr(jgrd) + t0Ptr(jgrd) + toff ! Eqn A1 + static correction
                  res = wt*(obs - est)                   ! Weighted residual
                  !objWork(igrd+jgrd-1) = objWork(igrd+jgrd-1) + res*res 
                  objPtr(jgrd) = objPtr(jgrd) + res*res 
   16          CONTINUE
               NULLIFY(ttPtr)
   15       CONTINUE
            NULLIFY(t0Ptr)
            NULLIFY(objPtr)
   11    CONTINUE
         !$OMP END PARALLEL DO
      ENDIF

!if (size(objWork) >= 474953) then
!do i=1,ngrd
! if (objWork(i) <= objWork(474953)) then != zero) then 
!   print *, i, objWork(i), t0Work(i)
!  !pause
! endif
!enddo
!endif
!print *, 'objWork:', objWork(474953), objWork(455965), objWork(474953+1), t0Work(474953)
 
!     DO 21 igrd=1,ngrd
!        t0 = t0Work(igrd) 
!        tt = ttptr(igrd)
!        objfn = 0.0
!        !$OMP SIMD ALIGNED(lskip, scptr, obsptr, wtptr: 64) REDUCTION(+:objfn)
!        DO 22 i=1,ntf
!           est = tt + t0 + scptr(i) ! travel time + origin time + static correction 
!           res = obsptr(i) - est
!           res2 = res*res
!           IF (lskip(i)) res2 = 0.0 
!           objfn = objfn + wtptr(i)*res2
!  22    CONTINUE 
!        objWork(igrd) = half*objfn
!  21 CONTINUE
      NULLIFY(obsPtr)
      NULLIFY(scPtr)
      NULLIFY(wtPtr)
      NULLIFY(skipPtr)
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
      !$OMP SIMD !aligned(tstatic: 64)
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
         !$OMP SIMD ALIGNED(obsptr, ttptr, wtsptr: 64) !tstatic, obsptr, ttptr, wtsptr: 64)
         DO i=1,nobs
            test = ttptr(i) + tori
            res = wtsptr(i)*(obsptr(i) - test)
            tstatic(i) = tstatic(i) + res
         ENDDO
      ENDDO
      RETURN
      END
!========================================================================================!
!                                 Begin the MPI modules                                  !
!========================================================================================!
#if defined(FTEIK_FORTRAN_USE_MPI)
!>    @brief Initializes the MPI based locator.
!>
!>    @param[in] root       Root process on communicator.
!>    @param[in] comm       MPI communicator handle.
!>    @param[in] nEventsIn  Number of input events.  This must be defined on the
!>                          root process.
!>    @param[in] ntfIn      Number of travel time fields.  This must be defined on
!>                          the root process.
!>    @param[in] ngrdIn     Total number of grid points in travel time field.
!>                          This must be defined on the root process.
!>
!>    @param[out] ierr      0 indicates success.
!>
!> @author Ben Baker
!>
!> @copyright MIT
!> 
      SUBROUTINE locate_initializeMPIF(root, comm,                     &
                                       nEventsIn, ntfIn, ngrdIn, ierr) &
      BIND(C, NAME='locate_initializeMPIF')
      USE MPI_F08
      USE ISO_C_BINDING
      TYPE(MPI_Comm), VALUE, INTENT(IN) :: comm
      INTEGER(C_INT), VALUE, INTENT(IN) :: root, nEventsIn, ntfIn, ngrdIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT), ALLOCATABLE :: grdPtr(:), lkeep(:)
      INTEGER(C_INT) dGrid, i, ierrLoc, nEventsTemp, &
                     ntfTemp, ngrdLoc, ngrdTemp
      TYPE(MPI_Group) :: group, keepGroup
      ierr = 0
      CALL MPI_Comm_rank(comm, mylocatorID, mpierr)
      CALL MPI_Comm_size(comm, nprocs,      mpierr)
      ! Broadcast anticipated sizes
      IF (mylocatorID == root) THEN
         nEventsTemp = nEventsIn
         ntfTemp     = ntfIn
         ngrdTemp    = ngrdIn
      ENDIF
      CALL MPI_Bcast(nEventsTemp, 1, MPI_INTEGER, root, comm, mpierr)
      CALL MPI_Bcast(ntfTemp,     1, MPI_INTEGER, root, comm, mpierr)
      CALL MPI_Bcast(ngrdTemp,    1, MPI_INTEGER, root, comm, mpierr)
      ! Create the group which will do the location s.t. the group size doesn't exceed
      ! the number of grid points. 
      ALLOCATE(lkeep(nprocs))
      lkeep(:) = 0
      DO i=1,MIN(nprocs, ngrdTemp)
         lkeep(i) = i - 1
      ENDDO
      CALL MPI_Comm_group(comm, group, mpierr)
      CALL MPI_Group_incl(group, MIN(ngrdTemp, nprocs), lkeep, keepGroup, mpierr)
      CALL MPI_Comm_create_group(comm, keepGroup, 99, locatorComm, mpierr)
      DEALLOCATE(lkeep)
      CALL MPI_Group_free(keepGroup, mpierr)
      CALL MPI_Group_free(group, mpierr)
      ! Initialize the group and set the grid points the process is response for
      linLocator = .FALSE.
      IF (mylocatorID < MIN(ngrdTemp, nprocs)) THEN
         ALLOCATE(grdPtr(nprocs+1))
         IF (nprocs == 1) THEN
            grdPtr(1) = 1
            grdPtr(2) = ngrdTemp + 1
         ELSE
            dGrid = MAX(1, INT(REAL(ngrdTemp)/REAL(nprocs) + 0.5))
            grdPtr(1) = 1
            DO i=1,nprocs
               grdPtr(i+1) = MIN(ngrdTemp+1, grdPtr(i) + dGrid)
            ENDDO
            grdPtr(nprocs+1) = ngrdTemp + 1
         ENDIF
         ngrdLoc = grdPtr(mylocatorID+2) - grdPtr(mylocatorID+1)
         CALL locate_initialize(nEventsTemp, ntfTemp, ngrdLoc, ierrLoc)
         IF (ierrLoc /= 0) THEN
            WRITE(*,900) mylocatorID
  900       FORMAT('locate_initializeMPIF: Error initializing on process ', I4)
            ierrLoc = 1
         ENDIF
         CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, locatorComm, mpierr)
         IF (ierr == 0) THEN
            CALL MPI_Comm_size(locatorComm, nprocs,      mpierr)
            CALL MPI_Comm_rank(locatorComm, mylocatorID, mpierr)
            linLocator = .TRUE.
            IF (ALLOCATED(sendCounts)) DEALLOCATE(sendCounts)
            IF (ALLOCATED(displs))     DEALLOCATE(displs)
            ALLOCATE(sendCounts(nprocs))
            ALLOCATE(displs(nprocs))
            DO i=1,nprocs
               sendCounts(i) = grdPtr(i+1) - grdPtr(i)
               displs(i) = grdPtr(i) - 1
            ENDDO
         ENDIF
         DEALLOCATE(grdPtr)
      ! This process process isn't in the group
      ELSE
         linit = .TRUE.
         nprocs = 0
      ENDIF 
      CALL MPI_Bcast(ierr, 1, MPI_INTEGER, root, comm, mpierr) 
      IF (ierr /= 0) luseMPI = .TRUE.
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gets the ranks of the process in the locator.
!>
!>    @param[out] rank    MPI rank of process in locator.  If negative then the process
!>                        is not in the locator.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_getRank(rank) BIND(C, NAME='locate_getRank')
      USE MPI_F08
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(OUT) :: rank
      INTEGER(C_INT) mpierr
      rank = 0
      IF (.NOT. luseMPI) RETURN
      IF (.NOT. linLocator) THEN
         rank =-1
         RETURN
      ENDIF
      CALL MPI_Comm_rank(locatorComm, rank, mpierr)
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the travel time fields on all processes.
!>
!>    @param[in] root       Root process ID on the locator communicator containing 
!>                          all inputs.
!>    @param[in] ngrdIn     Total number of grid points in model.  This must be defined
!>                          on the root process.
!>    @param[in] itfIn      Travel-time field number.  This must be defined on the root
!>                          process.
!>    @param[in] ttIn       The travel time field to distribute to all processes.  This
!>                          is a vector of dimension [ngrdIn] and must be defined on
!>                          the root process.
!>
!>    @param[out] ierr      0 indicates success. 
!>
!>    @result 0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_setTravelTimeField64fMPIF(root, ngrdIn, itfIn, &
                                                  ttIn, ierr)          &
      BIND(C, NAME='locate_setTravelTimeField64fMPIF')
      USE MPI_F08
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: root, ngrdIn, itfIn
      REAL(C_DOUBLE), INTENT(IN) :: ttIn(ngrdIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE), ALLOCATABLE :: work(:)
      INTEGER(C_INT) itf, mpierr
      ierr = 0
      ! Have the root send appropriate protions of travel time field to each process
      ALLOCATE(work(MAX(1, ngrd)))
      IF (mylocatorID == root) itf = itfIn
      CALL MPI_Bcast(itf, 1, MPI_INTEGER, root, locatorComm, mpierr)
      CALL MPI_Scatterv(ttIn, sendCounts, displs, MPI_DOUBLE_PRECISION, work,   &
                        ngrd, MPI_DOUBLE_PRECISION, root, locatorComm, mpierr)
      IF (mpierr /= MPI_SUCCESS) THEN
         WRITE(*,'("locate_setTravelTimeField64fMPIF: Error scattering travel times",A)')
         ierr = 1
         RETURN
      ENDIF
      ! And set it
      CALL locate_setTravelTimeField64f(ngrd, itf, work, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,905) mylocatorID
  905    FORMAT('locate_setTravelTimeField64fMPIF: Error on rank', I4)
      ENDIF
      DEALLOCATE(work)      
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Locates an event with MPI.
!>
!>    @param[in] evnmbr    Event number.
!>    @param[out] optIndx  Index in grid corresponding to the optimum.
!>    @param[out] t0Opt    The origin time (epochal seconds; UTC) at the optimal
!>                         grid point.
!>    @param[out] objOpt   Value of objective function at optimal grid point.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE locate_locateEventMPIF(evnmbr, optIndx, t0Opt, objOpt) &
      BIND(C, NAME='locate_locateEventMPIF')
      USE MPI_F08
      USE FTEIK_CONSTANTS32F, ONLY : zero
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: evnmbr
      REAL(C_DOUBLE), INTENT(OUT) :: t0Opt, objOpt
      INTEGER(C_INT), INTENT(OUT) :: optIndx
      REAL(C_DOUBLE)  objWork1(nprocs), objWork2(nprocs), &
                      t0Work1(nprocs), t0Work2(nprocs)
      INTEGER(C_INT) optIndx1(nprocs), optIndx2(nprocs)
      INTEGER(C_INT), PARAMETER :: master = 0
      ! Initialize
      optindx = 0 
      objWork1(:) = zero
      t0Work1(:) = zero
      optIndx1(:) = 0
      ! Locate the event on each processes's grid
      CALL locate_locateEvent(evnmbr, optIndx1(mylocatorID+1), &
                              t0Work1(mylocatorID+1), objWork1(mylocatorID+1))
      ! Reduce the results onto the master 
      CALL MPI_Reduce(objWork1, objWork2, nprocs, MPI_DOUBLE, MPI_SUM,  &
                      master, locatorComm, mpierr)
      CALL MPI_Reduce(t0Work1,  t0Work2,  nprocs, MPI_DOUBLE, MPI_SUM,  &
                      master, locatorComm, mpierr)
      CALL MPI_Reduce(optIndx1, optIndx2, nprocs, MPI_INTEGER, MPI_SUM,  &
                      master, locatorComm, mpierr) 
      ! Let the master compute the optimum
      IF (mylocatorID == master) THEN
print *, objWork2
         optIndx = MINLOC(objWork2, 1) ! Get the smallest objective function
         objOpt = objWork2(optIndx)    ! This is the optimum's value
         t0Opt  = t0Work2(optIndx)     ! This is the corresponding origin time
         optIndx = displs(optIndx) + optIndx2(optIndx) ! This is the optimum global index
      ENDIF
      CALL MPI_Bcast(objOpt,  1, MPI_DOUBLE,  master, locatorComm, mpierr)
      CALL MPI_Bcast(t0Opt,   1, MPI_DOUBLE,  master, locatorComm, mpierr)
      CALL MPI_Bcast(optIndx, 1, MPI_INTEGER, master, locatorComm, mpierr)
      print *, 'opt info:', optIndx, t0Opt, objOpt
      RETURN
      END

#endif
END MODULE
