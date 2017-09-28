      MODULE FTEIK_CONSTANTS64F
         USE ISO_C_BINDING
         IMPLICIT NONE 
         !> Default large value travel times will be set to.
         REAL(C_DOUBLE), PUBLIC, PARAMETER :: FTEIK_HUGE = 99999.d0
         !> Double preicision zero.
         REAL(C_DOUBLE), PUBLIC, PARAMETER :: zero = 0.d0 
         !> Double precision half
         REAL(C_DOUBLE), PUBLIC, PARAMETER :: half = 0.5d0
         !> Double precision one. 
         REAL(C_DOUBLE), PUBLIC, PARAMETER :: one = 1.d0 
         !> Double precision four.
         REAL(C_DOUBLE), PUBLIC, PARAMETER :: four = 4.d0
         !> Double precision nine.
         REAL(C_DOUBLE), PUBLIC, PARAMETER :: nine = 9.d0
         !> Smallest double precision number.
         REAL(C_DOUBLE), PUBLIC, PARAMETER :: DBL_EPSILON = EPSILON(1.d0)
         !> Source perturbation to avoid collocating to grid point.
         REAL(C_DOUBLE), PUBLIC, PARAMETER :: perturbSource = 0.0001d0
         !> Chunk size in level scheduling.
         INTEGER(C_INT), PUBLIC, PARAMETER :: chunkSize = 16
      END MODULE !FTEIK_CONSTANTS64
