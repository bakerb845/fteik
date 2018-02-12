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
         !> Double precision \f$ \sqrt{1/2} \f$.
         REAL(C_DOUBLE), PUBLIC, PARAMETER :: sqrtHalf = SQRT(half)
         !> Smallest double precision number.
         REAL(C_DOUBLE), PUBLIC, PARAMETER :: DBL_EPSILON = EPSILON(1.d0)
         !> Source perturbation to avoid collocating to grid point.
         REAL(C_DOUBLE), PUBLIC, PARAMETER :: perturbSource = 0.0001d0
         !> Chunk size in level scheduling.
         INTEGER(C_INT), PUBLIC, PARAMETER :: chunkSize = 32
         !> True
         LOGICAL(C_BOOL), PUBLIC, PARAMETER :: TRUE = .TRUE.
         !> False
         LOGICAL(C_BOOL), PUBLIC, PARAMETER :: FALSE = .FALSE.
         !> Natural solver column major ordering -> z, x, y in 3D or z, x  in 2D
         INTEGER(C_INT), PUBLIC, PARAMETER :: FTEIK_NATURAL_ORDERING = 0
         !> Natural solver column major ordering -> z, x, y
         INTEGER(C_INT), PUBLIC, PARAMETER :: FTEIK_ZXY_ORDERING = 0
         !> Natural solver column major ordering -> z, x in 2D 
         INTEGER(C_INT), PUBLIC, PARAMETER :: FTEIK_ZX_ORDERING = 0 
         !> Column major ordering -> x, y, z
         INTEGER(C_INT), PUBLIC, PARAMETER :: FTEIK_XYZ_ORDERING = 1
         !> Column major ordering -> x, z
         INTEGER(C_INT), PUBLIC, PARAMETER :: FTEIK_XZ_ORDERING = 1
         !> Column major ordering -> z, y, x  
         INTEGER(C_INT), PUBLIC, PARAMETER :: FTEIK_ZYX_ORDERING = 2
      END MODULE !FTEIK_CONSTANTS64

      MODULE FTEIK_CONSTANTS32F
         USE ISO_C_BINDING
         USE FTEIK_CONSTANTS64F, ONLY : FTEIK_NATURAL_ORDERING,  &
                                        FTEIK_ZXY_ORDERING,      &
                                        FTEIK_ZX_ORDERING,       &
                                        FTEIK_ZYX_ORDERING
         IMPLICIT NONE 
         !> Default large value travel times will be set to.
         REAL(C_FLOAT), PUBLIC, PARAMETER :: FTEIK_HUGE = 99999.0
         !> single preicision zero.
         REAL(C_FLOAT), PUBLIC, PARAMETER :: zero = 0.0 
         !> Single precision half
         REAL(C_FLOAT), PUBLIC, PARAMETER :: half = 0.50
         !> Single precision one. 
         REAL(C_FLOAT), PUBLIC, PARAMETER :: one = 1.0 
         !> Single  precision four.
         REAL(C_FLOAT), PUBLIC, PARAMETER :: four = 4.0
         !> Single precision nine.
         REAL(C_FLOAT), PUBLIC, PARAMETER :: nine = 9.0
         !> Single precision \f$ \sqrt{1/2} \f$.
         REAL(C_FLOAT), PUBLIC, PARAMETER :: sqrtHalf = SQRT(half)
         !> Smallest double precision number.
         REAL(C_FLOAT), PUBLIC, PARAMETER :: DBL_EPSILON = EPSILON(one)
         !> Source perturbation to avoid collocating to grid point.
         REAL(C_FLOAT), PUBLIC, PARAMETER :: perturbSource = 0.00010
         !> True
         LOGICAL(C_BOOL), PUBLIC, PARAMETER :: TRUE = .TRUE.
         !> False
         LOGICAL(C_BOOL), PUBLIC, PARAMETER :: FALSE = .FALSE.
      END MODULE
