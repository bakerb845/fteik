#define ERRORMSG(msg) WRITE(ERROR_UNIT,'("[ERROR] at ",I4," in file ",/,A,/,"Error message: ",A)') __LINE__,__FILE__,msg

!> @defgroup memory Memory
!> @brief Memory allocation/deallocation routines.
!> @copyright Ben Baker distributed under the MIT license.
MODULE FTEIK_MEMORY
   USE ISO_C_BINDING
   USE ISO_FORTRAN_ENV
   IMPLICIT NONE

   PUBLIC :: allocate64f
   PUBLIC :: allocate32f
   !PUBLIC :: free64
   !PUBLIC :: free32
   PUBLIC :: padLength64F
   PUBLIC :: padLength32F
   CONTAINS
   !-------------------------------------------------------------------------------------!
   !                                   Begin the Code                                    !
   !-------------------------------------------------------------------------------------!
   !> @brief Allocates aligned memory to a Fortran pointer.  For example, to 
   !>        allocate a 64 bit aligned vector with 250 elements one would call
   !>        CALL allocate64f(x, 64, 250).
   !>
   !> @param[in,out] x     On input this is a NULL Fortran pointer. \n
   !>                      On output contains a pointer to aligned memory.
   !>
   !> @param[in] align     Alignment value (e.g., 64 would be 64 bit alignment).
   !> @param[in] n         Number of elements in output double pointer. 
   !>
   !> @author Ben Baker
   !>
   !> @copyright MIT
   !>
   SUBROUTINE allocate64f(x, align, n)
   REAL(C_DOUBLE), POINTER, DIMENSION(:), INTENT(INOUT) :: x
   INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: align, n
   INTEGER(C_SIZE_T) alignment, ns
   TYPE(C_PTR) :: cptr = C_NULL_PTR
   INTERFACE
      ! void *aligned_alloc(size_t alignment, size_t size);
      TYPE(C_PTR) FUNCTION aligned_alloc(alignment, len) &
      BIND(C, NAME="aligned_alloc")
      USE ISO_C_BINDING, ONLY : C_PTR, C_SIZE_T
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: alignment
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: len
      END FUNCTION

      ! void posix_memalign(void **memptr, size_t alignment, size_t size);
      SUBROUTINE posix_memalign(memptr, alignment, sz) &
      BIND(C, NAME='posix_memalign')
      USE ISO_C_BINDING, ONLY : C_PTR, C_SIZE_T
      TYPE(C_PTR), INTENT(INOUT) :: memptr 
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: alignment
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: sz
      END SUBROUTINE
   END INTERFACE
   CALL free64f(x)
   alignment = align
   ns = SIZEOF(1.d0)*n
#if defined(USE_POSIX)
   CALL posix_memalign(cptr, alignment, ns)
#else
   cptr = aligned_alloc(alignment, ns)
#endif
   CALL C_F_POINTER(cptr, x, [n])
   x(:) = 0.d0
   cptr = C_NULL_PTR
   RETURN
   END SUBROUTINE
   !> @brief Allocates aligned memory to a Fortran pointer.  For example, to 
   !>        allocate a 64 bit aligned vector with 250 elements one would call
   !>        CALL allocate64f(x, 64, 250).
   !>  
   !> @param[in,out] x     On input this is a NULL Fortran pointer. \n
   !>                      On output contains a pointer to aligned memory.
   !>  
   !> @param[in] align     Alignment value (e.g., 64 would be 64 bit alignment).
   !> @param[in] n         Number of elements in output double pointer. 
   !> @ingroup memory
   SUBROUTINE allocate32f(x, align, n)
   REAL(C_FLOAT), POINTER, DIMENSION(:), INTENT(INOUT) :: x
   INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: align, n
   INTEGER(C_SIZE_T) alignment, ns
   TYPE(C_PTR) :: cptr = C_NULL_PTR
   INTERFACE
      ! void *aligned_alloc(size_t alignment, size_t size);
      TYPE(C_PTR) FUNCTION aligned_alloc(alignment, len) &
      BIND(C, NAME="aligned_alloc")
      USE ISO_C_BINDING, ONLY : C_PTR, C_SIZE_T
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: alignment
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: len
      END FUNCTION

      ! void posix_memalign(void **memptr, size_t alignment, size_t size);
      SUBROUTINE posix_memalign(memptr, alignment, sz) &
      BIND(C, NAME='posix_memalign')
      USE ISO_C_BINDING, ONLY : C_PTR, C_SIZE_T
      TYPE(C_PTR), INTENT(INOUT) :: memptr
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: alignment
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: sz
      END SUBROUTINE
   END INTERFACE
   CALL free32f(x)
   alignment = align
   ns = SIZEOF(1.0)*n
#if defined(USE_POSIX)
   CALL posix_memalign(cptr, alignment, ns)
#else
   cptr = aligned_alloc(alignment, ns)
#endif
   CALL C_F_POINTER(cptr, x, [n])
   x(:) = 0.0
   cptr = C_NULL_PTR
   RETURN
   END SUBROUTINE
   !                                                                                     !
   !=====================================================================================!
   !                                                                                     !
   !> @brief Frees a double precision Fortran pointer.
   !>
   !> @param[in,out] x   On input this is an associated Fortran pointer. \n
   !>                    On ouptut the pointer has been freed and set to nullified. 
   !> @ingroup memory
   SUBROUTINE free64f(x)
   REAL(C_DOUBLE), POINTER, INTENT(INOUT) :: x(:)
   TYPE(C_PTR), TARGET :: cptr = C_NULL_PTR
   INTERFACE
      ! void free(void *);
      SUBROUTINE freeF(ptr) BIND(C, NAME="free")
      USE ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: ptr
      END SUBROUTINE
   END INTERFACE
   IF (.NOT.ASSOCIATED(x)) RETURN
   cptr = C_LOC(x(1))
   CALL freeF(cptr)
   cptr = C_NULL_PTR
   NULLIFY(x)
   END SUBROUTINE
   !                                                                                     !
   !=====================================================================================!
   !                                                                                     !
   !> @brief Frees a double precision Fortran pointer.
   !> @param[in,out] x   On input this is an associated Fortran pointer. \n
   !>                    On ouptut the pointer has been freed and set to nullified. 
   !> @ingroup memory
   SUBROUTINE free32f(x)
   REAL(C_FLOAT), POINTER, INTENT(INOUT) :: x(:)
   TYPE(C_PTR), TARGET :: cptr = C_NULL_PTR
   INTERFACE
      ! void free(void *);
      SUBROUTINE freeF(ptr) BIND(C, NAME="free")
      USE ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: ptr
      END SUBROUTINE
   END INTERFACE
   IF (.NOT.ASSOCIATED(x)) RETURN
   cptr = C_LOC(x(1))
   CALL freeF(cptr)
   cptr = C_NULL_PTR
   NULLIFY(x)
   END SUBROUTINE
   !                                                                                     !
   !=====================================================================================!
   !                                                                                     !
   !> @brief Computes the requisite padding so that n + padLength will yield a double 
   !>        matrix that is aligned on the given bit boundary.
   !>
   !> @param[in] alignment   Bit alignment.  This must be a power of 2 (e.g., 64).
   !> @param[in] n           Number of elements for which to compute padding.
   !> @ingroup memory
   !>
   INTEGER(C_INT) FUNCTION padLength64F(alignment, n) &
   BIND(C, NAME='padLength64fF')
   INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: alignment
   INTEGER(C_INT), VALUE, INTENT(IN) :: n
   INTEGER(C_INT) xmod
   INTEGER, PARAMETER :: sizeof_double = 8
   padLength64F = 0
   xmod = MOD(n*sizeof_double, INT(alignment))
   IF (xmod /= 0) padLength64F = (INT(alignment) - xmod)/sizeof_double
   RETURN
   END FUNCTION
   !> @brief Computes the requisite padding so that n + padLength will yield a float
   !>        matrix that is aligned on the given bit boundary.
   !>  
   !> @param[in] alignment   Bit alignment.  This must be a power of 2 (e.g., 64).
   !> @param[in] n           Number of elements for which to compute padding.
   !>  
   !> @result The padding so that n + padLength is bit aligned row or column.
   !> @ingroup memory
   !>
   INTEGER(C_INT) FUNCTION padLength32F(alignment, n) &
   BIND(C, NAME='padLength32fF')
   INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: alignment
   INTEGER(C_INT), VALUE, INTENT(IN) :: n
   INTEGER(C_INT) xmod
   INTEGER, PARAMETER :: sizeof_float = 4
   padLength32F = 0 
   xmod = MOD(n*sizeof_float, INT(alignment))
   IF (xmod /= 0) padLength32F = (INT(alignment) - xmod)/sizeof_float
   RETURN
   END FUNCTION

   !veclen = alignment/sizeof(1.d0)
   !IF (veclen*sizeof(double) /= alignment) THEN
   !   ERRORMSG("Invalid alignment")
   !   ierr = 1
   !   RETURN
   !ENDIF
   !padLength64f = ((n + veclen - 1)/veclen)*veclen 
   !RETURN
   !END FUNCTION
    

END MODULE

!USE FTEIK_MEMORY
!use iso_c_binding
!real*8, pointer, dimension(:) :: x
!NULLIFY(x)
!call allocate64f(x, 64, 200)
!x(:) = 0.d0
!print *, size(x), maxval(x)
!call free64f(x)
!stop
!end
