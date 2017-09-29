#define ERRORMSG(msg) WRITE(0,'("[ERROR] at ",I4," in file ",/,A,/,"Error message: ",A)') __LINE__,__FILE__,msg

MODULE FTEIK_MEMORY
   USE ISO_C_BINDING
   CONTAINS
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
   USE ISO_C_BINDING
   REAL(C_DOUBLE), POINTER, DIMENSION(:), INTENT(INOUT) :: x
   INTEGER(C_INT), INTENT(IN) :: align, n
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
   !                                                                                     !
   !=====================================================================================!
   !                                                                                     !
   !> @brief Frees a double precision Fortran pointer.
   !>
   !> @param[in,out] x   On input this is an associated Fortran pointer. \n
   !>                    On ouptut the pointer has been freed and set to nullified. 
   !>
   !> @author Ben Baker
   !>
   !> @copyright MIT
   !>
   SUBROUTINE free64f(x)
   USE ISO_C_BINDING
   REAL(C_DOUBLE), POINTER, INTENT(INOUT) :: x(:)
   TYPE(C_PTR), TARGET :: cptr = C_NULL_PTR
   INTERFACE
      ! void free(void *);
      SUBROUTINE free(ptr) BIND(C, NAME="free")
      USE ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: ptr
      END SUBROUTINE
   END INTERFACE
   IF (.NOT.ASSOCIATED(x)) RETURN
   cptr = C_LOC(x(1))
   CALL free(cptr)
   cptr = C_NULL_PTR
   NULLIFY(x)
   END SUBROUTINE

   !> @brief Computes the requisite padding so that n + padLength will yield a double 
   !>        matrix that is aligned on the given bit boundary.
   !>
   !> @param[in] alignment   Bit alignment.  This must be a power of 2 (e.g., 64).
   !> @param[in] n           Number of elements for which to compute padding.
   !>
   !> @result The padding so that n + padLength is bit aligned row or column.
   !>
   !> @author Ben Baker
   !>
   !> @copyright MIT
   !>
   INTEGER(C_INT) FUNCTION padLength64F(alignment, n) &
   BIND(C, NAME='padLength64fF')
   USE ISO_C_BINDING
   IMPLICIT NONE
   INTEGER(C_INT), INTENT(IN) :: alignment, n
   INTEGER(C_INT) xmod
   INTEGER, PARAMETER :: sizeof_double = 8
   padLength64F = 0
   xmod = MOD(n*sizeof_double, alignment)
   IF (xmod /= 0) padLength64F = (alignment - xmod)/sizeof_double
   RETURN
   END FUNCTION

   !> @brief Computes the requisite padding so that n + padLength will yield a float
   !>        matrix that is aligned on the given bit boundary.
   !>  
   !> @param[in] alignment   Bit alignment.  This must be a power of 2 (e.g., 64).
   !> @param[in] n           Number of elements for which to compute padding.
   !>  
   !> @result The padding so that n + padLength is bit aligned row or column.
   !>  
   !> @author Ben Baker
   !>  
   !> @copyright MIT
   !>
   INTEGER(C_INT) FUNCTION padLength32F(alignment, n) &
   BIND(C, NAME='padLength32fF')
   USE ISO_C_BINDING
   IMPLICIT NONE
   INTEGER(C_INT), INTENT(IN) :: alignment, n
   INTEGER(C_INT) xmod
   INTEGER, PARAMETER :: sizeof_float = 4
   padLength32F = 0 
   xmod = MOD(n*sizeof_float, alignment)
   IF (xmod /= 0) padLength32F = (alignment - xmod)/sizeof_float
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
