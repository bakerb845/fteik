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
