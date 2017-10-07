MODULE IPPS_MODULE
#if defined(FTEIK_FORTRAN_USE_INTEL)
      INTERFACE
         !> @brief Converts a 64 bit float to a 32 bit float.
         !> 
         !> @param[in] pSrc   64 bit source array.  This is a vector of
         !>                   dimension [len].
         !> @param[out] pDst  32 bit destination array.  This is a vector of
         !>                   dimension [len].
         !> @param[in] len    The length of the arrays.
         !>
         INTEGER(C_INT) FUNCTION ippsConvert_64f32f(pSrc, pDst, len) &
         BIND(C, NAME='ippsConvert_64f32f')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT), VALUE, INTENT(IN) :: len 
         REAL(C_DOUBLE), INTENT(IN) :: pSrc(len)
         REAL(C_FLOAT), INTENT(OUT) :: pDst(len)
         END FUNCTION
      END INTERFACE
#endif
END MODULE
