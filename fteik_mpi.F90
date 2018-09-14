MODULE FTEIK_MPI
      USE ISO_FORTRAN_ENV
      USE ISO_C_BINDING
      USE MPI_F08

      CONTAINS
!========================================================================================!
!                                      Begin the Code                                    !
!========================================================================================!
!>    @brief Broadcasts the basic graph startup variables to the other processes
!>           on the communicator.
!>
!>    @param[in] root       The root process ID on the communicator.  This process
!>                          will have performed the graph initialization.
!>    @param[in] comm       MPI communicator on which to broadcast graph.
!>
!>    @param[out] mpierr    MPI_SUCCESS indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_mpi_broadcastGraph64fF(root, comm, mpierr) &
                 BIND(C, NAME='fteik_mpi_broadcastGraph64fF')
!     USE FTEIK_UTILS64F, ONLY : dx, dy, dz, x0, y0, z0
      USE FTEIK_UTILS64F, ONLY : ncell, ngrd, nx, ny, nz, nzx, nzm1, nzm1_nxm1
      USE FTEIK_UTILS64F, ONLY : levelPtr, ijkv1, ijkv2, ijkv3, ijkv4, &
                                 ijkv5, ijkv6, ijkv7, ijkv8, &
                                 nLevels, maxLevelSize
!     USE FTEIK_UTILS64F, ONLY : tt1
!     USE FTEIK_UTILS64F, ONLY : lupd1, lupd2, lupd3, lupd4, lupd5, lupd6, lupd7, lupd8
      USE FTEIK_UTILS64F, ONLY : lhaveGrid, lhaveGridSpacing
!     USE FTEIK_UTILS64F, ONLY : zero
      !USE MPI
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comm, root
      INTEGER, INTENT(OUT) :: mpierr
      INTEGER lreturn, myid
      ! who am i?
!     CALL MPI_Comm_rank(myid, comm, mpierr)
      ! is the root process ready to broadcast?
      IF (myid /= root) THEN
         lhaveGrid = .FALSE.
         lhaveGridSpacing = .FALSE.
      ELSE
         lreturn = 0
         IF (.NOT.lhaveGrid .OR. .NOT.lhaveGridSpacing) THEN
            WRITE(*,*) 'fteik_mpi_broadcastGraph: Graph not yet set'
            lreturn = 1
         ENDIF
      ENDIF
!     CALL MPI_Bcast(lreturn, 1, MPI_INTEGER, root, comm, mpierr)
      IF (lreturn == 1) THEN
         mpierr = MPI_ERR_ARG
         RETURN
      ENDIF
      ! grid sizes
!     CALL MPI_Bcast(ncell,     1, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(ngrd,      1, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(nz,        1, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(nx,        1, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(ny,        1, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(nzx,       1, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(nzm1,      1, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(nzm1_nxm1, 1, MPI_INT32_T, root, comm, mpierr)
      ! graph sizes
!     CALL MPI_Bcast(nLevels,      1, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(maxLevelSize, 1, MPI_INT32_T, root, comm, mpierr)
      ! grid spacing and origin
!     CALL MPI_Bcast(dz, 1, MPI_DOUBLE, root, comm, mpierr)
!     CALL MPI_Bcast(dx, 1, MPI_DOUBLE, root, comm, mpierr)
!     CALL MPI_Bcast(dy, 1, MPI_DOUBLE, root, comm, mpierr)
!     CALL MPI_Bcast(z0, 1, MPI_DOUBLE, root, comm, mpierr)
!     CALL MPI_Bcast(x0, 1, MPI_DOUBLE, root, comm, mpierr)
!     CALL MPI_Bcast(y0, 1, MPI_DOUBLE, root, comm, mpierr)
      ! number of Gauss-Seidel iterations
!     CALL fteik_mpi_broadcastNumberOfSweepsF(root, comm, mpierr)
      ! spherical to cartesian parameter
!     CALL fteik_mpi_broadcastEpsS2C64fF(root, comm, mpierr)
      ! set space
      IF (myid /= root) THEN
         IF (ALLOCATED(levelPtr)) DEALLOCATE(levelPtr)
         IF (ALLOCATED(ijkv1)) DEALLOCATE(ijkv1)
         IF (ALLOCATED(ijkv2)) DEALLOCATE(ijkv2)
         IF (ALLOCATED(ijkv3)) DEALLOCATE(ijkv3)
         IF (ALLOCATED(ijkv4)) DEALLOCATE(ijkv4)
         IF (ALLOCATED(ijkv5)) DEALLOCATE(ijkv5)
         IF (ALLOCATED(ijkv6)) DEALLOCATE(ijkv6)
         IF (ALLOCATED(ijkv7)) DEALLOCATE(ijkv7)
         IF (ALLOCATED(ijkv8)) DEALLOCATE(ijkv8)
         IF (ALLOCATED(lupd1)) DEALLOCATE(lupd1)
         IF (ALLOCATED(lupd2)) DEALLOCATE(lupd2)
         IF (ALLOCATED(lupd3)) DEALLOCATE(lupd3)
         IF (ALLOCATED(lupd4)) DEALLOCATE(lupd4)
         IF (ALLOCATED(lupd5)) DEALLOCATE(lupd5)
         IF (ALLOCATED(lupd6)) DEALLOCATE(lupd6)
         IF (ALLOCATED(lupd7)) DEALLOCATE(lupd7)
         IF (ALLOCATED(lupd8)) DEALLOCATE(lupd8)
         IF (ALLOCATED(tt1))   DEALLOCATE(tt1)
         ALLOCATE(levelPtr(nLevels+1))
         ALLOCATE(ijkv1(4*ngrd))
         ALLOCATE(ijkv2(4*ngrd))
         ALLOCATE(ijkv3(4*ngrd))
         ALLOCATE(ijkv4(4*ngrd))
         ALLOCATE(ijkv5(4*ngrd))
         ALLOCATE(ijkv6(4*ngrd))
         ALLOCATE(ijkv7(4*ngrd))
         ALLOCATE(ijkv8(4*ngrd))
         ALLOCATE(lupd1(ngrd))
         ALLOCATE(lupd2(ngrd))
         ALLOCATE(lupd3(ngrd))
         ALLOCATE(lupd4(ngrd))
         ALLOCATE(lupd5(ngrd))
         ALLOCATE(lupd6(ngrd))
         ALLOCATE(lupd7(ngrd))
         ALLOCATE(lupd8(ngrd))
         ALLOCATE(tt1(maxLevelSize))
         tt1(:) = zero
      ENDIF
!     CALL MPI_Bcast(levelPtr, nLevels+1, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(ijkv1, 4*ngrd, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(ijkv2, 4*ngrd, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(ijkv3, 4*ngrd, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(ijkv4, 4*ngrd, MPI_INT32_T, root, comm, mpierr) 
!     CALL MPI_Bcast(ijkv5, 4*ngrd, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(ijkv6, 4*ngrd, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(ijkv7, 4*ngrd, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(ijkv8, 4*ngrd, MPI_INT32_T, root, comm, mpierr)
!     CALL MPI_Bcast(lupd1, ngrd, MPI_C_BOOL, root, comm, mpierr)
!     CALL MPI_Bcast(lupd2, ngrd, MPI_C_BOOL, root, comm, mpierr)
!     CALL MPI_Bcast(lupd3, ngrd, MPI_C_BOOL, root, comm, mpierr)
!     CALL MPI_Bcast(lupd4, ngrd, MPI_C_BOOL, root, comm, mpierr)
!     CALL MPI_Bcast(lupd5, ngrd, MPI_C_BOOL, root, comm, mpierr)
!     CALL MPI_Bcast(lupd6, ngrd, MPI_C_BOOL, root, comm, mpierr)
!     CALL MPI_Bcast(lupd7, ngrd, MPI_C_BOOL, root, comm, mpierr)
!     CALL MPI_Bcast(lupd8, ngrd, MPI_C_BOOL, root, comm, mpierr)
      ! we're okay
!     CALL MPI_Bcast(lhaveGrid,        1, MPI_C_BOOL, root, comm, mpierr)
!     CALL MPI_Bcast(lhaveGridSpacing, 1, MPI_C_BOOL, root, comm, mpierr)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Broadcasts the spherical to cartesian parameter to the other processes
!>           on the communicator.
!>
!>    @param[in] root     Root process on which epsS2C is defined.
!>    @param[in] comm     MPI communicator on which to broadcast epsS2C.
!>
!>    @param[out] mpierr  MPI_SUCCESS indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_mpi_broadcastEpsS2C64fF(root, comm, mpierr)
      USE FTEIK_UTILS64F, ONLY : epsS2C
      USE MPI
      USE ISO_C_BINDING
      IMPLICIT NONE 
      INTEGER, INTENT(IN) :: root, comm
      INTEGER, INTENT(OUT) :: mpierr
      CALL MPI_Bcast(epsS2C, 1, MPI_DOUBLE, root, comm, mpierr)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Broadcasts the number of sweeps to perform in the iterative Gauss-Seidel
!>           method to the other processes on the communicator.
!>
!>
!>    @param[in] root     Root process on which nsweeps is defined.
!>    @param[in] comm     MPI communicator on which to broadcast nsweeps
!>
!>    @param[out] mpierr  MPI_SUCCESS indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_mpi_broadcastNumberOfSweepsF(root, comm, mpierr)
      USE FTEIK_UTILS64F, ONLY : nsweep
      USE MPI
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: root, comm
      INTEGER, INTENT(OUT) :: mpierr
      CALL MPI_Bcast(nsweep, 1, MPI_INT32_T, root, comm, mpierr)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE fteik_mpi_broadcastReceivers(root, comm, mpierr)
      USE FTEIK_UTILS64F, ONLY : nrec
      USE MPI
      USE ISO_C_BINDING
      !IMPLICIT NONE
      INTEGER, INTENT(IN) :: root, comm
      INTEGER, INTENT(OUT) :: mpierr
      INTEGER myid
      CALL MPI_Comm_Rank(myid, comm, mpierr)
      CALL MPI_Bcast(nrec, 1, MPI_INT32_T, root, comm, mpierr)
      IF (nrec < 1) THEN
         IF (myid == root) THEN
            WRITE(*,*) 'fteik_mpi_broadcastReceivers: No receivers to broadcast'
         ENDIF 
         RETURN
      ENDIF

      RETURN
      END
END MODULE
