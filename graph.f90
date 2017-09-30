      MODULE FTEIK_GRAPH
         USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_NULL_PTR
         IMPLICIT NONE
         PRIVATE
         TYPE GRAPH_TYPE
            PRIVATE
            TYPE(C_PTR) :: object = C_NULL_PTR
         END TYPE GRAPH_TYPE
         !-------------------------------------------------------------------------------!
         !                                C Wrappers                                     !
         !-------------------------------------------------------------------------------!
         INTERFACE
            FUNCTION C_GRAPH_INIT(nz, nx, ny, ierr) &
                     BIND(C, NAME='fteik_graph_initializeF') &
                     RESULT(this)
            USE ISO_C_BINDING
            IMPORT
            TYPE(C_PTR) :: this
            INTEGER(C_INT), INTENT(IN) :: nz, nx, ny
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END FUNCTION C_GRAPH_INIT

            FUNCTION C_GRAPH_GET_NLEVELS(this) &
            BIND(C, NAME='fteik_graph_getNumberOfLevels') &
            RESULT(nLevels)
            USE ISO_C_BINDING
            IMPORT
            IMPLICIT NONE
            TYPE(C_PTR), VALUE :: this
            INTEGER(C_INT) :: nLevels
            END FUNCTION C_GRAPH_GET_NLEVELS 

            FUNCTION C_GRAPH_GETIJKV(this, nwork, sweep, ijkv) &
                     BIND(C, NAME='fteik_graph_getIJKVF') &
                     RESULT(ierr)
            USE ISO_C_BINDING
            IMPORT
            TYPE(C_PTR), VALUE :: this
            INTEGER(C_INT), INTENT(IN) :: nwork, sweep
            INTEGER(C_INT), INTENT(OUT) :: ijkv(nwork)
            INTEGER(C_INT) :: ierr
            END FUNCTION C_GRAPH_GETIJKV

            FUNCTION C_GRAPH_GETLEVELPTR(this, nwork, sweep, levelPtr) &
            BIND(C, NAME='fteik_graph_getLevelPointerF') &
            RESULT(ierr)
            USE ISO_C_BINDING
            IMPORT
            TYPE(C_PTR), VALUE :: this
            INTEGER(C_INT), INTENT(IN) :: nwork, sweep
            INTEGER(C_INT), INTENT(OUT) :: levelPtr(nwork)
            INTEGER(C_INT) :: ierr
            END FUNCTION C_GRAPH_GETLEVELPTR


            FUNCTION C_GRAPH_GET_NGRDPTS(this) &
                     BIND(C, NAME='fteik_graph_getNumberOfGridPoints') &
                     RESULT(ngrd)
            USE ISO_C_BINDING
            IMPORT
            TYPE(C_PTR), VALUE :: this
            INTEGER(C_INT) :: ngrd
            END FUNCTION C_GRAPH_GET_NGRDPTS

            SUBROUTINE C_GRAPH_DELETE(this) &
                       BIND(C, NAME='fteik_graph_finalize')
            USE ISO_C_BINDING
            IMPORT
            TYPE(C_PTR), VALUE :: this
            END SUBROUTINE C_GRAPH_DELETE
         END INTERFACE
         !-------------------------------------------------------------------------------!
         !                                 Procedures                                    !
         !-------------------------------------------------------------------------------!
         PUBLIC :: INIT, NUMBER_OF_LEVELS, NUMBER_OF_GRDPTS, &
                   IJKV, LEVEL_PTR, DELETE, GRAPH_TYPE
         INTERFACE INIT
            MODULE PROCEDURE FTEIK_GRAPH_INIT 
         END INTERFACE INIT 
         INTERFACE DELETE
            MODULE PROCEDURE FTEIK_GRAPH_DELETE
         END INTERFACE DELETE
         INTERFACE IJKV
            MODULE PROCEDURE FTEIK_GRAPH_GETIJKV
         END INTERFACE IJKV 
         INTERFACE LEVEL_PTR
            MODULE PROCEDURE FTEIK_GRAPH_GETLEVELPTR
         END INTERFACE LEVEL_PTR
         INTERFACE NUMBER_OF_LEVELS
            MODULE PROCEDURE FTEIK_GRAPH_GET_NLEVELS
         END INTERFACE NUMBER_OF_LEVELS
         INTERFACE NUMBER_OF_GRDPTS
            MODULE PROCEDURE FTEIK_GRAPH_GET_NGRDPTS
         END INTERFACE NUMBER_OF_GRDPTS
!           FUNCTION TYPE(C_PTR) FTEIK_GRAPH_INITIALIZE(nz, nx, ny, ierr) &
!                                BIND(C, NAME='fteik_graphInitializeF')
!           USE ISO_C_BINDING
!           INTEGER(C_INT), INTENT(IN) :: nx, ny, nz
!           INTEGER(C_INT), INTENT(OUT) :: ierr
!           END FUNCTION
!        END INTERFACE
         !-------------------------------------------------------------------------------!
         !                             Actual subroutines                                !
         !-------------------------------------------------------------------------------!
         CONTAINS
         SUBROUTINE FTEIK_GRAPH_INIT(this, nz, nx, ny, ierr)
         USE ISO_C_BINDING
         IMPLICIT NONE
         TYPE(GRAPH_TYPE), INTENT(OUT) :: this
         INTEGER(C_INT), INTENT(IN) :: nz, nx, ny
         INTEGER(C_INT), INTENT(OUT) :: ierr
         this%object = C_GRAPH_INIT(nz, nx, ny, ierr)
         RETURN
         END SUBROUTINE FTEIK_GRAPH_INIT

         INTEGER(C_INT) FUNCTION FTEIK_GRAPH_GET_NLEVELS(this)
         USE ISO_C_BINDING
         IMPLICIT NONE
         TYPE(GRAPH_TYPE), INTENT(INOUT) :: this
         FTEIK_GRAPH_GET_NLEVELS = C_GRAPH_GET_NLEVELS(this%object)
         RETURN
         END FUNCTION FTEIK_GRAPH_GET_NLEVELS

         INTEGER(C_INT) FUNCTION FTEIK_GRAPH_GET_NGRDPTS(this)
         USE ISO_C_BINDING
         IMPLICIT NONE
         TYPE(GRAPH_TYPE), INTENT(INOUT) :: this 
         FTEIK_GRAPH_GET_NGRDPTS = C_GRAPH_GET_NGRDPTS(this%object)
         RETURN
         END FUNCTION FTEIK_GRAPH_GET_NGRDPTS 

         INTEGER(C_INT) FUNCTION FTEIK_GRAPH_GETIJKV(this, sweep, ijkv)
         USE ISO_C_BINDING
         IMPLICIT NONE
         TYPE(GRAPH_TYPE), INTENT(INOUT) :: this
         INTEGER(C_INT), ALLOCATABLE, INTENT(INOUT) :: ijkv(:)
         INTEGER(C_INT), INTENT(IN) :: sweep
         INTEGER(C_INT) ngrd, nwork
         ngrd = C_GRAPH_GET_NGRDPTS(this%object)
         FTEIK_GRAPH_GETIJKV = 0
         IF (ALLOCATED(ijkv)) DEALLOCATE(ijkv)
         IF (ngrd < 1) THEN
            WRITE(*,*) 'fteik_grah_getijkv: No grid points'
            FTEIK_GRAPH_GETIJKV =-1
            RETURN
         ENDIF
         nwork = 4*ngrd
         ALLOCATE(ijkv(nwork))
         FTEIK_GRAPH_GETIJKV = C_GRAPH_GETIJKV(this%object, nwork, sweep, ijkv)
         RETURN
         END FUNCTION FTEIK_GRAPH_GETIJKV

         INTEGER(C_INT) FUNCTION FTEIK_GRAPH_GETLEVELPTR(this, sweep, levelPtr)
         USE ISO_C_BINDING
         IMPLICIT NONE
         TYPE(GRAPH_TYPE), INTENT(INOUT) :: this
         INTEGER(C_INT), ALLOCATABLE, INTENT(INOUT) :: levelPtr(:)
         INTEGER(C_INT), INTENT(IN) :: sweep
         INTEGER(C_INT) nLevels, nwork
         nLevels = C_GRAPH_GET_NLEVELS(this%object)
         FTEIK_GRAPH_GETLEVELPTR = 0 
         IF (ALLOCATED(levelPtr)) DEALLOCATE(levelPtr)
         nwork = nLevels + 1 
         ALLOCATE(levelPtr(nwork))
         FTEIK_GRAPH_GETLEVELPTR &
            = C_GRAPH_GETLEVELPTR(this%object, nwork, sweep, levelPtr)
         RETURN
         END FUNCTION FTEIK_GRAPH_GETLEVELPTR

         SUBROUTINE FTEIK_GRAPH_DELETE(this)
         USE ISO_C_BINDING
         TYPE(GRAPH_TYPE), INTENT(INOUT) :: this
         CALL C_GRAPH_DELETE(this%object)
         this%object = C_NULL_PTR
         RETURN
         END SUBROUTINE FTEIK_GRAPH_DELETE

      END MODULE FTEIK_GRAPH
