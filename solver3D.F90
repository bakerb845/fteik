MODULE FTEIK_SOLVER64F
  USE FTEIK_CONSTANTS64F, ONLY : zero
  USE FTEIK_CONSTANTS64F, ONLY : FTEIK_NATURAL_ORDERING, &
                                 FTEIK_ZXY_ORDERING,     &
                                 FTEIK_XYZ_ORDERING,     &
                                 FTEIK_ZYX_ORDERING
  USE FTEIK_GRAPH3D, ONLY : ijkv1, ijkv2, ijkv3, ijkv4, &
                            ijkv5, ijkv6, ijkv7, ijkv8, &
                            levelPtr, nLevels, maxLevelSize
  USE ISO_C_BINDING
  IMPLICIT NONE
  !> Holds the travel-times (seconds).  This has dimension [ngrd].
  REAL(C_DOUBLE), PROTECTED, ALLOCATABLE, SAVE :: ttimes(:)
  !DIR$ ATTRIBUTES ALIGN: 64 :: ttimes
  !> Maps from the level'th level to the first node node in the level.
  !> This has dimension [nLevels+1].
! INTEGER(C_INT), PROTECTED, ALLOCATABLE, SAVE :: levelPtr(:)
! INTEGER(C_INT), PROTECTED, ALLOCATABLE, SAVE ::         &
!                 ijkv1(:), ijkv2(:), ijkv3(:), ijkv4(:), &
!                 ijkv5(:), ijkv6(:), ijkv7(:), ijkv8(:)
  LOGICAL(C_BOOL), PROTECTED, ALLOCATABLE, SAVE ::         &
                   lupd1(:), lupd2(:), lupd3(:), lupd4(:), &
                   lupd5(:), lupd6(:), lupd7(:), lupd8(:)
  LOGICAL(C_BOOL), PROTECTED, ALLOCATABLE, SAVE ::          &
                   lupdInit1(:), lupdInit2(:), &
                   lupdInit3(:), lupdInit4(:), &
                   lupdInit5(:), lupdInit6(:), &
                   lupdInit7(:), lupdInit8(:)
  !> Defines the transition from the spherical to the Cartesian solver solver
  !> during the initialization phase.  This has units of grid points.
  REAL(C_DOUBLE), PROTECTED, SAVE :: epsS2C = zero
  !> Defines the number of Gauss-Seidel iterations.
  INTEGER(C_INT), PROTECTED, SAVE :: nsweep = 0 
  !> Flag indicating whether or not the travel times were computed.
  LOGICAL(C_BOOL), PROTECTED, SAVE :: lhaveTimes = .FALSE. 
  !> Number of levels in level scheduling method.
! INTEGER(C_INT), PROTECTED, SAVE :: nLevels = 0 
  !> The number of nodes in the largest level.
! INTEGER(C_INT), PROTECTED, SAVE :: maxLevelSize = 0 
  !> This is for 3D 
  LOGICAL(C_BOOL), PARAMETER :: lis3d = .TRUE.
  !> Controls the verbosity
  INTEGER(C_INT), PROTECTED, SAVE :: verbose = 0
  ! Public/private subroutines
  PUBLIC :: fteik_solver3d_initialize64f
  PUBLIC :: fteik_solver3d_solveSourceLSM
  PUBLIC :: fteik_solver3d_free
  PUBLIC :: fteik_solver3d_setVelocityModel64f
  PUBLIC :: fteik_solver3d_setSources64f
  PUBLIC :: fteik_solver3d_getTravelTimeField64f
  PUBLIC :: fteik_solver3d_getTravelTimes64f
  PUBLIC :: fteik_solver3d_getNumberOfSources
  PUBLIC :: fteik_solver3d_getNumberOfReceivers
  PUBLIC :: fteik_solver3d_setNumberOfSweeps
  PUBLIC :: fteik_solver3d_setSphereToCartEpsilon
  PRIVATE :: fteik_solver3d_getSweepLimits
  PRIVATE :: fteik_solver3d_computeGraph
  PRIVATE :: fteik_solver3d_setInitialUpdateNodes
  PRIVATE :: fteik_solver3d_setUpdateNodes
  CONTAINS
!----------------------------------------------------------------------------------------!
!                                     Begin the Code                                     !
!----------------------------------------------------------------------------------------!
!
!>    @brief This program initializes the solver.  It will also perform the graph
!>           reordering.
!>
!>    @param[in] nzIn       Number of z grid points in model.  This must be at least 3.
!>    @param[in] nxIn       Number of x grid points in model.  This must be at least 3.
!>    @param[in] nyIn       Number of y grid points in model.  This must be at least 3.
!>    @param[in] z0In       z0 model origin (meters).
!>    @param[in] x0In       x0 model origin (meters).
!>    @param[in] y0In       y0 model origin (meters).
!>    @param[in] dzIn       Grid spacing (meters) in z.  This must be positive.
!>    @param[in] dxIn       Grid spacing (meters) in x.  This must be positive.
!>    @param[in] dyIn       Grid spacing (meters) in y.  This must be positive.
!>    @param[in] nsweepIn   This is the number of sweeps or iterations of the 
!>                          Gauss-Seidel method to be performed.  One is generally
!>                          sufficient.
!>    @param[in] epsIn      This is the radius, in number of grid points, around source
!>                          where the spherical approximation finite difference stencils
!>                          will be used.  This cannot be negative. 
!>    @param[in] verboseIn  Controsl the verbosity.  < 1 is quiet. 
!>
!>    @param[out] ierr      0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver3d_initialize64f(nzIn, nxIn, nyIn,    &
                                              z0In, x0In, y0In,    &
                                              dzIn, dxIn, dyIn,    &
                                              nsweepIn, epsIn,     &
                                              verboseIn, ierr)     &
                 BIND(C, NAME='fteik_solver3d_initialize64f')
      USE ISO_C_BINDING
      USE FTEIK_MODEL64F, ONLY : fteik_model_intializeGeometry, ngrd
      IMPLICIT NONE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x0In, y0In, z0In
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dxIn, dyIn, dzIn
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: epsIn
      INTEGER(C_INT), VALUE, INTENT(IN) :: nxIn, nyIn, nzIn
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsweepIn, verboseIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      WRITE(*,*) 'fteik_solver3d_initialize64f: Initializing...'
      CALL fteik_solver3d_free()
      CALL fteik_solver3d_setVerobosity(verboseIn)
      CALL fteik_model_intializeGeometry(lis3d,            &
                                         nzIn, nxIn, nyIn, &
                                         dzIn, dxIn, dyIn, &
                                         z0In, x0In, y0In, &
                                         ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver3d_initialize64f: Error setting grid geometry'
         RETURN
      ENDIF
      CALL fteik_solver3d_setNumberOfSweeps(nsweepIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver3d_initialize64f: Failed to set max number of sweeps'
         RETURN
      ENDIF
      CALL fteik_solver3d_setSphereToCartEpsilon(epsIn, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver3d_initialize64f: Failed to set epsilon'
         RETURN
      ENDIF
      CALL fteik_solver3d_computeGraph(ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver3d_initialize64f: Failed to compute graph'
         RETURN
      ENDIF
      ! Set space for the travel-times
      IF (.NOT.ALLOCATED(ttimes)) ALLOCATE(ttimes(ngrd))
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the verbosity on the module.
!>
!>    @param[in] verboseIn   Verbosity level to set.  Less than 1 is quiet.
!>
      SUBROUTINE fteik_solver3d_setVerobosity(verboseIn) &
      BIND(C, NAME='fteik_solver3d_setVerbosity')
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: verboseIn
      verbose = verboseIn
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the graph and pointers to be used by the level schedule solver.
!>
!>    @result 0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solver3d_computeGraph(ierr)
      USE ISO_C_BINDING
      USE FTEIK_MODEL64F, ONLY : ngrd, nx, ny, nz
      USE FTEIK_CONSTANTS64F, ONLY : FALSE, zero
      USE FTEIK_GRAPH3D, ONLY : fteik_graph3d_initialize, &
                                fteik_graph3d_makeLevelStructures
      !USE FTEIK_GRAPH
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: ierr
      !TYPE(GRAPH_TYPE) graph
      INTEGER(C_INT) ierrs(8)
      ierr = 0
      IF (nx < 1 .OR. ny < 1 .OR. nz < 1 .OR. ngrd < 1) THEN
         WRITE(*,*) 'fteik_solver3d_computeGraph: Error grid not yet set'
         ierr = 1
         RETURN
      ENDIF
      ! Compute the graph
      CALL fteik_graph3d_initialize(nz, nx, ny, verbose, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver3d_computeGraph: Error initializing graph'
         ierr = 1
         RETURN
      ENDIF
      CALL fteik_graph3d_makeLevelStructures(ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver3d_computeGraph: Error making level structures'
         ierr = 1
         RETURN
      ENDIF
!     CALL INIT(graph, nz, nx, ny, ierr)
!     IF (ierr /= 0) THEN
!        WRITE(*,*) 'fteik_solver3d_computeGraph: Error initializing graph'
!        ierr = 1
!        RETURN
!     ENDIF
!     nLevels = NUMBER_OF_LEVELS(graph)
!     IF (.NOT.ALLOCATED(ijkv1)) ALLOCATE(ijkv1(4*ngrd))
!     IF (.NOT.ALLOCATED(ijkv2)) ALLOCATE(ijkv2(4*ngrd))
!     IF (.NOT.ALLOCATED(ijkv3)) ALLOCATE(ijkv3(4*ngrd))
!     IF (.NOT.ALLOCATED(ijkv4)) ALLOCATE(ijkv4(4*ngrd))
!     IF (.NOT.ALLOCATED(ijkv5)) ALLOCATE(ijkv5(4*ngrd))
!     IF (.NOT.ALLOCATED(ijkv6)) ALLOCATE(ijkv6(4*ngrd))
!     IF (.NOT.ALLOCATED(ijkv7)) ALLOCATE(ijkv7(4*ngrd))
!     IF (.NOT.ALLOCATED(ijkv8)) ALLOCATE(ijkv8(4*ngrd))
!     ierr = IJKV(graph, 1, ijkv1)
!     ierr = IJKV(graph, 2, ijkv2) + ierr
!     ierr = IJKV(graph, 3, ijkv3) + ierr
!     ierr = IJKV(graph, 4, ijkv4) + ierr
!     ierr = IJKV(graph, 5, ijkv5) + ierr
!     ierr = IJKV(graph, 6, ijkv6) + ierr
!     ierr = IJKV(graph, 7, ijkv7) + ierr
!     ierr = IJKV(graph, 8, ijkv8) + ierr
!     IF (ierr /= 0) THEN
!        WRITE(*,*) 'fteik_solver3d_computeGraph: Errors encountered getting ijkv'
!        ierr = 1
!        RETURN
!     ENDIF
!     ierr = LEVEL_PTR(graph, 1, levelPtr)
!     IF (ierr /= 0) THEN
!        WRITE(*,*) 'fteik_solver3d_computeGraph: Error encountered getting levelPtr'
!        ierr = 1
!        RETURN
!     ENDIF
!     ! finished with the graph
!     CALL DELETE(graph)
!     maxLevelSize = 0
!     DO 1 i=1,nLevels
!        maxLevelSize = MAX(maxLevelSize, levelPtr(i+1) - levelPtr(i))
!   1 CONTINUE
!!    ALLOCATE(tt1(maxLevelSize))
!!    tt1(:) = zero
      ! Set the update grid.  This will mask nodes on the boundary.
      IF (.NOT.ALLOCATED(lupd1)) ALLOCATE(lupd1(ngrd))
      IF (.NOT.ALLOCATED(lupd2)) ALLOCATE(lupd2(ngrd))
      IF (.NOT.ALLOCATED(lupd3)) ALLOCATE(lupd3(ngrd))
      IF (.NOT.ALLOCATED(lupd4)) ALLOCATE(lupd4(ngrd))
      IF (.NOT.ALLOCATED(lupd5)) ALLOCATE(lupd5(ngrd))
      IF (.NOT.ALLOCATED(lupd6)) ALLOCATE(lupd6(ngrd))
      IF (.NOT.ALLOCATED(lupd7)) ALLOCATE(lupd7(ngrd))
      IF (.NOT.ALLOCATED(lupd8)) ALLOCATE(lupd8(ngrd))
      lupd1(:) = FALSE
      lupd2(:) = FALSE
      lupd3(:) = FALSE
      lupd4(:) = FALSE
      lupd5(:) = FALSE
      lupd6(:) = FALSE
      lupd7(:) = FALSE
      lupd8(:) = FALSE
      CALL fteik_solver3d_setUpdateNodes(1, nlevels, FALSE, 0, levelPtr, &
                                         ijkv1, lupd1, ierrs(1))
      CALL fteik_solver3d_setUpdateNodes(2, nlevels, FALSE, 0, levelPtr, &
                                         ijkv2, lupd2, ierrs(2))
      CALL fteik_solver3d_setUpdateNodes(3, nlevels, FALSE, 0, levelPtr, &
                                         ijkv3, lupd3, ierrs(3))
      CALL fteik_solver3d_setUpdateNodes(4, nlevels, FALSE, 0, levelPtr, &
                                         ijkv4, lupd4, ierrs(4))
      CALL fteik_solver3d_setUpdateNodes(5, nlevels, FALSE, 0, levelPtr, &
                                         ijkv5, lupd5, ierrs(5))
      CALL fteik_solver3d_setUpdateNodes(6, nlevels, FALSE, 0, levelPtr, &
                                         ijkv6, lupd6, ierrs(6))
      CALL fteik_solver3d_setUpdateNodes(7, nlevels, FALSE, 0, levelPtr, &
                                         ijkv7, lupd7, ierrs(7))
      CALL fteik_solver3d_setUpdateNodes(8, nlevels, FALSE, 0, levelPtr, &
                                         ijkv8, lupd8, ierrs(8))
      IF (MAXVAL(ABS(ierrs)) /= 0) THEN
         WRITE(*,*) 'fteik_solver3d_computeGraph: Error setting update nodes'
         ierr = 1
      ENDIF
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the update nodes for the initialization.  Because during initialization
!>           it may be possible that all nodes preceding the source location have no
!>           useful information to propagate forward they can be ignored.
!>
!>    @param[in] isrc     Source number.
!>
!>    @param[out] ierr    0 indicates.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solver3d_setInitialUpdateNodes(isrc, ierr)
      USE FTEIK_CONSTANTS64F, ONLY : TRUE, FALSE
      USE FTEIK_MODEL64F, ONLY : ngrd
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: isrc
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) ierrs(8)
      ! Set the initial update grid.  This will mask nodes on the boundary and
      ! points preceding the source location.
      IF (.NOT.ALLOCATED(lupdInit1)) ALLOCATE(lupdInit1(ngrd))
      IF (.NOT.ALLOCATED(lupdInit2)) ALLOCATE(lupdInit2(ngrd))
      IF (.NOT.ALLOCATED(lupdInit3)) ALLOCATE(lupdInit3(ngrd))
      IF (.NOT.ALLOCATED(lupdInit4)) ALLOCATE(lupdInit4(ngrd))
      IF (.NOT.ALLOCATED(lupdInit5)) ALLOCATE(lupdInit5(ngrd))
      IF (.NOT.ALLOCATED(lupdInit6)) ALLOCATE(lupdInit6(ngrd))
      IF (.NOT.ALLOCATED(lupdInit7)) ALLOCATE(lupdInit7(ngrd))
      IF (.NOT.ALLOCATED(lupdInit8)) ALLOCATE(lupdInit8(ngrd))
      lupdInit1(:) = FALSE
      lupdInit2(:) = FALSE
      lupdInit3(:) = FALSE
      lupdInit4(:) = FALSE
      lupdInit5(:) = FALSE
      lupdInit6(:) = FALSE
      lupdInit7(:) = FALSE
      lupdInit8(:) = FALSE
      CALL fteik_solver3d_setUpdateNodes(1, nlevels, TRUE, isrc, levelPtr, &
                                         ijkv1, lupdInit1, ierrs(1))
      CALL fteik_solver3d_setUpdateNodes(2, nlevels, TRUE, isrc, levelPtr, &
                                         ijkv2, lupdInit2, ierrs(2))
      CALL fteik_solver3d_setUpdateNodes(3, nlevels, TRUE, isrc, levelPtr, &
                                         ijkv3, lupdInit3, ierrs(3))
      CALL fteik_solver3d_setUpdateNodes(4, nlevels, TRUE, isrc, levelPtr, &
                                         ijkv4, lupdInit4, ierrs(4))
      CALL fteik_solver3d_setUpdateNodes(5, nlevels, TRUE, isrc, levelPtr, &
                                         ijkv5, lupdInit5, ierrs(5))
      CALL fteik_solver3d_setUpdateNodes(6, nlevels, TRUE, isrc, levelPtr, &
                                         ijkv6, lupdInit6, ierrs(6))
      CALL fteik_solver3d_setUpdateNodes(7, nlevels, TRUE, isrc, levelPtr, &
                                         ijkv7, lupdInit7, ierrs(7))
      CALL fteik_solver3d_setUpdateNodes(8, nlevels, TRUE, isrc, levelPtr, &
                                         ijkv8, lupdInit8, ierrs(8))
      IF (MAXVAL(ABS(ierrs)) /= 0) THEN
         WRITE(*,*) 'fteik_source3d_setInitialUpdateNodes: Error setting update nodes'
         ierr = 1 
         RETURN
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility routine for setting the nodes to update in the sweep.
!>
!>    @param[in] sweep      Sweep number.  This must be in the range [1,8].
!>    @param[in] nLevels    Number of levels in level set method.
!>    @param[in] linitk     If true then this sets the initialization nodes for the
!>                          given source. \n
!>                          Otherwise, this sets the initialization nodes for the
!>                          general solver.
!>    @param[in] isrc       Source number.
!>    @param[in] levelPtr   Maps from level'th level to start node in ijkv.
!>                          This is a vector of dimension [nLevels+1].
!>    @param[in] ijkv       Maps from in'th node in to the (iz,ix,iy,igrd) indices.
!>                          This is a vector of length [4*ngrd] where ngrd is the
!>                          number of grid points.
!>
!>    @param[in,out] lupd   On input this contains the nodes to update.  \n
!>                          On output this contains a new list of nodes to update based
!>                          on the sweep and potentially the source information.
!>
!>    @param[out] ierr      0 indicates success. 
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>    
      SUBROUTINE fteik_solver3d_setUpdateNodes(sweep, nLevels, linitk, &
                                               isrc, levelPtr, ijkv,   &
                                               lupd, ierr)
      USE ISO_C_BINDING
      USE FTEIK_MODEL64F, ONLY : nz, nx, ny
      USE FTEIK_SOURCE64F, ONLY : fteik_source_getSourceIndices32iF
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: sweep, nLevels
      LOGICAL(C_BOOL), VALUE, INTENT(IN) :: linitk
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      INTEGER(C_INT), DIMENSION(:), INTENT(IN) :: levelPtr, ijkv
      LOGICAL(C_BOOL), DIMENSION(:), INTENT(INOUT) :: lupd
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) i, ix, iy, iz, node, maxx, maxy, maxz, minx, miny, minz, &
                     x1, x2, xsi, y1, y2, ysi, z1, z2, zsi
      zsi =-1
      xsi =-1
      ysi =-1
      IF (linitk) THEN
         CALL fteik_source_getSourceIndices32iF(isrc, zsi, xsi, ysi, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'fteik_solver3d_setUpdateNodes: Error with source node init'
            ierr = 1
            RETURN
         ENDIF
      ENDIF
      CALL fteik_solver3d_getSweepLimits(sweep, linitk,          &
                                         nz, nx, ny,             &
                                         zsi, xsi, ysi,          &
                                         z1, z2, x1, x2, y1, y2, &
                                         ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver3d_setUpdateNodes: Invalid sweep', sweep
         ierr = 1
         RETURN
      ENDIF
      minz = MIN(z1, z2)
      maxz = MAX(z1, z2)
      minx = MIN(x1, x2)
      maxx = MAX(x1, x2)
      miny = MIN(y1, y2)
      maxy = MAX(y1, y2)
      DO 1 i=1,nLevels
         DO 2 node=levelPtr(i),levelPtr(i+1)-1 
            iz   = ijkv(4*(node - 1) + 1)
            ix   = ijkv(4*(node - 1) + 2)
            iy   = ijkv(4*(node - 1) + 3)
            IF (iz >= minz .AND. iz <= maxz .AND. &
                ix >= minx .AND. ix <= maxx .AND. &
                iy >= miny .AND. iy <= maxy) THEN
               lupd(node) = .TRUE.
            ENDIF 
    2    CONTINUE
    1 CONTINUE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Returns the sweep limits the given sweep.
!>
!>    @param[in] sweep    Sweep number.  This must be in the range [1,8].
!>    @param[in] linitk   If true then the sweep ranges are limited by the source
!>                        position.  This is useful on start-up when causality prevents
!>                        a chunk of nodes from being able to be updated.
!>                        Otherwise, these are the start and stop indices for the
!>                        entire grid.
!>    @param[in] nz       Number of z grid points in mesh.
!>    @param[in] nx       Number of x grid points in mesh.
!>    @param[in] ny       Number of y grid points in mesh.
!>    @param[in] zsi      z source grid point (Fortran numbered).
!>    @param[in] xsi      x source grid point (Fortran numbered).
!>    @param[in] ysi      y source grid point (Fortran numbered).
!>
!>    @param[out] z1      Start index in grid for loop on z (Fortran numbered).
!>    @param[out] z2      End index in grid for loop on z (Fortran numbered).
!>    @param[out] x1      Start index in grid for loop on x (Fortran numbered).
!>    @param[out] x2      End index in grid for loop on x (Fortran numbered).
!>    @param[out] y1      Start index in grid for loop on y (Fortran numbered).
!>    @param[out] y2      End index in grid for loop on y (Fortran numbered).
!>    @param[out] ierr    0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solver3d_getSweepLimits(sweep, linitk,          &
                                               nz, nx, ny,             &
                                               zsi, xsi, ysi,          &
                                               z1, z2, x1, x2, y1, y2, &
                                               ierr)
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: sweep, nz, nx, ny
      INTEGER(C_INT), VALUE, INTENT(IN) ::  zsi, xsi, ysi
      LOGICAL(C_BOOL), VALUE, INTENT(IN) :: linitk
      INTEGER(C_INT), INTENT(OUT) :: x1, x2, y1, y2, z1, z2, ierr
      ierr = 0
      IF (sweep == 1) THEN
         y1 = 2;  x1 = 2;  z1 = 2;
         y2 = ny; x2 = nx; z2 = nz;
         IF (linitk) THEN
            y1 = MAX(2, ysi); x1 = MAX(2, xsi); z1 = MAX(2, zsi);
            y2 = ny;          x2 = nx;          z2 = nz;
         ENDIF
      ELSE IF (sweep == 2) THEN
         y1 = 2;  x1 = nx - 1; z1 = 2; 
         y2 = ny; x2 = 1;      z2 = nz;
         IF (linitk) THEN
            y1 = MAX(2, ysi); x1 = xsi + 1;     z1 = MAX(2, zsi);
            y2 = ny;          x2 = 1;           z2 = nz;
         ENDIF
      ELSE IF (sweep == 3) THEN
         y1 = ny - 1; x1 = 2;  z1 = 2;
         y2 = 1;      x2 = nx; z2 = nz;
         IF (linitk) THEN
            y1 = ysi + 1;     x1 = MAX(2, xsi); z1 = MAX(2, zsi);
            y2 = 1;           x2 = nx;          z2 = nz;
         ENDIF
      ELSE IF (sweep == 4) THEN
         y1 = ny - 1; x1 = nx - 1; z1 = 2;
         y2 = 1;      x2 = 1;      z2 = nz;
         IF (linitk) THEN
            y1 = ysi + 1;     x1 = xsi + 1;     z1 = MAX(2, zsi);
            y2 = 1;           x2 = 1;           z2 = nz;
         ENDIF
      ELSE IF (sweep == 5) THEN
         y1 = 2;  x1 = 2;  z1 = nz - 1
         y2 = ny; x2 = nx; z2 = 1;
         IF (linitk) THEN 
            y1 = MAX(2, ysi); x1 = MAX(2, xsi); z1 = zsi + 1
            y2 = ny;          x2 = nx;          z2 = 1;
         ENDIF
      ELSE IF (sweep == 6) THEN
         y1 = 2;  x1 = nx - 1; z1 = nz - 1;
         y2 = ny; x2 = 1;      z2 = 1;
         IF (linitk) THEN 
            y1 = MAX(2, ysi); x1 = xsi + 1;     z1 = zsi + 1
            y2 = ny;          x2 = 1;           z2 = 1;
         ENDIF
      ELSE IF (sweep == 7) THEN
         y1 = ny - 1; x1 = 2;  z1 = nz - 1;
         y2 = 1;      x2 = nx; z2 = 1;
         IF (linitk) THEN 
            y1 = ysi + 1;     x1 = MAX(2, xsi); z1 = zsi + 1
            y2 = 1;           x2 = nx;          z2 = 1;
         ENDIF
      ELSE IF (sweep == 8) THEN
         y1 = ny - 1; x1 = nx - 1; z1 = nz - 1; 
         y2 = 1;      x2 = 1;      z2 = 1;
         IF (linitk) THEN 
            y1 = ysi + 1;     x1 = xsi + 1;     z1 = zsi + 1
            y2 = 1;           x2 = 1;           z2 = 1;
         ENDIF
      ELSE
         WRITE(*,*) 'fteik_solver3d_getSweepLimits: Invalid sweep', sweep
         ierr = 1
      ENDIF 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases all memory associated with the solver and resets all variables.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solver3d_free()             &
      BIND(C, NAME='fteik_solver3d_free')
      USE FTEIK_GRAPH3D, ONLY : fteik_graph3d_free
      USE FTEIK_RECEIVER64F, ONLY : fteik_receiver_free
      USE FTEIK_SOURCE64F, ONLY : fteik_source_free
      USE FTEIK_MODEL64F, ONLY : fteik_model_free
      USE FTEIK_CONSTANTS64F, ONLY : zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      lhaveTimes = .FALSE.
      epsS2C = zero
      nsweep = 0
      verbose = 0
      IF (ALLOCATED(ttimes)) DEALLOCATE(ttimes)
      IF (ALLOCATED(lupd1)) DEALLOCATE(lupd1)
      IF (ALLOCATED(lupd2)) DEALLOCATE(lupd2)
      IF (ALLOCATED(lupd3)) DEALLOCATE(lupd3)
      IF (ALLOCATED(lupd4)) DEALLOCATE(lupd4)
      IF (ALLOCATED(lupd5)) DEALLOCATE(lupd5)
      IF (ALLOCATED(lupd6)) DEALLOCATE(lupd6)
      IF (ALLOCATED(lupd7)) DEALLOCATE(lupd7)
      IF (ALLOCATED(lupd8)) DEALLOCATE(lupd8)
      IF (ALLOCATED(lupdInit1)) DEALLOCATE(lupdInit1)
      IF (ALLOCATED(lupdInit2)) DEALLOCATE(lupdInit2)
      IF (ALLOCATED(lupdInit3)) DEALLOCATE(lupdInit3)
      IF (ALLOCATED(lupdInit4)) DEALLOCATE(lupdInit4)
      IF (ALLOCATED(lupdInit5)) DEALLOCATE(lupdInit5)
      IF (ALLOCATED(lupdInit6)) DEALLOCATE(lupdInit6)
      IF (ALLOCATED(lupdInit7)) DEALLOCATE(lupdInit7)
      IF (ALLOCATED(lupdInit8)) DEALLOCATE(lupdInit8)
!     IF (ALLOCATED(ijkv1)) DEALLOCATE(ijkv1)
!     IF (ALLOCATED(ijkv2)) DEALLOCATE(ijkv2)
!     IF (ALLOCATED(ijkv3)) DEALLOCATE(ijkv3)
!     IF (ALLOCATED(ijkv4)) DEALLOCATE(ijkv4)
!     IF (ALLOCATED(ijkv5)) DEALLOCATE(ijkv5)
!     IF (ALLOCATED(ijkv6)) DEALLOCATE(ijkv6)
!     IF (ALLOCATED(ijkv7)) DEALLOCATE(ijkv7)
!     IF (ALLOCATED(ijkv8)) DEALLOCATE(ijkv8)
      CALL fteik_graph3d_free()
      CALL fteik_receiver_free()
      CALL fteik_source_free()
      CALL fteik_model_free()
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the spherical approximation transition tolerance around the source.
!>
!>    @param[in] epsIn   Radius in number of grid points around source where the
!>                       spherical approximation finite difference stencils will 
!>                       be used.  This is non-negative.
!>    @param[out] ierr   0 indicates success.
!>
!>    @author Mark Noble, Alexandrine Gesret, and Ben Baker
!>
!>    @version 2
!>
!>    @date July 2017
!>
!>    @copyright CeCILL-3
!>
!>    This is now a subroutine and incorporated with the fteik Fortran modules.
!>
      SUBROUTINE fteik_solver3d_setSphereToCartEpsilon(epsIn, ierr) &
      BIND(C, NAME='fteik_solver3d_setSphereToCartEpsilon')
      USE ISO_C_BINDING
      USE FTEIK_MODEL64F, ONLY : nz, nx, ny
      USE FTEIK_CONSTANTS64F, ONLY : zero
      IMPLICIT NONE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: epsIn
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      epsS2C = zero
      IF (nx < 1 .OR. ny < 1 .OR. nz < 1) THEN
         WRITE(*,*) 'fteik_solver3d_setSphereToCartEpsilon: Grid not yet set'
         ierr = 1
         RETURN
      ENDIF
      IF (INT(epsIn) > nz .OR. INT(epsIn) > nx .OR. INT(epsIn) > ny) THEN
         IF (INT(epsIn) > nz) THEN
            WRITE(*,*) 'fteik_solver3d_setSphereToCartEpsilon: eps bigger than nz', &
                       INT(epsIn), nz
         ENDIF
         IF (INT(epsIn) > nx) THEN
            WRITE(*,*) 'fteik_solver3d_setSphereToCartEpsilon: eps bigger than nx', &
                       INT(epsIn), nx
         ENDIF
         IF (INT(epsIn) > ny) THEN
            WRITE(*,*) 'fteik_solver3d_setSphereToCartEpsilon: eps bigger than ny', &
                       INT(epsIn), ny
         ENDIF
         ierr = 1
         RETURN
      ENDIF
      epsS2C = epsIn
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the number of Gauss-Seidel iterations  in the fast sweeping method.
!>
!>    @param[in] nsweepIn    Number of sweeps.  This cannot be negative and 1 is
!>                           usually sufficient.
!>
!>    @param[out] ierr       0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT 
!>
      SUBROUTINE fteik_solver3d_setNumberOfSweeps(nsweepIn, ierr) &
      BIND(C, NAME='fteik_solver3d_setNumberOfSweeps')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsweepIn 
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      nsweep = 0
      IF (nsweepIn < 0) THEN
         WRITE(*,*) 'fteik_solver3d_setNumberOfSweeps: nsweep must be positive', nsweep
         ierr = 1
         RETURN
      ENDIF
      nsweep = nsweepIn
      RETURN
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Initializes the source(s) on the solver.
!> 
!>    @param[in] nsrc     Number of sources.
!>    @param[in] zsrc     z locations (meters) of source.  This is a vector of dimension
!>                        [nsrc].
!>    @param[in] xsrc     x locations (meters) of source.  This is a vector of dimension
!>                        [nsrc].
!>    @param[in] ysrc     y locations (meters) of source.  This is a vector of dimension
!>                        [nsrc].
!>
!>    @param[out] ierr    0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solver3d_setSources64f(nsrc, zsrc, xsrc, ysrc, &
                                              ierr)                   &
      BIND(C, NAME='fteik_solver3d_setSources64f')
      USE FTEIK_SOURCE64F, ONLY : fteik_source_initialize64f
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsrc
      REAL(C_DOUBLE), INTENT(IN) :: zsrc(nsrc), xsrc(nsrc), ysrc(nsrc)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      CALL fteik_source_initialize64f(nsrc, zsrc, xsrc, ysrc, verbose, ierr)
      IF (ierr /= 0) WRITE(*,*) 'fteik_solver3d_setSources64f: Failed to set source'
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility routine to return the number of sources.
!>
!>    @param[out] nsrc   Number of sources.
!>    @param[out] ierr   0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver3d_getNumberOfSources(nsrc, ierr) & 
      BIND(C, NAME='fteik_solver3d_getNumberOfSources')
      USE FTEIK_SOURCE64F, ONLY : fteik_source_getNumberOfSources
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: nsrc, ierr
      CALL fteik_source_getNumberOfSources(nsrc, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver3d_getNumberOfSources: No sources initialized!'
         RETURN
      ENDIF 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Returns the travel time field.
!>
!>    @param[in] ngin    Number of input grid points.  This must equal [nz x nx x ny].
!>    @param[in] order   Desired ordering of output. \n
!>                       If order == FTEIK_NATURAL_ORDERING or FTEIK_ZXY_ORDERING then
!>                       ttimes will be [nz x nx ny] with first leading dimension nz
!>                       and second leading dimension nx. \n
!>                       If order == FTEIK_XYZ_ORDERING then ttimes will be [nx x ny x nz]
!>                       with first leading dimension nx and second leading
!>                       dimension ny. \n
!>                       If order == FTEIK_ZYX_ORDERING Then ttimes will be [nz x ny x nx]
!>                       with first leading dimension nz and second leading 
!>                       dimension ny. \n
!>
!>    @param[out] ttout  Travel time field at the grid points.  This is a [nz x nx x ny]
!>                       vector.  The ordering is defined by order.
!>    @param[out] ierr   0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver3d_getTravelTimeField64f(ngin, order, ttout, ierr) &
      BIND(C, NAME='fteik_solver3d_getTravelTimeField64f')
      USE FTEIK_MODEL64F, ONLY : ngrd, nx, ny, nz, nzx
      USE FTEIK_MODEL64F, ONLY : fteik_model_grid2indexF
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE, INTENT(IN) :: ngin, order
      REAL(C_DOUBLE), INTENT(OUT) :: ttout(ngin)
      INTEGER(C_INT), INTENT(OUT) :: ierr 
      INTEGER indx, ix, iy, iz, jndx
      ierr = 0
      IF (.NOT.lhaveTimes) THEN 
         WRITE(*,*) 'fteik_solver3d_getTravelTimeField64f: Travel times not yet computed'
         ierr = 1
         RETURN
      ENDIF
      IF (ngin /= ngrd) THEN 
         WRITE(*,*) 'fteik_solver3d_getTravelTimeField64f: No receivers'
         RETURN
      ENDIF
      IF (order == FTEIK_XYZ_ORDERING) THEN
         DO iz=1,nz
            DO iy=1,ny
               !$OMP SIMD
               DO ix=1,nx
                  indx = (iy - 1)*nzx + (ix - 1)*nz + iz 
                  !indx = fteik_model_grid2indexF(iz, ix, iy, nz, nzx)
                  jndx = (iz - 1)*nx*ny + (iy - 1)*nx + ix 
                  ttout(jndx) = ttimes(indx)
               ENDDO
            ENDDO
         ENDDO
      ELSEIF (order == FTEIK_ZYX_ORDERING) THEN
         DO iz=1,nz
            DO iy=1,ny
               !$OMP SIMD
               DO ix=1,nx
                  indx = (iy - 1)*nzx + (ix - 1)*nz + iz  
                  !indx = fteik_model_grid2indexF(iz, ix, iy, nz, nzx)
                  jndx = (ix - 1)*nz*ny + (iy - 1)*nz + iz
                  ttout(jndx) = ttimes(indx)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         IF (order /= FTEIK_NATURAL_ORDERING .AND. order /= FTEIK_ZXY_ORDERING) THEN
            WRITE(*,*) 'fteik_solver3d_getTravelTimeField64f: Defaulting to natural order'
         ENDIF
         ttout(:) = ttimes(:)
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Extracts the travel times at the receivers.
!>
!>    @param[in] nrec   Number of receivers.
!>
!>    @param[out] ttr   Travel times (seconds) at the receivers.  This is a vector of
!>                      dimension [nrec].
!>    @param[out] ierr  0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver3d_getTravelTimes64f(nrec, ttr, ierr) &
      BIND(C, NAME='fteik_solver3d_getTravelTimes64f')
      USE FTEIK_RECEIVER64F, ONLY : fteik_receiver_getTravelTimes64f
      USE FTEIK_MODEL64F, ONLY : ngrd
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nrec
      REAL(C_DOUBLE), INTENT(OUT) :: ttr(nrec)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      IF (nrec <= 0) THEN
         WRITE(*,*) 'fteik_solver3d_getTravelTimes64f: No receivers'
         RETURN
      ENDIF
      IF (.NOT.lhaveTimes) THEN
         WRITE(*,*) 'fteik_solver3d_getTravelTimes64f: Travel times not yet computed'
         ierr = 1
         RETURN
      ENDIF
      CALL fteik_receiver_getTravelTimes64f(nrec, ngrd, ttimes, ttr, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver3d_getTravelTimes64f: Error getting travel times'
         ierr = 1
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility routine to return the number of receivers.
!>
!>    @param[out] nrec   Number of receivers.
!>    @param[out] ierr   0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver3d_getNumberOfReceivers(nrec, ierr) &
      BIND(C, NAME='fteik_solver3d_getNumberOfReceivers')
      USE FTEIK_RECEIVER64F, ONLY : fteik_receiver_getNumberOfReceivers
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(OUT) :: nrec, ierr
      ierr = 0
      CALL fteik_receiver_getNumberOfReceivers(nrec, ierr)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the receivers in the model.
!>
!>    @param[in] nrec    Number of receivers to set.  
!>    @param[in] zrec    z locations (meters) of receivers.  This is a vector of
!>                       of dimension [nrec].
!>    @param[in] xrec    x locations (meters) of receivers.  This is a vector of
!>                       of dimension [nrec].
!>    @param[in] yrec    y locations (meters) of receivers.  This is a vector of
!>                       of dimension [nrec].
!>
!>    @param[out] ierr   0 indicates success.
!>
!>    @copyright Ben Baker distributed under the MIT license.
!>
      SUBROUTINE fteik_solver3d_setReceivers64f(nrec, zrec, xrec, yrec, &
                                                ierr)                   &
      BIND(C, NAME='fteik_solver3d_setReceivers64f')
      USE FTEIK_RECEIVER64F, ONLY : fteik_receiver_initialize64f
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nrec
      REAL(C_DOUBLE), INTENT(IN) :: zrec(nrec), xrec(nrec), yrec(nrec)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ! It's actually safe to have no receivers
      ierr = 0
      IF (nrec < 1) THEN
         WRITE(*,*) 'fteik_solver3d_setReceivers64f: No receivers to set'
         RETURN
      ENDIF
      CALL fteik_receiver_initialize64f(nrec, zrec, xrec, yrec, verbose, ierr)
      IF (ierr /= 0) WRITE(*,*) 'fteik_solver3d_setReceivers64f: Failed to set receivers'
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the velocity model on the solver.
!>
!>    @param[in] ncell    Number of cells in velocity model.  This should be 
!>                        (nz-1)*(nx-1)*(ny-1).
!>    @param[in] vel      Velocity model (meters/second) in model cells.  This is
!>                        a vector of dimension [ncell] whose fastest direction is
!>                        z and whose slowest direction is y.
!>
!>    @param[out] ierr    0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solver3d_setVelocityModel64f(ncell, vel, ierr) &
      BIND(C, NAME='fteik_solver3d_setVelocityModel64f')
      USE FTEIK_MODEL64F, ONLY : fteik_model_setVelocityModel64f
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: ncell
      REAL(C_DOUBLE), INTENT(IN) :: vel(ncell)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      CALL fteik_model_setVelocityModel64f(ncell, vel, ierr)
      IF (ierr /= 0) &
      WRITE(*,*) 'fteik_solver3d_setVelocityModel64f: Error setting velocity model'
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Solves the eikonal equation for the given source.
!>
!>    @param[in] isrc     Source number.  This must be in the range [1,nsrc].
!>
!>    @param[out] ierr    0 indicates success.
!>
!>    @author Ben Baker
!>
!>    @copyright MIT
!>
      SUBROUTINE fteik_solver3d_solveSourceLSM(isrc, ierr)   &
                 BIND(C, NAME='fteik_solver3d_solveSourceLSM')
      USE FTEIK_SOURCE64F, ONLY : nsrc
      USE FTEIK_MODEL64F, ONLY : lhaveModel
      USE FTEIK_LOCALSOLVER64F, ONLY : fteik_localSolver_initialize64fF
      USE FTEIK_CONSTANTS64F, ONLY : FTEIK_HUGE, TRUE, FALSE
      USE FTEIK_AUTOCODE, ONLY : fteik_evaluateSweep1LS64fF, &
                                 fteik_evaluateSweep2LS64fF, &
                                 fteik_evaluateSweep3LS64fF, &
                                 fteik_evaluateSweep4LS64fF, &
                                 fteik_evaluateSweep5LS64fF, &
                                 fteik_evaluateSweep6LS64fF, &
                                 fteik_evaluateSweep7LS64fF, &
                                 fteik_evaluateSweep8LS64fF
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: isrc
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_DOUBLE) ts8(8), t0, t1
      INTEGER(C_INT) dest(8), kiter
      ierr = 0
      lhaveTimes = FALSE
      IF (.NOT.lhaveModel) THEN
         WRITE(*,*) 'fteik_solver3d_solveSourceLSM: Model not yet set'
         ierr = 1
         GOTO 500
      ENDIF
      IF (isrc < 1 .OR. isrc > nsrc) THEN
         WRITE(*,*) 'fteik_solver3d_solveSourceLSM: Invalid source number:', isrc, 1, nsrc
         ierr = 1
         GOTO 500
      ENDIF
      ! Set the update nodes for the first iteration
      CALL fteik_solver3d_setInitialUpdateNodes(isrc, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver3d_solveSourceLSM: Error setting initial nodes'
         GOTO 500
      ENDIF
      ! Define the mesh constants and get the travel-times near the source
      CALL fteik_localSolver_initialize64fF(isrc, epsS2C, dest, ts8, ierr) 
      IF (ierr /= 0) THEN
         WRITE(*,*) 'fteik_solver3d_solveSourceLSM: Error setting mesh constants'
         GOTO 500
      ENDIF
      CALL CPU_TIME(t0)
      ttimes(:) = FTEIK_HUGE
      ttimes(dest(1:8)) = ts8(1:8) !+ t0(isrc) ! Add in initial time
      ! Solve the Eikonal Equation - first initialize
      CALL fteik_evaluateSweep1LS64fF(TRUE, ttimes, ierr)
      CALL fteik_evaluateSweep2LS64fF(TRUE, ttimes, ierr)
      CALL fteik_evaluateSweep3LS64fF(TRUE, ttimes, ierr)
      CALL fteik_evaluateSweep4LS64fF(TRUE, ttimes, ierr)
      CALL fteik_evaluateSweep5LS64fF(TRUE, ttimes, ierr)
      CALL fteik_evaluateSweep6LS64fF(TRUE, ttimes, ierr)
      CALL fteik_evaluateSweep7LS64fF(TRUE, ttimes, ierr)
      CALL fteik_evaluateSweep8LS64fF(TRUE, ttimes, ierr)
      ! Now perform the Gauss-Seidel iterations
      DO kiter=1,nsweep
         CALL fteik_evaluateSweep1LS64fF(FALSE, ttimes, ierr)
         CALL fteik_evaluateSweep2LS64fF(FALSE, ttimes, ierr)
         CALL fteik_evaluateSweep3LS64fF(FALSE, ttimes, ierr)
         CALL fteik_evaluateSweep4LS64fF(FALSE, ttimes, ierr)
         CALL fteik_evaluateSweep5LS64fF(FALSE, ttimes, ierr)
         CALL fteik_evaluateSweep6LS64fF(FALSE, ttimes, ierr)
         CALL fteik_evaluateSweep7LS64fF(FALSE, ttimes, ierr)
         CALL fteik_evaluateSweep8LS64fF(FALSE, ttimes, ierr)
      ENDDO 
      IF (verbose > 2) THEN
         CALL CPU_TIME(t1)
         WRITE(*,902) t1 - t0 
         !print *, minval(ttimes), maxval(ttimes)
      ENDIF
      lhaveTimes = .TRUE.
  500 CONTINUE
      IF (ierr /= 0) THEN
         ttimes(:) = FTEIK_HUGE
         lhaveTimes = .FALSE.
      ENDIF
  902 FORMAT(' fteik_solver3d_solveSourceLSM: Solver time in seconds=', F14.8)
      RETURN
      END
!----------------------------------------------------------------------------------------!
!                                       End the Code                                     !
!----------------------------------------------------------------------------------------!
END MODULE
