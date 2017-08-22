!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                          Futility Development Group                          !
!                             All rights reserved.                             !
!                                                                              !
! Futility is a jointly-maintained, open-source project between the University !
! of Michigan and Oak Ridge National Laboratory.  The copyright and license    !
! can be found in LICENSE.txt in the head directory of this repository.        !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module provides a linear system type and methods to solve systems
!> of equations via a multigrid method
!>
!> For valid reference lists
!> see @ref MatrixTypes::LinearSolverTypes_Declare_ValidParams
!> "LinearSolverTypes_Declare_ValidParams".
!>
!> Currently supported TPLs include:
!>  - PETSc (with interfaces to KSP)
!>
!> @par Module Dependencies
!>  - @ref IntrType "IntrType": @copybrief IntrType
!>  - @ref BLAS "BLAS": @copybrief BLAS
!>  - @ref Times "Times": @copybrief Times
!>  - @ref ExceptionHandler "ExceptionHandler": @copybrief ExceptionHandler
!>  - @ref Allocs "Allocs": @copybrief Allocs
!>  - @ref ParameterLists "ParameterLists": @copybrief ParameterLists
!>  - @ref ParallelEnv "ParallelEnv": @copybrief ParallelEnv
!>  - @ref VectorTypes "VectorTypes": @copybrief VectorTypes
!>  - @ref MatrixTypes "MatrixTypes": @copybrief MatrixTypes
!>  - @ref LinearSolverTypes "LinearSolverTypes": @copybrief LinearSolverTypes
!>
!> @par EXAMPLES
!> @code
!>
!> @endcode
!>
!> @author Ben C. Yee
!>   @date 08/22/2017
!>
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE LinearSolverTypes_Multigrid
  USE IntrType
  USE BLAS
  USE trilinos_interfaces
  USE Times
  USE ExceptionHandler
  USE Allocs
  USE ParameterLists
  USE ParallelEnv
  USE VectorTypes
  USE MatrixTypes
  USE PreconditionerTypes
  USE Strings
  USE IOUtil
  USE LinearSolverTypes
  IMPLICIT NONE

#ifdef FUTILITY_HAVE_PETSC
#include <finclude/petsc.h>
#include <petscversion.h>
!petscisdef.h defines the keyword IS, and it needs to be reset
#undef IS
#endif
!
! List of public members
  PUBLIC :: levelInfoForPetsc
  PUBLIC :: LinearSolverType_Multigrid
#ifdef FUTILITY_HAVE_PETSC
  PUBLIC :: interpolate_petsc
  PUBLIC :: restrict_petsc
#endif

  INTEGER(SIK),PARAMETER :: max_levels=8_SIK

  !> @brief The extended type for the Iterative Linear Solver
  TYPE,EXTENDS(LinearSolverType_Iterative) :: LinearSolverType_Multigrid
    !> Number of grids:
    INTEGER(SIK) :: nLevels=1_SIK
    !> Whether or not the restriciton, interpolation, and smoothing is ready:
    LOGICAL(SBK) :: isMultigridSetup=.FALSE.
    !> Size of each grid level_info(level,:) = (/num_eqns,nx,ny,nz/)
    INTEGER(SIK),ALLOCATABLE :: level_info(:,:)
#ifdef FUTILITY_HAVE_PETSC
    !> Array of pointers to petsc interpolation matrices
    Mat :: PetscInterpMats(max_levels)
#endif

    CONTAINS
      !> @copybrief TODO
      !> @copydetails TODO
      PROCEDURE,PASS :: init => init_LinearSolverType_Multigrid
      !> @copybrief TODO
      !> @copydetails TODO
      PROCEDURE,PASS :: setupInterpMats => setupInterpMats_LinearSolverType_Multigrid
      !> @copybrief TODO
      !> @copydetails TODO
      PROCEDURE,PASS :: clear => clear_LinearSolverType_Multigrid
  ENDTYPE LinearSolverType_Multigrid

  !> Logical flag to check whether the required and optional parameter lists
  !> have been created yet for the Linear Solver Type.
  LOGICAL(SBK),SAVE :: LinearSolverType_Paramsflag=.FALSE.

  !> Used by PETSc to get data from LinearSolverType_Multigrid
  INTEGER(SIK),ALLOCATABLE :: levelInfoForPetsc(:,:)

  !> Name of module
  CHARACTER(LEN=*),PARAMETER :: modName='LINEARSOLVERTYPES_MULTIGRID'
!
!===============================================================================
  CONTAINS
!
!-------------------------------------------------------------------------------
!> @brief Initializes the multigrid Linear Solver Type with a parameter list
!>
!> @param solver The linear solver to act on
!> @param Params the parameter list
!> @param A The A operator in Ax=b
!>
    SUBROUTINE init_LinearSolverType_Multigrid(solver,Params,A)
      CHARACTER(LEN=*),PARAMETER :: myName='init_LinearSolverType_Multigrid'
      CLASS(LinearSolverType_Multigrid),INTENT(INOUT) :: solver
      TYPE(ParamType),INTENT(IN) :: Params
      CLASS(MatrixType),POINTER,INTENT(INOUT),OPTIONAL :: A
      CLASS(ParamType),POINTER :: pListPtr
      TYPE(ParamType) :: validParams,matPList,vecxPList,vecbPList
      ! local variables
      INTEGER(SIK) :: n
      INTEGER(SIK) :: TPLType
      INTEGER(SIK) :: matType,matEngine
      INTEGER(SIK) :: MPI_Comm_ID,numberOMP
      CHARACTER(LEN=256) :: timerName
#ifdef FUTILITY_HAVE_PETSC
      KSP :: ksp_temp
      PetscErrorCode  :: iperr
#else
      CALL eLinearSolverType%raiseError(modName//"::"//myName//" - "// &
        "For now, LinearSolverType_Multigrid requires PETSc enabled.")
#endif
      !Check to set up required and optional param lists.
      IF(.NOT.LinearSolverType_Paramsflag) CALL LinearSolverType_Declare_ValidParams()

      !Validate against the reqParams and OptParams
      validParams=Params
      CALL validParams%validate(LinearSolverType_reqParams)

      !Pull LS data from the parameter list
      TPLType=-1
      MPI_Comm_ID=-1
      matType=-1
      matEngine=-1
      timerName=''
      numberOMP=-1
      CALL validParams%get('LinearSolverType->TPLType',TPLType)
      CALL validParams%get('LinearSolverType->MPI_Comm_ID',MPI_Comm_ID)
      CALL validParams%get('LinearSolverType->timerName',timerName)
      CALL validParams%get('LinearSolverType->matType',matType)
      CALL validParams%add('LinearSolverType->A->MatrixType->matType',matType)
      CALL validParams%get('LinearSolverType->numberOMP',numberOMP)
      ! pull data for matrix and vector parameter lists
      CALL validParams%get('LinearSolverType->A->MatrixType',pListPtr)
      matPList=pListPtr
      CALL validParams%get('LinearSolverType->x->VectorType',pListPtr)
      vecxPList=pListPtr
      CALL validParams%get('LinearSolverType->b->VectorType',pListPtr)
      vecbPList=pListPtr
      !add mpi communicator to parameter lists
      CALL matPList%add('MatrixType->MPI_Comm_ID',MPI_Comm_ID)
      CALL vecxPList%add('VectorType->MPI_Comm_ID',MPI_Comm_ID)
      CALL vecbPList%add('VectorType->MPI_Comm_ID',MPI_Comm_ID)
      !pull size from source vector
      CALL validParams%get('LinearSolverType->b->VectorType->n',n)

      CALL validParams%clear()

      !Initialize parallel environments based on input
      IF(MPI_Comm_ID /= -1) CALL solver%MPIparallelEnv%init(MPI_Comm_ID)
      IF(numberOMP > 0) CALL solver%OMPparallelEnv%init(numberOMP)

      IF(TPLType /= PETSC) THEN
        CALL eLinearSolverType%raiseError(modName//"::"//myName//" - "// &
          "For now, LinearSolverType_Multigrid only works with PETSC.")
      ENDIF

      IF(.NOT.solver%isInit) THEN
        solver%info=0
        IF(TPLType == PETSC) THEN
#ifdef FUTILITY_HAVE_PETSC
          solver%TPLType=PETSC
          matEngine=VM_PETSC

          !Should be irrelevant if using PETSc:
          solver%PCTypeName='NOPC'
          solver%pciters=0
          solver%pcsetup=0

          CALL matPList%add("MatrixType->engine",matEngine)
          ! allocate and initialize matrix (A)
          CALL MatrixFactory(solver%A, matPList)
          IF(PRESENT(A)) A=>solver%A

          CALL vecxPList%add('VectorType->engine',matEngine)
          CALL VectorFactory(solver%X, vecxPlist)
          CALL vecbPList%add('VectorType->engine',matEngine)
          CALL VectorFactory(solver%b, vecbPlist)

          solver%solverMethod=MULTIGRID

          CALL KSPCreate(solver%MPIparallelEnv%comm,solver%ksp,iperr)
          SELECTTYPE(A=>solver%A); TYPE IS(PETScMatrixType)
#if ((PETSC_VERSION_MAJOR>=3) && (PETSC_VERSION_MINOR>=5))
            CALL KSPSetOperators(solver%ksp,A%a,A%a,iperr)
#else
            CALL KSPSetOperators(solver%ksp,A%a,A%a, &
              DIFFERENT_NONZERO_PATTERN,iperr)
#endif
          ENDSELECT

          !KSPRICHARDSON+PCMG = Multigrid linear solver, not multigrid precon.
          CALL KSPSetType(solver%ksp,KSPRICHARDSON,iperr)
          CALL KSPGetPC(solver%ksp,solver%pc,iperr)
          CALL PCSetType(solver%pc,PCMG,iperr)

          !For now, only Galerkin coarse grid operators are supported.
          !  Galerkin means A_c = R*A*I
          CALL PCMGSetGalerkin(solver%pc,PETSC_TRUE,iperr)

          !ZZZZ why doesnt this work for multigrid?
          !CALL KSPSetInitialGuessNonzero(solver%ksp,PETSC_TRUE,iperr)
#endif
        ENDIF

        !assign values to solver
        CALL solver%SolveTime%setTimerName(timerName)
        solver%isInit=.TRUE.
        solver%isMultigridSetup=.FALSE.

      ELSE
        CALL eLinearSolverType%raiseError('Incorrect call to '// &
          modName//'::'//myName//' - LinearSolverType already initialized')
      ENDIF
      CALL vecbPList%clear()
      CALL vecxPList%clear()
      CALL matPList%clear()
    ENDSUBROUTINE init_LinearSolverType_Multigrid
!
!-------------------------------------------------------------------------------
!> @brief Setup the interpolation and restriction matrices.
!>
!> @param solver The linear solver to act on
!> @param Params the parameter list
!>
    SUBROUTINE setupInterpMats_LinearSolverType_Multigrid(solver,Params,weights)
      CHARACTER(LEN=*),PARAMETER :: myName='init_LinearSolverType_Multigrid'
      CLASS(LinearSolverType_Multigrid),INTENT(INOUT) :: solver
      TYPE(ParamType),INTENT(IN) :: Params
      INTEGER(SIK) :: n,n_old,i,iLevel,nnz,MPI_Comm_ID
      !Number of coupled equations (in neutronics, # of groups)
      INTEGER(SIK) :: num_eqns
      !Number of dimensions
      INTEGER(SIK) :: num_dims
      !  Right now, this only works for 3-,5-, and 7-point cnetered finite
      !    difference stencils
      INTEGER(SIK) :: nx,ny,nz,nx_old,ny_old,nz_old
      INTEGER(SIK) :: inx,iny,inz,inum_eqn

      REAL(SRK),INTENT(IN),OPTIONAL :: weights(:,:,:,:)
#ifdef FUTILITY_HAVE_PETSC
      KSP :: ksp_temp
      PC :: pc_temp
      PetscErrorCode  :: iperr
#endif

      IF(solver%isMultigridSetup) &
          CALL eLinearSolverType%raiseError(modName//"::"//myName//" - "// &
                 'Multigrid linear system is already setup!')

      !Determine the weight type:
      CALL Params%get('LinearSolverType->Multigrid->nx',nx)
      CALL Params%get('LinearSolverType->Multigrid->ny',ny)
      CALL Params%get('LinearSolverType->Multigrid->nz',nz)
      CALL Params%get('LinearSolverType->MPI_Comm_ID',MPI_Comm_ID)
      num_dims=1
      IF(ny > 1) num_dims=num_dims+1
      IF(nz > 1) num_dims=num_dims+1
      CALL Params%get('LinearSolverType->Multigrid->num_eqns',num_eqns)
      CALL Params%get('LinearSolverType->b->VectorType->n',n)
      IF(nx*ny*nz*num_eqns /= n) THEN
          CALL eLinearSolverType%raiseError(modName//"::"//myName//" - "// &
                 'number of unknowns (n) does not match provided '// &
                 'nx,ny,nz,num_eqns')
      ENDIF

      !Number of levels required to reduce down to ~5*num_eqns unknowns per processor:
      solver%nLevels=FLOOR(log(MAX(nx-1,ny-1,nz-1)/ &
            solver%MPIParallelEnv%nproc/5.0_SRK)/log(2.0_SRK))+1
      solver%nLevels=MIN(solver%nLevels,max_levels)

      IF(solver%nLevels < 2) &
          CALL eLinearSolverType%raiseDebug(modName//"::"//myName//" - "// &
                 'The grid is too small to coarsen, using multigrid with '// &
                 ' only 1 level!')

      IF(solver%TPLType == PETSC) THEN
#ifdef FUTILITY_HAVE_PETSC
        IF(ALLOCATED(levelInfoForPetsc)) &
          CALL eLinearSolverType%raiseError(modName//"::"//myName//" - "// &
                 'Only 1 PETSc LinearSolverType_Multigrid object permitted!')

        !Set # of levels:
        CALL PCMGSetLevels(solver%pc,solver%nLevels,PETSC_NULL_OBJECT,iperr) !TODO use some sort of mpi thing here?

        nx_old=nx
        ny_old=ny
        nz_old=nz
        ALLOCATE(solver%level_info(solver%nLevels,4))
        solver%level_info(solver%nLevels,:)=(/num_eqns,nx,ny,nz/)
        DO iLevel=solver%nLevels-1,1,-1
          !Set the smoother:
          CALL PCMGGetSmoother(solver%pc,iLevel,ksp_temp,iperr)
          CALL KSPSetType(ksp_temp,KSPPREONLY,iperr)
          CALL KSPGetPC(ksp_temp,pc_temp,iperr)
          !TODO use PCBJACOBI and set block size
          CALL PCSetType(pc_temp,PCSOR,iperr)

          !Create the interpolation operator:
          nx_old=nx
          ny_old=ny
          nz_old=nz
          n_old=n
          IF(nx > 1) nx=nx/2+1
          IF(ny > 1) ny=ny/2+1
          IF(nz > 1) nz=nz/2+1
          n=nx*ny*nz*num_eqns
          nnz=n
          IF(num_dims == 3) THEN
            nnz=nnz+((nx_old*ny_old*nz_old*num_eqns)-n)*6
          ELSEIF(num_dims == 2) THEN
            nnz=nnz+((nx_old*ny_old*nz_old*num_eqns)-n)*4 !ZZZZ what about the 2->1 part
          ELSE
            nnz=nnz+((nx_old*num_eqns)-n)*2 !ZZZZ what about the 2->1 part
          ENDIF

          IF(ny == 1 .AND. ny_old > 1) num_dims=num_dims-1
          IF(nz == 1 .AND. nz_old > 1) num_dims=num_dims-1

          !To be used by the interp/restrict functions:
          solver%level_info(iLevel,:)=(/num_eqns,nx,ny,nz/)
          IF(.NOT. PRESENT(weights)) THEN !uniform weights will be used
            !TODO determine nlocal for multiprocessor problems
            CALL MatCreateShell(solver%MPIparallelEnv%comm, &
                    n_old,n,n_old,n,PETSC_NULL_INTEGER, &
                    solver%PetscInterpMats(iLevel),iperr)
            CALL MatShellSetOperation(solver%PetscInterpMats(iLevel),MATOP_MULT, &
                    restrict_petsc,iperr);
            CALL MatShellSetOperation(solver%PetscInterpMats(iLevel), &
                    MATOP_MULT_TRANSPOSE_ADD,interpolate_petsc,iperr);
            CALL PCMGSetInterpolation(solver%pc,iLevel, &
                    solver%PetscInterpMats(iLevel),iperr);
            CALL PCMGSetRestriction(solver%pc,iLevel, &
                    solver%PetscInterpMats(iLevel),iperr);
          ELSE
            !TODO
            CALL eLinearSolverType%raiseError(modName//"::"//myName//" - "// &
              "Nonuniform multigrid weights not implemented yet.")
          ENDIF

        ENDDO
        !TODO determine coarsest smoother option
        !Set coarsest smoother options:
        CALL PCMGGetSmoother(solver%pc,0,ksp_temp,iperr)
        CALL KSPSetType(ksp_temp,KSPGMRES,iperr)
        CALL KSPSetInitialGuessNonzero(ksp_temp,PETSC_TRUE,iperr)

        ALLOCATE(levelInfoForPetsc(solver%nLevels,5))
        levelInfoForPetsc(:,1:4)=solver%level_info
        DO i=1,solver%nLevels
          levelInfoForPetsc(i,5)=PRODUCT(levelInfoForPetsc(i,1:4))
        ENDDO
#endif
      ENDIF
      solver%isMultigridSetup=.TRUE.

    ENDSUBROUTINE setupInterpMats_LinearSolverType_Multigrid
!
!-------------------------------------------------------------------------------
!> @brief Clears the Multigrid Linear Solver Type
!> @param solver The linear solver to act on
!>
!> This routine clears the data spaces for the iterative linear solver.
!>
    SUBROUTINE clear_LinearSolverType_Multigrid(solver)
      CLASS(LinearSolverType_Multigrid),INTENT(INOUT) :: solver

      INTEGER(SIK) :: iLevel

#ifdef FUTILITY_HAVE_PETSC
      PetscErrorCode :: iperr

      IF(solver%isMultigridSetup) THEN
        DO iLevel=1,solver%nLevels-1
          !CALL MatDestroy(solver%PetscInterpMats(iLevel),iperr) !ZZZZ
        ENDDO
      ENDIF
#endif

      solver%isMultigridSetup=.FALSE.
      IF(ALLOCATED(solver%level_info)) DEALLOCATE(solver%level_info)
      IF(ALLOCATED(levelInfoForPetsc)) DEALLOCATE(levelInfoForPetsc)
      solver%nLevels=1_SIK

      CALL solver%LinearSolverType_Iterative%clear()

    ENDSUBROUTINE clear_LinearSolverType_Multigrid

#ifdef FUTILITY_HAVE_PETSC
! TODO
    SUBROUTINE interpolate_petsc(mat,xx,yy,zz,iperr)
      Mat :: mat
      Vec :: xx
      Vec :: yy
      Vec :: zz
      PetscErrorCode :: iperr

      INTEGER(SIK) :: probsize
      INTEGER(SIK) :: iLevel,nLevels
      INTEGER(SIK) :: num_eqns,nx,ny,nz,nx_c,ny_c,nz_c
      INTEGER(SIK) :: inum_eqns,inx,iny,inz
      INTEGER(SIK) :: inx_f,iny_f,inz_f

      REAL(SRK),POINTER :: myxx(:),myyy(:),myzz(:)

      CALL VecGetSize(xx,probsize,iperr)

      nLevels=SIZE(levelInfoForPetsc,1)
      DO iLevel=1,nLevels
        IF(levelInfoForPetsc(iLevel,5)==probsize) EXIT
      ENDDO

      num_eqns=levelInfoForPetsc(iLevel,1)
      nx=levelInfoForPetsc(iLevel,2)
      ny=levelInfoForPetsc(iLevel,3)
      nz=levelInfoForPetsc(iLevel,4)
      nx_c=levelInfoForPetsc(iLevel-1,2)
      ny_c=levelInfoForPetsc(iLevel-1,3)
      nz_c=levelInfoForPetsc(iLevel-1,4)

      CALL VecGetArrayF90(xx,myxx,iperr)
      CALL VecGetArrayF90(yy,myyy,iperr)
      CALL VecGetArrayF90(zz,myzz,iperr)

      !TODO expand to higher num_eqns, ny, nz
      !TODO account for odd grids
      DO inx=1,nx
        IF(XOR(MOD(inx,2)==1, (inx>nx/2 .AND. MOD(nx,2) == 0))) THEN
          myzz(inx)=myxx((inx+1)/2-1)+myyy(inx)
        ELSE
          myzz(inx)=0.5*myxx(inx/2-1)+0.5*myxx(inx/2)+myyy(inx)
        ENDIF
      ENDDO

      CALL VecRestoreArrayF90(xx,myxx,iperr)
      CALL VecRestoreArrayF90(yy,myyy,iperr)
      CALL VecRestoreArrayF90(zz,myzz,iperr)

    ENDSUBROUTINE
! TODO
    SUBROUTINE restrict_petsc(mat,rr,bb,iperr)
      Mat :: mat
      Vec :: rr
      Vec :: bb
      PetscErrorCode :: iperr

      INTEGER(SIK) :: probsize
      INTEGER(SIK) :: iLevel,nLevels
      INTEGER(SIK) :: num_eqns,nx,ny,nz,nx_c,ny_c,nz_c
      INTEGER(SIK) :: inum_eqns,inx,iny,inz
      INTEGER(SIK) :: inx_f,iny_f,inz_f

      REAL(SRK),POINTER :: myrr(:),mybb(:)

      CALL VecGetSize(rr,probsize,iperr)

      nLevels=SIZE(levelInfoForPetsc,1)
      DO iLevel=1,nLevels
        IF(levelInfoForPetsc(iLevel,5)==probsize) EXIT
      ENDDO

      num_eqns=levelInfoForPetsc(iLevel,1)
      nx=levelInfoForPetsc(iLevel,2)
      ny=levelInfoForPetsc(iLevel,3)
      nz=levelInfoForPetsc(iLevel,4)
      nx_c=levelInfoForPetsc(iLevel-1,2)
      ny_c=levelInfoForPetsc(iLevel-1,3)
      nz_c=levelInfoForPetsc(iLevel-1,4)

      CALL VecGetArrayF90(rr,myrr,iperr)
      CALL VecGetArrayF90(bb,mybb,iperr)

      !TODO expand to higher num_eqns, ny, nz
      !TODO account for odd grids
      DO inx=1,nx_c
        inx_f=inx*2
        mybb(inx)=myrr(inx_f-1)
        IF(inx > 1) mybb(inx)=mybb(inx)+0.5*myrr(inx_f-2)
        IF(inx < nx_c) mybb(inx)=mybb(inx)+0.5*myrr(inx_f)
      ENDDO

      CALL VecRestoreArrayF90(rr,myrr,iperr)
      CALL VecRestoreArrayF90(bb,mybb,iperr)

      !TODO

    ENDSUBROUTINE
#endif

ENDMODULE LinearSolverTypes_Multigrid
