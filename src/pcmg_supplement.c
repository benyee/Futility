#include <petsc.h>
#include <petsc/private/pcmgimpl.h>
#include <petscsys.h>
#include <petscfix.h>
#include <petsc/private/fortranimpl.h>
#include <petscksp.h>

#ifdef PETSC_USE_POINTER_CONVERSION
extern void *PetscToPointer(void*);
extern int PetscFromPointer(void *); 
extern void PetscRmPointer(void*);
#else
#define PetscToPointer(a) (*(PetscFortranAddr *)(a))
#define PetscFromPointer(a) (PetscFortranAddr)(a)
#define PetscRmPointer(a)
#endif

#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pcmgsoftreset_ PCMGSOFTRESET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pcmgsoftreset_ pcmgsoftreset
#endif

#undef __FUNCT__
#define __FUNCT__ "PCMGSoftReset"

PetscErrorCode PCMGSoftReset(PC pc)
{
  PC_MG          *mg        = (PC_MG*)pc->data;
  PC_MG_Levels   **mglevels = mg->levels;
  PetscErrorCode ierr;
  PetscInt       i,n;

  PetscFunctionBegin;
  if (mglevels) {
    n = mglevels[0]->levels;
    for (i=0; i<n; i++) {
      ierr = MatDestroy(&mglevels[i]->A);CHKERRQ(ierr);
      if (mglevels[i]->smoothd != mglevels[i]->smoothu) {
        ierr = KSPReset(mglevels[i]->smoothd);CHKERRQ(ierr);
      }
      ierr = KSPReset(mglevels[i]->smoothu);CHKERRQ(ierr);
    }
  }
  //This may not work if DM's are used to define the MG structure.
  pc->setupcalled = PETSC_FALSE;
  PetscFunctionReturn(0);
}

PETSC_EXTERN void PETSC_STDCALL  pcmgsoftreset_(PC pc, int *__ierr ){
*__ierr = PCMGSoftReset(
  (PC)PetscToPointer((pc) ));
}
