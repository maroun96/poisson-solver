#include <petsc.h>

//Callback to initialize a vector
typedef PetscReal (*vecinit_fn)(PetscReal, PetscReal, PetscReal);

typedef enum
{
    DIRICHLET,
    NEUMANN
}boundary_condition_type;

typedef struct 
{
    boundary_condition_type bc_type;
    vecinit_fn bc_fn;
}boundary_condition;


typedef struct 
{
    PetscInt nx, ny, nz;
    PetscReal xmin, xmax;
    PetscReal ymin, ymax;
    PetscReal zmin, zmax;
    PetscReal dx, dy, dz;
    boundary_condition bc;
}Grid;

extern PetscErrorCode formExact(DM, Vec, Grid, vecinit_fn);
extern PetscErrorCode formMatrix(DM, Mat, Vec, Grid);
extern PetscErrorCode formRHS(DM, Vec, Grid, vecinit_fn);
extern PetscReal harmonicMean(PetscReal x, PetscReal y);