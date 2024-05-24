#include <stdio.h>
#include <petsc.h>
#include <petscviewerhdf5.h>
#include <poisson.h>


PetscReal exact_fn(PetscReal x, PetscReal y, PetscReal z){
    return cos(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);
}

PetscReal rhs_interior_fn(PetscReal x, PetscReal y, PetscReal z){
    return -12*pow(M_PI, 2)*cos(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);
}

PetscReal neumann_bc_fn(PetscReal x, PetscReal y, PetscReal z){
    return 0;
}

int main(int argc, char **argv)
{
    DM da;
    Mat A;
    Vec b, u, uexact, rho, beta;
    KSP ksp;
    PetscReal errnorm;
    DMDALocalInfo info;
    PetscViewer viewer;
    FILE *f;

    boundary_condition bc = {
        .bc_fn = &neumann_bc_fn,
        .bc_type = NEUMANN
    };


    Grid grid = {
        .nx=3,
        .ny=3,
        .nz=3,
        .xmin = -1.0,
        .xmax = 1.0,
        .ymin = -1.0,
        .ymax = 1.0,
        .zmin = -1.0,
        .zmax = 1.0,
        .bc = bc
    };

    PetscCall(PetscInitialize(&argc, &argv, PETSC_IGNORE, PETSC_IGNORE));

    PetscCall(DMDACreate3d(
        PETSC_COMM_WORLD,
        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
        DMDA_STENCIL_STAR, grid.nx, grid.ny, grid.nz,
        PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1,
        PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, &da
    ));

    PetscCall(DMSetFromOptions(da));
    PetscCall(DMSetUp(da));
    PetscCall(DMDASetUniformCoordinates(da, grid.xmin, grid.xmax, grid.ymin, grid.ymax, grid.zmin, grid.zmax));
    PetscCall(DMDAGetLocalInfo(da,&info));
    grid.dx = (grid.xmax-grid.xmin)/(info.mx-1);
    grid.dy = (grid.ymax-grid.ymin)/(info.my-1);
    grid.dz = (grid.zmax-grid.zmin)/(info.mz-1);
    

    PetscCall(DMCreateMatrix(da, &A));
    PetscCall(MatSetFromOptions(A));
    
    PetscCall(DMCreateGlobalVector(da, &b));
    PetscCall(DMCreateLocalVector(da, &rho)); PetscCall(VecDuplicate(rho, &beta));
    PetscCall(VecDuplicate(b, &u)); PetscCall(VecDuplicate(b, &uexact)); 
    PetscCall(PetscObjectSetName((PetscObject)uexact, "u_exact"));
    PetscCall(PetscObjectSetName((PetscObject)u, "u_sol"));
    PetscCall(formExact(da, uexact, grid, &exact_fn));
    PetscCall(VecSet(rho, 1)); PetscCall(VecSet(beta, 1));
    PetscCall(VecPointwiseDivide(beta, beta, rho));
    PetscCall(formMatrix(da, A, beta, grid));
    PetscCall(formRHS(da, b, grid, &rhs_interior_fn));

    if (grid.bc.bc_type == NEUMANN){
        MatNullSpace nullspace;
        PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace));
        PetscCall(MatSetNullSpace(A, nullspace));
        PetscCall(MatNullSpaceRemove(nullspace, b));
        PetscCall(MatNullSpaceDestroy(&nullspace));
    }
    
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, A, A));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPSolve(ksp, b, u));
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, "u.hdf5", FILE_MODE_WRITE, &viewer));
    PetscCall(VecView(u, viewer));

    PetscCall(VecAXPY(u,-1.0,uexact));    // u <- u + (-1.0) uxact
    PetscCall(VecNorm(u,NORM_INFINITY,&errnorm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                "on %d x %d x %d grid:  error |u-uexact|_inf = %g\n",
                info.mx,info.my,info.mz, errnorm));
    f = fopen("error_inf.dat", "a");
    PetscCall(PetscFPrintf(PETSC_COMM_WORLD, f, "%g\n", errnorm));
    fclose(f);

    PetscCall(DMDestroy(&da));
    PetscCall(MatDestroy(&A));
    PetscCall(KSPDestroy(&ksp));
    PetscCall(VecDestroy(&b)); PetscCall(VecDestroy(&u)); PetscCall(VecDestroy(&uexact)); PetscCall(VecDestroy(&beta));
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(PetscFinalize());
    return 0;
}