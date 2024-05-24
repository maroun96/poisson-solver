#include <petsc.h>
#include <poisson.h>

PetscReal harmonicMean(PetscReal x, PetscReal y){
    return (2*x*y)/(x+y);
}

static void formMatrixAtBoundary(Grid grid, DMDALocalInfo info, PetscReal *v, MatStencil *col, PetscInt *ncols, PetscInt i, PetscInt j, PetscInt k){
    switch (grid.bc.bc_type)
    {
    case DIRICHLET:
        col[*ncols].i = i; 
        col[*ncols].j = j; 
        col[*ncols].k = k; 
        v[(*ncols)++] = 1.0;
        break;
    case NEUMANN:
        if (i==0){
            col[*ncols].i = i; col[*ncols].j = j; col[*ncols].k = k; v[(*ncols)++] = -1;
            col[*ncols].i = i+1; col[*ncols].j = j; col[*ncols].k = k; v[(*ncols)++] = 1;
        }
        else if (i==info.mx-1){
            col[*ncols].i = i; col[*ncols].j = j; col[*ncols].k = k; v[(*ncols)++] = 1;
            col[*ncols].i = i-1; col[*ncols].j = j; col[*ncols].k = k; v[(*ncols)++] = -1;
        }
        else if (j==0){
            col[*ncols].i = i; col[*ncols].j = j; col[*ncols].k = k; v[(*ncols)++] = -1;
            col[*ncols].i = i; col[*ncols].j = j+1; col[*ncols].k = k; v[(*ncols)++] = 1;
        }
        else if (j==info.my-1){
            col[*ncols].i = i; col[*ncols].j = j; col[*ncols].k = k; v[(*ncols)++] = 1;
            col[*ncols].i = i; col[*ncols].j = j-1; col[*ncols].k = k; v[(*ncols)++] = -1;
        }
        else if (k==0){
            col[*ncols].i = i; col[*ncols].j = j; col[*ncols].k = k; v[(*ncols)++] = -1;
            col[*ncols].i = i; col[*ncols].j = j; col[*ncols].k = k+1; v[(*ncols)++] = 1;
        }
        else if (k==info.mz-1){
            col[*ncols].i = i; col[*ncols].j = j; col[*ncols].k = k; v[(*ncols)++] = 1;
            col[*ncols].i = i; col[*ncols].j = j; col[*ncols].k = k-1; v[(*ncols)++] = -1;
        }
        break;
    default:
        break;
    }
}

PetscErrorCode formMatrix(DM da, Mat A, Vec beta, Grid grid){
    DMDALocalInfo info;
    MatStencil row, col[7];
    PetscReal v[7];
    PetscReal beta_interp, ***abeta;
    PetscInt i,j,k,ncols;
    PetscReal xfact = 1/pow(grid.dx, 2);
    PetscReal yfact = 1/pow(grid.dy, 2);
    PetscReal zfact = 1/pow(grid.dz, 2);
    
    PetscCall(DMDAGetLocalInfo(da, &info));
    PetscCall(DMDAVecGetArray(da, beta, &abeta));

    for (k=info.zs; k<info.zs+info.zm; k++){
        for (j=info.ys; j<info.ys+info.ym; j++){
            for (i=info.xs; i<info.xs+info.xm; i++){
                row.i = i; row.j = j; row.k = k;
                ncols = 0;
                if (i==0 || i==info.mx-1 || j==0 || j==info.my-1 || k==0 || k==info.mz-1){
                    formMatrixAtBoundary(grid, info, v, col, &ncols, i, j, k);        
                }
                else{
                    col[ncols].k = k; col[ncols].j = j; col[ncols].i = i+1;
                    beta_interp = harmonicMean(abeta[k][j][i], abeta[k][j][i+1]);
                    v[ncols++] = xfact*beta_interp;

                    col[ncols].k = k; col[ncols].j = j; col[ncols].i = i-1;
                    beta_interp = harmonicMean(abeta[k][j][i-1], abeta[k][j][i]);
                    v[ncols++] = xfact*beta_interp;

                    col[ncols].k = k; col[ncols].j = j+1; col[ncols].i = i;
                    beta_interp = harmonicMean(abeta[k][j][i], abeta[k][j+1][i]);
                    v[ncols++] = yfact*beta_interp;

                    col[ncols].k = k; col[ncols].j = j-1; col[ncols].i = i;
                    beta_interp = harmonicMean(abeta[k][j-1][i], abeta[k][j][i]);
                    v[ncols++] = yfact*beta_interp;

                    col[ncols].k = k+1; col[ncols].j = j; col[ncols].i = i;
                    beta_interp = harmonicMean(abeta[k][j][i], abeta[k+1][j][i]);
                    v[ncols++] = zfact*beta_interp;

                    col[ncols].k = k-1; col[ncols].j = j; col[ncols].i = i;
                    beta_interp = harmonicMean(abeta[k-1][j][i], abeta[k][j][i]);
                    v[ncols++] = zfact*beta_interp;

                    col[ncols].k = k; col[ncols].j = j; col[ncols].i = i;
                    v[ncols] = 0;
                    for (int h=0; h<6; h++){
                        v[ncols] -= v[h];
                    }
                    ncols++;
                }
                PetscCall(MatSetValuesStencil(A, 1, &row, ncols, col, v, INSERT_VALUES));
            }
        }
    }

    PetscCall(DMDAVecRestoreArray(da, beta, &abeta));

    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

    return 0;
}

PetscErrorCode formRHS(DM da, Vec b, Grid grid, vecinit_fn rhs_interior_fn){
    PetscInt i,j,k;
    PetscReal x, y, z, ***ab;
    DMDALocalInfo info;

    PetscCall(DMDAGetLocalInfo(da, &info));
    PetscCall(DMDAVecGetArray(da, b, &ab));

    for (k=info.zs; k<info.zs+info.zm; k++){
        z = grid.zmin + k*grid.dz;
        for (j=info.ys; j<info.ys+info.ym; j++){
            y = grid.ymin + j*grid.dy;
            for (i=info.xs; i<info.xs+info.xm; i++){
                x = grid.xmin + i*grid.dx;
                if (i==0 || i==info.mx-1 || j==0 || j==info.my-1 || k==0 || k==info.mz-1){
                    ab[k][j][i] = grid.bc.bc_fn(x, y, z);       
                }
                else{
                    ab[k][j][i] = rhs_interior_fn(x, y, z);
                }
            }
        }
    }
    PetscCall(DMDAVecGetArray(da, b, &ab));
    return 0;

}

PetscErrorCode formExact(DM da, Vec uexact, Grid grid, vecinit_fn fn){
    PetscInt i,j,k;
    PetscReal x, y, z, ***auexact;
    DMDALocalInfo info;

    PetscCall(DMDAGetLocalInfo(da, &info));
    PetscCall(DMDAVecGetArray(da, uexact, &auexact));
    for (k=info.zs; k<info.zs+info.zm; k++){
        z = grid.zmin + k*grid.dz;
        for (j=info.ys; j<info.ys+info.ym; j++){
            y = grid.ymin + j*grid.dy;
            for (i=info.xs; i<info.xs+info.xm; i++){
                x = grid.xmin + i*grid.dx;
                auexact[k][j][i] = fn(x, y, z);
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(da, uexact, &auexact));
    return 0;
}