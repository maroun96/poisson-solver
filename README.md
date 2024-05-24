# PETSc Poisson equation solver
This repo shows how to use PETSc to solve a 3D poisson equation with dirichlet or neumann boundary conditions

## Petsc installation
```conda install petsc -c conda-forge```

## Running the test case in parallel
```make clean_project```

```make all```

```mpirun -np 2 ./bin/main -ksp_type bcgs -ksp_rtol 1e-12 -pc_type asm```

