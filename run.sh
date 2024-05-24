make clean_project;
make all;

for K in 0 1 2 3 4 5 6
do
    time mpirun -np 2 ./bin/main -da_refine $K -ksp_type bcgs -ksp_rtol 1e-12 -pc_type asm
done;