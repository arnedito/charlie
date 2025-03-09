#!/bin/bash
echo "Compiling mod_mesh.f90..."
mpif90 -c -I ../../../modules ../../../mesh/mod_mesh.f90 || exit 1

echo "Compiling test_read_mesh.f90..."
mpif90 -o test_read_mesh test_read_mesh.f90 mod_mesh.o -I ../../../modules || exit 1

echo "Running test..."
mpirun -np 1 ./test_read_mesh

