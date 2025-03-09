# Fortran Compiler
FC = mpif90
FLAGS = -O2 -J modules/ -I modules/

# Core Modules (Ensure correct order)
MODULES = modules/constants.f90 modules/utils.f90 \
          mesh/mod_mesh.f90 \
          parallel/mod_parallel.f90 parallel/mod_mpi.f90 \
          solver/mod_solver.f90

# Object files (store in modules/)
OBJ = $(patsubst %.f90, modules/%.o, $(notdir $(MODULES)))

# Test files
TEST_MPI = tests/parallel/test_mpi.f90
TEST_MESH_2D = tests/mesh/read_mesh_files/test_read_mesh.f90
TEST_MESH_3D = tests/mesh/read_mesh_3d/test_read_mesh_3d.f90
TEST_SOLVER = tests/solver/test_solver.f90
TEST_PARTITION = tests/parallel/test_mesh_partition.f90

# Executables
EXE_MPI = tests/parallel/test_mpi
EXE_MESH_2D = tests/mesh/read_mesh_files/test_read_mesh
EXE_MESH_3D = tests/mesh/read_mesh_3d/test_read_mesh_3d
EXE_SOLVER = tests/solver/test_solver
EXE_PARTITION = tests/parallel/test_mesh_partition

# Default rule: compile all tests
all: $(EXE_MPI) $(EXE_MESH_2D) $(EXE_MESH_3D) $(EXE_SOLVER) $(EXE_PARTITION)

# Compile core modules
modules/constants.o: modules/constants.f90
	$(FC) -c $< $(FLAGS) -o $@

modules/utils.o: modules/utils.f90
	$(FC) -c $< $(FLAGS) -o $@

# Compile mesh module first (needed by solver)
modules/mod_mesh.o: mesh/mod_mesh.f90
	$(FC) -c $< $(FLAGS) -o $@


modules/mod_parallel.o: parallel/mod_parallel.f90 modules/mod_mpi.o
	$(FC) -c $< $(FLAGS) -o $@


modules/mod_mpi.o: parallel/mod_mpi.f90 modules/mod_mesh.o
	$(FC) -c $< $(FLAGS) -o $@

# Compile solver module (depends on mod_mesh and MPI modules)
modules/mod_solver.o: solver/mod_solver.f90 modules/mod_mesh.o modules/mod_mpi.o modules/mod_parallel.o
	$(FC) -c $< $(FLAGS) -o $@

# Compile and link MPI test
$(EXE_MPI): $(OBJ) $(TEST_MPI)
	$(FC) $(FLAGS) -o $@ $(TEST_MPI) $(OBJ)

# Compile and link 2D Mesh test
$(EXE_MESH_2D): $(OBJ) $(TEST_MESH_2D)
	$(FC) $(FLAGS) -o $@ $(TEST_MESH_2D) $(OBJ)

# Compile and link 3D Mesh test
$(EXE_MESH_3D): $(OBJ) $(TEST_MESH_3D)
	$(FC) $(FLAGS) -o $@ $(TEST_MESH_3D) $(OBJ)

# Compile and link Solver test
$(EXE_SOLVER): $(OBJ) $(TEST_SOLVER)
	$(FC) $(FLAGS) -o $@ $(TEST_SOLVER) $(OBJ)

# Compile and link Mesh Partition test
$(EXE_PARTITION): $(OBJ) $(TEST_PARTITION)
	$(FC) $(FLAGS) -o $@ $(TEST_PARTITION) $(OBJ)

# Run tests
test_mpi: $(EXE_MPI)
	mpirun -np 4 $(EXE_MPI)

test_mesh_2d: $(EXE_MESH_2D)
	$(EXE_MESH_2D)

test_mesh_3d: $(EXE_MESH_3D)
	$(EXE_MESH_3D)

test_solver: $(EXE_SOLVER)
	$(EXE_SOLVER)

test_partition: $(EXE_PARTITION)
	mpirun -np 4 $(EXE_PARTITION)

# Clean build files
clean:
	rm -f modules/*.o modules/*.mod $(EXE_MPI) $(EXE_MESH_2D) $(EXE_MESH_3D) $(EXE_SOLVER) $(EXE_PARTITION)
