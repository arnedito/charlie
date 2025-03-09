# Fortran Compiler
FC = mpif90
FLAGS = -O2 -J modules/ -I modules/

# Core Modules (Ensure mod_mesh is compiled before mod_solver)
MODULES = modules/constants.f90 modules/utils.f90 \
          mesh/mod_mesh.f90 \
          parallel/mod_mpi.f90 parallel/mod_parallel.f90 \
          solver/mod_solver.f90

# Object files (store in modules/)
OBJ = $(patsubst %.f90, modules/%.o, $(notdir $(MODULES)))

# Test files
TEST_MPI = tests/parallel/test_mpi.f90
TEST_MESH_2D = tests/mesh/read_mesh_files/test_read_mesh.f90
TEST_MESH_3D = tests/mesh/read_mesh_3d/test_read_mesh_3d.f90  # NEW 3D TEST
TEST_SOLVER = tests/solver/test_solver.f90

# Executables
EXE_MPI = tests/parallel/test_mpi
EXE_MESH_2D = tests/mesh/read_mesh_files/test_read_mesh
EXE_MESH_3D = tests/mesh/read_mesh_3d/test_read_mesh_3d  # NEW 3D TEST EXECUTABLE
EXE_SOLVER = tests/solver/test_solver

# Default rule: compile all tests
all: $(EXE_MPI) $(EXE_MESH_2D) $(EXE_MESH_3D) $(EXE_SOLVER)

# Compile core modules
modules/constants.o: modules/constants.f90
	$(FC) -c $< $(FLAGS) -o $@

modules/utils.o: modules/utils.f90
	$(FC) -c $< $(FLAGS) -o $@

# Compile mesh module first (needed by solver)
modules/mod_mesh.o: mesh/mod_mesh.f90
	$(FC) -c $< $(FLAGS) -o $@

# Compile MPI modules
modules/mod_mpi.o: parallel/mod_mpi.f90
	$(FC) -c $< $(FLAGS) -o $@

modules/mod_parallel.o: parallel/mod_parallel.f90
	$(FC) -c $< $(FLAGS) -o $@

# Compile solver module (depends on mod_mesh)
modules/mod_solver.o: solver/mod_solver.f90 modules/mod_mesh.o
	$(FC) -c $< $(FLAGS) -o $@

# Compile and link MPI test
$(EXE_MPI): $(OBJ) $(TEST_MPI)
	$(FC) $(FLAGS) -o $@ $(TEST_MPI) $(OBJ)

# Compile and link 2D Mesh test
$(EXE_MESH_2D): $(OBJ) $(TEST_MESH_2D)
	$(FC) $(FLAGS) -o $@ $(TEST_MESH_2D) $(OBJ)

# Compile and link 3D Mesh test (NEW)
$(EXE_MESH_3D): $(OBJ) $(TEST_MESH_3D)
	$(FC) $(FLAGS) -o $@ $(TEST_MESH_3D) $(OBJ)

# Compile and link Solver test
$(EXE_SOLVER): $(OBJ) $(TEST_SOLVER)
	$(FC) $(FLAGS) -o $@ $(TEST_SOLVER) $(OBJ)

# Run tests
test_mpi: $(EXE_MPI)
	$(EXE_MPI)

test_mesh_2d: $(EXE_MESH_2D)
	$(EXE_MESH_2D)

test_mesh_3d: $(EXE_MESH_3D)
	$(EXE_MESH_3D)  # NEW 3D TEST RUN COMMAND

test_solver: $(EXE_SOLVER)
	$(EXE_SOLVER)

# Clean build files
clean:
	rm -f modules/*.o modules/*.mod $(EXE_MPI) $(EXE_MESH_2D) $(EXE_MESH_3D) $(EXE_SOLVER)
