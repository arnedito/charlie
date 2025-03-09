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

# Test files (Corrected mesh test name)
TEST_MPI = tests/parallel/test_mpi.f90
TEST_MESH = tests/mesh/read_mesh_files/test_read_mesh.f90
TEST_SOLVER = tests/solver/test_solver.f90

# Executables (Corrected mesh test executable name)
EXE_MPI = tests/parallel/test_mpi
EXE_MESH = tests/mesh/read_mesh_files/test_read_mesh
EXE_SOLVER = tests/solver/test_solver

# Default rule: compile all tests
all: $(EXE_MPI) $(EXE_MESH) $(EXE_SOLVER)

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

# Compile and link Mesh test (Corrected name)
$(EXE_MESH): $(OBJ) $(TEST_MESH)
	$(FC) $(FLAGS) -o $@ $(TEST_MESH) $(OBJ)

# Compile and link Solver test
$(EXE_SOLVER): $(OBJ) $(TEST_SOLVER)
	$(FC) $(FLAGS) -o $@ $(TEST_SOLVER) $(OBJ)

# Clean build files
clean:
	rm -f modules/*.o modules/*.mod $(EXE_MPI) $(EXE_MESH) $(EXE_SOLVER)
