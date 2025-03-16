# Fortran Compiler
FC = mpif90
FLAGS = -O2 -std=f2008 -Wall -Wextra -g -fcheck=all -fbounds-check -fbacktrace -J modules/ -I modules/
FLAGS += -Wimplicit-interface -Wsurprising -Wconversion -Warray-temporaries

# Core Modules (Correct Compilation Order)
MODULES = modules/constants.f90 modules/utils.f90 \
          mesh/mod_mesh.f90 \
          parallel/mod_partition.f90 \
          parallel/mod_mpi.f90 \
          parallel/mod_parallel.f90 \
          solver/mod_solver.f90

# Object files (store in modules/)
OBJ = $(patsubst %.f90, modules/%.o, $(notdir $(MODULES)))

# Test files
TEST_MPI               = tests/parallel/test_mpi.f90
TEST_MESH_2D           = tests/mesh/read_mesh_files/test_read_mesh.f90
TEST_MESH_3D           = tests/mesh/read_mesh_3d/test_read_mesh_3d.f90
TEST_SOLVER            = tests/solver/test_solver.f90
TEST_PARTITION_2D      = tests/parallel/test_mesh_partition.f90
TEST_PARTITION_3D      = tests/parallel/test_mesh_partition_3d.f90
TEST_ASSEMBLE_2D_MPI   = tests/solver/test_assemble_2d_mpi.f90
TEST_ASSEMBLE_3D_MPI   = tests/solver/test_assemble_3d_mpi.f90

# Executables
EXE_MPI                = tests/parallel/test_mpi
EXE_MESH_2D            = tests/mesh/read_mesh_files/test_read_mesh
EXE_MESH_3D            = tests/mesh/read_mesh_3d/test_read_mesh_3d
EXE_SOLVER             = tests/solver/test_solver
EXE_PARTITION_2D       = tests/parallel/test_mesh_partition
EXE_PARTITION_3D       = tests/parallel/test_mesh_partition_3d
EXE_ASSEMBLE_2D_MPI    = tests/solver/test_assemble_2d_mpi
EXE_ASSEMBLE_3D_MPI    = tests/solver/test_assemble_3d_mpi

# Default rule: compile all test executables.
all: $(EXE_MPI) $(EXE_MESH_2D) $(EXE_MESH_3D) $(EXE_SOLVER) $(EXE_PARTITION_2D) $(EXE_PARTITION_3D) \
     $(EXE_ASSEMBLE_2D_MPI) $(EXE_ASSEMBLE_3D_MPI)

# Compile core modules
modules/constants.o: modules/constants.f90
	$(FC) -c $< $(FLAGS) -o $@

modules/utils.o: modules/utils.f90
	$(FC) -c $< $(FLAGS) -o $@

modules/mod_mesh.o: mesh/mod_mesh.f90
	$(FC) -c $< $(FLAGS) -o $@

modules/mod_partition.o: parallel/mod_partition.f90 modules/mod_mesh.o
	$(FC) -c $< $(FLAGS) -o $@

modules/mod_mpi.o: parallel/mod_mpi.f90
	$(FC) -c $< $(FLAGS) -o $@

modules/mod_parallel.o: parallel/mod_parallel.f90 modules/mod_partition.o modules/mod_mpi.o
	$(FC) -c $< $(FLAGS) -o $@

modules/mod_solver.o: solver/mod_solver.f90 modules/mod_mesh.o modules/mod_mpi.o modules/mod_parallel.o
	$(FC) -c $< $(FLAGS) -o $@

# Link each test executable with the object files
$(EXE_MPI): $(OBJ) $(TEST_MPI)
	$(FC) $(FLAGS) -o $@ $(TEST_MPI) $(OBJ)

$(EXE_MESH_2D): $(OBJ) $(TEST_MESH_2D)
	$(FC) $(FLAGS) -o $@ $(TEST_MESH_2D) $(OBJ)

$(EXE_MESH_3D): $(OBJ) $(TEST_MESH_3D)
	$(FC) $(FLAGS) -o $@ $(TEST_MESH_3D) $(OBJ)

$(EXE_SOLVER): $(OBJ) $(TEST_SOLVER)
	$(FC) $(FLAGS) -o $@ $(TEST_SOLVER) $(OBJ)

$(EXE_PARTITION_2D): $(OBJ) $(TEST_PARTITION_2D)
	$(FC) $(FLAGS) -o $@ $(TEST_PARTITION_2D) $(OBJ)

$(EXE_PARTITION_3D): $(OBJ) $(TEST_PARTITION_3D)
	$(FC) $(FLAGS) -o $@ $(TEST_PARTITION_3D) $(OBJ)

# New linking rules for the 2D/3D assembly tests
$(EXE_ASSEMBLE_2D_MPI): $(OBJ) $(TEST_ASSEMBLE_2D_MPI)
	$(FC) $(FLAGS) -o $@ $(TEST_ASSEMBLE_2D_MPI) $(OBJ)

$(EXE_ASSEMBLE_3D_MPI): $(OBJ) $(TEST_ASSEMBLE_3D_MPI)
	$(FC) $(FLAGS) -o $@ $(TEST_ASSEMBLE_3D_MPI) $(OBJ)

# Targets for running the tests with MPI.
test_mpi:
	mpirun -np 4 $(EXE_MPI)

test_partition_2d:
	mpirun -np 4 $(EXE_PARTITION_2D)

test_partition_3d:
	mpirun -np 4 $(EXE_PARTITION_3D)

test_assemble_2d_mpi:
	mpirun -np 4 $(EXE_ASSEMBLE_2D_MPI)

test_assemble_3d_mpi:
	mpirun -np 4 $(EXE_ASSEMBLE_3D_MPI)

# Other test targets for mesh and solver
test_mesh_2d:
	./$(EXE_MESH_2D)

test_mesh_3d:
	./$(EXE_MESH_3D)

test_solver:
	./$(EXE_SOLVER)

# Clean rule to remove object files, modules, and executables.
clean:
	rm -f modules/*.o modules/*.mod \
	      $(EXE_MPI) $(EXE_MESH_2D) $(EXE_MESH_3D) $(EXE_SOLVER) \
	      $(EXE_PARTITION_2D) $(EXE_PARTITION_3D) $(EXE_ASSEMBLE_2D_MPI) $(EXE_ASSEMBLE_3D_MPI)
	rm -rf *.o *.mod *.dSYM tests/parallel/*.dSYM tests/mesh/*.dSYM tests/solver/*.dSYM
