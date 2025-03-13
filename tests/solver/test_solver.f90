!=============================================================================
! @file   test_solver.f90
! @brief  Unit Test for the Solver Module
!
! @details
!   This program tests the solver module's ability to:
!    - Read a mesh from test files.
!    - Allocate memory for system matrix (A), RHS (b), and a dummy solution vector.
!    - Call `assemble_system` to verify solver functionality.
!    - Print key results for validation.
!=============================================================================
program test_solver
   use mod_mesh,   only: mesh_type, read_domain
   use mod_solver, only: assemble_system
   use constants,  only: rp, ip
   implicit none

   !--------------------------------------------------------------
   ! Variable Declarations
   !--------------------------------------------------------------
   type(mesh_type)                  :: mesh
   real(rp), allocatable            :: solution(:)
   real(rp), allocatable            :: A(:,:), b(:)

   !--------------------------------------------------------------
   ! Step 1: Read Mesh from Test Case Files
   !--------------------------------------------------------------
   print*, "Reading mesh files..."
   call read_domain("tests/solver/case.dom.dat", &
                    "tests/solver/case.geo.dat", &
                    "tests/solver/case.fix.dat", &
                    mesh)

   print*, "Mesh successfully loaded!"
   print*, "  - Number of Nodes:     ", mesh%npoin
   print*, "  - Number of Elements:  ", mesh%nelem
   print*, "  - Number of Boundaries:", mesh%nboun

   !--------------------------------------------------------------
   ! Step 2: Allocate & Initialize Arrays (A, b, solution)
   !--------------------------------------------------------------
   allocate(A(mesh%npoin, mesh%npoin))
   allocate(b(mesh%npoin))
   allocate(solution(mesh%npoin))

   A = 0.0_rp
   b = 0.0_rp
   solution = 0.0_rp

   print*, "Allocated system arrays:"
   print*, "  - A shape: ", shape(A)
   print*, "  - b shape: ", shape(b)
   print*, "  - solution size:", size(solution)

   !--------------------------------------------------------------
   ! Step 3: Assemble System (Dummy Test)
   !--------------------------------------------------------------
   print*, "Assembling system..."
   call assemble_system(mesh, A, b)

   ! Optionally, you could call solve_newton(mesh, A, b, solution) here.

   !--------------------------------------------------------------
   ! Step 4: Print Validation Output
   !--------------------------------------------------------------
   print*, "Solver test completed successfully!"
   print*, "  - A(1,:) =", A(1,:)
   print*, "  - b(1:5) =", b(1:min(5, size(b)))

   !--------------------------------------------------------------
   ! Step 5: Clean Up Allocated Memory
   !--------------------------------------------------------------
   deallocate(A, b, solution)

end program test_solver
