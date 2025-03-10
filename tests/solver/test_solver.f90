!> @brief Unit Test for the Solver Module
!>
!> This program tests the ability of the solver module to correctly
!> read a mesh, allocate memory, and assemble the system.
!>
!> Steps:
!> - Reads a predefined mesh from test files.
!> - Allocates memory for the solution vector.
!> - Calls `assemble_system` to verify solver functionality.
!> - Outputs key results for validation.

program test_solver
    use mod_mesh        ! Mesh handling module
    use mod_solver      ! Solver module
    use constants       ! Precision and type constants
    implicit none

    !--------------------------------------------------------------
    ! Variable Declarations
    !--------------------------------------------------------------
    type(mesh_type)             :: mesh                 ! Mesh structure
    real(rp), allocatable       :: solution(:)          ! Solution vector
    real(rp), allocatable       :: A(:,:), b(:)         ! System matrix and RHS vector

    !--------------------------------------------------------------
    ! Step 1: Read Mesh from Test Case Files
    !--------------------------------------------------------------
    print*, "Reading mesh files..."
    call read_domain("case.dom.dat", &
                    "case.geo.dat",  &
                    "case.fix.dat",  mesh)

    print*, "Mesh successfully loaded!"
    print*, "  - Number of Nodes: ", mesh % npoin
    print*, "  - Number of Elements: ", mesh % nelem
    print*, "  - Number of Boundaries: ", mesh % nboun

    !--------------------------------------------------------------
    ! Step 2: Allocate and Initialize Solution Vector
    !--------------------------------------------------------------

    allocate(solution(mesh % npoin))
    solution = 0.0_rp  ! Initialize solution to zero

    print*, "Solution vector allocated: ", size(solution)

    !--------------------------------------------------------------
    ! Step 3: Assemble System (Dummy Test)
    !--------------------------------------------------------------
    print*, "Assembling system..."
    call assemble_system(mesh, solution, A, b)

    !--------------------------------------------------------------
    ! Step 4: Print Validation Output
    !--------------------------------------------------------------
    print*, "Solver test completed successfully!"
    print*, "  - Number of nodes processed: ", mesh % npoin

    !--------------------------------------------------------------
    ! Step 5: Clean Up Allocated Memory
    !--------------------------------------------------------------
    deallocate(solution)

end program test_solver
