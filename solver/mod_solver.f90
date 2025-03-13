!=============================================================================
! @file   mod_solver.f90
! @brief  Provides local assembly and a basic solver for a finite element problem.
!
! @details
!   - assemble_system(mesh, A_local, b_local):
!       Loops over local elements (mesh%connectivity) to compute element-level
!       matrices/vectors, then assembles them into A_local and b_local.
!   - solve_newton(mesh, A_local, b_local, x_local):
!       A placeholder routine demonstrating how you'd proceed with a solver.
!   - compute_element_matrices(mesh, iel, ke, fe):
!       Dummy element routine that sets diagonal values and a simple force vector.
!
!=============================================================================
module mod_solver
   use constants,    only: ip, rp
   use mod_mesh,     only: mesh_type
   implicit none

   public :: assemble_system
   public :: solve_newton

contains

    !---------------------------------------------------------------------------
    !> @brief Assembles the local system matrix (A_local) and RHS vector (b_local)
    !>        for each MPI rank based on the local mesh data.
    !
    !! - mesh: Partitioned mesh structure, containing local or full elements
    !! - A_local: The system matrix, dimension (npoin x npoin)
    !! - b_local: The RHS vector, dimension (npoin)
    !---------------------------------------------------------------------------
    subroutine assemble_system(mesh, A_local, b_local)
        type(mesh_type), intent(in)               :: mesh
        real(rp), allocatable, intent(inout)      :: A_local(:,:), b_local(:)
        integer(ip)                               :: iel, i, j
        integer(ip)                               :: ndof
        integer(ip), allocatable                  :: e_dof(:)
        real(rp), allocatable                     :: ke(:,:), fe(:)

        ! Determine number of nodes per element (2D tri -> 3, 3D tet -> 4)
        ndof = mesh%ndim + 1

        ! Temporary arrays for local element computations
        allocate(e_dof(ndof))
        allocate(ke(ndof, ndof))
        allocate(fe(ndof))

        ! Loop over local elements
        do iel = 1, size(mesh%connectivity, 1)

            ! e_dof: the global DOF (node) indices for this element
            e_dof = mesh%connectivity(iel, :)

            ! Compute local ke, fe
            call compute_element_matrices(mesh, iel, ke, fe)

            ! Assemble ke, fe into A_local, b_local
            do i = 1, ndof
            do j = 1, ndof
                A_local(e_dof(i), e_dof(j)) = A_local(e_dof(i), e_dof(j)) + ke(i,j)
            end do
            b_local(e_dof(i)) = b_local(e_dof(i)) + fe(i)
            end do

        end do

        deallocate(e_dof, ke, fe)
    end subroutine assemble_system

    !---------------------------------------------------------------------------
    !> @brief A placeholder routine demonstrating Newton's method structure.
    !>
    !! - mesh: The mesh data (not strictly used here).
    !! - A_local: The local system matrix built from assemble_system().
    !! - b_local: The local RHS vector.
    !! - x_local: The solution vector (allocated to npoin).
    !---------------------------------------------------------------------------
    subroutine solve_newton(mesh, A_local, b_local, x_local)
        type(mesh_type), intent(in)             :: mesh
        real(rp), allocatable, intent(inout)    :: A_local(:,:), b_local(:), x_local(:)

        ! Just a debug print for now
        print *, "solve_newton: First row of A_local => ", A_local(1, :)
        print *, "solve_newton: First entries of b_local =>", b_local(1: min(5, size(b_local)))
        print *, "solve_newton: (Placeholder) Setting x_local = 0"
        x_local = 0.0_rp

        ! Real code will:
        ! 1) Apply boundary conditions
        ! 2) Solve A_local * x_local = b_local (iterative or direct)
        ! 3) Update x_local
    end subroutine solve_newton

    !---------------------------------------------------------------------------
    !> @brief Computes element-level stiffness (ke) and force (fe).
    !
    !! - mesh : local or global mesh data if needed for shape function integrals
    !! - iel  : which element index in mesh%connectivity
    !! - ke   : ndof x ndof matrix (3x3 for tri, 4x4 for tet)
    !! - fe   : ndof array
    !---------------------------------------------------------------------------
    subroutine compute_element_matrices(mesh, iel, ke, fe)
        type(mesh_type), intent(in)             :: mesh
        integer(ip), intent(in)                 :: iel
        real(rp), intent(out)                   :: ke(:,:), fe(:)
        integer(ip)                             :: ndof, i, j

        ndof = size(ke, 1)  ! 3 for tri, 4 for tet (linear)

        ! Zero them first
        do i = 1, ndof
            fe(i) = 0.0_rp
        end do

        do i = 1, ndof
            do j = 1, ndof
            ke(i, j) = 0.0_rp
            end do
        end do

        ! Example diagonal
        do i = 1, ndof
            ke(i, i) = 1.0_rp
            fe(i)    = real(i, rp)
        end do

    end subroutine compute_element_matrices

end module mod_solver
