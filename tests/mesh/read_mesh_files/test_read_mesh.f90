!=============================================================================
! @file   test_mesh.f90
! @brief  Tests reading and loading a 2D mesh using mod_mesh.
!
! @details
!   This program demonstrates reading a mesh from input files using read_domain().
!   It outputs basic mesh information to confirm successful loading, and
!   prints the first few node coordinates for visual inspection.
!
! @note
!   - Boundary data may not be displayed if read_fix_dat is not yet active.
!=============================================================================
program test_mesh
   use mod_mesh,     only: mesh_type, read_domain
   use constants,    only: ip, rp
   implicit none

   type(mesh_type) :: mesh

   !-------------------------------------------------------------------------
   !> @brief Load the mesh from test files and print basic stats.
   !-------------------------------------------------------------------------
   call read_domain( &
        "tests/mesh/read_mesh_files/case.dom.dat", &
        "tests/mesh/read_mesh_files/case.geo.dat", &
        "tests/mesh/read_mesh_files/case.fix.dat", &
        mesh )

   print*, "------Mesh successfully loaded.------"
   print*, "Number of nodes:     ", mesh%npoin
   print*, "Number of elements:  ", mesh%nelem
   print*, "Number of boundaries:", mesh%nboun
   print*, "Spatial dimensions:  ", mesh%ndim

   ! Display the first 5 node coordinates
   print*, "First 5 node coordinates:"
   call print_matrix(mesh%coords, min(5, mesh%npoin), mesh%ndim)

contains

   !-------------------------------------------------------------------------
   !> @brief Prints the first nrows x ncols portion of a 2D matrix.
   !>
   !! @param[in] mat   (real(rp)) The matrix to print, dimensioned at least (nrows, ncols).
   !! @param[in] nrows (integer)  Number of rows to print.
   !! @param[in] ncols (integer)  Number of columns to print.
   !-------------------------------------------------------------------------
   subroutine print_matrix(mat, nrows, ncols)
      real(rp), intent(in) :: mat(:,:)
      integer, intent(in)  :: nrows, ncols
      integer              :: i, j

      do i = 1, nrows
         write(*, '(100(F8.3,1x))') (mat(i,j), j=1, ncols)
      end do
   end subroutine print_matrix

end program test_mesh
