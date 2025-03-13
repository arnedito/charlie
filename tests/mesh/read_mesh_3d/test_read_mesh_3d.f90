!=============================================================================
! @file   test_mesh_3d.f90
! @brief  Tests reading and loading a 3D mesh using mod_mesh.
!
! @details
!   This program uses the subroutine read_domain() from mod_mesh to load a
!   3D mesh from the specified test files. It then prints basic mesh statistics
!   (nodes, elements, boundaries, dimensions) and displays the first few node
!   coordinates to verify correctness.
!
! @note
!   - The boundary data may not be fully utilized if read_fix_dat is not yet
!     implemented in read_domain().
!=============================================================================
program test_mesh_3d
   use mod_mesh,     only: mesh_type, read_domain
   use constants,    only: ip, rp
   implicit none

   type(mesh_type) :: mesh

   !-------------------------------------------------------------------------
   !> @brief Load the 3D mesh from test files and print its statistics.
   !-------------------------------------------------------------------------
   call read_domain( &
        "tests/mesh/read_mesh_3d/case.dom.dat", &
        "tests/mesh/read_mesh_3d/case.geo.dat", &
        "tests/mesh/read_mesh_3d/case.fix.dat", &
        mesh )

   print*, "------3D Mesh successfully loaded.------"
   print*, "Number of nodes:     ", mesh%npoin
   print*, "Number of elements:  ", mesh%nelem
   print*, "Number of boundaries:", mesh%nboun
   print*, "Spatial dimensions:  ", mesh%ndim

   ! Display the first 10 node coordinates
   print*, "First 10 node coordinates:"
   call print_matrix(mesh%coords, min(10, mesh%npoin), mesh%ndim)

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
         write(*, '(100(F10.4,1x))') (mat(i,j), j=1, ncols)
      end do
   end subroutine print_matrix

end program test_mesh_3d
