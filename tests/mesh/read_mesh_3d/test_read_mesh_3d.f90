    program test_mesh_3d
        use mod_mesh
        use constants
        implicit none

        type(mesh_type) :: mesh

        ! Load the 3D mesh from test files
        call read_domain("tests/mesh/read_mesh_3d/case.dom.dat", &
                        "tests/mesh/read_mesh_3d/case.geo.dat", &
                        "tests/mesh/read_mesh_3d/case.fix.dat", mesh)

        ! Output mesh information
        print*, "✅ 3D Mesh successfully loaded."
        print*, "Number of nodes: ", mesh % npoin
        print*, "Number of elements: ", mesh % nelem
        print*, "Number of boundaries: ", mesh % nboun
        print*, "Spatial dimensions: ", mesh % ndim

        ! Display the first few node coordinates
        print*, "First 10 node coordinates:"
        call print_matrix(mesh % coords, min(10, mesh % npoin), mesh % ndim)

    contains

        ! Prints a matrix with the specified number of rows and columns
        subroutine print_matrix(mat, nrows, ncols)
            real(rp), intent(in)    ::    mat(:,:)
            integer, intent(in)     ::    nrows, ncols
            integer :: i, j
            do i = 1, nrows
                write(*, '(100(F10.4,1x))') (mat(i,j), j=1, ncols)
            end do
        end subroutine print_matrix

    end program test_mesh_3d
