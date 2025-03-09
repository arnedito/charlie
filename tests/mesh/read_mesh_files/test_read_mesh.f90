    program test_mesh
        use mod_mesh        ! Mesh module for handling mesh data
        use constants       ! Constants for data types
        implicit none

        type(mesh_type) :: mesh

        ! Load the mesh from test files
        call read_domain("tests/mesh/read_mesh_files/case.dom.dat", &
                        "tests/mesh/read_mesh_files/case.geo.dat", &
                        "tests/mesh/read_mesh_files/case.fix.dat", mesh)

        ! Output mesh information
        print*, "✅ Mesh successfully loaded."
        print*, "Number of nodes: ", mesh % npoin
        print*, "Number of elements: ", mesh % nelem
        print*, "Number of boundaries: ", mesh % nboun

        ! Display the first few node coordinates
        print*, "First 5 node coordinates:"
        call print_matrix(mesh % coords, min(5, mesh % npoin), mesh % ndim)

    contains

        ! Prints a matrix with the specified number of rows and columns
        subroutine print_matrix(mat, nrows, ncols)
            real(rp), intent(in)    ::    mat(:,:)
            integer, intent(in)     ::    nrows, ncols
            integer :: i, j
            do i = 1, nrows
                write(*, '(100(F8.3,1x))') (mat(i,j), j=1, ncols)
            end do
        end subroutine print_matrix

    end program test_mesh
