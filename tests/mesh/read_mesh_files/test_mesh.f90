    program test_mesh
        use mod_mesh        ! Uses our module with mesh_type and read_domain
        use constants       ! For ip, rp
        implicit none

        type(mesh_type) :: mesh

        ! Read a sample mesh from test files (replace with your actual test files)
        call read_domain("case.dom.dat", "case.geo.dat", "case.fix.dat", mesh)

        ! Print basic information about the mesh
        print*, "Number of nodes: ", mesh % npoin
        print*, "Number of elements: ",mesh % nelem
        print*, "Number of boundaries: ", mesh % nboun

        ! Optionally, print first few node coordinates
        print*, "First 5 node coordinates:"
        call print_matrix(mesh % coords, min(5, mesh % npoin), mesh % ndim)

    contains

        subroutine print_matrix(mat, nrows, ncols)
            real(rp), intent(in)    ::    mat(:,:)
            integer, intent(in)     ::    nrows, ncols
            integer :: i, j
            do i = 1, nrows
                write(*, '(100(F8.3,1x))') (mat(i,j), j=1, ncols)
            end do
        end subroutine print_matrix

    end program test_mesh
