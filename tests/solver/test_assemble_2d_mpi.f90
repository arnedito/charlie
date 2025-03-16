program test_assemble_2d_mpi
    use mpi
    use mod_mpi,      only: mpi_init_wrapper, mpi_finalize_wrapper, mpirank, mpisize, mpi_comm_global
    use mod_mesh,     only: mesh_type, read_domain
    use mod_solver,   only: assemble_system
    use mod_partition,only: compute_partition
    use mod_parallel, only: distribute_mesh
    use constants,    only: ip, rp
    implicit none

    type(mesh_type)                 :: mesh
    real(rp), allocatable           :: A(:,:), b(:)
    integer(ip)                     :: rank, size, ierr
    integer(ip)                     :: i, j, nprint, start_elem, end_elem
    integer(ip)                     :: n_rows, n_cols, ncoords, nboundary

    call mpi_init_wrapper()
    rank = mpirank
    size = mpisize

    if (rank == 0_ip) then
        print *, "=== MPI 2D Assembly Test Running on", size, "ranks ==="
    end if

    if (rank == 0_ip) then
        call read_domain("tests/mesh/read_mesh_files/case.dom.dat", &
                         "tests/mesh/read_mesh_files/case.geo.dat", &
                         "tests/mesh/read_mesh_files/case.fix.dat", &
                         mesh)
    end if

    call MPI_Bcast(mesh%npoin, 1, MPI_INTEGER, 0, mpi_comm_global, ierr)
    call MPI_Bcast(mesh%nelem, 1, MPI_INTEGER, 0, mpi_comm_global, ierr)
    call MPI_Bcast(mesh%nboun, 1, MPI_INTEGER, 0, mpi_comm_global, ierr)
    call MPI_Bcast(mesh%ndim,  1, MPI_INTEGER, 0, mpi_comm_global, ierr)

    if (rank /= 0_ip) then
        allocate(mesh%coords(mesh%npoin, mesh%ndim))
        allocate(mesh%boundary(mesh%nboun))
    end if

    ncoords   = mesh%npoin * mesh%ndim
    nboundary = mesh%nboun

    call MPI_Bcast(mesh%coords,   ncoords,   MPI_DOUBLE_PRECISION, 0, mpi_comm_global, ierr)
    call MPI_Bcast(mesh%boundary, nboundary, MPI_INTEGER,          0, mpi_comm_global, ierr)

    call compute_partition(mesh, rank, size, start_elem, end_elem)
    print *, "Rank", rank, "owns elements:", start_elem, "to", end_elem

    call distribute_mesh(mesh, rank, size)
    call MPI_Barrier(mpi_comm_global, ierr)

    allocate(A(mesh%npoin, mesh%npoin))
    allocate(b(mesh%npoin))
    A = 0.0_rp
    b = 0.0_rp

    call assemble_system(mesh, A, b)

    call MPI_Barrier(mpi_comm_global, ierr)

    if (rank == 0_ip) then
        print *, "======= 2D Assembly in MPI (Global Summation) ======="
        print *, "Global matrix A shape:", shape(A)

        n_rows = 5
        n_cols = 10
        do i = 1, n_rows
            if (i > mesh%npoin) exit
            print *, "Row", i, ":"
            do j = 1, n_cols
                if (j > mesh%npoin) exit
                write(*, '(F12.6)', advance='no') A(i,j)
            end do
            print *, "" ! New line after each row
        end do
    end if

    call mpi_finalize_wrapper()

end program test_assemble_2d_mpi
