!=============================================================================
! @file   test_assemble_3d_mpi.f90
! @brief  MPI test for assembling a 3D mesh system (tetrahedral elements).
!
! @details
!   Steps:
!    1) MPI init
!    2) Rank 0 reads the 3D mesh files
!    3) Broadcast mesh metadata to all ranks
!    4) Each rank allocates the full mesh arrays
!    5) Broadcast coords, boundary
!    6) Allocate global system matrix A and vector b
!    7) Call assemble_system
!    8) Print partial results
!    9) MPI finalize
!
!=============================================================================

program test_assemble_3d_mpi
    use mpi
    use mod_mpi,      only: mpi_init_wrapper, mpi_finalize_wrapper, mpirank, mpisize, mpi_comm_global
    use mod_mesh,     only: mesh_type, read_domain
    use mod_solver,   only: assemble_system
    use constants,    only: ip, rp
    implicit none

    type(mesh_type)                 :: mesh
    real(rp), allocatable           :: A(:,:), b(:)
    integer(ip)                     :: rank, size, ierr
    integer                         :: ncoords, nboundary, i, nprint

    !---------------------------------------------------------------------
    ! 1) MPI init
    !---------------------------------------------------------------------
    call mpi_init_wrapper()
    rank = mpirank
    size = mpisize

    if (rank == 0_ip) then
        print *, "=== MPI 3D Assembly Test: Reading mesh on Rank 0 ==="
    end if

    !---------------------------------------------------------------------
    ! 2) Rank 0 reads the 3D mesh from test files
    !---------------------------------------------------------------------
    if (rank == 0_ip) then
        call read_domain("tests/mesh/read_mesh_3d/case.dom.dat", &
                        "tests/mesh/read_mesh_3d/case.geo.dat", &
                        "tests/mesh/read_mesh_3d/case.fix.dat", &
                        mesh)
    end if

    !---------------------------------------------------------------------
    ! 3) Broadcast mesh metadata
    !---------------------------------------------------------------------
    call MPI_Bcast(mesh%npoin, 1, MPI_INTEGER,          0, mpi_comm_global, ierr)
    call MPI_Bcast(mesh%nelem, 1, MPI_INTEGER,          0, mpi_comm_global, ierr)
    call MPI_Bcast(mesh%nboun, 1, MPI_INTEGER,          0, mpi_comm_global, ierr)
    call MPI_Bcast(mesh%ndim,  1, MPI_INTEGER,          0, mpi_comm_global, ierr)

    !---------------------------------------------------------------------
    ! 4) Each rank allocates mesh arrays
    !---------------------------------------------------------------------
    if (.not. allocated(mesh%coords)) then
        allocate(mesh%coords(mesh%npoin, mesh%ndim))
    end if
    if (.not. allocated(mesh%boundary)) then
        allocate(mesh%boundary(mesh%nboun))
    end if

    !---------------------------------------------------------------------
    ! 5) Broadcast coords, boundary
    !---------------------------------------------------------------------
    ncoords   = mesh%npoin * mesh%ndim
    nboundary = mesh%nboun

    call MPI_Bcast(mesh%coords,   ncoords,   MPI_DOUBLE_PRECISION, 0, mpi_comm_global, ierr)
    call MPI_Bcast(mesh%boundary, nboundary, MPI_INTEGER,          0, mpi_comm_global, ierr)

    if (rank == 0_ip) then
        print *, "==> 3D mesh broadcast done"
        print *, "   npoin =", mesh%npoin, ", nelem =", mesh%nelem
    end if

    !---------------------------------------------------------------------
    ! 6) Allocate global system A, b
    !---------------------------------------------------------------------
    allocate(A(mesh%npoin, mesh%npoin))
    allocate(b(mesh%npoin))
    A = 0.0_rp
    b = 0.0_rp

    !---------------------------------------------------------------------
    ! 7) Call assemble_system
    !---------------------------------------------------------------------
    if (rank == 0_ip) then
        call assemble_system(mesh, A, b)
    end if

    !---------------------------------------------------------------------
    ! 8) Print partial results on rank 0
    !---------------------------------------------------------------------
    call MPI_Barrier(mpi_comm_global, ierr)
    if (rank == 0_ip) then
        print *, "=== 3D Assembly in MPI (Rank 0 final output) ==="
        print *, "Global matrix A shape: ", shape(A)

        ! Print first 10 columns of row 1
        print *, "Sample row of A (1st row, 10 columns): "
        do i = 1, min(10, mesh%npoin)
            write(*, '(F6.2, 1x)', advance='no') A(1, i)
        end do
        write(*,*)

        nprint = min(5, mesh%npoin)
        print *, "Sample of b: "
        do i = 1, nprint
            write(*, '(F6.2, 1x)', advance='no') b(i)
        end do
        write(*,*)
    end if

    !---------------------------------------------------------------------
    ! 9) MPI finalize
    !---------------------------------------------------------------------
    call mpi_finalize_wrapper()

end program test_assemble_3d_mpi
