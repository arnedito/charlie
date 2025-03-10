program test_mesh_partition
use mpi
use mod_mpi,      only: mpi_init_wrapper, mpi_finalize_wrapper, mpi_comm_global, mpirank, mpisize
use mod_parallel, only: distribute_mesh
use mod_mesh,     only: mesh_type, read_domain, local_elements
implicit none

type(mesh_type)      :: mesh
integer              :: rank, size, total_local_elems, expected_total, ierr
integer, allocatable :: recv_counts(:)

! Initialize MPI environment.
call mpi_init_wrapper()
rank = mpirank
size = mpisize

! On rank 0, read the mesh from files.
if (rank == 0) then
    print*, "Rank 0: Reading mesh..."
    call read_domain("tests/mesh/read_mesh_3d/case.dom.dat", &
                        "tests/mesh/read_mesh_3d/case.geo.dat", &
                        "tests/mesh/read_mesh_3d/case.fix.dat", mesh)
end if

! Distribute mesh partitions across all MPI processes.
call distribute_mesh(mesh, rank, size)

! Synchronize all processes before checking the results.
call MPI_Barrier(mpi_comm_global, ierr)

! Use the helper function to get the number of local elements.
total_local_elems = local_elements(mesh)
if (total_local_elems > 0) then
    print*, "Rank ", rank, " received ", total_local_elems, " elements."
else
    print*, "Warning: Rank ", rank, " received no elements!"
end if

! Gather the number of elements from all processes on rank 0.
allocate(recv_counts(size))
call MPI_Gather(total_local_elems, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, 0, mpi_comm_global, ierr)

! On rank 0, verify that the sum of all elements matches the original mesh size.
if (rank == 0) then
    expected_total = sum(recv_counts)
    if (expected_total == mesh%nelem) then
        print*, "Mesh successfully partitioned: Total elements =", expected_total
    else
        print*, "Error: Distributed elements do not match original!"
        print*, "   Expected:", mesh%nelem, "but got:", expected_total
    end if
end if

! Optionally, print the first element from each partition for quick verification.
if (total_local_elems > 0) then
    print*, "Rank ", rank, " first element:", mesh%connectivity(1, :)
end if

! Clean up allocated memory.
if (allocated(recv_counts)) deallocate(recv_counts)

! Finalize the MPI environment.
call mpi_finalize_wrapper()
end program test_mesh_partition
