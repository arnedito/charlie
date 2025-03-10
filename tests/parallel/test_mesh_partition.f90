    program test_mesh_partition
        use mpi
        use mod_mpi,        only: mpi_init_wrapper, mpi_finalize_wrapper, mpi_comm_global, distribute_mesh, mpirank, mpisize
        use mod_mesh,       only: mesh_type, read_domain
        use mod_parallel,   only: master_worker
        implicit none

        type(mesh_type)                 :: mesh
        integer                         :: rank, size, total_local_elems, expected_total, ierr
        integer, allocatable            :: recv_counts(:)

        ! Initialize MPI
        call mpi_init_wrapper()
        rank = mpirank
        size = mpisize

        ! Rank 0 loads the mesh, then distributes it
        if (rank == 0) then
            print*, "🔍 Rank 0: Reading mesh..."
            call read_domain("tests/mesh/read_mesh_3d/case.dom.dat", &
                            "tests/mesh/read_mesh_3d/case.geo.dat", &
                            "tests/mesh/read_mesh_3d/case.fix.dat", mesh)
        end if

        ! Distribute mesh partitions to all ranks
        call distribute_mesh(mesh)

        ! Ensure all processes reach this point before checking results
        call MPI_Barrier(mpi_comm_global, ierr)

        ! Check partitioning results
        total_local_elems = size(mesh % connectivity, 1)

        if (total_local_elems > 0) then
            print*, "✅ Rank ", rank, " received ", total_local_elems, " elements."
        else
            print*, "⚠ Warning: Rank ", rank, " received no elements!"
        end if

        ! Gather total elements across all ranks
        allocate(recv_counts(size))
        call MPI_Gather(total_local_elems, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, 0, mpi_comm_global, ierr)

        ! Verify correct partitioning
        if (rank == 0) then
            expected_total = sum(recv_counts)
            if (expected_total == mesh % nelem) then
                print*, "✅ Mesh successfully partitioned: Total elements =", expected_total
            else
                print*, "❌ Error: Distributed elements do not match original!"
                print*, "   Expected:", mesh % nelem, "but got:", expected_total
            end if
        end if

        ! Print first element of each rank for verification
        if (total_local_elems > 0) then
            print*, "🔍 Rank ", rank, " first element:", mesh % connectivity(1, :)
        end if

        ! Free memory
        if (allocated(recv_counts)) deallocate(recv_counts)

        ! Finalize MPI
        call mpi_finalize_wrapper()
    end program test_mesh_partition
