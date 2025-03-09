    program test_mesh_partition
        use mod_mpi
        use mod_mesh
        use mod_parallel
        implicit none

        type(mesh_type)                 :: mesh
        integer                         :: rank, size, total_local_elems, expected_total
        integer, allocatable            :: recv_counts(:)

        ! Initialize MPI
        call mpi_init_wrapper()
        rank = mpirank
        size = mpisize

        ! Rank 0 loads the mesh, then distributes it
        if (rank == 0) then
            print*, "🔍 Rank 0: Reading mesh..."
            call read_domain("tests/mesh/read_mesh_3d/case_3d.dom.dat", &
                            "tests/mesh/read_mesh_3d/case_3d.geo.dat", &
                            "tests/mesh/read_mesh_3d/case_3d.fix.dat", mesh)
        end if

        ! Use master-worker structure for partitioning
        call master_worker(mesh)

        ! Check partitioning results
        total_local_elems = size(mesh % connectivity, 1)

        ! Verify each rank received some elements
        if (total_local_elems > 0) then
            print*, "✅ Rank ", rank, " received ", total_local_elems, " elements."
        else
            print*, "⚠ Warning: Rank ", rank, " received no elements!"
        end if

        ! Gather total elements across all ranks
        allocate(recv_counts(size))
        call MPI_Gather(total_local_elems, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, 0, mpi_comm_global, mpierr)

        if (rank == 0) then
            expected_total = sum(recv_counts)
            if (expected_total == mesh % nelem) then
                print*, "✅ Mesh successfully partitioned: Total elements =", expected_total
            else
                print*, "❌ Error: Distributed elements do not match original!"
                print*, "   Expected:", mesh % nelem, "but got:", expected_total
            end if
        end if

        !  Print first element of each rank for verification
        if (total_local_elems > 0) then
            print*, "🔍 Rank ", rank, " first element:", mesh % connectivity(1, :)
        end if

        ! Finalize MPI
        call mpi_finalize_wrapper()
    end program test_mesh_partition
