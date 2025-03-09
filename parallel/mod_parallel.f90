    !-----------------------------------------------------------------
    ! Module: mod_parallel
    ! Purpose: Provides a basic master-worker structure utilizing MPI.
    !          This module leverages the MPI functionalities defined in
    !          mod_mpi to implement parallel tasks distribution.
    !          The master process (rank 0) is responsible for distributing
    !          tasks (such as mesh distribution), while worker processes wait
    !          for tasks.
    !-----------------------------------------------------------------

    module mod_parallel
        use mpi
        use mod_mpi,        only: mpirank, mpisize, mpi_comm_global, mpi_barrier_wrapper
        use mod_mesh,       only: mesh_type
        implicit none

        public :: master_worker, partition_mesh

    contains

        !-----------------------------------------------------------------
        subroutine master_worker()
            !-------------------------------------------------------------
            ! Purpose: Implements the master-worker paradigm.
            !          - If the current process is the master (rank 0), it will
            !            distribute tasks (e.g., mesh data).
            !          - Worker processes (non-zero rank) will wait for tasks.
            !          - A barrier is used to synchronize all processes before
            !            proceeding further.
            !-------------------------------------------------------------
            if (mpirank == 0) then
                print *, 'Master process (rank ', mpirank, '): distributing tasks...'
                ! Future implementation: Call subroutine to distribute the mesh or tasks.
            else
                print *, 'Worker process (rank ', mpirank, '): waiting for tasks...'
                ! Future implementation: Call subroutine to receive and process tasks.
            end if

            ! Synchronize all processes to ensure all have reached this point.
            call mpi_barrier_wrapper()
        end subroutine master_worker
        !-----------------------------------------------------------------

        !-----------------------------------------------------------------
        subroutine partition_mesh(mesh, rank, nprocs)
            !-------------------------------------------------------------
            ! Purpose: Distributes elements among MPI processes.
            ! - Rank 0 partitions the mesh and sends data to workers.
            ! - Other ranks receive their subdomain.
            !-------------------------------------------------------------
            type(mesh_type), intent(inout)      :: mesh
            integer, intent(in)                 :: rank, nprocs
            integer                             :: i, start_elem, end_elem, num_local_elems
            integer                             :: ierr, sendcount

            ! Compute partitioning
            num_local_elems = mesh % nelem / nprocs
            start_elem = rank * num_local_elems + 1
            end_elem = start_elem + num_local_elems - 1

            if (rank == nprocs - 1) then
                ! Last process gets any remainder
                end_elem = mesh % nelem
            end if

            ! Allocate memory for local elements
            allocate(mesh % connectivity(num_local_elems, size(mesh % connectivity,2)))

            ! Rank 0 sends partitions to other ranks
            if (rank == 0) then
                do i = 1, nprocs-1
                    start_elem = i * num_local_elems + 1
                    end_elem = start_elem + num_local_elems - 1
                    if (i == nprocs - 1) end_elem = mesh % nelem

                    sendcount = (end_elem - start_elem + 1) * size(mesh % connectivity,2)
                    call MPI_Send(mesh % connectivity(start_elem:end_elem, :), sendcount, MPI_INTEGER, i, 0, mpi_comm_global, ierr)
                end do
            else
                call receive_partition(mesh, rank, nprocs)
            end if
        end subroutine partition_mesh
        !-----------------------------------------------------------------

        !-----------------------------------------------------------------
        subroutine receive_partition(mesh, rank, nprocs)
            !-------------------------------------------------------------
            ! Purpose: Worker processes (ranks > 0) receive their mesh subdomain.
            !-------------------------------------------------------------
            type(mesh_type), intent(inout)      :: mesh
            integer, intent(in)                 :: rank, nprocs
            integer                             :: num_local_elems, ierr

            ! Compute local partition size
            num_local_elems = mesh%nelem / nprocs
            if (rank == nprocs - 1) then
                num_local_elems = mesh % nelem - (num_local_elems * (nprocs-1))
            end if

            ! Allocate memory for received data
            allocate(mesh % connectivity(num_local_elems, size(mesh % connectivity, 2)))

            ! Receive data from Rank 0
            call MPI_Recv(mesh % connectivity, &
                          num_local_elems * size(mesh % connectivity, 2), &
                          MPI_INTEGER, 0, 0, mpi_comm_global, MPI_STATUS_IGNORE, ierr)

        end subroutine receive_partition
        !-----------------------------------------------------------------

    end module mod_parallel
