    !**********************************************************************
    ! Module: mod_parallel
    !
    ! Purpose:
    !   This module handles all MPI communications related to mesh partitioning.
    !   It calls the partitioning routine from mod_partition and takes care of
    !   sending and receiving data between processes.
    !
    ! Usage:
    !   Use 'distribute_mesh' to partition the mesh and distribute it among MPI
    !   processes. The rank 0 process sends the relevant partitions to other ranks.
    !**********************************************************************
    module mod_parallel
        use mpi
        use mod_partition,      only: compute_partition
        use mod_mesh,           only: mesh_type
        implicit none

        public :: distribute_mesh

    contains

        subroutine distribute_mesh(mesh, rank, nprocs)
            !------------------------------------------------------------------
            ! Inputs/Outputs:
            !   mesh   - The mesh structure. On rank 0, this contains the full mesh.
            !            On other ranks, it will hold only the local partition.
            !   rank   - The MPI rank of the calling process.
            !   nprocs - Total number of MPI processes.
            !------------------------------------------------------------------
            type(mesh_type), intent(inout)          :: mesh
            integer, intent(in)                     :: rank, nprocs
            integer                                 :: start_elem, end_elem, num_local_elems, ierr
            integer                                 :: sendcount, i

            ! Get the partition indices for the current process.
            call compute_partition(mesh, rank, nprocs, start_elem, end_elem)
            num_local_elems = end_elem - start_elem + 1

            ! Allocate space for the local connectivity data.
            allocate(mesh % connectivity(num_local_elems, size(mesh % connectivity, 2)))

            ! If this is the root process, send out partitions to the others.
            if (rank == 0) then
                do i = 1, nprocs-1
                    call compute_partition(mesh, i, nprocs, start_elem, end_elem)
                    sendcount = (end_elem - start_elem + 1) * size(mesh % connectivity, 2)
                    call MPI_Send(mesh % connectivity(start_elem:end_elem, :), sendcount, MPI_INTEGER, i, 0, MPI_COMM_WORLD, ierr)
                end do
            else
                ! Non-root processes receive their partition.
                call receive_partition(mesh, rank, nprocs)
            end if
        end subroutine distribute_mesh

        subroutine receive_partition(mesh, rank, nprocs)
            !------------------------------------------------------------------
            ! Inputs/Outputs:
            !   mesh   - The mesh structure which will be allocated for the local partition.
            !   rank   - The MPI rank of the calling process.
            !   nprocs - Total number of MPI processes.
            !
            ! This subroutine receives the partitioned mesh data from the root.
            !------------------------------------------------------------------
            type(mesh_type), intent(inout)          :: mesh
            integer, intent(in)                     :: rank, nprocs
            integer                                 :: num_local_elems, ierr

            ! Here, we assume the partition size is computed similarly.
            call compute_partition(mesh, rank, nprocs, num_local_elems, num_local_elems)

            ! Allocate memory for the local mesh partition.
            allocate(mesh % connectivity(num_local_elems, size(mesh % connectivity, 2)))

            ! Receive the mesh partition data from process 0.
            call MPI_Recv(mesh % connectivity, num_local_elems * size(mesh % connectivity, 2), &
                        MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        end subroutine receive_partition

    end module mod_parallel
