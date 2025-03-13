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
    use constants,          only: rp, ip
    use mod_partition,      only: compute_partition, compute_local_nodes
    use mod_mesh,           only: mesh_type
    implicit none

    public :: distribute_mesh

contains

    subroutine distribute_mesh(mesh, rank, nprocs)
        type(mesh_type), intent(inout)          :: mesh
        integer(ip), intent(in)                 :: rank, nprocs
        integer(ip)                             :: start_elem, end_elem, num_local_elems
        integer(ip)                             :: sendcount, i, ierr
        integer(ip)                             :: nnodes_per_elem, num_local_nodes
        integer(ip), allocatable                :: temp_connectivity(:,:), temp_send(:,:)

        ! Compute partition indices for the current rank
        call compute_partition(mesh, rank, nprocs, start_elem, end_elem)
        num_local_elems = end_elem - start_elem + 1

        ! Broadcast nnodes_per_elem to all ranks
        if (rank == 0) then
            nnodes_per_elem = size(mesh % connectivity, 2)
        end if

        call MPI_Bcast(nnodes_per_elem, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        ! Allocate local mesh storage
        allocate(mesh % connectivity(num_local_elems, nnodes_per_elem))

        call compute_local_nodes(mesh, rank, nprocs, num_local_nodes)

        allocate(mesh % coords(num_local_nodes, mesh % ndim))
        allocate(mesh % boundary(mesh % nboun))

        print*, "Rank ", rank, ": Allocated local connectivity, coordinates, and boundary."

        ! Root rank prepares connectivity data
        if (rank == 0) then
            allocate(temp_connectivity(size(mesh % connectivity, 1), size(mesh % connectivity, 2)))
            temp_connectivity = mesh % connectivity
        end if

        ! Rank 0 sends partitions
        if (rank == 0) then
            do i = 1, nprocs - 1
                call compute_partition(mesh, i, nprocs, start_elem, end_elem)
                sendcount = (end_elem - start_elem + 1) * nnodes_per_elem
                print*, "Rank 0 sending ", sendcount, " elements to rank ", i

                allocate(temp_send(end_elem - start_elem + 1, nnodes_per_elem))
                temp_send = temp_connectivity(start_elem:end_elem, :)

                call MPI_Send(temp_send, sendcount, MPI_INTEGER, i, 0, MPI_COMM_WORLD, ierr)

                deallocate(temp_send)
            end do

            deallocate(temp_connectivity)
        else
            ! Non-root ranks receive their partitions
            call receive_partition(mesh, rank, nprocs)
        end if

        call MPI_Barrier(MPI_COMM_WORLD, ierr)

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
        integer(ip), intent(in)                 :: rank, nprocs
        integer(ip)                             :: num_local_elems, nnodes_per_elem, ierr
        integer(ip)                             :: start_elem, end_elem, num_local_nodes

        call compute_partition(mesh, rank, nprocs, start_elem, end_elem)
        num_local_elems = end_elem - start_elem + 1

        call MPI_Bcast(nnodes_per_elem, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        allocate(mesh % connectivity(num_local_elems, nnodes_per_elem))

        call compute_local_nodes(mesh, rank, nprocs, num_local_nodes)
        allocate(mesh % coords(num_local_nodes, mesh % ndim))

        call MPI_Recv(mesh % connectivity, num_local_elems * nnodes_per_elem, &
                      MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

        print*, "Rank ", rank, ": Received partition from Rank 0"
    end subroutine receive_partition

end module mod_parallel
