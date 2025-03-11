!**********************************************************************
! Module: mod_partition
!
! Purpose:
!   This module computes the partition indices for a mesh. It divides
!   the mesh elements among the available processes. Note that this
!   module does not perform any MPI communication to avoid circular
!   dependencies.
!
! Usage:
!   Call the subroutine 'compute_partition' with the mesh and MPI rank
!   details to get the start and end element indices for the process.
!**********************************************************************

module mod_partition
    use mod_mesh, only: mesh_type
    implicit none

    public :: compute_partition, compute_local_nodes

contains

    subroutine compute_partition(mesh, rank, nprocs, start_elem, end_elem)
        !------------------------------------------------------------------
        ! Inputs:
        !   mesh   - The mesh data structure containing the total number of elements.
        !   rank   - The current process's rank.
        !   nprocs - Total number of processes.
        !
        ! Outputs:
        !   start_elem - The starting element index for this process.
        !   end_elem   - The ending element index for this process.
        !------------------------------------------------------------------
        type(mesh_type), intent(in)         :: mesh
        integer, intent(in)                 :: rank, nprocs
        integer, intent(out)                :: start_elem, end_elem
        integer                             :: num_local_elems

        ! Calculate the number of elements each process should get.
        num_local_elems = mesh % nelem / nprocs
        start_elem = rank * num_local_elems + 1
        end_elem = start_elem + num_local_elems - 1

        ! If this is the last process, assign any remaining elements.
        if (rank == nprocs - 1) then
            end_elem = mesh % nelem
        end if
    end subroutine compute_partition

    !----------------------------------------------------------------------
    ! Subroutine: compute_local_nodes
    ! Purpose   : Computes the number of unique local nodes used in a partition.
    !----------------------------------------------------------------------
    subroutine compute_local_nodes(mesh, rank, nprocs, num_local_nodes)
        type(mesh_type), intent(in)         :: mesh
        integer, intent(in)                 :: rank, nprocs
        integer, intent(out)                :: num_local_nodes
        integer, allocatable                :: node_marker(:)
        integer                             :: i, j, node, start_elem, end_elem

        call compute_partition(mesh, rank, nprocs, start_elem, end_elem)

        allocate(node_marker(mesh % npoin))
        node_marker = 0

        do i = start_elem, end_elem
            do j = 1, size(mesh % connectivity, 2)
                node = mesh % connectivity(i, j)
                node_marker(node) = 1
            end do
        end do

        num_local_nodes = sum(node_marker)

        deallocate(node_marker)
    end subroutine compute_local_nodes


end module mod_partition
