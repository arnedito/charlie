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
    public :: compute_partition

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

end module mod_partition
