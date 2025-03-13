!=============================================================================
! @file   mod_partition.f90
! @brief  Provides routines to partition mesh elements among MPI processes.
!
! @details
!   This module divides the total number of elements in a mesh among multiple
!   processes for parallel computations. The partition is currently performed
!   by a simple block distribution: each rank gets a contiguous block of elements.
!   It also offers a routine to compute how many unique nodal points belong to
!   a given sub-partition (without doing any direct MPI communication).
!
!   @note
!   - This module does not perform MPI communications directly to avoid circular
!     dependencies. MPI rank and size are passed in as arguments.
!   - By default, leftover elements (if any) are all assigned to the last rank.
!     More advanced partitioning strategies could distribute leftover elements
!     more evenly.
!
!=============================================================================
module mod_partition
   use constants, only: ip, rp
   use mod_mesh, only: mesh_type
   implicit none

   !--------------------------------------------------------------------------
   !> Public subroutines
   !--------------------------------------------------------------------------
   public :: compute_partition
   public :: compute_local_nodes
   public :: local_elements

contains

   !---------------------------------------------------------------------------
   !> @brief Compute the block-partition indices for a mesh among multiple ranks.
   !>
   !! Given the total number of elements (mesh%nelem) and the current rank
   !! plus total process count, this subroutine calculates a contiguous block
   !! of elements [start_elem, end_elem] for the local process.
   !!
   !! @param[in]  mesh       (mesh_type)  A mesh structure (read-only).
   !! @param[in]  rank       (integer(ip)) Rank of the calling process.
   !! @param[in]  nprocs     (integer(ip)) Total number of processes.
   !! @param[out] start_elem (integer(ip)) The first element index assigned to this rank.
   !! @param[out] end_elem   (integer(ip)) The last element index assigned to this rank.
   !
   !! @note
   !!  - The last rank gets all leftover elements if nelem is not perfectly
   !!    divisible by nprocs.
   !---------------------------------------------------------------------------
    subroutine compute_partition(mesh, rank, nprocs, start_elem, end_elem)
      type(mesh_type), intent(in)  :: mesh
      integer(ip), intent(in)      :: rank, nprocs
      integer(ip), intent(out)     :: start_elem, end_elem
      integer(ip)                  :: num_local_elems

      ! Number of elements per process (integer division)
      num_local_elems = mesh%nelem / nprocs

      start_elem = rank * num_local_elems + 1_ip
      end_elem   = start_elem + num_local_elems - 1_ip

      ! If this is the last process, take the remaining elements
      if (rank == nprocs - 1_ip) then
         end_elem = mesh%nelem
      end if

      ! Possible Improvement:
      !   - If you want a more balanced approach for leftover elements, you
      !     could distribute them among the first 'remainder' ranks instead
      !     of dumping all on the last rank. For example:
      !
      !       remainder = mod(mesh%nelem, nprocs)
      !       if (rank < remainder) then
      !         ! each rank < remainder gets an extra element
      !       end if
      !
   end subroutine compute_partition

   !---------------------------------------------------------------------------
   !> @brief Determines the number of unique local nodes for a given sub-partition.
   !>
   !! This subroutine first calls compute_partition to find the element range
   !! [start_elem, end_elem] for a specific rank, then iterates through the
   !! connectivity of those elements. A simple marker array (node_marker) is
   !! used to count unique nodes.
   !!
   !! @param[in]  mesh            (mesh_type)  The global mesh (read-only).
   !! @param[in]  rank            (integer(ip)) Rank of the calling process.
   !! @param[in]  nprocs          (integer(ip)) Total number of processes.
   !! @param[out] num_local_nodes (integer(ip)) The count of unique nodes in
   !!                                            this rank's sub-partition.
   !
   !! @note
   !!  - Because no actual MPI communication is done here, the "local nodes"
   !!    concept is purely based on the block partition; real parallel codes
   !!    would also need to gather or communicate node ownership across ranks.
   !---------------------------------------------------------------------------
   subroutine compute_local_nodes(mesh, rank, nprocs, num_local_nodes)
      type(mesh_type), intent(in)   :: mesh
      integer(ip), intent(in)       :: rank, nprocs
      integer(ip), intent(out)      :: num_local_nodes

      integer(ip), allocatable      :: node_marker(:)
      integer(ip)                   :: i, j
      integer(ip)                   :: node
      integer(ip)                   :: start_elem, end_elem

      ! Determine which elements this rank owns
      call compute_partition(mesh, rank, nprocs, start_elem, end_elem)

      ! Mark array used to track which nodes are present
      allocate(node_marker(mesh%npoin))
      node_marker = 0_ip

      do i = start_elem, end_elem
         do j = 1_ip, size(mesh%connectivity, 2)
            node = mesh%connectivity(i, j)
            node_marker(node) = 1_ip
         end do
      end do

      ! The sum of the marker array gives the count of unique nodes
      num_local_nodes = sum(node_marker)

      deallocate(node_marker)
   end subroutine compute_local_nodes

    !---------------------------------------------------------------------
    !> @brief Returns the number of elements in mesh%connectivity
    !!       (i.e., how many local elements are allocated).
    !---------------------------------------------------------------------
    function local_elements(mesh) result(nlocal)
        type(mesh_type), intent(in) :: mesh
        integer(ip)                 :: nlocal

        if (.not. allocated(mesh%connectivity)) then
            nlocal = 0_ip
        else
            nlocal = size(mesh%connectivity, 1)
        end if
    end function local_elements

end module mod_partition
