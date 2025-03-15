!=============================================================================
! @file   mod_parallel.f90
! @brief  Contains subroutines that handle MPI-based data distribution for a mesh.
!
! @details
!   - The main routine here is `distribute_mesh`, which uses a simple "block"
!     partitioning and physically redistributes the connectivity array so each
!     rank only holds the subset of elements assigned to it.
!   - Rank 0 holds the original global connectivity, then sends each sub-block
!     to the appropriate rank. Rank 0 also re-allocates a smaller array for its
!     own local portion.
!
! @note
!   - This example only distributes the connectivity. If you want partial
!     coordinates, you'd do a similar approach for `coords` – but that
!     requires re-indexing node IDs for each rank.
!   - Make sure you do not deallocate the global connectivity on Rank 0 before
!     you finish sending. A separate array or a careful copy is needed.
!
!=============================================================================
module mod_parallel
   use mpi
   use constants,     only: ip, rp
   use mod_mesh,      only: mesh_type
   use mod_partition, only: compute_partition
   implicit none

   public :: distribute_mesh

contains

   !---------------------------------------------------------------------------
   !> @brief Distributes the global mesh connectivity from Rank 0 to other ranks.
   !>
   !! This routine:
   !! 1) Rank 0 has a global connectivity array (stored in mesh%connectivity).
   !! 2) For each rank, compute a block partition of elements: [start_elem, end_elem].
   !! 3) Send the relevant sub-block of connectivity to that rank (via MPI_Send).
   !! 4) The receiving rank (including rank 0) deallocates any existing connectivity,
   !!    allocates a local connectivity array sized to (num_local_elems, ndim+1),
   !!    and receives (or copies) the data.
   !!
   !! @param[inout] mesh (mesh_type)  The mesh structure to redistribute.
   !! @param[in]    rank (integer(ip)) The current MPI rank.
   !! @param[in]    size (integer(ip)) Total number of MPI ranks.
   !
   !! @warning
   !! - If `mesh%connectivity` on Rank 0 is deallocated before the Send calls
   !!   complete, you'll lose the data. We handle that carefully below.
   !---------------------------------------------------------------------------
   subroutine distribute_mesh(mesh, rank, size)
      use mpi
      type(mesh_type), intent(inout)        :: mesh
      integer(ip), intent(in)               :: rank, size
      integer(ip)                           :: start_elem, end_elem, num_local_elems
      integer(ip)                           :: i, j
      integer                               :: ierr
      integer(ip)                           :: r, s_elem, e_elem, n_elems
      integer(ip)                           :: send_count, recv_count
      integer(ip), allocatable              :: tmp(:,:)

      !-----------------------------------------------------------------------
      ! 1) Determine this rank's element range with compute_partition.
      !    We'll do the same for all ranks from 0..size-1 on Rank 0.
      !-----------------------------------------------------------------------
      call compute_partition(mesh, rank, size, start_elem, end_elem)
      num_local_elems = end_elem - start_elem + 1_ip

      !-----------------------------------------------------------------------
      ! 2) If rank = 0, loop over all ranks to send sub-blocks of connectivity.
      !    Otherwise, each rank receives its sub-block.
      !-----------------------------------------------------------------------
      if (rank == 0_ip) then
         ! On Rank 0, we still have the full mesh%connectivity from read_domain.
         ! We'll iterate over all ranks and send their portion.

         do r = 1_ip, size-1_ip
            call compute_partition(mesh, r, size, s_elem, e_elem)
            n_elems = e_elem - s_elem + 1_ip

            if (n_elems > 0_ip) then
               send_count = n_elems * (mesh%ndim+1_ip)
               ! Temporarily allocate a buffer to hold that chunk
               allocate(tmp(n_elems, mesh%ndim+1_ip))

               ! Copy the chunk from the big global connectivity
               do i = 1, n_elems
                  tmp(i,:) = mesh%lnods(:, s_elem + i - 1_ip)
               end do

               ! Send to rank r
               call MPI_Send(tmp, send_count, MPI_INTEGER, r, 0, MPI_COMM_WORLD, ierr)

               deallocate(tmp)
            else
               ! If n_elems=0, we don't send anything
            end if
         end do

         !---------------------------------------------------------------------
         ! Now handle rank 0's local portion
         ! We'll copy that sub-block for ourselves.
         !  - First, store sub-block in a temp array
         !  - Deallocate the giant connectivity
         !  - Reallocate a smaller local array
         !  - Copy from temp into the new local connectivity
         !---------------------------------------------------------------------
         if (num_local_elems > 0_ip) then
            allocate(tmp(num_local_elems, mesh%ndim+1_ip))

            do i = 1_ip, num_local_elems
               tmp(i,:) = mesh%lnods(:, start_elem + i - 1_ip)
            end do

            ! Deallocate the original big array
            if (allocated(mesh%lnods)) deallocate(mesh%lnods)

            ! Reallocate for local portion
            allocate(mesh%lnods(mesh%ndim+1_ip, num_local_elems))

            ! Copy from tmp
            mesh%lnods = tmp
            deallocate(tmp)
         else
            ! If no local elems, just deallocate everything
            if (allocated(mesh%lnods)) deallocate(mesh%lnods)
         end if

      else
         !-------------------------------------------------------------
         ! RANK != 0
         !  - Deallocate any existing connectivity (if it was allocated).
         !  - Allocate local connectivity for our share.
         !  - Receive the sub-block from Rank 0.
         !-------------------------------------------------------------
         if (allocated(mesh%lnods)) deallocate(mesh%lnods)

         if (num_local_elems > 0_ip) then
            allocate(mesh%lnods(mesh%ndim+1_ip, num_local_elems))
            recv_count = num_local_elems * (mesh%ndim+1_ip)
            call MPI_Recv(mesh%lnods, recv_count, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         end if
      end if

end subroutine distribute_mesh

end module mod_parallel
