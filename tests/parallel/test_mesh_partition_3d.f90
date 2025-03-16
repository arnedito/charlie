!=============================================================================
! @file   test_mesh_partition_3d.f90
! @brief  Tests parallel mesh partitioning for a 3D tetrahedral mesh.
!
! @details
!   Steps:
!   1) Initializes MPI.
!   2) Rank 0 reads the 3D mesh (full connectivity).
!   3) Broadcasts metadata and (optionally) full coords/boundaries to all ranks.
!   4) Calls distribute_mesh() to assign a local sub-block of the connectivity array.
!   5) Verifies partitioning by checking element distribution across ranks.
!
!=============================================================================
program test_mesh_partition_3d
   use mpi
   use mod_mpi,       only: mpi_init_wrapper, mpi_finalize_wrapper, mpi_comm_global, mpirank, mpisize
   use mod_partition, only: local_elements
   use mod_parallel,  only: distribute_mesh
   use mod_mesh,      only: mesh_type, read_domain
   use constants,     only: ip, rp
   implicit none

   type(mesh_type)      :: mesh
   integer(ip)          :: rank, size
   integer(ip)          :: total_local_elems, expected_total
   integer              :: ierr
   integer, allocatable :: recv_counts(:)

   !-------------------------------------------------------------------------
   !> Initialize MPI environment
   !-------------------------------------------------------------------------
   call mpi_init_wrapper()
   rank = mpirank
   size = mpisize

   !-------------------------------------------------------------------------
   !> On rank 0, read the full 3D mesh
   !-------------------------------------------------------------------------
   if (rank == 0_ip) then
      print *, "Rank 0: Reading 3D mesh..."
      call read_domain("tests/mesh/read_mesh_3d/case.dom.dat", &
                       "tests/mesh/read_mesh_3d/case.geo.dat", &
                       "tests/mesh/read_mesh_3d/case.fix.dat", &
                       mesh)
      ! At this point, rank 0 has a full mesh%lnods of size (4, nelem) for tetrahedral elements.
   end if

   !-------------------------------------------------------------------------
   !> Broadcast basic mesh metadata
   !-------------------------------------------------------------------------
   call MPI_Bcast(mesh%npoin, 1, MPI_INTEGER, 0, mpi_comm_global, ierr)
   call MPI_Bcast(mesh%nelem, 1, MPI_INTEGER, 0, mpi_comm_global, ierr)
   call MPI_Bcast(mesh%nboun, 1, MPI_INTEGER, 0, mpi_comm_global, ierr)
   call MPI_Bcast(mesh%ndim,  1, MPI_INTEGER, 0, mpi_comm_global, ierr)

   !-------------------------------------------------------------------------
   !> Optionally, broadcast coords & boundary to all ranks
   !-------------------------------------------------------------------------
   if (rank /= 0_ip) then
      allocate(mesh%coords(mesh%npoin, mesh%ndim))
      allocate(mesh%boundary(mesh%nboun))
   end if
   call MPI_Bcast(mesh%coords,   product(shape(mesh%coords)),   MPI_REAL,    0, mpi_comm_global, ierr)
   call MPI_Bcast(mesh%boundary, product(shape(mesh%boundary)), MPI_INTEGER, 0, mpi_comm_global, ierr)

   !-------------------------------------------------------------------------
   !> Distribute connectivity so each rank gets a sub-block
   !-------------------------------------------------------------------------
   call distribute_mesh(mesh, rank, size)

   ! Synchronize
   call MPI_Barrier(mpi_comm_global, ierr)

   !-------------------------------------------------------------------------
   !> Each rank now has a local connectivity array. Let's see how many we have.
   !-------------------------------------------------------------------------
   if (allocated(mesh%lnods)) then
      total_local_elems = local_elements(mesh)
   else
      total_local_elems = 0_ip
   end if

   if (total_local_elems > 0_ip) then
      print *, 'Rank', rank, 'received', total_local_elems, 'elements.'
   else
      print *, 'Warning: Rank', rank, 'received no elements!'
   end if

   !-------------------------------------------------------------------------
   !> Gather local element counts on rank 0
   !-------------------------------------------------------------------------
   allocate(recv_counts(size))
   call MPI_Gather(total_local_elems, 1, MPI_INTEGER, &
                   recv_counts, 1, MPI_INTEGER, 0, mpi_comm_global, ierr)

   if (rank == 0_ip) then
      expected_total = sum(recv_counts)
      if (expected_total == mesh%nelem) then
         print *, 'Mesh successfully partitioned: total elements =', expected_total
      else
         print *, 'Error: Distributed elements do not match original!'
         print *, '   Expected:', mesh%nelem, 'but got:', expected_total
      end if
   end if

   !-------------------------------------------------------------------------
   !> Optionally, print first element's connectivity on each rank
   !-------------------------------------------------------------------------
   if (total_local_elems > 0_ip .and. allocated(mesh%lnods)) then
      print *, 'Rank', rank, 'first local element:', mesh%lnods(:, 1)
   else
      print *, 'Rank', rank, ': No elements to print.'
   end if

   !-------------------------------------------------------------------------
   !> Clean up and finalize
   !-------------------------------------------------------------------------
   if (allocated(recv_counts)) deallocate(recv_counts)
   call mpi_finalize_wrapper()

end program test_mesh_partition_3d
