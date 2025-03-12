!=============================================================================
! @file   mod_mpi.f90
! @brief  Provides basic MPI initialization, finalization, and utility wrappers.
!
! @details
!   This module isolates MPI-specific calls from the rest of the application.
!   It manages:
!   - Global communicator (mpi_comm_global)
!   - Process rank (mpirank)
!   - Total number of processes (mpisize)
!
!   Routines are provided to initialize and finalize MPI, and to perform
!   barrier synchronizations. You can extend this module with further wrappers
!   (e.g. send/recv, broadcast) to reduce direct MPI calls throughout your code.
!
!=============================================================================
module mod_mpi
   use mpi
   implicit none

   public :: mpi_init_wrapper
   public :: mpi_finalize_wrapper
   public :: mpi_barrier_wrapper
   public :: get_mpi_comm

   public :: mpirank, mpisize, mpi_comm_global

   !-------------------------------------------------------------------------
   !> Global MPI communicator and process info
   !-------------------------------------------------------------------------
   integer :: mpi_comm_global = MPI_COMM_WORLD
   integer :: mpirank         = -1
   integer :: mpisize         = -1

   integer :: mpierr          = 0  ! for storing MPI error codes

contains

   !-------------------------------------------------------------------------
   !> @brief Initializes the MPI environment.
   !>
   !! This subroutine wraps the MPI_Init call, sets up the global communicator,
   !! and retrieves the current process rank and total size.
   !-------------------------------------------------------------------------
   subroutine mpi_init_wrapper()
      ! Initialize MPI and retrieve rank and size
      call MPI_Init(mpierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpierr)
      call MPI_Comm_size(MPI_COMM_WORLD, mpisize, mpierr)

      ! Possible Improvement:
      !  - Check mpierr for errors and handle them (e.g., print/log message or
      !    abort). Currently, we do not validate mpierr.
   end subroutine mpi_init_wrapper

   !-------------------------------------------------------------------------
   !> @brief Finalizes the MPI environment.
   !>
   !! This subroutine calls MPI_Finalize, indicating that no further MPI calls
   !! will be made in the application.
   !-------------------------------------------------------------------------
   subroutine mpi_finalize_wrapper()
      call MPI_Finalize(mpierr)
      ! Possible Improvement:
      !  - Same as above, you could check mpierr for potential errors.
   end subroutine mpi_finalize_wrapper

   !-------------------------------------------------------------------------
   !> @brief Performs a barrier synchronization among all processes.
   !-------------------------------------------------------------------------
   subroutine mpi_barrier_wrapper()
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      ! This ensures that all processes reach this point before any can proceed.
   end subroutine mpi_barrier_wrapper

   !-------------------------------------------------------------------------
   !> @brief Returns the global MPI communicator in use by this application.
   !>
   !! @return comm (integer)
   !!   The handle to the global MPI communicator (MPI_COMM_WORLD by default).
   !-------------------------------------------------------------------------
   function get_mpi_comm() result(comm)
      integer :: comm
      comm = mpi_comm_global
   end function get_mpi_comm

end module mod_mpi
