!**********************************************************************
! Module: mod_mpi
!
! Purpose:
!   This module provides a wrapper around basic MPI functionality such as
!   initialization, finalization, and barriers. It isolates MPI-specific
!   calls from the rest of the application.
!
! Usage:
!   Call 'mpi_init_wrapper' at the start and 'mpi_finalize_wrapper' at
!   the end of your program. Other MPI utilities are also available.
!**********************************************************************
module mod_mpi
    use mpi
    implicit none
    public :: mpi_init_wrapper, mpi_finalize_wrapper, mpi_barrier_wrapper
    public :: get_mpi_comm, mpirank, mpisize, mpi_comm_global

    integer :: mpi_comm_global = MPI_COMM_WORLD
    integer :: mpirank = -1
    integer :: mpisize = -1
    integer :: mpierr

contains

    subroutine mpi_init_wrapper()
        ! Initialize MPI and retrieve rank and size.
        call MPI_Init(mpierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpierr)
        call MPI_Comm_size(MPI_COMM_WORLD, mpisize, mpierr)
    end subroutine mpi_init_wrapper

    subroutine mpi_finalize_wrapper()
        ! Finalize the MPI environment.
        call MPI_Finalize(mpierr)
    end subroutine mpi_finalize_wrapper

    subroutine mpi_barrier_wrapper()
        ! Perform a barrier synchronization among all processes.
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    end subroutine mpi_barrier_wrapper

    function get_mpi_comm() result(comm)
        ! Return the global MPI communicator.
        integer :: comm
        comm = mpi_comm_global
    end function get_mpi_comm

end module mod_mpi
