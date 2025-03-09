    !-----------------------------------------------------------------
    ! Module: mod_mpi
    ! Purpose: Provides a wrapper for MPI initialization, finalization,
    !          and basic synchronization. It defines global variables
    !          to store the process rank, the number of processes, and
    !          the global communicator.
    !-----------------------------------------------------------------
    module mod_mpi
        use mpi
        implicit none

        ! Public entities accessible from other modules
        public :: mpi_init_wrapper, mpi_finalize_wrapper, mpi_barrier_wrapper
        public :: get_mpi_comm, mpirank, mpisize, mpi_comm_global

        ! Global communicator for all MPI processes (default: MPI_COMM_WORLD)
        integer :: mpi_comm_global = MPI_COMM_WORLD

        ! Variables to store the rank (ID) of the current process and the total number of processes
        integer :: mpirank = -1
        integer :: mpisize = -1

        ! Variable to capture MPI error codes for error checking
        integer :: mpierr

    contains

        !-----------------------------------------------------------------
        subroutine mpi_init_wrapper()
            !-------------------------------------------------------------
            ! Purpose: Initializes the MPI environment.
            !          - Calls MPI_Init to set up the MPI library.
            !          - Retrieves the rank of the current process.
            !          - Retrieves the total number of processes.
            !          - Checks for errors during initialization and stops
            !            execution if any error occurs.
            !-------------------------------------------------------------
            call MPI_Init(mpierr)
            if (mpierr /= MPI_SUCCESS) then
                print *, "Error initializing MPI!"
                stop
            end if
            call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpierr)
            call MPI_Comm_size(MPI_COMM_WORLD, mpisize, mpierr)

            ! Informative output to confirm successful initialization
            print *, "Process ", mpirank, " of ", mpisize, " initialized."
        end subroutine mpi_init_wrapper
        !-----------------------------------------------------------------

        !-----------------------------------------------------------------
        subroutine mpi_finalize_wrapper()
            !-------------------------------------------------------------
            ! Purpose: Finalizes the MPI environment.
            !          - Calls MPI_Finalize to properly shut down MPI.
            !          - Checks for errors during finalization.
            !-------------------------------------------------------------
            call MPI_Finalize(mpierr)
            if (mpierr /= MPI_SUCCESS) then
                print *, "Error finalizing MPI!"
            end if
        end subroutine mpi_finalize_wrapper
        !-----------------------------------------------------------------

        !-----------------------------------------------------------------
        subroutine mpi_barrier_wrapper()
            !-------------------------------------------------------------
            ! Purpose: Implements an MPI Barrier to synchronize processes.
            !          - Ensures that all processes reach this point before
            !            any of them can proceed.
            !-------------------------------------------------------------
            call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        end subroutine mpi_barrier_wrapper
        !-----------------------------------------------------------------

        !-----------------------------------------------------------------
        function get_mpi_comm() result(comm)
            !-------------------------------------------------------------
            ! Purpose: Returns the global MPI communicator.
            ! Returns:
            !    comm - The global MPI communicator (typically MPI_COMM_WORLD)
            !-------------------------------------------------------------
            integer :: comm
            comm = mpi_comm_global
        end function get_mpi_comm
        !-----------------------------------------------------------------

    end module mod_mpi
