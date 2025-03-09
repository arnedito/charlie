    !-----------------------------------------------------------------
    ! Module: mod_parallel
    ! Purpose: Provides a basic master-worker structure utilizing MPI.
    !          This module leverages the MPI functionalities defined in
    !          mod_mpi to implement parallel tasks distribution.
    !          The master process (rank 0) is responsible for distributing
    !          tasks (such as mesh distribution), while worker processes wait
    !          for tasks.
    !-----------------------------------------------------------------

    module mod_parallel
        use mod_mpi  ! Import MPI functionalities from mod_mpi
        implicit none
        public :: master_worker  ! Make the master_worker subroutine accessible

    contains

        !-----------------------------------------------------------------
        subroutine master_worker()
            !-------------------------------------------------------------
            ! Purpose: Implements the master-worker paradigm.
            !          - If the current process is the master (rank 0), it will
            !            distribute tasks (e.g., mesh data).
            !          - Worker processes (non-zero rank) will wait for tasks.
            !          - A barrier is used to synchronize all processes before
            !            proceeding further.
            !-------------------------------------------------------------
            if (mpirank == 0) then
                print *, 'Master process (rank ', mpirank, '): distributing tasks...'
                ! Future implementation: Call subroutine to distribute the mesh or tasks.
            else
                print *, 'Worker process (rank ', mpirank, '): waiting for tasks...'
                ! Future implementation: Call subroutine to receive and process tasks.
            end if

            ! Synchronize all processes to ensure all have reached this point.
            call mpi_barrier_wrapper()
        end subroutine master_worker
        !-----------------------------------------------------------------

    end module mod_parallel
