program test_mpi
    use mpi
    use mod_mpi
    implicit none

    integer :: test_error = 0, ierr
    integer :: msg, status(MPI_STATUS_SIZE)
    integer :: mpi_err

    call mpi_init_wrapper()

    ! ------------------- MPI Process Check -------------------
    print *, "Process ", mpirank, " of ", mpisize, " initialized."

    if (mpisize < 2) then
        print *, "Error: This test requires at least 2 MPI processes."
        test_error = 1
    end if

    ! ------------------- MPI Barrier Test -------------------
    print *, "Process ", mpirank, " reaching barrier..."
    call mpi_barrier_wrapper()
    print *, "Process ", mpirank, " passed the barrier."

    ! ------------------- MPI Message Passing Test -------------------
    if (mpirank == 0) then
        ! Master process sends message to all workers
        msg = 42
        do ierr = 1, mpisize - 1
            call MPI_Send(msg, 1, MPI_INTEGER, ierr, 0, MPI_COMM_WORLD, mpi_err)
        end do
        print *, "Master process sent message:", msg
    else
        ! Workers receive message from master
        call MPI_Recv(msg, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, status, mpi_err)
        if (msg /= 42) then
            print *, "Error: Process ", mpirank, " received incorrect message!"
            test_error = 1
        else
            print *, "Process ", mpirank, " received message:", msg
        end if
    end if

    ! ------------------- Finalize MPI -------------------
    if (test_error == 0) then
        print *, "MPI Test passed successfully on process ", mpirank
    else
        print *, "MPI Test failed on process ", mpirank
    end if

    call mpi_finalize_wrapper()

end program test_mpi
