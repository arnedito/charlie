    module mod_mpi
        use constants
        use mpi
        implicit none
        public :: mpi_init_wrapper, mpi_finalize_wrapper, get_mpi_comm

        integer :: mpi_comm_global = MPI_COMM_WORLD

    contains

        subroutine mpi_init_wrapper(rank, nprocs)
            integer, intent(out) :: rank, nprocs
            integer :: ierr
            call MPI_Init(ierr)
            call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
            call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
        end subroutine mpi_init_wrapper

        subroutine mpi_finalize_wrapper()
            integer :: ierr
            call MPI_Finalize(ierr)
        end subroutine mpi_finalize_wrapper

        function get_mpi_comm() result(comm)
            integer :: comm
            comm = MPI_COMM_WORLD
        end function get_mpi_comm

    end module mpi_module
