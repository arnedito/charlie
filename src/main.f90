    program main
        use constants
        use mod_mpi
        use mod_mesh
        use mod_solver
        use mod_post
        implicit none

        integer :: rank, nprocs, ierr
        type(Mesh) :: mesh
        real(rp), allocatable :: sol(:)
        real(rp) :: dt
        integer :: nsteps

        ! Initialize MPI.
        call mpi_init_wrapper(rank, nprocs)
        if (rank == 0) then
            print*, 'MPI initialized with ', nprocs, ' processes.'
        end if

        ! Read the mesh (for demonstration, using the GMSH reader on rank 0).
        if (rank == 0) then
            call read_gmsh(mesh, 'mesh.msh')
        end if

        ! In a complete code, partition and broadcast the mesh to all MPI ranks.

        ! Initialize the state vector (e.g., velocity/pressure field).
        allocate(sol(mesh % npoin))
        sol = 0.0_rp

        ! Set time step and number of time steps.
        dt = 0.01_rp
        nsteps = 5_ip

        ! Perform time integration (calls Newton solver internally).
        call time_integration(mesh, state, dt, nsteps)

        ! Write output for post-processing (each rank writes its own file).
        call write_output(mesh, state, nsteps, rank)

        ! Finalize MPI.
        call mpi_finalize_wrapper()
    end program main
