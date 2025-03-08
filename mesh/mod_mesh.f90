    module mod_mesh
        use constants
        implicit none
        public :: Mesh, read_gmsh, read_fensap

        ! mesh data structure (supports 1D, 2D, and 3D meshes)
        type :: Mesh
            integer(ip) :: ndim = 0                         ! dimension: 1, 2, or 3
            integer(ip) :: npoin = 0                        ! total number of nodes (gpoin)
            integer(ip) :: nelem = 0                        ! total number of elements (nelem)
            integer(ip) :: nnode_per_elem = 0               ! nodes per element (e.g., 3 for triangles)
            real(rp), allocatable :: coords(:,:)            ! node coordinates: coords(npoin, ndim)
            integer(ip), allocatable    ::  connect(:,:)    ! element connectivity: connect(nelem, nnode_per_elem)
        end type Mesh

    contains

        subroutine read_gmsh(mesh, filename)
            type(Mesh), intent(inout) :: mesh
            character(len=*), intent(in) :: filename
            ! Dummy implementation: Replace with actual GMSH file parsing.
            print*, 'Reading GMSH mesh from file: ', trim(filename)
            ! For demonstration, create a simple 2D mesh of 4 nodes and 2 triangular elements.
            mesh % ndim = 2
            mesh % npoin = 4
            mesh % nelem = 2
            mesh % nnode_per_elem = 3  ! Triangular elements
            allocate(mesh % coords(mesh % npoin, mesh % ndim))
            allocate(mesh % connect(mesh % nelem, mesh % nnode_per_elem))
            ! Define dummy node coordinates
            mesh % coords = reshape([0.0_rp, 0.0_rp, &
                                1.0_rp, 0.0_rp, &
                                1.0_rp, 1.0_rp, &
                                0.0_rp, 1.0_rp], shape(mesh % coords))

            ! Define dummy element connectivity (using 1-based indexing)
            mesh % connect = reshape([1, 2, 3, &
                                1, 3, 4], shape(mesh%connect))
        end subroutine read_gmsh

        subroutine read_fensap(mesh, filename)
            type(Mesh), intent(inout) :: mesh
            character(len=*), intent(in) :: filename
            ! Dummy implementation: Replace with actual FENSAP file parsing.
            print*, 'Reading FENSAP mesh from file: ', trim(filename)
            ! For demonstration, call read_gmsh (assuming similar structure).
            call read_gmsh(mesh, filename)
        end subroutine read_fensap

    end module mesh_module
