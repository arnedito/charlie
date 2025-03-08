    module mod_mesh
        use constants, only: ip, rp
        implicit none
        public :: Mesh, read_domain

        !>----------------------------------------------------------------------
        !> Unified Mesh Type
        !>----------------------------------------------------------------------
        type :: mesh_type
            integer(ip)                 :: ndim = 0             ! Spatial dimension (e.g., 2 or 3)
            integer(ip)                 :: npoin = 0            ! Total number of nodes
            integer(ip)                 :: nelem = 0            ! Total number of elements
            integer(ip)                 :: nboun = 0            ! Total number of boundaries
            real(rp), allocatable       :: coords(:,:)          ! Node coordinates: (npoin x ndim)
            integer(ip), allocatable    :: connectivity(:,:)    ! Element connectivity: (nelem x nodes_per_element)
            integer(ip), allocatable    :: boundary(:)          ! Boundary condition codes (size: nboun)
        end type mesh_type

    contains

        !>-----------------------------------------------------------------------
        !> read_domain: Reads the unified mesh format using case.dom.dat, case.geo.dat,
        !> and case.fix.dat. It extracts metadata, geometry, and boundary conditions.
        !>-----------------------------------------------------------------------

        subroutine read_domain(dom_filename, geo_filename, fix_filename, mesh)
            character(len=*), intent(in) :: dom_filename, geo_filename, fix_filename
            type(mesh_type), intent(out) :: mesh
            character(len=256) :: line, dummy
            integer(ip) :: ios

            ! Initialize temporary values
            mesh % npoin = 0
            mesh % nelem = 0
            mesh % nboun = 0
            mesh % ndim = 0

            !--- Read domain file: extract metadata ---
            open(unit=10, file=trim(dom_filename), status="old", action="read", iostat=ios)
            if (ios /= 0) then
                print*, "Error opening ", trim(dom_filename)
                stop
            end if

            do
                read(10, '(A)', iostat=ios) line
                if (ios /= 0) exit
                if (index(line, "NODAL_POINTS") > 0) then
                    read(line, *) dummy, mesh % npoin
                else if (index(line, "ELEMENTS") > 0) then
                    read(line, *) dummy, mesh % nelem
                else if (index(line, "BOUNDARIES") > 0) then
                    read(line, *) dummy, mesh % nboun
                else if (index(line, "SPACE_DIMENSIONS") > 0) then
                    read(line, *) dummy, mesh % ndim
                else if (index(line, "END_DIMENSIONS") > 0) then
                    exit
                end if
            end do

            close(10)

            !--- Read geometry file: node coordinates and connectivity ---
            call read_geo_dat(geo_filename, mesh)

            !--- Read fix file: boundary conditions ---
            call read_fix_dat(fix_filename, mesh)
        end subroutine read_domain

        !-----------------------------------------------------------------------
        ! read_geo_dat: Reads the COORDINATES and ELEMENTS sections from the geo file.
        !-----------------------------------------------------------------------
        subroutine read_geo_dat(geo_filename, mesh)
            character(len=*), intent(in) :: geo_filename
            type(Mesh), intent(inout) :: mesh
            character(len=256) :: line
            integer(ip) :: ios, i
            logical :: inside_coords, inside_elems

            inside_coords = .false.
            inside_elems = .false.

            open(unit=20, file=trim(geo_filename), status="old", action="read", iostat=ios)
            if (ios /= 0) then
                print*, "Error opening ", trim(geo_filename)
                stop
            end if

            ! Allocate arrays based on domain file values
            allocate(mesh % coords(mesh % npoin, mesh % ndim))
            ! Assuming triangular elements: 3 nodes per element
            allocate(mesh % connectivity(mesh % nelem, 3))

            do
                read(20, '(A)', iostat=ios) line
                if (ios /= 0) exit

                select case (trim(line))
                case ("COORDINATES")
                    inside_coords = .true.
                case ("END_COORDINATES")
                    inside_coords = .false.
                case ("ELEMENTS")
                    inside_elems = .true.
                case ("END_ELEMENTS")
                    inside_elems = .false.
                case default
                    if (inside_coords) then
                        ! Expected format: node_id x y [z]
                        integer(ip) :: node_id
                        real(rp) :: x, y, z  ! Declare all possible needed variables

                        if (mesh % ndim == 2) then
                            read(line, *) node_id, x, y
                            mesh % coords(node_id, 1) = x
                            mesh % coords(node_id, 2) = y
                        else if (mesh % ndim == 3) then
                            read(line, *) node_id, x, y, z
                            mesh % coords(node_id, 1) = x
                            mesh % coords(node_id, 2) = y
                            mesh % coords(node_id, 3) = z
                        end if

                    else if (inside_elems) then
                        ! Expected format: elem_id node1 node2 node3
                        integer(ip) :: elem_id, n1, n2, n3
                        read(line, *) elem_id, n1, n2, n3
                        mesh % connectivity(elem_id, :) = (/ n1, n2, n3 /)
                    end if
                end select
            end do

            close(20)
        end subroutine read_geo_dat

        !>-----------------------------------------------------------------------
        !> read_fix_dat: Reads the boundary conditions from the case.fix.dat file.
        !>-----------------------------------------------------------------------
        subroutine read_fix_dat(fix_filename, mesh)
            character(len=*), intent(in) :: fix_filename
            type(Mesh), intent(inout) :: mesh
            character(len=256) :: line
            integer(ip) :: ios, i
            ! For simplicity, assume fix file has a header line followed by:
            ! element_or_node_id, boundary_code
            open(unit=30, file=trim(fix_filename), status="old", action="read", iostat=ios)
            if (ios /= 0) then
                print*, "Error opening ", trim(fix_filename)
                stop
            end if

            ! Allocate boundary condition array (size: nboun)
            allocate(mesh%boundary(mesh%nboun))

            ! Skip header (e.g., "ON_BOUNDARIES")
            read(30, '(A)', iostat=ios) line

            do i = 1, mesh % nboun
                integer(ip) :: id, code
                read(30, *) id, code
                ! Here, you may wish to map id to a particular node/element.
                mesh % boundary(i) = code
            end do

            close(30)

        end subroutine read_fix_dat

end module mesh_mod_mesh
