module mod_mesh
    use iso_fortran_env
    use constants, only: ip, rp
    implicit none

    public :: mesh_type, read_domain, local_elements

    !>----------------------------------------------------------------------
    !> Data structure representing a computational mesh.
    !>----------------------------------------------------------------------
    type :: mesh_type
        integer(ip)                 :: ndim                 ! Number of spatial dimensions
        integer(ip)                 :: npoin                ! Number of nodal points
        integer(ip)                 :: nelem                ! Number of elements
        integer(ip)                 :: nboun                ! Number of boundary elements
        real(rp), allocatable       :: coords(:,:)          ! Nodal coordinates
        integer(ip), allocatable    :: connectivity(:,:)    ! Element connectivity
        integer(ip), allocatable    :: boundary(:)          ! Boundary conditions
    end type mesh_type

contains

    !>-----------------------------------------------------------------------
    !> Reads mesh data from user-provided input files.
    !>-----------------------------------------------------------------------
    subroutine read_domain(dom_filename, geo_filename, fix_filename, mesh)
        character(len=*), intent(in)            :: dom_filename, geo_filename, fix_filename
        type(mesh_type), intent(out)            :: mesh
        character(len=256)                      :: line, dummy
        integer(ip)                             :: ios

        print *, "Opening domain file: ", trim(dom_filename)
        print *, "Opening geometry file: ", trim(geo_filename)
        print *, "Opening fix file: ", trim(fix_filename)

        ! Initialize mesh values
        mesh % npoin = 0
        mesh % nelem = 0
        mesh % nboun = 0
        mesh % ndim = 0

        ! Read domain metadata
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

        ! Allocate memory for mesh components
        if (.not. allocated(mesh % coords)) then
            allocate(mesh % coords(mesh % npoin, mesh % ndim))
            print*, "Allocated mesh coordinates: ", mesh % npoin, " x ", mesh % ndim
        endif

        if (.not. allocated(mesh % connectivity)) then
            allocate(mesh % connectivity(mesh % nelem, 3)) ! Assuming triangular elements
            print*, "Allocated mesh connectivity: ", mesh % nelem, " x ", 3
        endif

        if (.not. allocated(mesh % boundary)) then
            allocate(mesh % boundary(mesh % nboun))
            print*, "Allocated mesh boundary: ", mesh % nboun
        endif

        call read_geo_dat(geo_filename, mesh)
        call read_fix_dat(fix_filename, mesh)
    end subroutine read_domain

    !>-----------------------------------------------------------------------
    !> Reads nodal coordinates and element connectivity.
    !>-----------------------------------------------------------------------
    subroutine read_geo_dat(geo_filename, mesh)
        character(len=*), intent(in)   :: geo_filename
        type(mesh_type), intent(inout) :: mesh
        character(len=256) :: line
        integer(ip) :: ios, node_id, elem_id, n1, n2, n3
        logical :: inside_coords, inside_elems
        real(rp) :: x, y, z

        inside_coords = .false.
        inside_elems = .false.

        open(unit=20, file=trim(geo_filename), status="old", action="read", iostat=ios)
        if (ios /= 0) then
            print*, "Error opening ", trim(geo_filename)
            stop
        end if

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
                    read(line, *) elem_id, n1, n2, n3
                    mesh % connectivity(elem_id, :) = (/ n1, n2, n3 /)
                end if
            end select
        end do

        close(20)
    end subroutine read_geo_dat

end module mod_mesh
