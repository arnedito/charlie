module mod_mesh
    use iso_fortran_env
    use constants, only: ip, rp
    implicit none
    public :: mesh_type, read_domain, local_elements

    !>----------------------------------------------------------------------
    !> Unified Mesh Type
    !>----------------------------------------------------------------------
    type :: mesh_type
        integer(ip)                 :: ndim = 0
        integer(ip)                 :: npoin = 0
        integer(ip)                 :: nelem = 0
        integer(ip)                 :: nboun = 0
        real(rp), allocatable       :: coords(:,:)
        integer(ip), allocatable    :: connectivity(:,:)
        integer(ip), allocatable    :: boundary(:)
    end type mesh_type

contains

    !>-----------------------------------------------------------------------
    !> read_domain: Reads the mesh using user-provided file paths
    !>-----------------------------------------------------------------------
    subroutine read_domain(dom_filename, geo_filename, fix_filename, mesh)
        character(len=*), intent(in) :: dom_filename, geo_filename, fix_filename
        type(mesh_type), intent(out) :: mesh
        character(len=256) :: line, dummy
        integer(ip) :: ios

        ! Debugging: Print the file paths
        print *, "Trying to open domain file: ", trim(dom_filename)
        print *, "Trying to open geometry file: ", trim(geo_filename)
        print *, "Trying to open fix file: ", trim(fix_filename)

        ! Initialize mesh values
        mesh % npoin = 0
        mesh % nelem = 0
        mesh % nboun = 0
        mesh % ndim = 0

        ! Open domain file
        open(unit=10, file=trim(dom_filename), status="old", action="read", iostat=ios)
        if (ios /= 0) then
            print*, "❌ Error opening ", trim(dom_filename)
            stop
        end if

        ! Read domain metadata
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

        ! Read geometry and fix files
        call read_geo_dat(geo_filename, mesh)
        call read_fix_dat(fix_filename, mesh)
    end subroutine read_domain

    !>-----------------------------------------------------------------------
    !> read_geo_dat: Reads node coordinates and element connectivity.
    !>-----------------------------------------------------------------------
    subroutine read_geo_dat(geo_filename, mesh)
        character(len=*), intent(in) :: geo_filename
        type(mesh_type), intent(inout) :: mesh
        character(len=256) :: line
        integer(ip) :: ios, node_id, elem_id, n1, n2, n3
        logical :: inside_coords, inside_elems
        real(rp) :: x, y, z

        inside_coords = .false.
        inside_elems = .false.

        open(unit=20, file=trim(geo_filename), status="old", action="read", iostat=ios)
        if (ios /= 0) then
            print*, "❌ Error opening ", trim(geo_filename)
            stop
        end if

        allocate(mesh % coords(mesh % npoin, mesh % ndim))
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

    !>-----------------------------------------------------------------------
    !> read_fix_dat: Reads boundary conditions from fix file.
    !>-----------------------------------------------------------------------
    subroutine read_fix_dat(fix_filename, mesh)
        character(len=*), intent(in) :: fix_filename
        type(mesh_type), intent(inout) :: mesh
        character(len=256) :: line
        integer(ip) :: ios, i, id, code

        open(unit=30, file=trim(fix_filename), status="old", action="read", iostat=ios)
        if (ios /= 0) then
            print*, "❌ Error opening ", trim(fix_filename)
            stop
        end if

        allocate(mesh % boundary(mesh % nboun))

        ! Skip header
        read(30, '(A)', iostat=ios) line

        do i = 1, mesh % nboun
            read(30, *) id, code
            mesh % boundary(i) = code
        end do

        close(30)
    end subroutine read_fix_dat

    !----------------------------------------------------------------------
    ! Function: local_elements
    ! Purpose : Returns the number of local elements by computing the size
    !           of the first dimension of the connectivity array.
    !----------------------------------------------------------------------
    function local_elements(m) result(n)
        type(mesh_type), intent(in) :: m
        integer :: n
        if (allocated(m % connectivity)) then
            n = size(m % connectivity, 1)
        else
            n = 0
        end if
    end function local_elements

end module mod_mesh
