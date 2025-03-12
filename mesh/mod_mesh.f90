!=============================================================================
! @file   mod_mesh.f90
! @brief  Manages the mesh data structure and associated I/O routines.
!
! @details
!   This module defines the mesh_type, which holds information about nodes,
!   elements, and boundary data. It also provides routines to read mesh data
!   from various input files. The code uses integer(ip) and real(rp) kinds for
!   consistency and potential future precision changes.
!
!=============================================================================
module mod_mesh
   use iso_fortran_env
   use constants, only: ip, rp
   implicit none

   !---------------------------------------------------------------------------
   !> @brief Public type that represents a computational mesh.
   !>
   !! - ndim  : Number of spatial dimensions (2D or 3D).
   !! - npoin : Number of nodal points in the domain.
   !! - nelem : Number of elements in the domain.
   !! - nboun : Number of boundary entities (e.g. edges or faces).
   !! - coords: Allocatable array holding nodal coordinates of size (npoin, ndim).
   !! - connectivity: Allocatable array of element connectivity of size (nelem, ?).
   !! - boundary    : Allocatable array describing boundary conditions or boundary
   !!                 element references.
   !---------------------------------------------------------------------------
   type :: mesh_type
      integer(ip)                 :: ndim                 ! Number of spatial dimensions
      integer(ip)                 :: npoin                ! Number of nodal points
      integer(ip)                 :: nelem                ! Number of elements
      integer(ip)                 :: nboun                ! Number of boundary elements
      real(rp), allocatable       :: coords(:,:)          ! Nodal coordinates
      integer(ip), allocatable    :: connectivity(:,:)    ! Element connectivity
      integer(ip), allocatable    :: boundary(:)          ! Boundary conditions
   end type mesh_type

   public :: mesh_type, read_domain, local_elements   ! local_elements presumably defined elsewhere?

contains

   !---------------------------------------------------------------------------
   !> @brief Reads mesh data from user-provided input files.
   !>
   !! This routine opens and reads domain metadata (number of points, elements,
   !! boundaries, etc.) from dom_filename. It then allocates the necessary arrays
   !! in the mesh_type structure before delegating actual geometry and boundary
   !! data reading to read_geo_dat() and read_fix_dat() (not shown here).
   !!
   !! @param[in]  dom_filename  (character(*)) File with domain metadata (NODAL_POINTS, ELEMENTS, etc.).
   !! @param[in]  geo_filename  (character(*)) File containing actual node and element data.
   !! @param[in]  fix_filename  (character(*)) File with boundary or "fix" data.
   !! @param[out] mesh          (mesh_type)    The mesh structure to be populated.
   !
   !! @note
   !!  - The subroutine stops if it fails to open a file (see improvements).
   !!  - By default, triangular(2d) and tetrahedral(3d) connectivity is assumed (3/4 nodes/element).
   !---------------------------------------------------------------------------
   subroutine read_domain(dom_filename, geo_filename, fix_filename, mesh)
      character(len=*), intent(in)            :: dom_filename, geo_filename, fix_filename
      type(mesh_type), intent(out)            :: mesh
      character(len=256)                      :: line, dummy
      integer(ip)                             :: ios

      print *, "Opening domain file: ", trim(dom_filename)
      print *, "Opening geometry file: ", trim(geo_filename)
      print *, "Opening fix file: ", trim(fix_filename)

      ! Initialize mesh values
      mesh%npoin = 0_ip
      mesh%nelem = 0_ip
      mesh%nboun = 0_ip
      mesh%ndim  = 0_ip

      !-----------------------------------------------------------------------
      ! Possible Improvement:
      !   - Use a robust error-handling mechanism (e.g., check return codes,
      !     handle file missing scenarios gracefully, etc.).
      !   - You could use an error flag or condition instead of `stop`.
      !-----------------------------------------------------------------------
      open(unit=10, file=trim(dom_filename), status="old", action="read", iostat=ios)
      if (ios /= 0_ip) then
         print*, "Error opening ", trim(dom_filename)
         stop
      end if

      do
         read(10, '(A)', iostat=ios) line
         if (ios /= 0_ip) exit
         if (index(line, "NODAL_POINTS") > 0) then
            read(line, *) dummy, mesh%npoin
         else if (index(line, "ELEMENTS") > 0) then
            read(line, *) dummy, mesh%nelem
         else if (index(line, "BOUNDARIES") > 0) then
            read(line, *) dummy, mesh%nboun
         else if (index(line, "SPACE_DIMENSIONS") > 0) then
            read(line, *) dummy, mesh%ndim
         else if (index(line, "END_DIMENSIONS") > 0) then
            exit
         end if
      end do

      close(10)

      !-----------------------------------------------------------------------
      ! Allocate memory for mesh components
      !-----------------------------------------------------------------------
      if (.not. allocated(mesh%coords)) then
         allocate(mesh%coords(mesh%npoin, mesh%ndim))
         print*, "Allocated mesh coordinates: ", mesh%npoin, " x ", mesh%ndim
      end if

      !-----------------------------------------------------------------------
      ! Possible Improvement:
      !   - Instead of forcing triangular or tetra connectivity (3/4 nodes per element),
      !     parse the element type from the input file or domain file to
      !     allocate connectivity appropriately (e.g., 4 for quads, 8 for hex,
      !     etc.).
      !-----------------------------------------------------------------------
      if (.not. allocated(mesh%connectivity)) then
         allocate(mesh%connectivity(mesh%nelem, mesh%ndim + 1))
         print*, "Allocated mesh connectivity: ", mesh%nelem, " x ", mesh%ndim + 1
      end if

      if (.not. allocated(mesh%boundary)) then
         allocate(mesh%boundary(mesh%nboun))
         print*, "Allocated mesh boundary: ", mesh%nboun
      end if

      !-----------------------------------------------------------------------
      ! Reading geometry and boundary data from files
      !-----------------------------------------------------------------------
      call read_geo_dat(geo_filename, mesh)
      call read_fix_dat(fix_filename, mesh)  ! <--- Assuming you have a read_fix_dat subroutine
    end subroutine read_domain

    !---------------------------------------------------------------------------
    !> @brief Reads nodal coordinates and element connectivity from geo_filename.
    !>
    !! The subroutine looks for markers COORDINATES and ELEMENTS in the file. It
    !! then populates the mesh%coords array with node coordinates, and
    !! mesh%connectivity with element connectivity data.
    !!
    !! @param[in]    geo_filename (character(*)) Input file with node coords and elements.
    !! @param[inout] mesh         (mesh_type)    The mesh structure (updated in-place).
    !
    !! @note
    !!  - The read logic toggles between reading coordinates and elements based
    !!    on specific keywords in the file (COORDINATES, ELEMENTS, etc.).
    !!  - The dimensionality (mesh%ndim) determines whether to read (x,y) or
    !!    (x,y,z).
    !---------------------------------------------------------------------------
    subroutine read_geo_dat(geo_filename, mesh)
       character(len=*), intent(in)     :: geo_filename
       type(mesh_type), intent(inout)   :: mesh
       character(len=256)               :: line
       integer(ip)                      :: ios, node_id, elem_id, n1, n2, n3
       logical                          :: inside_coords, inside_elems
       real(rp)                         :: x, y, z

       inside_coords = .false.
       inside_elems  = .false.

       open(unit=20, file=trim(geo_filename), status="old", action="read", iostat=ios)
       if (ios /= 0_ip) then
          print*, "Error opening ", trim(geo_filename)
          stop
       end if

       do
          read(20, '(A)', iostat=ios) line
          if (ios /= 0_ip) exit

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
                if (mesh%ndim == 2_ip) then
                   read(line, *) node_id, x, y
                   mesh%coords(node_id, 1) = x
                   mesh%coords(node_id, 2) = y
                else if (mesh%ndim == 3_ip) then
                   read(line, *) node_id, x, y, z
                   mesh%coords(node_id, 1) = x
                   mesh%coords(node_id, 2) = y
                   mesh%coords(node_id, 3) = z
                end if
             else if (inside_elems) then
                read(line, *) elem_id, n1, n2, n3
                mesh%connectivity(elem_id, :) = (/ n1, n2, n3 /)
             end if
          end select
       end do

       close(20)
    end subroutine read_geo_dat

    !---------------------------------------------------------------------------
    ! Possible additional subroutines, e.g., read_fix_dat, local_elements, etc.
    !---------------------------------------------------------------------------

end module mod_mesh
