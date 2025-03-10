module mod_post
    use constants
    use mod_mesh
    implicit none
    public :: write_output


contains

! > Write the solution to a binary file (mpio.bin) for later conversion.
    subroutine write_output(mesh, sol, step, rank)
        type(mesh_type), intent(in) :: mesh
        real(rp), intent(in) :: sol(:)
        integer, intent(in) :: step, rank
        character(len=100) :: filename
        integer :: ios, n

        ! Create a filename that includes the time step and MPI rank.
        write(filename, '(A,I0,A,I0,A)') 'sol_step_', step, '_rank_', rank, '.mpio.bin'
        open(unit=10, file=filename, form='unformatted', status='replace', iostat=ios)
        if (ios /= 0) then
            print*, 'Error opening file ', trim(filename)
        else
            n = size(state)
            write(10) state
            close(10)
            print*, 'Output written to ', trim(filename)
        end if
    end subroutine write_output

end module post_module
