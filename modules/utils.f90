module utils
    use constants
    implicit none
    public :: norm

contains

    ! Compute the Euclidean norm of a vector of size n.
    function norm(vec, n) result(res)
        real(rp), intent(in) :: vec(:)
        integer, intent(in) :: n
        real(rp) :: res
        res = sqrt(sum(vec(1:n)**2))
    end function norm

end module utils
