module interpolation_module
    implicit none
contains
    ! Subroutine to perform 2D linear interpolation
    subroutine linear_interp_2d(x_input, f_input, x, f_output, n)
        real(8), intent(in) :: x_input(2, :)      ! Coordinates of known points
        real(8), intent(in) :: f_input(:)          ! Function values at known points
        real(8), intent(in) :: x(2)                ! Coordinates of the interpolation point
        real(8), intent(out) :: f_output           ! Interpolated function value
        integer, intent(in) :: n                   ! Number of known points

        integer :: i, nearest_idx
        real(8) :: min_dist, dist

        ! Initialize minimum distance to a large value
        min_dist = 1.0d20
        nearest_idx = 1

        ! Find the nearest neighbor
        do i = 1, n
            dist = sqrt((x(1) - x_input(1, i))**2 + (x(2) - x_input(2, i))**2)
            if (dist < min_dist) then
                min_dist = dist
                nearest_idx = i
            end if
        end do

        ! Use the function value at the nearest neighbor as the interpolated value
        f_output = f_input(nearest_idx)
    end subroutine linear_interp_2d
    
end module interpolation_module
