subroutine centroidBdrySegment(x_centroid, x_v, dim)
    implicit none
    integer(4), intent(in) :: dim
    real(8), intent(in) :: x_v(dim, dim)
    real(8), intent(out) :: x_centroid(dim)
    integer(4) :: d

    x_centroid = 0.0d0
    do d = 1, dim
        x_centroid = x_centroid + x_v(:, d)
    end do
    x_centroid = x_centroid / dble(dim)

end subroutine