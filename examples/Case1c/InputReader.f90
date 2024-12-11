module InputReader
    implicit none
    real(8),DIMENSION(:,:),ALLOCATABLE :: x_input
    real(8),DIMENSION(:,:),ALLOCATABLE :: vx_input
    real(8),DIMENSION(:),ALLOCATABLE:: p_input
contains
    subroutine read_csv_with_header(file_path, x, vx, p, num_lines)
        implicit none
        character(len=*), intent(in) :: file_path
        integer, intent(out) :: num_lines
        real(8), allocatable, intent(out) :: x(:,:), vx(:,:), p(:)

        ! Temporary variables
        real(8) :: not_needed_value1, not_needed_value2, not_needed_value3
        character(len=100) :: header_line
        integer :: i, n
        integer :: unit_id = 10
        integer :: io_status

        ! Open the file and initialize line counter
        open(unit=unit_id, file=file_path, status='old', action='read')
        num_lines = 0

        ! Skip the header line
        read(unit_id, '(A)') header_line

        ! Count the number of lines to determine the array size
        do
            read(unit_id, *, iostat=io_status)
            if (io_status /= 0) exit
            num_lines = num_lines + 1
        end do

        ! Allocate arrays based on the number of data lines
        allocate(x(2, num_lines), vx(2, num_lines), p(num_lines))

        ! Rewind and skip the header line again
        rewind(unit_id)
        read(unit_id, '(A)') header_line

        ! Read the data into the arrays
        do n = 1, num_lines
            read(unit_id, *) p(n), vx(1,n), vx(2,n), not_needed_value1, not_needed_value2, &
                             x(1,n), x(2,n), not_needed_value3
        end do

        ! Close the file
        close(unit_id)
    end subroutine read_csv_with_header
end module csv_reader
