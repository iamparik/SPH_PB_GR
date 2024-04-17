subroutine createOutputFolder(directory_name)
    implicit none
    
    character(len= *), intent (in) :: directory_name
    character(len=100) :: command
    integer :: stat_result
    integer :: status

    ! Check if the directory exists
    inquire(DIRECTORY=directory_name, exist=stat_result)
    
    write(*,*) ' value of stat_result is :' , stat_result
    if (stat_result == 0) then
        ! Directory does not exist, create it
        command = 'mkdir ' // trim(directory_name)
        call execute_command_line(command, wait=.true., exitstat=status)
    
        if (status == 0) then
            print *, 'Directory "', trim(directory_name), '" created successfully.'
        else
            print *, 'Error creating directory "', trim(directory_name), '".'
        endif
    else
        ! Directory already exists
        print *, 'Directory "', trim(directory_name), '" already exists.'
    endif
    
end subroutine