subroutine printTimeStep(itimestep,print_step)
    implicit none
    integer(4), intent(in) :: itimestep,print_step
!   [itimestep mod print step] like 6mod4=2 is measured to check if
!   itimestep is a multiple of printstep, if yes then value of modulo
!   will be zero and the below if loop is executed    
    if (mod(itimestep,print_step).eq.0) then 
        write(*,*)'______________________________________________'
        write(*,*)'  current    number     of     time     step =',     itimestep
        write(*,*)'______________________________________________'
    endif  
    
end subroutine