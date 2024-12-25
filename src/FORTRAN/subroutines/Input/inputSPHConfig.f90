subroutine inputSPHConfig
    use config_parameter
    implicit none
    
    integer(4) file_type, run_packed_simulation
    
    ! Read parameters to configure the SPH simulation from config_data.dat
    call read_config_file

    
    !A series of subroutines are now run to initialize the particles, and obtain input file
    write(*,*)'  ***************************************************'
    write(*,*) 'Are you going to input external input file or back up file?'
    write(*,*) '(0=No Input file, 1=new external input file, 2=continue from backup file)'
    write(*,*)'  ***************************************************'
    read (*,*) file_type
    
    ! All particles are created, labeled, and retrieved to start calculations
    if(file_type.eq.0) then    
        !call input
    elseif(file_type.eq.1) then
        ! run packing algorithm according to user input
        call inputExt
        
        if(packagingIterations) then
            write(*,*) "BIPI algorithm implemented, press enter to now use BIPI config as input"
            !set external CAD input to type 5 
           ! ExtInputMeshType = 5
            packagingIterations = .false.
            
            call inputExt
            
        endif
        
    elseif(file_type.eq.2) then
        !call backupInput(itimestep)
    endif
    
end subroutine