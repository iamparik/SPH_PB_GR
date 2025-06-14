!****************************************************************************
!
!  SUBROUTINE: IICinput
!
!  PURPOSE:  Subroutine that allows imposing different types of   
!            Initial conditions for different particle types (itypes)
!
!****************************************************************************
subroutine ICinput(itype, tempType)
    implicit none
    integer(2), intent(out) :: itype
    integer(4), intent(in) :: tempType

    if(tempType .eq. 1) itype = 1
        
    if(tempType .eq. 2) itype = 2
    
    itype = tempType

end subroutine

    
    
subroutine ICinputValue(initalVal_particle,num_var, itype)
    use config_parameter , only : rho_init, c_sound, hsml_const, mu_const
    implicit none
    integer(2), intent(in) :: itype
    integer(4), intent(in) ::  num_var
    real(8), dimension(num_var), intent(inout) :: initalVal_particle
    integer(4) :: d
    
    
    if(itype .eq. 3) then
        
        initalVal_particle(3)= 0.D0 ! vx
        initalVal_particle(4)= 0.D0 ! vy
        initalVal_particle(5)= rho_init*((7.D0*initalVal_particle(6)/(rho_init*c_sound**2)+1.D0)**(1.D0/7.D0))
        !initalVal_particle(6)= (rho_init*c_sound**2/7.D0)*((initalVal_particle(5)/rho_init)**7.D0-1.D0) !pressure
        initalVal_particle(7)= hsml_const ! smoothing length
        initalVal_particle(8)= mu_const !viscosity mu
        
    endif
    
    
    end subroutine
    
    
    subroutine external_file_IC
    use config_parameter , only : DataConfigPath
    use particle_data, only: p
    integer(4) :: k
    real(8) :: tempx,tempy
    
        open(1,file= DataConfigPath // '/input_particles_param_value.dat',status='old')
        k=0
        do while (.not.eof(1))
            k=k+1
            read(1,*) tempx,tempy, p(k)   
        enddo

    end subroutine
    
