subroutine ParticlePressureEOS(p, rho, itype, itype_virtual) 
    !he pressure value is calculated from e.o.s, which relates pressure and density
    use config_parameter, only: rho_init,c_sound, P_EOS_extra
    implicit none
    
    
    integer(2) , intent(in):: itype, itype_virtual
    real(8), intent(in):: rho
    real(8), intent(inout):: p
    real(8) :: ecf
    
    
    if( mod(itype,itype_virtual) .eq. 2) then
        p=rho*c_sound**2
        !rho= p/c_sound**2
        
    elseif ( mod(itype,itype_virtual) .eq. 3) then
        ecf=7.D0
        p=(rho_init*c_sound**2/ecf)*((rho/rho_init)**ecf-1.D0) + P_EOS_extra
        !rho=rho_init*(p/(rho_init*c_sound**2/ecf)+1.D0)**(1.D0/ecf)
        ! there might be subtractive cancelling above
        
    elseif ( mod(itype,itype_virtual) .eq. 4) then
        p=0.D0
            
    endif
    
    
    end
   