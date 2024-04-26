!****************************************************************************
!
!  SUBROUTINE: pressureEOS
!                                                
!  PURPOSE:The pressure value is calculated from e.o.s
!           
!   CREATED:        08/22/2021       by  PARIKSHIT BOREGOWDA
!   Last Modified:  02/06/2021       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine pressureEOS

    use config_parameter, only: itype_virtual, &
        & rho_init, c_sound
    use particle_data, only: rho, p, ntotal, itype
    
    implicit none
    
    integer(4):: a
    real(8):: ecf
    p=0.D0
    do a=1,ntotal
        
        if( mod(itype(a),itype_virtual) .eq. 2) then
            ecf=7.D0
            !p(a)=(rho_init*c_sound**2/ecf)*((rho(a)/rho_init)**ecf-1.D0) 
            p(a)=rho(a)*c_sound**2
            ! there might be subtractive cancelling above
        endif
        
        if( mod(itype(a),itype_virtual) .eq. 3) then
            ecf=7.D0
            p(a)=(rho_init*c_sound**2/ecf)*((rho(a)/rho_init)**ecf-1.D0) 
            ! there might be subtractive cancelling above
        endif
        
        if( mod(itype(a),itype_virtual) .eq. 1) then
            ! this can be anything, as needed to be output for ghost points.
            !p(a)=(rho_init*c_sound**2/ecf)*((rho(a)/rho_init)**ecf-1.D0)  
            p(a)=0.D0
            
        endif
        
    enddo
    
    
end
    