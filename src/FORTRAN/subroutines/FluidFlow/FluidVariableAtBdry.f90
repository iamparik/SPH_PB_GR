
    subroutine FluidVariableAtBdry(num_bdry_var_seg,ID_prsrBdryType)
    use particle_data, only: etotal, rho_s,prsr_bdry_val, epair_a, epair_s, eniac, &
            & hsml, mid_pt_for_edge, vol, rho,x, vx, itype,bdryVal_seg, mass
    use config_parameter, only: SPH_dim
    implicit none
    integer(4), intent(in) :: num_bdry_var_seg, ID_prsrBdryType
    real(8), DIMENSION(:), allocatable :: gamma_discrt_s
    integer(4) :: a,s,k,d
    real(8) ::temp_scalar,temp_vector(SPH_dim), driac, w_temp, dxiac(SPH_dim),Sca_Bdry_val

    
    allocate( prsr_bdry_val(etotal),gamma_discrt_s(etotal) ,rho_s(etotal))
        prsr_bdry_val=0.D0
        gamma_discrt_s=0.D0
        rho_s=0.D0
            
            
        do k= 1,eniac
            a= epair_a(k)
            s= epair_s(k)
                
            driac=0.D0
            do d=1,SPH_dim
                dxiac(d) = x(d,a) - mid_pt_for_edge(d,s)
                driac    = driac + dxiac(d)*dxiac(d)
            enddo                
            driac=sqrt(driac)
            call kernel(driac,dxiac,hsml(a),w_temp,temp_vector)                
            rho_s(s) = rho_s(s)+ mass(a)*w_temp
                
            gamma_discrt_s(s) = gamma_discrt_s(s) + vol(a)*w_temp
        enddo  
            
        do s=1,etotal
            if( rho_s(s) .lt. 1.D-10) then
                rho_s(s) = 0.D0
            else
                rho_s(s)=rho_s(s)/gamma_discrt_s(s)
            endif
        enddo
            
        do k= 1,eniac
            a= epair_a(k)
            s= epair_s(k)
                
            driac=0.D0
            do d=1,SPH_dim
                dxiac(d) = x(d,a) - mid_pt_for_edge(d,s)
                driac    = driac + dxiac(d)*dxiac(d)
            enddo                
            driac=sqrt(driac)
            call kernel(driac,dxiac,hsml(a),w_temp,temp_vector)  
                
            call PressureBdryValue(Sca_Bdry_val,rho(a),x(:,a), vx(:,a), itype(a),bdryVal_seg(:,s), num_bdry_var_seg, a,s, ID_prsrBdryType)
            prsr_bdry_val(s)=prsr_bdry_val(s) + Sca_Bdry_val*vol(a)*w_temp
        enddo    
            
        do s=1,etotal
            if( prsr_bdry_val(s) .lt. 1.D-10) then
                prsr_bdry_val(s) = 0.D0
            else
                prsr_bdry_val(s)=prsr_bdry_val(s)/gamma_discrt_s(s)
            endif
        enddo
            
        deallocate(gamma_discrt_s)
    
end subroutine