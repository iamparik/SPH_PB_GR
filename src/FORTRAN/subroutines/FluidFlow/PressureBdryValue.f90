subroutine PressureBdryValue(prsr_,rho_, x_, vx_, itype_,bdryVal_seg_, num_seg, s, ID)
    use particle_data ,only: mid_pt_for_edge, surf_norm
    use config_parameter, only : dx_r, itype_virtual,c_sound, SPH_dim
    
implicit none
  
integer(4), intent(in) :: ID, s, num_seg
integer(2), intent(in) :: itype_
real(8), intent(in) :: bdryVal_seg_(num_seg), x_(SPH_dim), vx_(SPH_dim)
real(8), intent(inout) :: prsr_,rho_
real(8) :: rho_comp 

if(ID .eq. 1) then
    call ParticlePressureEOS(prsr_, rho_, itype_, itype_virtual)    
     prsr_ = prsr_ - rho_*c_sound*dot_product(vx_-bdryVal_seg_(1:SPH_dim), surf_norm(:,s))
     
elseif(ID .eq. 2) then
    call ParticlePressureEOS(prsr_, rho_comp, itype_, itype_virtual)    
    prsr_ = prsr_ + rho_*norm2(vx_-bdryVal_seg_(1:SPH_dim))**2
    
elseif(ID .eq. 3) then    
    rho_comp=max(rho_,(rho_/(dx_r/4.D0))*(2.D0*(dx_r/4.D0)-norm2(x_-mid_pt_for_edge(:,s))))    
    call ParticlePressureEOS(prsr_, rho_comp, itype_, itype_virtual)    
    prsr_ = prsr_ - rho_*c_sound*dot_product(vx_-bdryVal_seg_(1:SPH_dim), surf_norm(:,s)) 
    
elseif(ID .eq. 4) then    
    rho_comp=max(rho_,(rho_/(dx_r/2.D0))*(2.D0*(dx_r/2.D0)-norm2(x_-mid_pt_for_edge(:,s))))    
    call ParticlePressureEOS(prsr_, rho_comp, itype_, itype_virtual)    
    prsr_ = prsr_ - rho_*c_sound*dot_product(vx_-bdryVal_seg_(1:SPH_dim), surf_norm(:,s)) 
    
elseif(ID .eq. 5) then
    !rho_wall_compress=max(rho(a),rho(a)*(gamma_cont(a) + 0.5D0 - 2.D0*gamma_wall_cutoff)/(0.5D0-gamma_wall_cutoff))
    rho_comp=max(rho_,(rho_/(dx_r/4.D0))*(2.D0*(dx_r/4.D0)-norm2(x_-mid_pt_for_edge(:,s))))
    
    call ParticlePressureEOS(prsr_, rho_comp, itype_, itype_virtual)    
    prsr_ = prsr_ - rho_*c_sound*dot_product(vx_-bdryVal_seg_(1:SPH_dim), surf_norm(:,s))  &
                            !& - rho(a)*dot_product(F_ext, x(:,a)- mid_pt_for_edge(:,s)) &
                            !& + rho(a)*c_sound*norm2(vx(:,a)-bdryVal_seg(1:SPH_dim,s))*max(0.6D0-gamma_cont(a), 0.D0) &
                            & + 0.D0
else 
    call ParticlePressureEOS(prsr_, rho_, itype_, itype_virtual) 
    
endif

end