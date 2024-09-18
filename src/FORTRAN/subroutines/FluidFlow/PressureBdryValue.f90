subroutine PressureBdryValue(prsr_s,rho_, x_, vx_, itype_,bdryVal_seg_, num_seg,a, s, ID)
    use particle_data ,only: mid_pt_for_edge, surf_norm, rho_s, p, gamma_dens_cut_off, gamma_cont
    use config_parameter, only : dx_r, itype_virtual,c_sound, SPH_dim, F_ext, hsml_const, pi
    
implicit none
  
integer(4), intent(in) :: ID, s, num_seg,a
integer(2), intent(in) :: itype_
real(8), intent(in) :: bdryVal_seg_(num_seg), x_(SPH_dim), vx_(SPH_dim)
real(8), intent(inout) :: prsr_s,rho_
real(8) :: rho_comp 
real(8) ::xEdgeTemp(SPH_dim,SPH_dim), xrefPoint(SPH_dim), xEdge_surfNorm(SPH_dim),scale_k

if(ID .eq. 1) then
     prsr_s = p(a) - rho_*c_sound*dot_product(vx_-bdryVal_seg_(1:SPH_dim), surf_norm(:,s))
     
elseif(ID .eq. 2) then    
    prsr_s = p(a) + 0.5D0*rho_*(norm2(vx_)**2 - norm2(bdryVal_seg_(1:SPH_dim))**2) - rho_*dot_product(F_ext, x_-mid_pt_for_edge(:,s))
    
elseif(ID .eq. 3) then    
    if( .not. allocated(gamma_dens_cut_off) ) then
        gamma_dens_cut_off=0.D0
        !This needs to happen only once per particle close to the boundary
        call sml_mult_factor(scale_k)
        xEdgeTemp(:,:)= 0.D0
        xrefPoint(:)=0.D0
        xEdge_surfNorm(:)=0.D0
        xEdgeTemp(1,1)= -scale_k*hsml_const 
        xEdgeTemp(1,2)= scale_k*hsml_const 
        xrefPoint(2) = -dx_r/2.D0
        xEdge_surfNorm(2)= -1.D0
        
        call gamma_analytical_leroy_2D(xrefPoint,xEdgeTemp,xEdge_surfNorm,gamma_dens_cut_off,hsml_const,pi)
        gamma_dens_cut_off=1-gamma_dens_cut_off
    endif
    rho_comp=max(rho_,(rho_/(0.5D0-gamma_dens_cut_off))*(gamma_cont(a) - 2.D0*gamma_dens_cut_off + 0.5D0))    
    call ParticlePressureEOS(prsr_s, rho_comp, itype_, itype_virtual)    
    prsr_s = prsr_s - rho_*c_sound*dot_product(vx_-bdryVal_seg_(1:SPH_dim), surf_norm(:,s)) 
    
elseif(ID .eq. 4) then    
    rho_comp=max(rho_,(rho_/(dx_r/2.D0))*(2.D0*(dx_r/2.D0)-norm2(x_-mid_pt_for_edge(:,s))))    
    call ParticlePressureEOS(prsr_s, rho_comp, itype_, itype_virtual)    
    prsr_s = prsr_s - rho_*c_sound*dot_product(vx_-bdryVal_seg_(1:SPH_dim), surf_norm(:,s)) 
    
elseif(ID .eq. 5) then
    prsr_s = p(a) + rho_*0.5D0*norm2(vx_-bdryVal_seg_(1:SPH_dim))**2- rho_*dot_product(F_ext, x_-mid_pt_for_edge(:,s))

elseif(ID .eq. 6) then
    call ParticlePressureEOS(prsr_s, rho_, itype_, itype_virtual)    
    prsr_s = (prsr_s/rho_ + 0.5D0*(norm2(vx_)**2 - norm2(bdryVal_seg_(1:SPH_dim))**2) - dot_product(F_ext, x_-mid_pt_for_edge(:,s)))*rho_s(s)
else
    prsr_s = p(a) 
endif

end subroutine
    
    
    