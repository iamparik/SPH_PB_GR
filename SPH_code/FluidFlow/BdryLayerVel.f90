subroutine BdryLayerVel(SPH_dim,bdryValWallVel)
    use config_geometry, only: dx_r
    use particle_data, only: epair_a,epair_s, eniac, &
        & surf_norm,nedge_rel_edge, &
        & mu, x, vx
    
    implicit none
    
    integer(4), intent(in) :: SPH_dim
	real(8):: bdryValWallVel(SPH_dim,eniac) 
    integer(4) :: k, a, s
    real(8) ras(SPH_dim), zas, vas_t(SPH_dim)
    
       
    do k =1, eniac
    
        a = epair_a(k)
        s = epair_s(k)
        
        ras(:)=x(:,a) - x(:, nedge_rel_edge(s))
        
        zas = max(dot_product(ras,surf_norm(:,s)), dx_r)
        
        vas_t(:) = vx(:,a) - dot_product(vx(:,a),surf_norm(:,s)) * surf_norm(:,s)
        
        bdryValWallVel(:,k) = bdryValWallVel(:,k) + mu(a)*(vas_t(:)/zas)
    
    
    enddo
    
end