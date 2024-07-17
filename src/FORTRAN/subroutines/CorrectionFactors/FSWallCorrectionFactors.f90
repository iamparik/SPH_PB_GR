subroutine FSWallCorrectionFactors
    use particle_data, only: nreal, gamma_discrt, del_gamma, gamma_cont, xi1_mat, beta_mat, gamma_mat, &
                            & gamma_mat_FSWall1_inv, gamma_mat_FSWall2_inv
    use config_paramter, only: SPH_dim
    implicit none
    integer(4) :: a,d
    
    
    allocate(gamma_mat_FSWall1_inv(SPH_dim,SPH_dim,ntotal), gamma_mat_FSWall2_inv(SPH_dim,SPH_dim,ntotal))
    gamma_mat_FSWall1_inv = 0.D0
    gamma_mat_FSWall2_inv = 0.D0
    
    do a = 1,nreal
        if( (free_surf_particle(a) .eq. 1) .and. (gamma_cont(a) .lt. 1.D0)) then
            gamma_mat_FSWall1_inv=0.D0
            gamma_mat_FSWall2_inv=0.D0
            do d =1,SPH_dim
                gamma_mat_FSWall1_inv(d,d,a) = 1.D0/gamma_cont(a)
            enddo
            gamma_mat_FSWall2_inv(:,:,a) = xi1_mat(:,:,a)
        endif
         
    enddo    
    
end
    