!****************************************************************************
!
!  SUBROUTINE: thermalBoundaryConditionFerrandInspired
!
!  PURPOSE:  Subroutine to apply thermal boudnary condition for a   
!         for thermal diffusion problem inspired from FErrands homogeneuous neumann bc 
!
!   CREATED:        09/14/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  09/20/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine thermalBoundaryConditionFerrandInspired

    use config_parameter,   only: SPH_dim, etype_thermal_dirichlet, etype_thermal_neumann, dx_r
    use particle_data,      only: ntotal,etotal,etype,&
        & epair_a, epair_s, eniac, temp, hsml, mass, &
        & rho, nedge_rel_edge, x,surf_norm !, bdryVal_temp
    
    implicit none
      
    real(8) :: tdwdx(SPH_dim), dxiac(SPH_dim), driac, w ,fn, mv_in, scale_k , bdryVal_temp(etotal)    
    real(8),DIMENSION(:),ALLOCATABLE :: alpha, T   
    integer(4) :: a,b,s,k,d      

    ! assuming each edge has one edge point, we initialize alpha and T
    ! for more than one edge point per edge initialize(eppe) as 
    ! alpha(etotal, eppe), T(etotal, eppe)
    allocate(alpha(etotal), T(etotal))
    alpha=0.D0
    T=0.D0
    
    mv_in=dx_r*2.D0 !dx_r*1/20.D0 !dx_r*2.D0 ! 0.D0

    !call sml_mult_factor(scale_k)

    do k=1,eniac
        b=epair_a(k)
        s=epair_s(k)
    
        a = nedge_rel_edge(s)
        
        ! Apply neumann boundary condition for wall with dT/dn=0
        ! Here calcualte T_si=∑_b(T_b*W_sb*m_b/ρ_b)  & 
        ! α_s=∑_b(W_sb*m_b/ρ_b)
        if ( etype(s) .eq. etype_thermal_neumann ) then
        
            driac = 0.D0 
            
            do d=1,SPH_dim
                    !since surf_norm is inward we displace inwards by adding to x(d,a)
                    dxiac(d) = x(d,b) - (x(d,a) + surf_norm(d,s)*mv_in)
                    driac    = driac + dxiac(d)*dxiac(d)
            enddo
            driac=sqrt(driac)
        
            call kernel(driac,dxiac,hsml(b),w,tdwdx)
        
            call fncnApproxOperator(T(s),temp(b),mass(b),rho(b),w)
       
            call fncnApproxOperator(alpha(s),1.D0,mass(b),rho(b),w)
        endif
    
    enddo


    do s=1,etotal
        
         a = nedge_rel_edge(s)
        
        ! Apply neumann boundary condition for wall with dT/dn= r
        ! Here T_si=  ∑_b(T_b*W_sb*m_b/ρ_b)/α_s is calculated to
        ! determine T_s = T_si- dT/dn * (mv_in)
        if ( etype(s) .eq. etype_thermal_neumann) then
            
            temp(a)=T(s)/alpha(s) - bdryVal_temp(s)*mv_in

        endif
        
        !Apply Dirichlet boundary condition for wall with T= T_boundary
        
        if ( etype(s) .eq. etype_thermal_dirichlet) then
            
            temp(a)=bdryVal_temp(s)

        endif
    
    enddo
    
    deallocate(alpha,T)

end
    
    