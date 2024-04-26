!****************************************************************************
!
!  SUBROUTINE: BoundaryPressureUpdate
!
!  PURPOSE:  Subroutine to apply pressure boudnary condition for a   
!         for a flow problem inspired from ferrands homogeneuous neumann bc 
!           𝑝_𝑠/𝜌_𝑠   =(∑1_(𝑏≠𝑠)▒〖(𝑝_𝑏/𝜌_𝑏    −  𝒈∙(𝒙_𝑏−𝒙_𝑠 )  +  (|𝒗_𝑏 |^2−|𝒗_𝑠 |^2)/2)  𝑚_𝑏/𝜌_𝑏   𝑊_(𝑠𝑖,𝑏) 〗)/(∑1_(𝑏≠𝑠)▒𝑚_𝑏/𝜌_𝑏   𝑊_(𝑠𝑖,𝑏) )    
!
!   CREATED:        10/03/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/13/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine boundaryPressureUpdate

    use config_parameter,   only: SPH_dim, etype_SolidWall1, etype_SolidWall2, etype_FreeSurface1, &
         & pi, dx_r,g_const,hsml_const, c_sound, prsrBdryType
    use particle_data,      only: ntotal,etotal,etype,&
        & epair_a, epair_s, eniac, p, hsml, mass, &
        & rho, nedge_rel_edge, x,vx, bdryVal_prs,surf_norm, &
        & p_counter, gamma_cont
    
    
    implicit none
      
    real(8) ::  tdwdx(SPH_dim), dxiac(SPH_dim), g(SPH_dim), driac, w ,Psb, mv_in, reduced_hsml
    real(8),DIMENSION(:),ALLOCATABLE :: alpha, psi  
    real(8),DIMENSION(:,:),ALLOCATABLE :: xsi
    integer(4) :: a,b,s,k,d, max_a      
    real(8) :: r_V, a_V, theta, V_ratio, max_V_ratio 
    
    max_V_ratio=1.D0
    max_a=1
    
    if(.not. allocated(P_counter) ) allocate(P_counter(ntotal))
    P_counter=0.D0

    ! assuming each edge has one edge point, we initialize alpha and psi
    ! for more than one edge point per edge initialize(eppe) as 
    ! alpha(etotal, eppe), psi(etotal, eppe)
    allocate(alpha(etotal), psi(etotal),xsi(SPH_dim,etotal))
    alpha=0.D0
    psi=0.D0
    xsi=0.D0
    
    ! Use gravity with correct direction
    g(:)=0.D0
    g(2)= - abs(g_const)
    
    
    ! when mv_in is zero, we get Ferrand's homogeneous pressure neumann bc
    mv_in=dx_r*0.D0

    !call sml_mult_factor(scale_k)
      
    reduced_hsml= hsml_const!/2.D0  
    
    
    do s=1,etotal
        a = nedge_rel_edge(s)
        do d=1,SPH_dim
            xsi(d,s) = x(d,a) + surf_norm(d,s)*mv_in   
        enddo
    enddo
    

    do k=1,eniac
        b=epair_a(k)
        s=epair_s(k)
    
        a = nedge_rel_edge(s)
        
        ! Apply neumann boundary condition for wall with dPd/dn=0
        ! Here calcualte T_si=∑_b(T_b*W_sb*m_b/ρ_b)  & 
        ! α_s=∑_b(W_sb*m_b/ρ_b)
        if ( (etype(s) .eq. etype_SolidWall1) .or. (etype(s) .eq. etype_SolidWall2)) then
        
            driac = 0.D0 
            
            do d=1,SPH_dim
                !since surf_norm is inward we displace inwards by adding to x(d,a)
                dxiac(d) = x(d,b) - xsi(d,s)
                driac    = driac + dxiac(d)*dxiac(d)
            enddo
            driac=sqrt(driac)
        
            !call kernel(driac,dxiac,hsml(b),w,tdwdx)
            call kernel(driac,dxiac,reduced_hsml,w,tdwdx)
            
            
            if (prsrBdryType .eq. 1) then
                Psb= (p(b)/rho(b) &
                    & - dot_product(g, (x(:,b)-x(:,a))) &
                    & + 0.5D0*(dot_product(vx(:,b),vx(:,b))-dot_product(vx(:,a),vx(:,a)))  &
                    !& + 0.5D0*(dot_product(vx(:,b)-vx(:,a), surf_norm(:,s))**2)  &
                     &  )
            
            elseif(prsrBdryType .eq. 2) then
                Psb= (p(b)/rho(b) &
                    & - dot_product(g, (x(:,b)-x(:,a))) &
                    & - c_sound*dot_product(vx(:,b)-vx(:,a), surf_norm(:,s))  &
                     &  )
            endif
        
            call fncnApproxOperator(psi(s),Psb,mass(b),rho(b),w)

            ! alpha_s =∑1_(𝑏≠𝑠)𝑚_𝑏/𝜌_𝑏 𝑊_si,b  
            call fncnApproxOperator(alpha(s),1.D0,mass(b),rho(b),w)
            
        endif
    
    enddo


    do s=1,etotal
        
         a = nedge_rel_edge(s)
        
        ! Apply neumann boundary condition for wall with dP/dn= r
        ! Here P_si=  ∑_b(P_b*W_sb*m_b/ρ_b)/α_s is calculated to
        ! determine P_s = P_si- dP/dn * (mv_in)
        if (etype(s) .eq. etype_SolidWall1)  then
            p(a)= rho(a)* psi(s)/alpha(s)
            
            if(alpha(s) .eq. 0.D0)   p(a)= rho(a)* psi(s)
           ! p(a)=rho(a)*c_sound**2 
            
        endif
        
        !Apply Dirichlet boundary condition for wall with psi= P_boundary        
        if ( etype(s) .eq. etype_FreeSurface1) then                    
            p(a)=bdryVal_prs(s)
        
        endif
        
        
    
    enddo
    
    deallocate(alpha, psi,xsi)

end
    
    