!****************************************************************************
!
!  SUBROUTINE: BoundaryDensityUpdate
!
!  PURPOSE:  Subroutine to apply density boundary condition in weakly Compressible or Commpressible Schemse  
!
!   CREATED:        10/03/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/13/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine boundaryDensityUpdate

    use config_parameter,   only: SPH_dim, etype_SolidWall1, etype_FreeSurface1, &
        & dx_r,g_const,hsml_const, rho_init
    use particle_data,      only: ntotal,etotal,etype,&
        & epair_a, epair_s, eniac, p, hsml, mass, &
        & rho, nedge_rel_edge, x,vx !, bdryVal_rho
    
    implicit none
      
    real(8) ::  tdwdx(SPH_dim), dxiac(SPH_dim), driac, w, reduced_hsml!scale_k  
    real(8),DIMENSION(:),ALLOCATABLE :: alpha, rhosi  
    integer(4) :: a,b,s,k,d      

    ! assuming each edge has one edge point, we initialize alpha and psi
    ! for more than one edge point per edge initialize(eppe) as 
    ! alpha(etotal, eppe), psi(etotal, eppe)
    allocate(alpha(etotal), rhosi(etotal))
    alpha=0.D0
    rhosi=0.D0

    
   !call sml_mult_factor(scale_k)
      
    reduced_hsml= hsml_const!/2.D0
   

    do k=1,eniac
        b=epair_a(k)
        s=epair_s(k)
    
        a = nedge_rel_edge(s)
        
        !Calculate  𝜌_𝑠=(∑1_(𝑏≠𝑠)〖(𝜌_𝑏 )  𝑚_𝑏/𝜌_𝑏   𝑊〗)/(∑1_(𝑏≠𝑠)𝑚_𝑏/𝜌_𝑏 𝑊)
        
        if(etype(s) .eq. etype_SolidWall1) then
        
            driac = 0.D0 
            
            do d=1,SPH_dim
                dxiac(d) = x(d,b) - x(d,a)
                driac    = driac + dxiac(d)*dxiac(d)
            enddo
            driac=sqrt(driac)
        
            !call kernel(driac,dxiac,hsml(b),w,tdwdx)
            call kernel(driac,dxiac,reduced_hsml,w,tdwdx)               
        
            call fncnApproxOperator(rhosi(s),rho(b),mass(b),rho(b),w)
            
            ! alpha_s =∑1_(𝑏≠𝑠)𝑚_𝑏/𝜌_𝑏 𝑊            
            call fncnApproxOperator(alpha(s),1.D0,mass(b),rho(b),w)
        endif
    
    enddo


    do s=1,etotal
        
         a = nedge_rel_edge(s)

        ! Apply neumann boundary condition for wall with dP/dn= r
        ! Here 𝜌_si=  ∑_b(𝜌_b*W_sb*m_b/ρ_b)/α_s is calculated to
        if (etype(s) .eq. etype_SolidWall1)  then
            
            rho(a)= 1000.D0!rhosi(s)/alpha(s)
            
        endif
        
        !Apply Dirichlet boundary condition for wall with psi= P_boundary
        
        if ( etype(s) .eq. etype_FreeSurface1) then                    
            ! change below to
            rho(a)= rho_init !bdryVal_rho(s)
        
        endif
        
        ! The below is added to avoided NAN
        
        if ((rho(a) .le. 0.D0) .or. isNAN(rho(a))) then
            
            rho(a) = rho_init
        endif
        
    
    enddo
    
    deallocate(alpha,rhosi)

end