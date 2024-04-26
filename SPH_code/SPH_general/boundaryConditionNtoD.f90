!****************************************************************************
!
!  SUBROUTINE: boundaryConditionNtoD
!
!  PURPOSE:  Subroutine to apply Neumann BC approximately    
!         for certain type of formulations, as Dirichlet 
!
!   CREATED:        03/28/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  03/28/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine boundaryConditionNtoD(fncn_s,bdryVal,fncn)
    
    use config_parameter,   only: SPH_dim, etype_SolidWall1, etype_SolidWall2, &
                            & dx_r, hsml_const
    use particle_data,      only: ntotal,etotal,etype,&
        & epair_a, epair_s, eniac, hsml, mass, &
        & rho, nedge_rel_edge, x, surf_norm
    
    implicit none
      
    !real(8), allocatable, intent(inout) :: fncn_s(:)
    !real(8), allocatable, intent(in) :: bdryVal(:)
    real(8) :: fncn_s(etotal),bdryVal(etotal), fncn(ntotal)
    real(8) :: tdwdx(SPH_dim), dxiac(SPH_dim), driac, w , mv_in
    real(8),DIMENSION(:),ALLOCATABLE :: alpha
    integer(4) :: a,b,s,k,d      

    ! assuming each edge has one edge point, we initialize alpha 
    ! for more than one edge point per edge initialize(eppe) as 
    ! alpha(etotal, eppe),
    allocate(alpha(etotal))
    alpha=0.D0
    
    fncn_s=0.D0
    
    mv_in=dx_r*1/20.D0!2.D0*hsml_const!dx_r*2.D0 !dx_r*1/20.D0 !dx_r*2.D0 ! 0.D0


    do k=1,eniac
        b=epair_a(k)
        s=epair_s(k)
    
        a = nedge_rel_edge(s)
        
        ! Apply neumann boundary condition for wall with dT/dn=0
        ! Here calcualte T_si=∑_b(T_b*W_sb*m_b/ρ_b)  & 
        ! α_s=∑_b(W_sb*m_b/ρ_b)
        if ( etype(s) .eq. etype_SolidWall2 ) then
        
            driac = 0.D0 
            
            do d=1,SPH_dim
                    !since surf_norm is inward we displace inwards by adding to x(d,a)
                    dxiac(d) = x(d,b) - (x(d,a) + surf_norm(d,s)*mv_in)
                    driac    = driac + dxiac(d)*dxiac(d)
            enddo
            driac=sqrt(driac)
        
            call kernel(driac,dxiac,hsml(b),w,tdwdx)
        
            call fncnApproxOperator(fncn_s(s),fncn(b),mass(b),rho(b),w)
       
            call fncnApproxOperator(alpha(s),1.D0,mass(b),rho(b),w)
        endif
    
    enddo


    do s=1,etotal
        
         a = nedge_rel_edge(s)
        
        ! Apply neumann boundary condition for wall with df/dn= r
        ! Here f_si=  ∑_b(f_b*W_sb*m_b/ρ_b)/α_s is calculated to
        ! determine f_s = f_si- df/dn * (mv_in)
        if ( etype(s) .eq. etype_SolidWall2) then
            
            fncn_s(s)=fncn_s(s)/alpha(s) - bdryVal(s)*mv_in

        endif
        
        !Apply Dirichlet boundary condition for wall         
        if ( etype(s) .eq. etype_SolidWall1) then
            
            fncn_s(s)=bdryVal(s)

        endif
    
    enddo
    
    deallocate(alpha)

    
end