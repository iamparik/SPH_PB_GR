!
!  SUBROUTINE: GradientKulasegaram
!			    This calculates the gradient as:
!               ⟨𝜵𝑃⟩_𝑎≅𝜌_𝑎 ∑_b〖m_b (𝑃_𝑎/(γ_a 𝜌_𝑎^2 )+𝑃_𝑏/(γ_b 𝜌_𝑏^2)) 𝜵_a W_𝑎𝑏 〗−(𝑃_𝑎/γ_a)∑_s 𝜵γ_as 
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        06/08/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/25/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine GradientKulasegaram(d1f,f0, fs0, Gma, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, dgmas, mass, rho, rhos, itype, ntotal, num_edges)
    use particle_data, only: vx    
    use config_geometry,    only: dx_r,g_const

    implicit none
    integer(4), intent(in):: ntotal, SPH_dim, niac, eniac, num_edges
    integer(2), intent(in):: itype(ntotal)
    real(8), intent(in):: f0(ntotal), Gma(ntotal), mass(ntotal), rho(ntotal), rhos(num_edges), fs0(num_edges)
    real(8), intent(in):: dwdx(SPH_dim,niac), dgmas(SPH_dim,eniac) 
    integer(4), intent(in):: pair_i(niac), pair_j(niac) 
    integer(4), intent(in):: epair_a(eniac), epair_s(eniac)
    real(8), intent(inout):: d1f(SPH_dim,ntotal)
    integer(4)::k,a,b,s, d
    real(8) ::  delf,temp_dwdx(SPH_dim),CdwdxA(SPH_dim),CdwdxB(SPH_dim)
    real(8) ::  temp_dgmas(SPH_dim),Cdgmas(SPH_dim)
    
          
    ! Call the interpolation operation for first order derivative of a function
    ! (ρ_a/𝛾_𝑎)∑_b〖(f_a/ρ_b**2 + f_b/ρ_b**2 ) ∂_j W_ab m_b 〗
    
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)

        temp_dwdx(:)= dwdx(:,k)
        CdwdxA= temp_dwdx!/Gma(a)
        
        temp_dwdx(:)= -dwdx(:,k)
        CdwdxB= temp_dwdx!/Gma(b)
        
        do d= 1,SPH_dim
            delf= rho(a)*(f0(b)/(Gma(b)*rho(b)**2)+f0(a)/(Gma(a)*rho(a)**2))
            call fncnApproxOperator(d1f(d,a),delf,mass(b),1.D0,CdwdxA(d))
            delf= rho(b)*(f0(b)/(Gma(b)*rho(b)**2)+f0(a)/(Gma(a)*rho(a)**2))
            call fncnApproxOperator(d1f(d,b),delf,mass(a),1.D0,CdwdxB(d))
        enddo
        
        
     enddo
     
   
    ! call the boundary itnegral terms of the first order derivative  of function interpolations
    ! [γ_a^-1 ]_ij ∑_s∫_(∂(ω ∩ ω_w )_s〖(f_a-f_(s^' ) ) w_(as' ) n_(s_j ) ds' 〗
    do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)
        
        temp_dgmas(:)= dgmas(:,k)
        cdgmas= temp_dgmas/gma(a)
        
         do d=1,sph_dim
             delf= -1.d0*f0(a)
            call fncnapproxoperator(d1f(d,a),delf,1.d0,1.d0,cdgmas(d))
        enddo
    
    enddo
     
   
    end
    
    
    
    
    