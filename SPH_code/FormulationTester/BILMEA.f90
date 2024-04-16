!
!  SUBROUTINE: BILUSAW
!               ⟨⟨𝜵∙ 𝜵𝑓_𝑎 ⟩⟩BIL=  2/𝛾_𝑎 ∑_𝑠〖(𝛁𝑓_𝑎)∙𝛁𝛾_𝑎𝑠〗 
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        09/02/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/09/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine BILMEA(lapf, f0, fs0, x, xs, gma, SPH_dim, eniac, epair_a, epair_s, &
            & dgmas, etype, etotal, ntotal ,etype_periodic, Neumann, Dirichlet, &
            & bdryVal, bdryValDim, surf_norm, kk, kk_s)

    use config_geometry, only: dx_r,WallBoundaryLayer
    implicit none
    integer(4), intent(in):: SPH_dim, eniac, etotal, ntotal, &
                & Neumann, Dirichlet, etype_periodic, bdryValDim
    real(8), intent(in)::gma(ntotal), x(SPH_dim,ntotal), xs(SPH_dim,etotal)
    real(8), intent(in):: f0(ntotal) ,fs0(etotal), surf_norm(SPH_dim, etotal)
    real(8), intent(in):: dgmas(SPH_dim,eniac),kk(ntotal),kk_s(etotal),bdryVal(bdryValDim)
    integer(4), intent(in):: epair_a(eniac), epair_s(eniac), etype(etotal)
    real(8), intent(inout):: lapf(ntotal)
    integer(4):: k,a,b,s,d
    real(8) :: lap_coeff, r_as
  
    ! Call the boundary itnegral terms of the Laplacian of function interpolations
    !  2/𝛾_𝑎 ∑_𝑠〖( 𝛁𝑓_𝑎)∙𝛁𝛾_𝑎𝑠〗 
    do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)
        
        if (mod(etype(s),etype_periodic) .eq. Neumann) then
            !df/dn=lap_coeff
            lap_coeff=(kk(a)+kk_s(s))*bdryVal(s)*norm2(dgmas(:,k))
            call fncnApproxOperator(lapf(a),lap_coeff,1.D0,1.D0,-1.D0/gma(a))
      
        elseif(mod(etype(s),etype_periodic) .eq. Dirichlet) then
            if(WallBoundaryLayer) then
                lap_coeff=(kk(a)+kk_s(s))*bdryVal(k)*norm2(dgmas(:,k))
                call fncnApproxOperator(lapf(a),lap_coeff,1.D0,1.D0,-1.D0/gma(a))
            else
                
                r_as=dot_product(x(:,a)-xs(:,s), surf_norm(:,s))
                
        
                if(r_as .ne.0.D0) then
                    
                    !if( abs(r_as) .lt. dx_r/2.D0) r_as= sign(dx_r/2.D0, r_as)
                    !the laplacian is only valid for particles close to boundary and not boundary points themselves
                    lap_coeff=norm2(dgmas(:,k))*(kk(a)+kk_s(s))*(f0(a)-fs0(s))/r_as
                    call fncnApproxOperator(lapf(a),lap_coeff,1.D0,1.D0,-1.D0/gma(a))

                endif
                
            endif
            
		endif
            
    enddo
    
     
end



