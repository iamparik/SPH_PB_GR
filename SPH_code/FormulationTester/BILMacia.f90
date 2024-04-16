!
!  SUBROUTINE: BILMacia
!               ⟨⟨𝜵∙kk 𝜵𝑓_s ⟩⟩BILMacia=  kk_s/𝛾_𝑎 ∑_𝑠〖(𝛁𝑓_s)∙𝛁𝛾_𝑎𝑠〗 
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        09/02/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/09/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine BILMacia(lapf, f0, fs0, x, xs, gma, SPH_dim, eniac, epair_a, epair_s, dgmas, etype, &
                & etotal, ntotal, kk_s)
    use config_geometry, only: WallBoundaryLayer
    implicit none
    integer(4), intent(in):: SPH_dim, eniac, etotal, ntotal
    real(8), intent(in)::gma(ntotal), x(SPH_dim,ntotal) ,xs(SPH_dim,etotal)
    real(8), intent(in):: f0(ntotal) ,fs0(etotal)
    real(8), intent(in):: dgmas(SPH_dim,eniac), kk_s(etotal) 
    integer(4), intent(in):: epair_a(eniac), epair_s(eniac), etype(etotal)
    real(8), intent(inout):: lapf(ntotal)
    integer(4):: k,a,b,s,d
    real(8) :: lap_coeff, r_as
    

   
    ! Call the boundary itnegral terms of the Laplacian of function interpolations
    !  2/𝛾_𝑎 ∑_𝑠〖( 𝛁𝑓_𝑎)∙𝛁𝛾_𝑎𝑠〗 
     do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)
        
        r_as=norm2(x(:,a)-xs(:,s))
        
        if(r_as .ne.0.D0) then
            !the laplacian is only valid for particles close to boundary and not boundary points themselves
            do d=1,SPH_dim
                lap_coeff=2.D0*kk_s(s)*(f0(a)-fs0(s))*(x(d,a)-xs(d,s))/(r_as**2)
                call fncnApproxOperator(lapf(a),lap_coeff,1.D0,1.D0,-dgmas(d,k)/gma(a))
            enddo
        endif
        
    enddo
     
end



