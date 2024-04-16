!
!  SUBROUTINE: BIL_Laplacian
!               ⟨⟨𝜵∙ 𝜵𝑓_𝑎 ⟩⟩BIL1=  2/𝛾_𝑎 ∑_𝑠〖(𝛁𝑓_𝑎)∙𝛁𝛾_𝑎𝑠〗 
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        09/02/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/09/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine BIL1(lapf, df0, gma, SPH_dim, eniac, epair_a, epair_s, dgmas, etype, etotal, ntotal,etype_periodic, Neumann, Dirichlet, bdryVal, coeff)
    implicit none
    integer(4), intent(in):: SPH_dim, eniac, etotal, ntotal, Neumann, Dirichlet, etype_periodic
    real(8), intent(in)::gma(ntotal), bdryVal(etotal), coeff
    real(8), intent(in):: df0(SPH_dim,ntotal)
    real(8), intent(in):: dgmas(SPH_dim,eniac) 
    integer(4), intent(in):: epair_a(eniac), epair_s(eniac), etype(etotal)
    real(8), intent(inout):: lapf(ntotal)
    integer(4):: k,a,b,s,d
    real(8) :: lap_coeff
    

   
    ! Call the boundary itnegral terms of the Laplacian of function interpolations
    !  2/𝛾_𝑎 ∑_𝑠〖( 𝛁𝑓_𝑎)∙𝛁𝛾_𝑎𝑠〗 
     do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)
        
        if (mod(etype(s),etype_periodic) .eq. Neumann) then
            !df/dn=lap_coeff
            lap_coeff=coeff*bdryVal(s)*norm2(dgmas(:,k))
            call fncnApproxOperator(lapf(a),lap_coeff,1.D0,1.D0,-1.D0/gma(a))
      
        elseif(mod(etype(s),etype_periodic) .eq. Dirichlet) then
            do d=1,SPH_dim
                lap_coeff=coeff*df0(d,a)*dgmas(d,k)
                call fncnApproxOperator(lapf(a),lap_coeff,1.D0,1.D0,-1.D0/gma(a))
            enddo
        endif       
   
    enddo
     
end



