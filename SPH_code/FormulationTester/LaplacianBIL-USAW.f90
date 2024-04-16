!
!  SUBROUTINE: Ferrand_Laplacian
!               ⟨⟨𝜵∙𝑘_a 𝜵𝑓_𝑎 ⟩⟩= −1/𝛾_𝑎  ∑_𝑏〖(𝑘_𝑎+𝑘_𝑏)(𝑓_𝑎 − 𝑓_𝑏)((𝒙_𝑎 − 𝒙_b ))/(𝑟_𝑎𝑏^2 )∙𝛁𝑊 𝑚_𝑏/𝜌_𝑏 〗
!                                                   + 1/𝛾_𝑎 ∑_𝑠〖(𝑘_𝑎 𝛁𝑓_𝑎+𝑘_𝑠 𝛁𝑓_𝑠)∙𝛁𝛾_𝑎𝑠〗 
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        09/02/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/06/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine LaplacianBI_USAW(lapf, f0, df0, dfs0, x, kk, kk_s, gma, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, dgmas, mass, rho, itype, ntotal, etotal)
    implicit none
    integer(4), intent(in):: ntotal, SPH_dim, niac, eniac, etotal
    integer(2), intent(in):: itype(ntotal)
    real(8), intent(in):: f0(ntotal), gma(ntotal), mass(ntotal), rho(ntotal), kk(ntotal), kk_s(etotal)
    real(8), intent(in):: df0(SPH_dim,ntotal), dfs0(SPH_dim,etotal), x(SPH_dim,ntotal)
    real(8), intent(in):: dwdx(SPH_dim,niac), dgmas(SPH_dim,eniac) 
    integer(4), intent(in):: pair_i(niac), pair_j(niac)
    integer(4), intent(in):: epair_a(eniac), epair_s(eniac)
    real(8), intent(inout):: lapf(ntotal)
    integer(4)::k,a,b,s,d
    real(8) ::  delf, r_ab, lap_coeff
    
          
    ! Call the interpolation operation for finding the Laplacian of a function
    ! −1/𝛾_𝑎  ∑_𝑏〖(𝑘_𝑎+𝑘_𝑏)(𝑓_𝑎 − 𝑓_𝑏)((𝒙_𝑎 − 𝒙_b ))/(𝑟_𝑎𝑏^2 )∙𝛁𝑊 𝑚_𝑏/𝜌_𝑏 〗
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)
        
        r_ab=norm2(x(:,a)-x(:,b))
        if ( .not. (isNAN(1.D0/r_ab))) then
            do d=1,SPH_dim
            
                lap_coeff=(kk(a)+kk(b))*(f0(a)-f0(b))*(x(d,a)-x(d,b))/(r_ab**2)
                call fncnApproxOperator(lapf(a),lap_coeff,mass(b),rho(b),dwdx(d,k)/gma(a))
                
                lap_coeff=(kk(b)+kk(a))*(f0(b)-f0(a))*(x(d,b)-x(d,a))/(r_ab**2)
                call fncnApproxOperator(lapf(b),lap_coeff,mass(a),rho(a),-dwdx(d,k)/gma(b))
                
            enddo
        endif

     enddo
     
   
    ! Call the boundary itnegral terms of the Laplacian of function interpolations
    !  1/𝛾_𝑎 ∑_𝑠〖(𝑘_𝑎 𝛁𝑓_𝑎+𝑘_𝑠 𝛁𝑓_𝑠)∙𝛁𝛾_𝑎𝑠〗 
     do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)
        
        do d=1,SPH_dim
            lap_coeff=(kk(a)*df0(d,a)+kk_s(s)*dfs0(d,s))
           call fncnApproxOperator(lapf(a),lap_coeff,1.D0,1.D0,-dgmas(d,k)/gma(a))
        enddo
   
    enddo
     
end



