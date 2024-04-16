!
!  SUBROUTINE: BIL1_Laplacian
!               ⟨⟨𝜵∙ 𝜵𝑓_𝑎 ⟩⟩= −2/𝛾_𝑎  ∑_𝑏〖(𝑓_𝑎 − 𝑓_𝑏)((𝒙_𝑎 − 𝒙_b ))/(𝑟_𝑎𝑏^2 )∙𝛁𝑊 𝑚_𝑏/𝜌_𝑏 〗
!                                                   + 2/𝛾_𝑎 ∑_𝑠〖𝛁𝑓_𝑠∙𝛁𝛾_𝑎𝑠〗 
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        09/02/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/06/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine LaplacianBIL_Macia(lapf, f0, fs0, x, xs, gma, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, dgmas, mass, rho, itype, ntotal, etotal)
    implicit none
    integer(4), intent(in):: ntotal, SPH_dim, niac, eniac, etotal
    integer(2), intent(in):: itype(ntotal)
    real(8), intent(in):: f0(ntotal),fs0(etotal), gma(ntotal), mass(ntotal), rho(ntotal)
    real(8), intent(in)::   x(SPH_dim,ntotal), xs(SPH_dim,etotal)
    real(8), intent(in):: dwdx(SPH_dim,niac), dgmas(SPH_dim,eniac) 
    integer(4), intent(in):: pair_i(niac), pair_j(niac)
    integer(4), intent(in):: epair_a(eniac), epair_s(eniac)
    real(8), intent(inout):: lapf(ntotal)
    integer(4)::k,a,b,s,d
    real(8) ::  delf, r_ab, lap_coeff
    
          
    ! Call the interpolation operation for finding the Laplacian of a function
    ! −2/𝛾_𝑎  ∑_𝑏〖(𝑓_𝑎 − 𝑓_𝑏)((𝒙_𝑎 − 𝒙_b ))/(𝑟_𝑎𝑏^2 )∙𝛁𝑊 𝑚_𝑏/𝜌_𝑏 〗
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)
        
        r_ab=norm2(x(:,a)-x(:,b))
        if ( .not. (isNAN(1.D0/r_ab))) then
            do d=1,SPH_dim
            
                lap_coeff=2.D0*(f0(a)-f0(b))*(x(d,a)-x(d,b))/(r_ab**2)
                call fncnApproxOperator(lapf(a),lap_coeff,mass(b),rho(b),dwdx(d,k)/gma(a))
                
                lap_coeff=2.D0*(f0(b)-f0(a))*(x(d,b)-x(d,a))/(r_ab**2)
                call fncnApproxOperator(lapf(b),lap_coeff,mass(a),rho(a),-dwdx(d,k)/gma(b))
                
            enddo
        endif

     enddo
     
   
    ! Call the boundary itnegral terms of the Laplacian of function interpolations
    !  2/𝛾_𝑎  ∑_s▒〖(f_a-f_s )  ((x_a-x_s )  )/|x_a-x_s |^2 ∙n_s  |∇𝛾_𝑎s| 〗
     do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)
        
        r_ab=norm2(x(:,a)-xs(:,s))
        
        if(r_ab .ne.0.D0) then
            !the laplacian is only valid for particles close to boundary and not boundary points themselves
            do d=1,SPH_dim
                lap_coeff=2.D0*(f0(a)-fs0(s))*(x(d,a)-xs(d,s))/(r_ab**2)
                call fncnApproxOperator(lapf(a),lap_coeff,1.D0,1.D0,-dgmas(d,k)/gma(a))
            enddo
        endif

    enddo
     
end



