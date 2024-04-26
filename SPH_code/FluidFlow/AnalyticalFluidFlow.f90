!****************************************************************************
!
!  SUBROUTINE: Transient_Fluid_Flow
!
!  PURPOSE:  Subroutine that calculates analytical value of  
!         various fluid flow
!
!   CREATED:        01/16/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/16/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine AnalyticalFluidFlow(itimestep,dt,vel,bdryVals)

    use config_parameter, only: SPH_dim,pi, itype_real_min, itype_real_max, &
        & NumericalSimCase,F_ext, rho_init, mu_const
    use particle_data, only: ntotal, etotal, itype,x, rho, mu
    
    implicit none
    
    real(8), intent(in):: dt
    real(8) :: vel(ntotal)
    integer(4), intent(in)::itimestep
    logical, intent(in)::bdryVals
    integer(4) :: a, n, n_max 
    real(8) nu, Bn, BBn, lbd, time, width !Bn,lbd, K, T1, T2, f2, T0, Tss, BBn , time
    
    vel=0.D0
    Bn=0.D0
    BBn=0.D0
    
    width= 0.D0 ! Input width here 
    ! Case1 value is 1, and case 2 value is 2

    if (NumericalSimCase .eq. 3) then
        
    
        n_max= 2000

        time=dt*itimestep
    
        do a=1,ntotal
            if(((itype(a) .gt. itype_real_min) .and. (itype(a) .le. itype_real_max)) .or. bdryVals ) then
                
                  
                BBn=0.D0   
                if ( (rho(a) .eq. 0.D0) .or. (mu(a) .eq. 0.D0)) then
                    rho(a)= rho_init
                    mu(a)=mu_const
                endif
                
                nu=mu(a)/rho(a)
                
                vel(a)=-F_ext(1)*x(2,a)*(x(2,a)-width)/(2.D0*nu)
                
                do n=0,n_max
                    
                    lbd=4.D0*F_ext(1)*(width**2)/((nu*pi**3)*(2.D0*n+1.D0)**3)
                    Bn= dsin(pi*x(2,a)*(2.D0*n+1.D0)/width)*dexp(-((2.D0*n+1.D0)**2)*(pi**2)*nu*time/(width**2))        
                    BBn=BBn - Bn*lbd     
                    !write(*,*) n
                enddo        
                vel(a)= (vel(a) + BBn)
                !write(*,*) vx(1,a)
            endif
                
        enddo
        
    elseif (NumericalSimCase .eq. 4) then
        
        
    endif
    
        
end
    
    