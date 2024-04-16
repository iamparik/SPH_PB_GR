!****************************************************************************
!
!  SUBROUTINE: Transient_Thermal_Diffusion_Analytical
!
!  PURPOSE:  Subroutine that calculates analytical value of  
!         Temperature distribution  
!
!   CREATED:        08/31/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/23/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine AnalyticalThermalDiffusion(itimestep,dt,temp,bdryVals)

    use config_parameter, only: pi, itype_real_min, itype_real_max
    use config_geometry, only: NumericalSimCase
    use particle_data, only: ntotal, etotal, itype,x
    
    implicit none
    
    real(8), intent(in):: dt
    real(8) :: temp(ntotal)
    integer(4), intent(in)::itimestep
    logical, intent(in)::bdryVals
    integer(4) :: a, n, n_max 
    real(8) Bn,lbd, K, T1, T2, f2, T0, Tss, BBn , time, width
    
    temp=0.D0
    
    width= 0.D0 ! Input width here 

    
    ! Case1 value is 1, and case 2 value is 2

    if (NumericalSimCase .eq. 1) then
        T1=1000.D0+273.D0
        T2=500.D0+273.D0
        T0=750.D0+273.D0
        K=1.D-1
    
        n_max= 2000

        time=dt*itimestep
    
        do a=1,ntotal
            if(((itype(a) .gt. itype_real_min) .and. (itype(a) .le. itype_real_max)) .or. bdryVals ) then
                
                Tss=T1-(T1-T2)*x(2,a)/width  
                BBn=0.D0        
                do n=1,n_max
                    lbd=(n*pi/width)**2
                    Bn= 2.D0*((T2-T0)*((-1)**n) +T0-T1)/(width*sqrt(lbd))        
                    BBn=BBn + Bn*exp(-time*lbd*K)*sin(sqrt(lbd)*x(2,a))            
                enddo        
                temp(a)= Tss + BBn
            endif
                
        enddo
        
    else if (NumericalSimCase .eq. 2) then
        
        T1=1000.D0+273.D0
        f2=-(500.D0+273.D0)
        T0=750.D0+273.D0
        K=1.D-1
    
        n_max= 2000

        time=dt*itimestep
    
        do a=1,ntotal
            if(((itype(a) .gt. itype_real_min) .and. (itype(a) .le. itype_real_max)) .or. bdryVals ) then
        
                Tss=T1+f2*x(2,a)  
                BBn=0.D0        
                do n=1,n_max
                    lbd=((2*n-1)*pi/(2*width))**2
                    Bn=2.D0*(-(f2/sqrt(lbd))*(-1)**(n-1) +T0-T1)/(width*sqrt(lbd))
                    BBn=BBn + Bn*exp(-time*lbd*K)*sin(sqrt(lbd)*x(2,a))
                enddo        
                temp(a)= Tss + BBn
            
            endif
                
        enddo

    endif
    
        
end
    
    