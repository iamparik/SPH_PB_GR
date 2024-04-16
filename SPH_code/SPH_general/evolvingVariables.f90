!****************************************************************************
!
!  SUBROUTINE: evolvingVariables
!
!  PURPOSE:  Variables denoting systems energy, or L2norms of variables of itnerest are calculated for each time step
!
!   CREATED:        12/09/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/09/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine evolvingVariables(itimestep,max_step)
    use config_parameter, only: itype_real_max, itype_real_min
    use config_geometry, only:  g_const
    use particle_data, only: ntotal, itype, delC, mass,rho, vx,x, gamma_cont, gamma_discrt, &
        & max_vel, KE, PE, TE, delC, delCAvg, delCMax, delCL2, xStart, KE_prev, &
        & ga_Max,ga_Avg,ga_L2, g_a_min, g_a_max,FirstCutoffStep,PP_variable
    implicit none
    
    integer(4), intent(in):: itimestep, max_step
    integer(4) :: a, b, c,d, energyStepT
       
    if(itimestep .eq. 1) KE_prev =0.D0   

    if (.not.(Allocated(TE))) then
        Allocate(KE(max_step),PE(max_step),TE(max_step),max_vel(max_step),&
            & delcAvg(max_step), delCMax(max_step),delCL2(max_step), &
            & ga_Max(max_step),ga_Avg(max_step),ga_L2(max_step), g_a_min(max_step), g_a_max(max_step))
        KE=0.D0
        PE=0.D0
        TE=0.D0
        max_vel=0.D0
        delcAvg=0.D0
        delcMax=0.D0
        delCL2=0.D0
        ga_Avg=0.D0
        ga_Max=0.D0
        ga_L2=0.D0
        g_a_min=1.D0
        g_a_max=1.D0
        !KE_change=0.D0
    endif
    
    if(mod(itimestep,max_step).eq. 0) then
        energyStepT=max_step
    else
        energyStepT=mod(itimestep,max_step)
    endif
    
    
    b=0
    c=0
    do a=1,ntotal
        
        if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
            KE(energyStepT)=KE(energyStepT) + 0.5D0* mass(a)*norm2(vx(:,a))**2 
            max_vel(energyStepT)= max(norm2(vx(:,a)), max_vel(energyStepT))           
            !PE(energyStepT)=PE(energyStepT) + norm2(xStart(:,a)-x(:,a))**2 !mass(a)*g_const*x(2,a) 
            delcAvg(energyStepT)=delcAvg(energyStepT) + norm2(delC(:,a))
            !delcMax(energyStepT)=max(norm2(delC(:,a)), delcMax(energyStepT)) 
            delCL2(energyStepT)=delCL2(energyStepT) + (mass(a)/rho(a))*norm2(delC(:,a))**2
            ga_Avg(energyStepT)=ga_Avg(energyStepT) + abs(gamma_cont(a)-gamma_discrt(a))
            ga_Max(energyStepT)=max(abs(gamma_cont(a)-gamma_discrt(a)), ga_Max(energyStepT)) 
            g_a_min(energyStepT)=min(gamma_cont(a)/gamma_discrt(a),g_a_min(energyStepT))
            g_a_max(energyStepT)= max(gamma_cont(a)/gamma_discrt(a),g_a_max(energyStepT))
            ga_L2(energyStepT)=ga_L2(energyStepT) + (mass(a)/rho(a))*abs(gamma_cont(a)-gamma_discrt(a))**2
            !TE(energyStepT)= max(norm2(xStart(:,a)-x(:,a)),TE(energyStepT))
            TE(energyStepT) = TE(energyStepT) + norm2(xStart(:,a)-x(:,a))
            if(gamma_cont(a) .lt. 1.D0) then
                !c=c+1
                !TE(energyStepT) = TE(energyStepT) + norm2(xStart(:,a)-x(:,a))!TE(energyStepT) + 0.5D0* mass(a)*norm2(vx(:,a))**2
                !delcAvg(energyStepT)=delcAvg(energyStepT) + norm2(delC(:,a))
                !delcMax(energyStepT)=max(norm2(delC(:,a)), delcMax(energyStepT)) 
            else
                d=d+1
                PE(energyStepT) = PE(energyStepT) + norm2(xStart(:,a)-x(:,a))
                delcMax(energyStepT)=max(norm2(delC(:,a)), delcMax(energyStepT)) 
                !KE(energyStepT)=KE(energyStepT) + 0.5D0* mass(a)*norm2(vx(:,a))**2 
                !max_vel(energyStepT)= max(norm2(vx(:,a)), max_vel(energyStepT)) 
            endif
            
             
            b=b+1
        endif      
        
    enddo      
    
    !PE(energyStepT)=sqrt(PE(energyStepT))
    TE(energyStepT)=TE(energyStepT)/b
    PE(energyStepT)=PE(energyStepT)/d
    !TE(energyStepT)= KE(energyStepT)-KE_prev 
    !KE_prev=KE(energyStepT)
    
    if(b .gt. 0) then
        delcAvg(energyStepT)=delcAvg(energyStepT)/b
        delCL2(energyStepT) = sqrt(delCL2(energyStepT))
        ga_Avg(energyStepT)=ga_Avg(energyStepT)/b
        ga_L2(energyStepT) = sqrt(ga_L2(energyStepT))
    endif
    
    !TE(energyStepT)= abs(delcAvg(energyStepT)-KE_prev)
    KE_prev=delcAvg(energyStepT)
    !TE(energyStepT)=KE(energyStepT)+PE(energyStepT)
    
    !if(mod(itimestep,1000) .eq. 0) then
    !    if(FirstCutoffStep) then
    !        PP_Variable= TE(energyStepT)
    !    else
    !        PP_variable= delcMax(energyStepT)
    !    endif
    !endif
    
    
end