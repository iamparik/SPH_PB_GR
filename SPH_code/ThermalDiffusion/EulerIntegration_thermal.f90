!****************************************************************************
!
!  SUBROUTINE: EulerIntegration_thermal
!
!  PURPOSE:  Subroutine to march in time, using a backward marching scheme  
!         for thermal diffusion problem  
!
!   CREATED:        08/20/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  09/25/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine EulerIntegration_thermal(itimestep,dt)

    use config_parameter, only: SPH_dim,itype_real_max, itype_real_min , &
        & BILtype, save_step
    use particle_data, only: ntotal, etotal, itype,x, temp, dtemp, temp_prev,  &
        &  rho, gamma_mat_inv, pBC_edges
    
    implicit none
    
    
    real(8), intent(in):: dt
    integer(4), intent(in)::itimestep
    integer(4) :: a, i, j 
    logical:: AnalyticalSimulation= .false.
    
    
    !Either run the Analytical Simulation/SPHSimulation
    if(AnalyticalSimulation) then
        !if runnign Analytical Simulation, make sure the solution for the given case exists and 
        ! is included in the file 
        
        ! The last input defines the analytical case icnluded in the code
        !call AnalyticalThermalDiffusion(itimestep,dt,temp,AnalyticalSimulation)
        call AnalyticalThermalDiffusion(itimestep,dt,temp(1:ntotal),AnalyticalSimulation)
        
    else
    
        ! if there was periodic Bcapplied, for periodic particles and points update variables used 
        if (Allocated(pBC_edges)) then
            call PeriodicParameter(temp)
        endif

        ! Apply Thermal Boundary Condition 
        !call thermalBoundaryConditionFerrand
        call thermalBoundaryConditionFerrandInspired
    
        ! if there was periodic Bcapplied, for periodic particles and points update the other boundary condition applied 
        if (Allocated(pBC_edges)) then
            call PeriodicParameter(temp)
        endif

    
        ! Thermal Diffusion Equation is now solved
        allocate(dtemp(ntotal),temp_prev(ntotal))
        dtemp=0.D0
        temp_prev(1:ntotal)=temp(1:ntotal)

        !The below needs to be changed according to the formulation used for laplacian
        ! Use 2 for PCBI, 8 for BIL-USAW, 9 for BIL1
        call thermalDiffusionLaplacianOperator(temp_prev,dtemp,BILtype)

        do a=1,ntotal
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min))   then
                temp(a)=temp_prev(a) + (1.D-1)*dt*dtemp(a)
            endif
            if(itype(a) .eq. itype_real_min)  temp(a)=0.D0
        enddo
    
        ! if there was periodic Bcapplied, for periodic particles and points update temperature after tranisent update 
        if (Allocated(pBC_edges)) then
            call PeriodicParameter(temp)
        endif
    
        deallocate(dtemp,temp_prev)
    endif
    
    if ((mod(itimestep,save_step).eq.0) .or. (itimestep.eq.1)) call output_time_thermal(itimestep, dt)
    

end
     

    