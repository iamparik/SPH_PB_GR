!****************************************************************************
!
!  SUBROUTINE: WCSPH_minTimeStep
!
!  PURPOSE:  Subroutine to calculate min time step to be used in simulations
!
!   CREATED:        10/03/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  10/03/2023         by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine WCSPH_minTimeStep(dt)
    
    use config_geometry, only: SPH_dim, time_step, hsml_const, rho_init, F_ext, g_const, mu_const, c_sound
    
    implicit none
    
    real Fext_tot(SPH_dim), dt_min
    integer yesNo
    Fext_tot(:)=F_ext
    Fext_tot(2)=Fext_tot(2) - g_const
    
    dt_min=minval((/ 0.25D0*hsml_const/c_sound, 0.25D0*(hsml_const/norm2(Fext_tot))**0.5D0,0.125D0*hsml_const**2*rho_init/mu_const /))
    write(*,*) ' -------------------time step -----------------------------'
    write(*,*) 'the minimum permissible value of time step is : ', dt_min , 'for c_sound', c_sound
    write(*,*) ' time step used in simulation is : ', time_step
    write(*,*) ' enter 0 to use default time step, and 1 for suggested time step :'
    write(*,*) yesNo
    write(*,*) ' -------------------time step -----------------------------'
    
    if(yesNo .eq. 0) then
        dt = time_step
    else
        dt = dt_min
    endif
    
    
    
    
    end
    