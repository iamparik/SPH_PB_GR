!****************************************************************************
!
!  SUBROUTINE: maxTimeStep
!
!  PURPOSE:  Subroutine to calculate max time step to be used in simulations
!
!   CREATED:        10/03/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/29/2024         by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine minTimeStep
    
    use config_parameter, only: SPH_dim, time_step, hsml_const, &
        & rho_init, F_ext, g_const, mu_const, c_sound
    
    implicit none
    
    real Fext_tot(SPH_dim), dt_max
    
    Fext_tot(:)=F_ext
    Fext_tot(2)=Fext_tot(2) - g_const
    
    ! Maximum possible time step, is determined as the minimum value of the various theoretical criteria
    dt_max=minval((/ 0.25D0*hsml_const/c_sound, 0.25D0*(hsml_const/norm2(Fext_tot))**0.5D0,0.125D0*hsml_const**2*rho_init/mu_const /))
    write(*,*) ' -------------------time step -----------------------------'
    write(*,*) 'the minimum permissible value of time step is : ', dt_max , 'for c_sound', c_sound
    write(*,*) ' tiem step used in simulation is : ', time_Step
    write(*,*) ' -------------------time step -----------------------------'
    
    end
    