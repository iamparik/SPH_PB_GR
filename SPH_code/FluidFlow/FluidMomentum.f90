!****************************************************************************
!
!  SUBROUTINE: FluidMomentum
!
!  PURPOSE:  Subroutine to calculate accelaration of particles subjected to forces
!
!   CREATED:        08/20/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  06/29/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine FluidMomentum(dstress, x, vx, rho, p, itype, ntotal)
    use config_parameter, only: SPH_dim
    use config_geometry, only: BILtype, PrsrGradtype, artViscType
    use particle_data, only: maxn, pBC_edges

    implicit none
    real(8), intent(inout):: dstress(SPH_dim,ntotal)
    real(8), intent(in):: x(SPH_dim,maxn)
    real(8), intent(in):: vx(SPH_dim,maxn)
    real(8), intent(in):: rho(maxn)
    real(8), intent(in):: p(maxn)
    integer(2), intent(in):: itype(maxn)
    integer(4), intent(in):: ntotal
    
!--------------Calculate Deviatoric Stress --------------
    !Calculate accelaration due to Viscouse Stress
    call ViscousStressOperator(dstress, vx, rho, BILtype)
        
    !--------------Calculate acceleration due to External Force --------------        
    !Calculate accelaration due to external force
    call ExternalForceAcceleration(dstress, rho, itype(1:ntotal), ntotal)

        
    !--------------Calculate Pressure Gradient--------------
        
    ! Calculate Pressure using EOS
    call pressureEOS
        
    ! For periodic condition, update the calculated pressure of all periodic reals
    if (Allocated(pBC_edges)) then
        call PeriodicParameter(p)
    endif
        
    call boundaryPressureUpdate
        
    ! For periodic condition, update the calculated pressure of all periodic reals
    if (Allocated(pBC_edges)) then
        call PeriodicParameter(p)
    endif
        
        
    ! Calculate Accelaration due to Pressure field using Pressure gradient
    call PressureGradientOperator(dstress,p, rho, PrsrGradtype)
        
    !----------------Calculate Artificial Viscosity
        
    call artificialViscosityOperator(dstress, vx, x, rho, artViscType)
    
end