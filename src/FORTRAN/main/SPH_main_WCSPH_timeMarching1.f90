﻿program SPH3D
 
!-------------------------------------------------------
! Program for running Boundary Integral SPH
! Author- Parikshit B            
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/17/2022       by  PARIKSHIT BOREGOWDA 
!-------------------------------------------------------

use particle_data 
use config_parameter

implicit none

!--------------------------------------
!input_timeStep : An integer Paramter that decides the number of time step iterations 
!                to be run for the simulation
!maztimestep : An integer Paramter that logs the maximum time step iterations
!yesorno :   An integer paramter,that implies no for 0
!            and yes for 1
!ic1,ic2 :   Counters to note clock time 
!crate1, cmax1  :   Itneger kind 8 coutns time in Microseconds
!                    and is used for calling function SYSTEM_CLOCK
!--------------------------------------

real(8) dt
integer(4) input_timeStep,max_timeSteps, yesorno, input_file_type, current_ts
integer(8) :: ic1, crate1, cmax1, ic2
logical ::  runSPH = .true.
integer(4) :: k, a, b , d, s, i, j, Scalar0Matrix1
integer(4) :: correction_types, CF_density, ID_density, CF_pressure, ID_pressure, &
    & CF_BIL_visc, ID_BIL_visc,dirich0Neum1, num_bdry_var_seg, num_bdry_var_ve, &
    & CF_densDiff, ID_densDiff, ID_prsrBdryType
real(8) :: scalar_factor, Sca_Bdry_val, current_time, prsr_wall_compress,gamma_wall_cutoff, rho_wall_compress
real(8), DIMENSION(:), allocatable  :: F_a, F_b,DF_a, DF_b, Cdwdx_a, Cdwdx_b, Cdgmas,temp_vector, dxiac
real(8), DIMENSION(:,:,:), allocatable :: grad_vel
real(8), DIMENSION(:,:), allocatable :: matrix_factor, stress,bdryVal_ve, grad_rho ,grad_vel_s_temp, x_ve_temp
real(8), DIMENSION(:), allocatable :: div_vel, delx_ab, dens_diffusion, rho_temp
logical ::prsr_bdry_preCalc = .false.
real(8) ::temp_scalar, driac, w_temp


correction_types=10
   
!start system clock to evealuate time taken to run
! System_clock returns number of seconds from 00:00 CUT on 1 Jan 1970
    call system_clock(count=ic1, count_rate=crate1, count_max=cmax1)

! input particle configuration data
    call inputSPHConfig
    
    allocate(temp_vector(SPH_dim),dxiac(SPH_dim))
    
! Allocate variables important for physics simulation
    ALLOCATE(vx(SPH_dim,maxn), mass(maxn), rho(maxn), p(maxn), hsml(maxn), mu(maxn))
    vx=0.D0
    mass=0.D0
    rho= rho_init
    p=0.D0
    hsml = hsml_const
    mu =0.D0
    
    if(external_input_InitialCondition) call external_file_IC
    
! the size of the vector is the same as size of all variables assosciated 
! with particle position (including particle position)
    allocate(F_a(8)) 
    
! Initialize the variables
    do a=1,nreal
        
        F_a(1:SPH_dim) = x(:,a)
        F_a(SPH_dim+1:SPH_dim*2)= vx(:,a) 
        F_a(5)=rho(a)
        F_a(SPH_dim*2+2) = p(a) 
        F_a(SPH_dim*2+3) = hsml(a)
        F_a(SPH_dim*2+4) = mu(a) 
        
        call ICinputValue(F_a,8, itype(a))      
        
        vx(:,a) = F_a(SPH_dim+1:SPH_dim*2)
        rho(a) = F_a(5)
        p(a) = F_a(SPH_dim*2+2)
        hsml(a) = F_a(SPH_dim*2+3) 
        mu(a) = F_a(SPH_dim*2+4) 
        
        ! mass is calculated with volume from geoemtry file and density from initial condition file
        mass(a)= rho(a)*vol(a)
    enddo
    
! deallocate F_a so it can be sued later again
    deallocate(F_a)
    
! Allocate other variables used in time itnegration steps
    allocate(F_a(SPH_dim), F_b(SPH_dim), Cdwdx_a(SPH_dim), Cdwdx_b(SPH_dim), Cdgmas(SPH_dim), &
    & matrix_factor(SPH_dim,SPH_dim),  delx_ab(SPH_dim), DF_a(SPH_dim), DF_b(SPH_dim))

! theoretical maximum for time step.
    call maxTimeStep
    dt = time_step


! Maximum number of timesteps, and ith time step variables are initiatilized as 0 
    max_timeSteps=0
    itimestep=0

!The below requires the user to input on the terminal the maximum time step to run the simulation
    write(*,*)'  ***************************************************'
    write(*,*)'          Please input the additional number of time steps to run simulations'
    write(*,*)'  ***************************************************'
    read(*,*) max_timeSteps

!
    max_timeSteps=max_timeSteps+itimestep
    Allocate(xStart(SPH_dim,ntotal))
    xStart=x
    
    current_ts= itimestep
    
    !define the numbe of boundary variables that are used in the problem
    !this is velocity, density,pressure, temperature etc
    !num_bdry_var_seg = (SPH_dim)*vector_vaiables + scalar variables
    num_bdry_var_seg = (SPH_dim)*1 + 1
    
    !dimilarly define number of boudnary variables consdiered at vertices
    ! vertices contain an additional spatial cordinate values, that is,
    !  location, velocity, density, pressure, temp, etc
    num_bdry_var_ve = (SPH_dim)*(1 + 1) + 1
    
    
    CF_density=mod( ConDivtype, correction_types)
    ID_density=int( ConDivtype/correction_types)
    CF_densDiff= mod(densDiffType,correction_types)
    ID_densDiff= int(densDiffType/correction_types)

    CF_pressure=mod( PrsrGradtype, correction_types)
    ID_pressure=int( PrsrGradtype/correction_types)        
    CF_BIL_visc=mod( BILtype, correction_types)
    ID_BIL_visc=int( BILtype/correction_types)
    
    
     
    if(int(prsrBdryType/correction_types) .ge. 1) prsr_bdry_preCalc = .true.
    ID_prsrBdryType = mod(prsrBdryType, correction_types)
    
    allocate(bdryVal_seg(num_bdry_var_seg,maxedge), bdryVal_ve(num_bdry_var_ve, SPH_dim), vx_ve(SPH_dim,maxnv))
    
    do itimestep = current_ts+1, max_timesteps
        
        !initialize velocity to 0
        vx_ve=0.D0
        
        current_time=dble(itimestep)*dt
        
        ! Add time dependent boundary conditions to the boundary:
        do s =1, etotal
            bdryVal_ve = 0.D0
            do d=1,SPH_dim
                a= edge(d,s)
                bdryVal_ve(1:SPH_dim,d)=x_ve(:,a)
            enddo
            
            call BCinputValue(bdryVal_seg(:,s),num_bdry_var_seg,bdryVal_ve,num_bdry_var_ve,etype(s),SPH_dim,current_time)

            
            do d=1,SPH_dim
                a= edge(d,s)
                vx_ve(:,a)=bdryVal_ve(SPH_dim+1:SPH_dim*2,d)
            enddo
            
            !if(etype(s) .eq. 4) pause
        enddo
        
        call printTimeStep(itimestep,print_step)
        
        ! Calculate all variables and correction factors taht depend on particle position
        call PositionDependentFactors
        
        !if (itimestep .eq. 1) call output_initial_config
        
        ! if there was periodic BC applied, for periodic particles and points update variables used 
        if (Allocated(pBC_edges)) then
            call PeriodicParameterVector(vx,SPH_dim)
            call PeriodicParameterScalar(rho,SPH_dim)
            call PeriodicParameterScalar(p,SPH_dim)
            call PeriodicParameterScalar(mass,SPH_dim)
            call PeriodicParameterScalar(vol,SPH_dim)
        endif
        
        ! Allocate variables necessary
        ! Density updated related paramters
        !allocate(drho(ntotal),rho_prev(ntotal))
        !drho=0.D0
        
        allocate(stress(SPH_dim, ntotal),grad_vel(SPH_dim, SPH_dim, ntotal),  &
            &   grad_vel_s_temp(SPH_dim,SPH_dim), free_surf_particle(ntotal), free_surf_val(ntotal)) ! this can be reduced by accoutnign for nreal and nedge correctly
        
        grad_vel =0.D0
        stress =0.D0
        free_surf_particle=0
        free_surf_val=0.D0
        
        
        
        
        !Find freeSurfaceValues
        !Particle-particle itneraction for finding freeSurfaceValues
        do k= 1,niac
            a= pair_i(k)
            b= pair_j(k)  
            
            !------------------- Find divergence of position (to determine free_surf_particle) -------------------------!
            call CorrectedVecDivPtoP(free_surf_val(a),free_surf_val(b),a,b,x(:,a),x(:,b),dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                    & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                    & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                    & 0,0, SPH_dim, 1, 2) ! SPH_dim, correctionFactorID, divType
            ! -------------------------------------------------------------------------------------------------------------!
        enddo
        
        !Particle-boundary itneraction for finding freeSurfaceValues
        do k= 1,eniac
            a= epair_a(k)
            s= epair_s(k)
            
            !------ Find divergence to calculate freeSurfaceValues-------------!
            F_a(:) = x(:,a)
            F_b(:) = mid_pt_for_edge(:,s)

             call CorrectedVecDivPtoB(free_surf_val(a),a,s,F_a,F_b,del_gamma_as(:,k),  &
                    & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                    & 1,SPH_dim, 1, 2) ! SPH_dim, correctionFactorID, divType
            ! -----------------------------------------------------------------------!
        enddo
        
        
        ! Update freesurface values, as 1 - and the non free surface values are 0
        do a= 1,nreal
            if(free_surf_val(a) .lt. FScutoff) free_surf_particle(a) = 1
        enddo
        
        if(summationDensity) then
            
            allocate(rho_temp(ntotal))
            rho_temp=0.D0
            
            if( .not. allocated(dgrho_prev)) then
                allocate(dgrho_prev(ntotal))
                dgrho_prev=0.D0
            endif
            
            call sumDens(rho_temp)
            
            
             
            if(itimestep .eq. 1) then
                do a = 1,nreal
                    call summationDensityOperatorPtoP( rho_temp(a), a, .true., SumDenstype)
                    
                    rho(a) = rho_temp(a)
                    
                    !Update Volume, since density is updated
                    vol(a) = mass(a)/rho(a)
            
                    ! Update Pressure as it depends on density for WCSPH
                    call ParticlePressureEOS(p(a), rho(a), itype(a), itype_virtual)    
                enddo
            else
                do a = 1,nreal
                    call summationDensityOperatorPtoP( rho_temp(a), a, .false., SumDenstype)
                    
                    rho(a) = rho_temp(a)
                    
                    !Update Volume, since density is updated
                    vol(a) = mass(a)/rho(a)
                    
            
                    ! Update Pressure as it depends on density for WCSPH
                    call ParticlePressureEOS(p(a), rho(a), itype(a), itype_virtual)   
                enddo
            endif
            
            deallocate(rho_temp)
            
        else
            
            allocate(grad_rho(SPH_dim, ntotal),dens_diffusion(ntotal),div_vel(ntotal)) 
            grad_rho=0.D0
            dens_diffusion=0.D0
            div_vel = 0.D0
            
            ! Use all particle-particle interaction to find non boundary terms
            do k= 1,niac
                a= pair_i(k)
                b= pair_j(k)  

                !------------------- Find divergence of velocity (to be used in continuity equation) -------------------------!
                call CorrectedVecDivPtoP(div_vel(a),div_vel(b),a,b,vx(:,a),vx(:,b),dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                        & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                        & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                        & free_surf_particle(a),free_surf_particle(b),SPH_dim, CF_density, ID_density) ! SPH_dim, correctionFactorID, divType
                ! -------------------------------------------------------------------------------------------------------------!
            
                !-------------- Find density Gradient term (to be used in density diffusion equation) --------------!
                call CorrectedScaGradPtoP(grad_rho(:,a), grad_rho(:,b),a,b,rho(a),rho(b),dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                        & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                        & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                        & free_surf_particle(a),free_surf_particle(b),SPH_dim, CF_densDiff, ID_densDiff) ! SPH_dim, correctionFactorID, grad_type
                !-------------------------------------------------------------------------------------------------------------!
                do d=1,SPH_dim
                    if(isNAN(grad_rho(d,a))) then
                        write(*,*) "error calculating in PtoP grad_rho(",d,",",a,")=", grad_rho(d,a)
                        write(*,*) " ++ gamma_mat_inv(:,:,",a,")=", gamma_mat_inv(:,:,a)
                        write(*,*) " ++ gamma_mat(:,:,",a,")=", gamma_mat(:,:,a)
                        write(*,*) " ++ xi1_mat_inv(:,:,",a,")=", xi1_mat_inv(:,:,a)
                        write(*,*) " ++ xi1_mat(:,:,",a,")=", xi1_mat(:,:,a)
                        write(*,*) " ++ dwdx(:,",k,")=", dwdx(:,k)
                
                        pause    
                    endif
                enddo
            
            enddo
        
            ! Use all particle-edge interactions to find boundary terms
            do k= 1,eniac
                a= epair_a(k)
                s= epair_s(k)
            
                !------ Find divergence of velocity (to be used in continuity equation)-------------!
                F_a(:) = vx(:,a)
                F_b(:) = bdryVal_seg(1:SPH_dim,s)

                 call CorrectedVecDivPtoB(div_vel(a),a,s,F_a,F_b,del_gamma_as(:,k),  &
                        & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                        & free_surf_particle(a),SPH_dim, CF_density, ID_density) ! SPH_dim, correctionFactorID, divType
                ! -----------------------------------------------------------------------!
             
                !------ Find density Gradient term (to be used in density diffusion equation) -------------!         
                Sca_Bdry_val = rho(a)
            
                call CorrectedScaGradPtoB(grad_rho(:,a),a,s,rho(a),Sca_Bdry_val,del_gamma_as(:,k),  &
                        & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                        & free_surf_particle(a),SPH_dim, CF_densDiff, ID_densDiff) ! SPH_dim, correctionFactorID, grad_type
                ! -----------------------------------------------------------------------!
                do d=1,SPH_dim
                    if(isNAN(grad_rho(d,a))) then
                        write(*,*) "error calculating in PtoB grad_rho(",d,",",a,")=", grad_rho(d,a)
                        write(*,*) " ++ gamma_mat_inv(:,:,",a,")=", gamma_mat_inv(:,:,a)
                        write(*,*) " ++ gamma_mat(:,:,",a,")=", gamma_mat(:,:,a)
                        write(*,*) " ++ xi1_mat_inv(:,:,",a,")=", xi1_mat_inv(:,:,a)
                        write(*,*) " ++ xi1_mat(:,:,",a,")=", xi1_mat(:,:,a)
                        write(*,*) " ++ del_gamma_as(:,",k,")=", del_gamma_as(:,k)
                
                
                        pause    
                    endif
                enddo
            enddo
        
        
             ! Use all particle-particle interaction to find non boundary terms
            do k= 1,niac
                a= pair_i(k)
                b= pair_j(k)  
            
                delx_ab(:)= x(:,a)- x(:,b)
            
                call densDiffOperatorPtoP(dens_diffusion(a),dens_diffusion(b), &
                    & grad_rho(:,a),grad_rho(:,b),dwdx(:,k),mass(a), mass(b), rho(a), rho(b), &
                    & SPH_dim, hsml_const, c_sound, delx_ab, delta_SPH, ID_densDiff, a,b )  
            
            enddo

            deallocate(grad_rho)
        
            if(itimestep .gt. 1) then
                ! Update variables for the next time step
                do a =1, nreal
                    ! Calcualate density as (𝐷𝜌_𝑎)/𝐷𝑡=− 𝜌_𝑎  ∇∙𝑣_𝑎
                    rho(a) = rho(a) - dt* rho(a) * div_vel(a) + dt*dens_diffusion(a)
            
                    ! Use Hughes density correction if necessary
                    if ((HG_density_correction) .and. (rho(a) .le. rho_init)) rho(a)=rho_init
            
                    ! for free surface impose 𝜌_𝑎= rho_free_surface
                    if(FS_density_correction) rho(a)=dble(1-free_surf_particle(a))*rho(a)+dble(free_surf_particle(a))*rho_init
            
                    !Update Volume, since density is updated
                    vol(a) = mass(a)/rho(a)
            
                    ! Update Pressure as it depends on density for WCSPH
                    call ParticlePressureEOS(p(a), rho(a), itype(a), itype_virtual)    
            
                enddo
            endif
        
            deallocate(dens_diffusion, div_vel)
        endif
        
        ! The varioables need to be updated for periodic particles
        if (Allocated(pBC_edges)) then
            call PeriodicParameterScalar(rho,SPH_dim)
            call PeriodicParameterScalar(p,SPH_dim)
            call PeriodicParameterScalar(vol,SPH_dim)
            call PeriodicParameterScalar(mass,SPH_dim)
        endif
        
        

        !Calcualte density and pressure boundary   
        if (prsr_bdry_preCalc) then
            call FluidVariableAtBdry(num_bdry_var_seg,ID_prsrBdryType)
        else
            call FluidVariableAtBdry(num_bdry_var_seg,0)
        endif
        

        
        
        ! Use all particle-particle interaction to find non boundary terms
        do k= 1,niac
            a= pair_i(k)
            b= pair_j(k)
            
            !-------------- Find Pressure Gradient term (to be used in momentum equation) --------------!
            DF_a=0.D0
            DF_b=0.D0
            call CorrectedScaGradPtoP( DF_a, DF_b,a,b,P(a),P(b),dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                    & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                    & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                    & free_surf_particle(a),free_surf_particle(b),SPH_dim, CF_pressure, ID_pressure) ! SPH_dim, correctionFactorID, grad_type
            stress(:,a) = stress(:,a)+ DF_a(:)
            stress(:,b) = stress(:,b)+ DF_b(:)
            !-------------------------------------------------------------------------------------------------------------!
            
            !-------------- Find velocity gradient term (to be used to find viscous stress in momentum equation) --------------!
            do d =1, SPH_dim
                call CorrectedScaGradPtoP(grad_vel(d,:,a),grad_vel(d,:,b),a,b,vx(d,a),vx(d,b),dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                        & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                        & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                        & free_surf_particle(a),free_surf_particle(b),SPH_dim, CF_BIL_visc, 2) ! SPH_dim, correctionFactorID, grad_type
            enddo
            !-------------------------------------------------------------------------------------------------------------!

                
        enddo
        
        
        
        ! Use all particle-edge interactions to find boundary terms
        do k= 1,eniac
            a= epair_a(k)
            s= epair_s(k)
            
            !------ Find Pressure Gradient term (to be used in momentum equation) -------------!
            DF_a=0.D0
            !gamma_wall_cutoff=0.6D0
            
            if(prsr_bdry_preCalc) then
                Sca_Bdry_val = prsr_bdry_val(s)
            else
                call PressureBdryValue(Sca_Bdry_val,rho(a),x(:,a), vx(:,a), itype(a),bdryVal_seg(:,s), num_bdry_var_seg, a,s, ID_prsrBdryType)
            endif
            
            call CorrectedScaGradPtoB(DF_a,a,s,P(a),Sca_Bdry_val,del_gamma_as(:,k),  &
                    & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                    & free_surf_particle(a),SPH_dim, CF_pressure, ID_pressure) ! SPH_dim, correctionFactorID, grad_type
            
            stress(:,a) = stress(:,a)+ DF_a(:)
            ! -----------------------------------------------------------------------!
            
             !------ Find velocity gradient term (to be used to find viscous stress in momentum equation) ------------
            F_b(:) = bdryVal_seg(1:SPH_dim,s)
            do d = 1, SPH_dim
                call CorrectedScaGradPtoB(grad_vel(d,:,a),a,s,vx(d,a),F_b(d),del_gamma_as(:,k),  &
                        & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                        & free_surf_particle(a),SPH_dim, CF_BIL_visc, 2) ! SPH_dim, correctionFactorID, grad_type
            enddo
            ! -----------------------------------------------------------------------!
            
        enddo
        
        
        
        ! The varioables need to be updated for periodic particles
        if (Allocated(pBC_edges)) then
           call PeriodicParameterTensor(grad_vel,SPH_dim)            
        endif
        
        
        ! Use all particle-particle interaction to find non boundary terms
        do k= 1,niac
            a= pair_i(k)
            b= pair_j(k)
            
            delx_ab(:)= x(:,a)- x(:,b)
            
            !-------------- Find viscous stress term (to be used in momentum equation) --------------!
            DF_a=0.D0
            DF_b=0.D0
            do d= 1,SPH_dim
                call CorrectedBILapPtoP(DF_a(d),DF_b(d),vx(d,a),vx(d,b),dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                        & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                        & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                        & free_surf_particle(a),free_surf_particle(b),SPH_dim, delx_ab(:),ID_BIL_visc)
            enddo
            
            stress(:,a) = stress(:,a)- mu(a)*DF_a(:)
            stress(:,b) = stress(:,b)- mu(b)*DF_b(:)
            
            !-------------------------------------------------------------------------------------------------------------!

            !------------------ Artificial Viscosity term used if applicable -----------------------!
            DF_a=0.D0
            DF_b=0.D0
            call artViscOperatorPtoP(DF_a,DF_b,vx(:,a),vx(:,b),dwdx(:,k),mass(a), mass(b), rho(a), rho(b), &
                    & SPH_dim, hsml_const, c_sound, delx_ab,artViscType )  
            stress(:,a) = stress(:,a)+ DF_a(:)
            stress(:,b) = stress(:,b)+ DF_b(:)
            !-------------------------------------------------------------------------------------------------------------!

            
        enddo
        
        
        ! Use all particle-edge interactions to find boundary terms
        do k= 1,eniac
            a= epair_a(k)
            s= epair_s(k)
            
            ! The below lines find distance between particle to edge mid point            
            delx_ab(:)= x(:,a)- mid_pt_for_edge(:,s) ! this is form is important mostly for Macia formulation
            ! otherwise dot_product(x(:,a), surf_norm(:,s)) is enough mostly
            
            !------ Find viscous stress term (to be used in momentum equation) -------------!           
            F_b(:) = bdryVal_seg(1:SPH_dim,s)
            grad_vel_s_temp(:,:)=0.D0
            
            DF_a=0.D0          
            
            do d= 1,SPH_dim
                
                call BILViscousBdry(F_a,Sca_Bdry_val, dirich0Neum1, grad_vel(d,:,a), grad_vel_s_temp(d,:), vx(d,a), F_b, &
                    & delx_ab, dx_r, surf_norm(:,s), d, SPH_dim, ID_BIL_visc) ! for BIL type using ∇v_s input, v_s needs to be first defined
                
                call CorrectedBILapPtoB(DF_a(d), F_a, Sca_Bdry_val, del_gamma_as(:,k), &
                    & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                    & free_surf_particle(a),SPH_dim, dirich0Neum1,1) ! SPH_dim, correctionFactorID, BIL_type(1 = BIL-PCg, Macia, 2=BIL-NTG), dirich0Neum1
            enddo
            
            stress(:,a) = stress(:,a)- mu(a)*DF_a(:)
            
            ! -----------------------------------------------------------------------!


            
        enddo
        
        
        ! Update real particle position and velocity for the next time step
        do a =1, nreal
            !Update Velocity
            vx(:,a) = vx(:,a) + dt* (-stress(:,a)/rho(a) + F_ext(:))
            
            !Update Position
            x(:,a) = x(:,a) + dt* vx(:,a)
                
            
        enddo
        
        if (Allocated(pBC_edges)) then
            call PeriodicParameterVector(vx,SPH_dim)
        endif
        
        allocate(x_ve_temp(SPH_dim,SPH_dim))
        ! Update edges for next step
        do s = 1,etotal
            if((etype(s) .le. etype_real_max) .and. (etype(s) .gt. etype_real_min)) then
                
                 !Update the vertex positions
                do d=1,SPH_dim
                    a= edge(d,s)
                    !Update position
                    x_ve(:,a) = x_ve(:,a) + dt* vx_ve(:,a)
                    x_ve_temp(:,d)= x_ve(:,a)
                enddo
                
                
                ! now update the integration points and midpoints, that will be used again in the simulation
        
                ! Determine the location of edge mid point
                call centroidBdrySegment(mid_pt_for_edge(:, s), x_ve_temp, SPH_dim)
                ! currently integration poitns are not updated, as those are tnot used in the code for any boundary itnegration,
                ! if it is used,that needs to be updated here.
                
            endif
        enddo
        deallocate(x_ve_temp)
        
        deallocate(grad_vel, stress, grad_vel_s_temp)
        
        
        !---------------------- free surface detection and PST algorithm -------------------------------------!
        !call FreeSurfaceDetection
        call ParticleShiftingTechnique(PSTtype,PSTcoeff, PST_Step, itimestep)
        
        allocate(grad_rho(SPH_dim, ntotal),grad_vel(SPH_dim, SPH_dim, ntotal)) 
        grad_rho =0.D0
        grad_vel =0.D0
        if ((mod(itimestep,PST_step).eq.0) .and. PSTtype .ge. 1) then
            do k= 1,niac
                a= pair_i(k)
                b= pair_j(k)
            
                !-------------- Find density Gradient term --------------!
                    call CorrectedScaGradPtoP(grad_rho(:,a), grad_rho(:,b),a,b,rho(a),rho(b),dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                            & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                            & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                            & free_surf_particle(a),free_surf_particle(b),SPH_dim, CF_densDiff, ID_densDiff) ! SPH_dim, correctionFactorID, grad_type
                
                !-------------- Find velocity gradient term --------------!
                do d =1, SPH_dim
                    call CorrectedScaGradPtoP(grad_vel(d,:,a),grad_vel(d,:,b),a,b,vx(d,a),vx(d,b),dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                            & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                            & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                            & free_surf_particle(a),free_surf_particle(b),SPH_dim, CF_BIL_visc, 2) ! SPH_dim, correctionFactorID, grad_type
                enddo
                !-------------------------------------------------------------------------------------------------------------!
            enddo
        
             do k= 1, eniac
                a=epair_a(k)
                s=epair_s(k)

            
                !------ Find velocity gradient term  ------------
                F_b(:) = bdryVal_seg(1:SPH_dim,s)
                do d = 1, SPH_dim
                    call CorrectedScaGradPtoB(grad_vel(d,:,a),a,s,vx(d,a),F_b(d),del_gamma_as(:,k),  &
                            & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                            & free_surf_particle(a),SPH_dim, 3, 2) ! SPH_dim, correctionFactorID, grad_type
                enddo
            
                !------ Find density Gradient term  -------------!         
                Sca_Bdry_val = rho_s(s)
            
                call CorrectedScaGradPtoB(grad_rho(:,a),a,s,rho(a),Sca_Bdry_val,del_gamma_as(:,k),  &
                        & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                        & free_surf_particle(a),SPH_dim, 5, 2) ! SPH_dim, correctionFactorID, grad_type
            
            enddo  
        
            do a=1,nreal
                rho(a) = rho(a) + dot_product(grad_rho(:,a),x(:,a) - x_prev(:,a))
                do d=1,SPH_dim
                    vx(d,a)= vx(d,a) + dot_product(grad_vel(d,:,a),x(:,a) - x_prev(:,a))
                enddo            
            enddo
            
            deallocate(x_prev)
        endif
        
        deallocate(grad_rho,grad_vel)
        
        !------------------------ export parameter values as output -----------------!
               
        if(time_ev_par_op)call caseBasedOutput(itimestep,dt)
        
        if ((mod(itimestep,save_step).eq.0) .or. (itimestep.eq.1)) call output_flow_simplified(itimestep,dt)   
        
        if( allocated(prsr_bdry_val)) deallocate(prsr_bdry_val)
        if(allocated(rho_s)) deallocate(rho_s)
        
        !-----------------------------------------------------------
        deallocate(delC)
        !deallocate(free_surf_particle)
        
       !call EulerIntegration_fluid(itimestep,dt)
        
        if (Allocated(pBC_edges)) call PeriodicBCreset
        
        ! the below needs to be deallocated here because periodic bc can sometimes increase the total number of particles in next loop
        deallocate(w_aa, gamma_discrt,del_gamma_as, del_gamma, xi1_mat, beta_mat,gamma_mat, &
                   & gamma_mat_inv,xi1_mat_inv , free_surf_particle,free_surf_val)
        
        
        ! if (mod(itimestep,backup_step).eq.0)  call backupOutput(itimestep)
        
    enddo
        
    
deallocate(bdryVal_seg, bdryVal_ve, vx_ve)

if(allocated(dgrho_prev)) deallocate(dgrho_prev)


!Now net time is calculated
call system_clock(count=ic2)
write (*,*)'        Elapsed time = ', (ic2-ic1)/real(crate1), 'sec'

!All allocated variables need to be deallocated so that the memory is freed

if (Allocated(pBC_edges)) deallocate(pBC_epair_a, pBC_epair_s)

DEALLOCATE(x,vx,mass, rho, p, hsml, itype, mu, xStart) 
Deallocate(surf_norm, edge, etype )

if(Allocated(dgrho_prev)) DEALLOCATE(dgrho_prev)
if(Allocated(drho)) DEALLOCATE(drho)
if(Allocated(rho_prev)) DEALLOCATE(rho_prev)
if(Allocated(gamma_dens_cut_off)) Deallocate(gamma_dens_cut_off)

deallocate(temp_vector,dxiac)

!Pause is used so that the terminal is closed only after hitting enter/return on keyboard
write(*,*) 'The code has finished executing, press return to exit'
pause 
end program