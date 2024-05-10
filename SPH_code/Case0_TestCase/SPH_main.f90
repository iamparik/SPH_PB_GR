program SPH3D
 
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
real(8) :: scalar_factor, Sca_Bdry_val
real(8), DIMENSION(:), allocatable  :: F_a, F_b, Cdwdx_a, Cdwdx_b, Cdgmas
real(8), DIMENSION(:,:,:), allocatable :: grad_vel
real(8), DIMENSION(:,:), allocatable :: matrix_factor, grad_P, visc_stress, x_ve_temp
real(8), DIMENSION(:), allocatable :: div_vel, delx_ab


   
!start system clock to evealuate time taken to run
! System_clock returns number of seconds from 00:00 CUT on 1 Jan 1970
    call system_clock(count=ic1, count_rate=crate1, count_max=cmax1)

! input particle configuration data
    call inputSPHConfig
    
    allocate(F_a(SPH_dim), F_b(SPH_dim), Cdwdx_a(SPH_dim), Cdwdx_b(SPH_dim), Cdgmas(SPH_dim), &
    & matrix_factor(SPH_dim,SPH_dim),  delx_ab(SPH_dim))

! Add boundary conditions to the boundary:
    allocate( bdryVal_vel(SPH_dim,maxedge),bdryVal_prs(maxedge), bdryVal_rho(maxedge), bdryVal_temp(maxedge))
    do s =1, etotal
        call BCinputValue(etype(s),bdryVal_vel(:,s), bdryVal_prs(s), bdryVal_rho(s), bdryVal_temp(s))
    enddo

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
    
    do itimestep = current_ts+1, max_timesteps
        
        
        call printTimeStep(itimestep,print_step)
        
        call PositionDependentFactors
        
        if (itimestep .eq. 1) call output_initial_config
        
        ! if there was periodic BC applied, for periodic particles and points update variables used 
        if (Allocated(pBC_edges)) then
            call PeriodicParameter(vx(1,:))
            call PeriodicParameter(vx(2,:))
            call PeriodicParameter(rho)
            call PeriodicParameter(p)
        endif
        
        ! Allocate variables necessary
        ! Density updated related paramters
        !allocate(drho(ntotal),rho_prev(ntotal))
        !drho=0.D0
        
        allocate(div_vel(ntotal), grad_P(SPH_dim, ntotal), visc_stress(SPH_dim, ntotal), &
            & grad_vel(SPH_dim, SPH_dim, ntotal), x_ve_temp(SPH_dim,SPH_dim)) ! this can be reduced by accoutnign for nreal and nedge correctly
        div_vel=0.D0
        grad_P=0.D0
        grad_vel =0.D0
        visc_stress =0.D0
        x_ve_temp=0.D0
        
        
        CF_density=mod( ConDivtype, correction_types)
        ID_density=int( ConDivtype/correction_types)
        
        ! Use all particle-particle interaction to find non boundary terms
        do k= 1,niac
            a= pair_i(k)
            b= pair_j(k)  

           
            !------------------- Find divergence of velocity (to be used in continuity equation) -------------------------!
            call CorrectedVecDivPtoP(div_vel(a),div_vel(b),vx(:,a),vx(:,b),dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                    & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                    & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                    & SPH_dim, CF_div, ID_div) ! SPH_dim, correctionFactorID, divType
            ! -------------------------------------------------------------------------------------------------------------!
            
        enddo
        
        ! Use all particle-edge interactions to find boundary terms
        do k= 1,eniac
            a= epair_a(k)
            s= epair_s(k)
            
            
            !------ Find divergence of velocity for calculating density -------------!
            F_a(:) = vx(:,a)
            F_b(:) = bdryVal_vel(:,s)

             call CorrectedVecDivPtoB(div_vel(a),F_a,F_b,del_gamma_as(:,k),  &
                    & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                    & SPH_dim, CF_div, ID_div) ! SPH_dim, correctionFactorID, divType
            ! -----------------------------------------------------------------------!
        enddo
        
        
        ! Update variables for the next time step
        do a =1, ntotal
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                
                ! Calcualate density as (𝐷𝜌_𝑎)/𝐷𝑡=− 𝜌_𝑎  ∇∙𝑣_𝑎
                rho(a) = rho(a) - dt* rho(a) * div_vel(a)
                ! Use Hughes density correction if necessary
                if ((HG_density_correction) .and. (rho(a) .le. rho_init)) rho(a)=rho_init
            
                ! Update Pressure as it depends on density for WCSPH
                call ParticlePressureEOS(p(a), rho(a), itype(a), itype_virtual)               
                
            
            endif
        enddo
        
        ! The varioables need to be updated for periodic particles
        if (Allocated(pBC_edges)) then
            call PeriodicParameter(rho)
            call PeriodicParameter(p)
        endif
        
        
        CF_pressure=mod( PrsrGradtype, correction_types)
        ID_pressure=int( PrsrGradtype/correction_types)
        
        CF_BIL_visc=mod( BILtype, correction_types)
        ID_BIL_visc=int( BILtype/correction_types)
        
        ! Use all particle-particle interaction to find non boundary terms
        do k= 1,niac
            a= pair_i(k)
            b= pair_j(k)
            
            !-------------- Find Pressure Gradient term (to be used in momentum equation) --------------!
            
            call CorrectedScaGradPtoP(grad_P(:,a),grad_P(:,b),P(a),P(b),dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                    & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                    & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                    & SPH_dim, CF_pressure, ID_pressure) ! SPH_dim, correctionFactorID, grad_type
            !-------------------------------------------------------------------------------------------------------------!
            
            !-------------- Find velocity gradient term (to be used to find viscous stress in momentum equation) --------------!
            do d =1, SPH_dim
                call CorrectedScaGradPtoP(grad_vel(d,:,a),grad_vel(d,:,b),vx(d,a),vx(d,b),dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                        & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                        & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                        & SPH_dim, CF_BIL_visc, 2) ! SPH_dim, correctionFactorID, grad_type
            enddo
            !-------------------------------------------------------------------------------------------------------------!

        enddo
        
        
        ! Use all particle-edge interactions to find boundary terms
        do k= 1,eniac
            a= epair_a(k)
            s= epair_s(k)
            
            !------ Find Pressure Gradient term (to be used in momentum equation) -------------!
            Sca_Bdry_val = P(a) -rho(a)*c_sound*dot_product(vx(:,a)-bdryVal_vel(:,s), surf_norm(:,s))
            
            call CorrectedScaGradPtoB(grad_P(:,a),P(a),Sca_Bdry_val,del_gamma_as(:,k),  &
                    & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                    & SPH_dim, CF_pressure, ID_pressure) ! SPH_dim, correctionFactorID, grad_type
            ! -----------------------------------------------------------------------!
            
             !------ Find velocity gradient term (to be used to find viscous stress in momentum equation) ------------
            F_b(:) = bdryVal_vel(:,s)
            do d = 1, SPH_dim
                call CorrectedScaGradPtoB(grad_vel(d,:,a),vx(d,a),F_b(d),del_gamma_as(:,k),  &
                        & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                        & SPH_dim, CF_BIL_visc, 2) ! SPH_dim, correctionFactorID, grad_type
            enddo
            ! -----------------------------------------------------------------------!
            
        enddo
        
        ! The varioables need to be updated for periodic particles
        if (Allocated(pBC_edges)) then
            do i =1, SPH_dim
                do j =1,SPH_dim
                    call PeriodicParameter(grad_vel(i,j,:))
                enddo
            enddo
            
        endif
        
        
        ! Use all particle-particle interaction to find non boundary terms
        do k= 1,niac
            a= pair_i(k)
            b= pair_j(k)
            
            !-------------- Find viscous stress term (to be used in momentum equation) --------------!
            delx_ab(:)= x(:,a)- x(:,b)
            do d= 1,SPH_dim
                call CorrectedBILapPtoP(visc_stress(d,a),visc_stress(d,b),vx(d,a),vx(d,b),dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                        & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                        & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                        & SPH_dim, delx_ab(:))
            enddo
            !-------------------------------------------------------------------------------------------------------------!


        enddo
        
        
        ! Use all particle-edge interactions to find boundary terms
        do k= 1,eniac
            a= epair_a(k)
            s= epair_s(k)
            b= nedge_rel_edge(s)
            
            !------ Find viscous stress term (to be used in momentum equation) -------------!
            delx_ab(:)= x(:,a)- x(:,b)
            
            do d= 1,SPH_dim
                
                F_a(:) = 2.D0*grad_vel(d,:,a) !+ grad_vel(d,:,b) !0.D0
                !Sca_Bdry_val=2.D0*(vx(d,a)-bdryVal_vel(d,s))/dot_product(delx_ab, surf_norm(:,s)) !0.D0
                call BILViscousBdry(F_a,Sca_Bdry_val, grad_vel(d,:,a), grad_vel(d,:,b), vx(d,a), bdryVal_vel(d,s), delx_ab, &
                    & surf_norm(:,s), ID_BIL_visc)
                
                call CorrectedBILapPtoB(visc_stress(d,a), F_a, Sca_Bdry_val, del_gamma_as(:,k), &
                    & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                    & SPH_dim) ! SPH_dim, correctionFactorID, BIL_type(1 = BIL-PCg, Macia, 2=BIL-NTG), dirich0Neum1
            enddo
            ! -----------------------------------------------------------------------!


            
        enddo
        
        
        ! Update variables for the next time step
        do a =1, ntotal
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                

                !Update Velocity
                vx(:,a) = vx(:,a) + dt* (-grad_P(:,a)/rho(a) + (mu(a)/rho(a))*visc_stress(:,a) + F_ext(:))
            

                !Update position
                x(:,a) = x(:,a) + dt* vx(:,a)
                
            
            endif
        enddo
        
        
        ! Update edges for next step
        do s = 1,etotal
            if((etype(s) .le. etype_real_max) .and. (etype(s) .gt. etype_real_min)) then
                
                 !Update the vertex positions
                do d=1,SPH_dim
                    a= edge(d,s)
                    !Update position
                    x_ve(:,a) = x_ve(:,a) + dt* bdryVal_vel(:,s)
                enddo
                
                ! now update the midpoint of edge, or poitns,
                ! which are used in numerical integrations
                k=nedge_rel_edge(s)        
                do d =1,SPH_dim
                    x_ve_temp(:,d)=x_ve(:,edge(d,s))
                enddo        
                call centroidBdrySegment(x(:,k), x_ve_temp, SPH_dim)

            
            endif
        enddo
        
         deallocate(div_vel,grad_P, grad_vel, visc_stress, x_ve_temp)
        
        
        !---------------------- free surface detection and PST algorithm -------------------------------------!
        call FreeSurfaceDetection
        if(PSTtype .gt. 0) then
            call ParticleShiftingTechnique
        endif
        
        if( .not. allocated(delC)) call ConcGradient  
        
        
        if ((mod(itimestep,save_step).eq.0) .or. (itimestep.eq.1)) call output_flow_simplified(itimestep,dt)   
        deallocate(delC)
        deallocate(FreeSurfaceVar)
        
       !call EulerIntegration_fluid(itimestep,dt)
        
        if (Allocated(pBC_edges)) call PeriodicBCreset
        
        ! the below needs to be deallocated here because periodic bc can sometimes increase the total number of particles in next loop
        deallocate(w_aa, gamma_discrt,del_gamma_as, del_gamma, xi1_mat, beta_mat,gamma_mat,xi_cont_mat, &
                   & gamma_mat_inv,xi1_mat_inv,xi_cont_mat_inv,b_MLS )
        
        
        ! if (mod(itimestep,backup_step).eq.0)  call backupOutput(itimestep)
        
    enddo
        




!Now net time is calculated
call system_clock(count=ic2)
write (*,*)'        Elapsed time = ', (ic2-ic1)/real(crate1), 'sec'

!All allocated variables need to be deallocated so that the memory is freed

if (Allocated(pBC_edges)) deallocate(pBC_epair_a, pBC_epair_s)

DEALLOCATE(x,vx,mass, rho, p, hsml, itype, mu, temp, xStart) 
Deallocate(surf_norm, edge,nedge_rel_edge, etype )

if(Allocated(dgrho_prev)) DEALLOCATE(dgrho_prev)
if(Allocated(drho)) DEALLOCATE(drho)
if(Allocated(rho_prev)) DEALLOCATE(rho_prev)

!Pause is used so that the terminal is closed only after hitting enter/return on keyboard
pause 
end program