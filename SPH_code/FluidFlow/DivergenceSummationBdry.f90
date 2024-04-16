!
!  SUBROUTINE: DivergenceSummationBdry
!			    This calculates divergence of velocity using continuity eqn away from wall  
!               and uses summation eqn near wall
!               This only works with the Euler scheme currentyl
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        06/08/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/04/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine DivergenceSummationBdry( d1f, f0, rho_in)
    use config_parameter,   only: SPH_dim, itype_real_max, itype_real_min
    use config_geometry, only:SumDenstype
    use particle_data, only: pair_i,pair_j,niac,ntotal, w, w_aa, &
        & dwdx, &
        & mass, gamma_cont, gamma_discrt, itype, maxn, dgrho_prev

    implicit none
    real(8), intent(in):: f0(SPH_dim,ntotal)
    real(8), intent(inout):: d1f(ntotal), rho_in(ntotal)
    integer(4)::k,a,b, d
    real(8) ::  delf,CdwdxA(SPH_dim),CdwdxB(SPH_dim)
    real(8) ::  rho_temp, Kbt, bt
    logical:: firstStep
    real(8),DIMENSION(:),ALLOCATABLE :: rho_sum   
    
    
    firstStep= .false.
    allocate(rho_sum(ntotal))    
    rho_sum=0.D0
    
    Kbt=30000.D0
          
    ! Use continuity equation to update drho for all particles away from wall    
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)

        CdwdxA= dwdx(:,k)
        CdwdxB= -dwdx(:,k)
        
        
        if(gamma_cont(a) .eq. 1.D0) then
            do d= 1,SPH_dim
                delf= f0(d,b)-f0(d,a)
                call fncnApproxOperator(d1f(a),delf,mass(b),rho_in(b),CdwdxA(d))
            enddo
        endif
        
        if(gamma_cont(b) .eq. 1.D0) then
            do d= 1,SPH_dim
                delf= f0(d,a)-f0(d,b)
                call fncnApproxOperator(d1f(b),delf,mass(a),rho_in(a),CdwdxB(d))
            enddo                    
        endif
       
        call fncnApproxOperator(rho_sum(a),1.D0,mass(b),1.D0,w(k))
        call fncnApproxOperator(rho_sum(b),1.D0,mass(a),1.D0,w(k))
        
     enddo
     
    if( .not. Allocated(dgrho_prev)) then
        ALLOCATE(dgrho_prev(maxn))                
        dgrho_prev=0.D0
        firstStep=.true.
    endif
     
    
    do a=1,ntotal        
        
        if((gamma_cont(a) .lt. 1.D0) .and. (gamma_cont(a) .gt. 0.D0)) then
            
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                
                ! Add 𝑚_a 𝑊_aa to 𝜌_a for particles close to wall
                rho_sum(a)=rho_sum(a)+w_aa(a)*mass(a)     
                  
                if(SumDenstype .eq. 1) then
                ! Correct SPH density with gamma_discrt
                    rho_sum(a)=rho_sum(a)/gamma_discrt(a)
        
                elseif(SumDenstype .eq. 2) then
                ! Correct density at the free surface + wall boudanry 
                ! differently form just wall boundary
                    bt=exp(-Kbt*(min(gamma_discrt(a)/gamma_cont(a),1.D0)-1.D0)**2)            
                    rho_sum(a)=rho_sum(a)/(gamma_cont(a)*bt+(1-bt)*gamma_discrt(a))    

                elseif(SumDenstype .eq. 3) then
        
                ! This initializes appropriate density and dgrho_prev for next time step
                ! Here density is initialized as recommended by Ferrand while accounting for correction at free surface
                
                    if(firstStep) then
                        ! Correct density at the free surface and around wall boundary
                        ! with correction factor
                        rho_temp= rho_sum(a)
                        bt=exp(-Kbt*(min(gamma_discrt(a)/gamma_cont(a),1.D0)-1.D0)**2)            
                        rho_sum(a)=rho_sum(a)/(gamma_cont(a)*bt+(1-bt)*gamma_discrt(a))  
                        dgrho_prev(a)=rho_sum(a)*gamma_cont(a) - rho_temp

                    else
                        rho_temp= rho_sum(a)
                        rho_sum(a)= (rho_sum(a) + dgrho_prev(a))/ gamma_cont(a)              
                        dgrho_prev(a)=rho_sum(a)*gamma_cont(a) - rho_temp            
                    endif
                    
                endif 
                
                !Assign summation density value, with correction as density and assign drho as 0
                rho_in(a)=rho_sum(a)
                d1f(a)=0.D0
        
            endif
            
        elseif(gamma_cont(a) .eq. 1.D0) then
            if(SumDenstype .eq. 3) then
                ! Add 𝑚_a 𝑊_aa to 𝜌_a for all particles to account for dghro_prev of particles 
                ! coming from free surface to near wall
                rho_sum(a)=rho_sum(a)+w_aa(a)*mass(a)  
                
                rho_temp= rho_sum(a)
                bt=exp(-Kbt*(min(gamma_discrt(a)/gamma_cont(a),1.D0)-1.D0)**2)            
                rho_sum(a)=rho_sum(a)/(gamma_cont(a)*bt+(1-bt)*gamma_discrt(a))  
                dgrho_prev(a)=rho_sum(a)*gamma_cont(a) - rho_temp
                
            endif
            
        endif
    enddo
    
    

     
    end
    
    
    
    
    