!****************************************************************************
!
!  SUBROUTINE: sum_density
!
!  PURPOSE:  This subroutine determines density using summation formula   
!
!   CREATED:       -       by  GR lab
!   Last Modified:  01/17/2021       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine summationDensityOperator(rho,OprtType)
    use config_parameter,   only: SPH_dim, itype_real_max, itype_real_min
    use particle_data, only: pair_i,pair_j,niac,ntotal, w, w_aa, &
        & mass, gamma_cont, gamma_discrt, itype, maxn, dgrho_prev, gamma_density_cont,FreeSurfaceVar
    
    implicit none
    integer(4)  a,b,k
    integer(4), intent(in)::OprtType
    real(8) Kbt, bt, rho(maxn) , rho_temp
!    real(8),DIMENSION(:),ALLOCATABLE :: drho   
!    allocate(drho(ntotal))
!    drho=0.D0
    rho=0.D0
    
    Kbt=30000.D0
    
! We calculate 𝜌_a=∑_b 𝑚_b 𝑊_ab = ∑_b≠a 𝑚_b 𝑊_ab + 𝑚_a 𝑊_aa
    
! Add ∑_b≠a 𝑚_b 𝑊_ab to 𝜌_a
    do k =1, niac
        a=pair_i(k)
        b=pair_j(k)        
        call fncnApproxOperator(rho(a),1.D0,mass(b),1.D0,w(k))
        call fncnApproxOperator(rho(b),1.D0,mass(a),1.D0,w(k))        
    enddo
    
! Add 𝑚_a 𝑊_aa to 𝜌_a
    do a=1,ntotal        
        if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
            rho(a)=rho(a)+w_aa(a)*mass(a)     
        endif
    enddo
    

    if(OprtType .eq. 1) then
  ! Correct SPH density everywhere  
        do a=1,ntotal        
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                rho(a)=rho(a)/gamma_discrt(a)    
            endif
        enddo
        
    elseif(OprtType .eq. 2) then
 ! Correct density at the free surface and around wall boundary
 ! with correction factor
        do a=1,ntotal        
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                bt=exp(-Kbt*(min(gamma_discrt(a)/gamma_cont(a),1.D0)-1.D0)**2)            
                rho(a)=rho(a)/(gamma_cont(a)*bt+(1-bt)*gamma_discrt(a))    
            endif
        enddo
        
    elseif(OprtType .eq. 3) then
        
        ! This initializes appropriate density and dgrho_prev for next time step
        ! Here density is initialized as recommended by Ferrand while accounting for correction at free surface
        if( .not. Allocated(dgrho_prev)) then
            ALLOCATE(dgrho_prev(maxn))
            dgrho_prev=0.D0
            ! Correct density at the free surface and around wall boundary
            ! with correction factor
            
            do a=1,ntotal     
                rho_temp= rho(a)
                if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                    bt=exp(-Kbt*(min(gamma_discrt(a)/gamma_cont(a),1.D0)-1.D0)**2)            
                    rho(a)=rho(a)/(gamma_cont(a)*bt+(1-bt)*gamma_discrt(a))    
                endif
                dgrho_prev(a)=rho(a)*gamma_cont(a) - rho_temp
            enddo
        else
            do a=1, ntotal
                rho_temp= rho(a)
                rho(a)= (rho(a) + dgrho_prev(a))/ gamma_cont(a)              
                dgrho_prev(a)=rho(a)*gamma_cont(a) - rho_temp
            enddo
            
        endif
        
    elseif(OprtType .eq. 4) then
        
        ! This initializes appropriate density and dgrho_prev for next time step
        ! Here density is initialized as recommended by Ferrand while accounting for correction at free surface
        if( .not. Allocated(dgrho_prev)) then
            ALLOCATE(dgrho_prev(maxn))
            dgrho_prev=0.D0
            ! Correct density at the free surface and around wall boundary
            ! with correction factor
            
            do a=1,ntotal     
                !rho_temp= rho(a)
                if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                    bt=exp(-Kbt*(min(gamma_discrt(a)/gamma_cont(a),1.D0)-1.D0)**2)            
                    rho_temp=rho(a)/(gamma_cont(a)*bt+(1-bt)*gamma_discrt(a))    
                endif
                dgrho_prev(a)=rho_temp*gamma_cont(a) - rho(a)                
            enddo
        else
            do a=1, ntotal
                rho_temp= rho(a)
                rho(a)= (rho(a) + dgrho_prev(a))/ gamma_cont(a)              
                dgrho_prev(a)=rho(a)*gamma_cont(a) - rho_temp
            enddo
            
        endif
        
    elseif(OprtType .eq. 5) then
        call gamma_density_continuous_leroy
        call FreeSurfaceDetection
        do a=1,ntotal        
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                rho(a)=rho(a)/(gamma_discrt(a)*gamma_density_cont(a))    
                 if( FreeSurfaceVar(a) .lt. 0.D0) then 
                    rho(a)=1000.D0
                endif
            endif
        enddo
        
        deallocate(FreeSurfaceVar)
        
    elseif(OprtType .eq. 6) then
        call gamma_density_continuous_leroy
        call FreeSurfaceDetection
        do a=1,ntotal        
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then  
                rho(a)=rho(a)/(gamma_cont(a))
                 if( FreeSurfaceVar(a) .lt. 0.D0) then 
                    rho(a)=max(1000.D0, rho(a)/(gamma_discrt(a)))
                 endif
                 rho(a)=rho(a)/gamma_density_cont(a)
            endif
        enddo
        
        deallocate(FreeSurfaceVar)
    endif

    
end