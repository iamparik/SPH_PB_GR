!****************************************************************************
!
!  SUBROUTINE: gamma_continuous_Ferrand
!
!  PURPOSE: This defines the continuous wall renormalization factor gamma in FErrand et.al paper
!
!   CREATED:        05/24/2021       by  PARIKSHIT BOREGOWDA
!   Last Modified:  10/14/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine gamma_continuous_Ferrand(dt)

    use config_parameter, only: SPH_dim, pi, itype_real_max, itype_real_min
    use particle_data, only:gamma_cont, eniac, epair_a, epair_s, &
                & x,x_prev,vx, surf_norm, gamma_cont_prev, hsml, del_gamma, &
                & del_gamma_as_prev, del_gamma_as,  &
                & ntotal, edge, etype, etotal, nedge_rel_edge, &
                & maxn, itype
    implicit none

    integer(4) k,a,s, steps_moved,v, nn 
    real(8) xa(SPH_dim), xs_v(SPH_dim,size(edge,1)), dx_r(SPH_dim)
    real(8) length_moved, scale_k, dt, khsml
    real(8), dimension(:,:),allocatable::dgamma_as_prev, dgamma_as
    logical :: interact_edge_particle
    

 

    call sml_mult_factor(scale_k)
    
    ! If gamma_cont is not calculated before, then it needs to be calculated through an iteration
    ! and is executed via the code in the first if statement below,
    ! in timesteps after 1st, where gamma_cont is calculated, the else statement is executed 
    if( .NOT. allocated(gamma_cont)) then
        allocate(gamma_cont(maxn)) 
        gamma_cont=0.D0
    
        !Here previous gamma is initiatied as 1, as in Ferrand's
        !algorithm the particles close to boundary are moved back into the domain from where 
        ! it is iterated to its originial position
        if( .NOT. allocated(gamma_cont_prev)) allocate(gamma_cont_prev(maxn))
        gamma_cont_prev=1.D0
    
        allocate(dgamma_as(SPH_dim,eniac), dgamma_as_prev(SPH_dim,eniac))
        dgamma_as=0.D0
        dgamma_as_prev=0.D0
    
        ! The particle is moved from region gamma is 1, to a region where gamma
        ! is to be determined (real, current position) in many steps, 
        ! the # of steps is defined by steps_moved
        steps_moved= 2000
        
        ! calcualte 𝛾_𝑎^(𝑛+1)=𝛾_𝑎^𝑛+ 1/2 ∑1_𝑠▒〖(∇^𝑛 𝛾_𝑎𝑠+∇^(𝑛+1) 𝛾_𝑎𝑠 )∙(𝑟_𝑎^𝑛+𝑟_𝑎^(𝑛+1) ) 〗
        !
    
        do nn=1,steps_moved
            gamma_cont=0.D0 
            do k=1,eniac
                a=epair_a(k)
                s=epair_s(k)
                
                
                !The real particle next to edge is moved back by twice the radius of the 
                !compact support
                length_moved= 2.D0*hsml(a)*scale_k 
                
                khsml= hsml(a)*scale_k
        
                dx_r(:)= length_moved*(del_gamma(:,a)/norm2(del_gamma(:,a)))/steps_moved
        
                xa(:)= x(:,a)+ (steps_moved-nn)*dx_r(:)
        
                if (SPH_dim .eq. 2) then
			        do v=1,size(edge,1)
                        xs_v(:,v)=x(:,edge(v,s))    
                    enddo
        
                    call edge_particle_pair_2D(khsml,xa,xs_v(:,1),xs_v(:,2),surf_norm(:,s),interact_edge_particle, dble(5.D-6))
                    if(interact_edge_particle) call dgamma_analytical_Ferrand_2D(xa,xs_v,surf_norm(:,s), dgamma_as(:,k),hsml(a), pi, scale_k)                                                
                    
                else
                
                endif
				
                ! (𝑑𝛾_𝑎)/𝑑𝑡= 1/2 ∑1_𝑠▒〖(∇^𝑛 𝛾_𝑎𝑠+∇^(𝑛+1) 𝛾_𝑎𝑠 )∙(𝑟_𝑎^𝑛+𝑟_𝑎^(𝑛+1) ) 〗
                gamma_cont(a)= gamma_cont(a) + dot_product(0.5D0*(dgamma_as(:,k)+dgamma_as_prev(:,k)), -dx_r(:))
                
                dgamma_as_prev(:,k)=dgamma_as(:,k)
                    
                
        
            enddo
            
            ! 𝛾_𝑎^(𝑛+1)=(𝑑𝛾_𝑎)/𝑑𝑡 + 𝛾_𝑎^𝑛
            gamma_cont(:)=gamma_cont(:) + gamma_cont_prev(:)
            
            !  store 𝛾_𝑎^(𝑛+1) as 𝛾_𝑎^n for next step
            gamma_cont_prev(:)= gamma_cont(:)
            
        enddo
        
        deallocate(dgamma_as, dgamma_as_prev)
        
    
    else 
    !If timestep is not 1, then gamma value from previous step can be directly used,
    ! to calculate gamma, defined as 
    ! ɣ(n+1)=ɣ(n) + dt/2* Σ(∇ɣ_as(n+1) + ∇ɣ_as(n)). u_a_Rs(n+1)
        do k=1,eniac
            a=epair_a(k)
            s=epair_s(k)
            
            
            !  ∇ɣ_as(n) needs to be calcualted, and sincek k (pair) values
            !  might be different for previous loop and this loop 
            if (SPH_dim .eq. 2) then
			    do v=1,size(edge,1)
                    xs_v(:,v)=x(:,edge(v,s))    
                enddo
                xa(:)=x_prev(:,a)
                call edge_particle_pair_2D(khsml,xa,xs_v(:,1),xs_v(:,2),surf_norm(:,s),interact_edge_particle, dble(5.D-6))
                if(interact_edge_particle) call dgamma_analytical_Ferrand_2D(xa,xs_v,surf_norm(:,s), del_gamma_as_prev(:),hsml(a), pi, scale_k)                                                
                    
            else
                
            endif
            
            if((itype(a) .ge. itype_real_min) .and. (itype(a) .le. itype_real_max)) then
            !Here u_a_Rs(n+1) is just used as vx(:,a), but it needs to be (vx(:,a) -vx(:,nedge_rel_edge(s)))
                gamma_cont_prev(a)=gamma_cont_prev(a)+ &
                    & dt*dot_product(del_gamma_as(:,k)+del_gamma_as_prev(:),vx(:,a))/2.D0    
            endif
            
        enddo
        
        do a=1,ntotal
            gamma_cont(a)=gamma_cont_prev(a)
        enddo                  
        
        !! Now assign all  gamma_cont away from boundary  as 1.DO
        !! and for those under boudnary's influence, update with gamma_cont_prev
        !gamma_cont(:)=1.D0
        !
        !do s=1,etotal
        !    a=nedge_rel_edge(s)
        !    gamma_cont(a)=gamma_cont_prev(a)        
        !enddo
        !
        !! now assign current gamma to gamma_prev to be used in next loop
        !gamma_cont_prev=gamma_cont
        
    endif
  
  ! ensure  ∇ɣ_as(n) is allocated, to ensure, x_prev is allocated and updated in the 
  ! time step
  if( .NOT. allocated(del_gamma_as_prev)) allocate(del_gamma_as_prev(SPH_dim))
   
    
end
    
    
    
    
