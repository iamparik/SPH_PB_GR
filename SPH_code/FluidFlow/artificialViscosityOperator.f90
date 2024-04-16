!****************************************************************************
!
!  SUBROUTINE: artificialViscosityOperator
!
!  PURPOSE:  Subroutine to apply artificial viscosity to flow   
!
!   CREATED:        06/26/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  06/26/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine artificialViscosityOperator(dstress, vx, x, rho, oprtrType)

    use config_parameter,   only:SPH_dim
    use config_geometry,    only: hsml_const,c_sound
    use particle_data,      only: pair_i,pair_j,niac,ntotal,itype, dwdx, &
        & mass
    
    
    implicit none
      
    real(8) ::   dstress(SPH_dim,ntotal), vx(SPH_dim,ntotal), x(SPH_dim,ntotal), rho(ntotal)
    integer(4) :: oprtrType
    !real(8),DIMENSION(:),ALLOCATABLE :: 
    real(8),DIMENSION(:,:),ALLOCATABLE :: dfncn
    integer(4) :: a,b,d, k     
    real(8) :: MU_ab, PI_ab,eta, alpha, beta
     real(8) :: CdwdxA(SPH_dim),CdwdxB(SPH_dim)
    
    ALLOCATE(dfncn(SPH_dim,ntotal))
    dfncn=0.D0

    if(oprtrType .eq. 1) then
        !---------------------Monoghan artificial viscosity term -------------------------
        
        eta=0.1*hsml_const
        
        alpha = 1.D0
        beta = 2.D0
        
        do k = 1,niac
            a=pair_i(k)
            b=pair_j(k)
            
            if((mass(a).gt. 0.D0) .and. mass(b) .gt. 0.D0) then
            
                MU_ab= hsml_const*dot_product((vx(:,a)-vx(:,b)),(x(:,a)-x(:,b)))/(dot_product((x(:,a)-x(:,b)),(x(:,a)-x(:,b)))+eta**2)
            
                if(MU_ab .lt. 0.D0) then
            
                    PI_ab= (-alpha*c_sound*MU_ab+ beta* MU_ab**2)/((rho(a)+rho(b))/2.D0)

                    CdwdxA= dwdx(:,k)
    
                    CdwdxB= -dwdx(:,k)
        
                    do d= 1,SPH_dim
                        call fncnApproxOperator(dfncn(d,a),PI_ab,mass(b),1.D0,CdwdxA(d)*rho(a))
                        
                        call fncnApproxOperator(dfncn(d,b),PI_ab,mass(a),1.D0,CdwdxB(d)*rho(b))
                        
                    enddo     
                endif
                
            endif
        
        enddo
    
    endif
    
    do d=1,SPH_dim
        dstress(d,1:ntotal)= dstress(d,1:ntotal) - dfncn(d,1:ntotal)
    enddo
    
    
    
end
    