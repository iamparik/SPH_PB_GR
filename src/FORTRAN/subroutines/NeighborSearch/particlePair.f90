!****************************************************************************
!
!  SUBROUTINE: particlePair 
!
!  PURPOSE: This subroutine checks if a pair of particle interact, if they do,
!           functions for the pair are calculated.
!
!   CREATED:        11/27/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  11/27/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************  
    
    subroutine particlePair(a,b,scale_k,mhsml)
    use config_parameter, only:SPH_dim, itype_virtual
    use particle_data ,   only: niac,max_interaction, &     
        & pair_i,pair_j,w,dwdx, x      
     
    implicit none
    integer(4), intent(in):: a,b
    real(8), intent(in):: scale_k, mhsml
    integer(4) d
    real(8)  driac, dxiac(SPH_dim), tdwdx(SPH_dim)

    ! We intiialize the paramter to determine the distance between any two particles
    driac = 0.D0 
            
    do d=1,SPH_dim
    ! Here dx is stored as x(i)-x(j), 
    ! changing this order will change the sign of radial kernel derivatives
        dxiac(d) = x(d,a) - x(d,b)
        driac    = driac + dxiac(d)*dxiac(d)
    enddo
    
    driac=sqrt(driac)
            
    if (driac.le.(scale_k*mhsml))  then     
        if (niac.lt.max_interaction) then
            !     Neighboring pair list, and totalinteraction number and
            !     the interaction number for each particle 
            niac = niac + 1
            pair_i(niac) = a
            pair_j(niac) = b
                    
            !     Kernel and derivations of kernel
            call kernel(driac,dxiac,mhsml,w(niac),tdwdx)
            do d=1,SPH_dim
                dwdx(d,niac) = tdwdx(d)
            enddo
        else
            write(*,*) ' >>> ERROR <<< : Too many interactions'
            write(*,*) ' max_interaction needs to be more than ', max_interaction
            pause
        endif
    endif
    
    end subroutine