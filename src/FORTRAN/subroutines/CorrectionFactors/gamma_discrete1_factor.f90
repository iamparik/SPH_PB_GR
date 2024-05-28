!****************************************************************************
!
!  SUBROUTINE: gamma_discrete1_factor
!
!  PURPOSE: This defines the discrete renormalization factor gamma, usually
!           referred to as Shephards factor in literature   W_ab*m_b/rho_b
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/18/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine gamma_discrete1_factor
    use config_parameter, only: itype_real_max, itype_real_min, itype_virtual
    use particle_data, only: gamma_discrt, niac, pair_i, pair_j, ntotal, &
                & w,itype,mass, rho, w_aa,x, maxn
    implicit none

    integer(4) k, a, b

    if( .NOT. allocated(gamma_discrt)) allocate(gamma_discrt(maxn)) 

   gamma_discrt=0.D0

    do k=1,niac
        a=pair_i(k)
        b=pair_j(k)
        
        ! The below can be accomodated better by changing the algorithm slighty while calculating gamma_discrt  
        gamma_discrt(a) = gamma_discrt(a) + mass(b)*w(k)/rho(b)
        gamma_discrt(b) = gamma_discrt(b) + mass(a)*w(k)/rho(a)
        
    
    enddo

    do a= 1, ntotal
        ! The below can be accomodated better by changing the algorithm slighty while calculating gamma_discrt  
        if((mod(itype(a),itype_virtual) .le. itype_real_max) .and. (mod(itype(a),itype_virtual) .ge. itype_real_min)) &
                & gamma_discrt(a)= gamma_discrt(a)+ mass(a)*w_aa(a)/rho(a) 
    enddo


end
    





