!
!  SUBROUTINE: CorrectionFactorMatrixInversion
!			   Gives the inverse of matrices to be used as correctionfactors
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        05/02/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  05/05/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine CorrectionFactorMatrixInversion
    use config_parameter, only: SPH_dim,itype_real_max,itype_real_min, itype_virtual
    use particle_data, only: xi1_mat,xi1_mat_inv,gamma_mat, gamma_mat_inv, ntotal, itype

    implicit none        
   
    integer a
    
    if( .NOT. allocated(gamma_mat_inv)) allocate(gamma_mat_inv(SPH_dim,SPH_dim,ntotal)) 
    if( .NOT. allocated(xi1_mat_inv)) allocate(xi1_mat_inv(SPH_dim,SPH_dim,ntotal))
    xi1_mat_inv=0.D0
    gamma_mat_inv=0.D0
    
    if(SPH_dim .eq. 2) then
        do a=1,ntotal
            !if((mod(itype(a),itype_virtual) .le. itype_real_max) .and. (mod(itype(a),itype_virtual) .ge. itype_real_min)) then
            if((itype(a) .le. itype_real_max) .and. (itype(a) .ge. itype_real_min)) then
                 call inv2d(xi1_mat_inv(:,:,a),xi1_mat(:,:,a))
                 call inv2d(gamma_mat_inv(:,:,a), gamma_mat(:,:,a))
             endif
        enddo
    else
        write(*,*) ' inversion only defined for 2 X 2 matrix'
    endif
    
    
    
    
end
    
    