!****************************************************************************
!
!  SUBROUTINE: MLSnearBoundary
!
!  PURPOSE: This defines the MLS correction factor close to the boundary 
!            
!
!   CREATED:        07/19/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  07/20/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine MLSnearBoundary
   use config_parameter, only: SPH_dim, itype_real_min, &
        & itype_real_max, itype_virtual
    use particle_data, only: x, mass, rho, w, niac,w_aa, itype,  &
        & b_MLS, ntotal, pair_i, pair_j, gamma_cont, gamma_discrt
        
    implicit none
    
    integer(4):: k,a,b 
    real(8),DIMENSION(:,:,:),ALLOCATABLE :: AA_MLS
    real(8),DIMENSION(:,:,:),ALLOCATABLE :: AA_MLS_inv
    

    allocate(AA_MLS(SPH_dim+1, SPH_dim+1, ntotal),AA_MLS_inv(SPH_dim+1, SPH_dim+1, ntotal))
    AA_MLS=0.D0
    AA_MLS_inv=0.D0
    
    
    
    ! Determine the matrix AA_MLS for all particles close to the wall boundary
    do k=1,niac
        
        a=pair_i(k)
        b=pair_j(k)        

        if( gamma_cont(a) .lt. 1.D0) then
            if((mod(itype(a),itype_virtual) .le. itype_real_max) .and. (mod(itype(a),itype_virtual) .ge. itype_real_min)) &
                & call MLS_AA_matrix_linear(AA_MLS(:,:,a),x(:,a),x(:,b),w(k),mass(b),rho(b),SPH_dim)
        endif
        
        if( gamma_cont(b) .lt. 1.D0) then
            if((mod(itype(b),itype_virtual) .le. itype_real_max) .and. (mod(itype(b),itype_virtual) .ge. itype_real_min)) &
                & call MLS_AA_matrix_linear(AA_MLS(:,:,b),x(:,b),x(:,a),w(k),mass(a),rho(a),SPH_dim)        
        endif
        
    enddo
    
    ! UpDATE aa_mls for particles effecr on itself, only Aa(1,1) is effected, again close to bournday only
    do a=1,ntotal
        if( gamma_cont(a) .lt. 1.D0) then
            if((mod(itype(a),itype_virtual) .le. itype_real_max) .and. (mod(itype(a),itype_virtual) .ge. itype_real_min)) &
                & AA_MLS(1,1,a)=AA_MLS(1,1,a) + mass(a)*w_aa(a)/rho(a) 
        endif
        
    enddo
    
    
    allocate(b_MLS(SPH_dim+1,ntotal))
    
    b_MLS=0.D0
    b_MLS(1,:)=1.D0
    
    ! Inverse the matrix AA_MLS for particles near boundary, and update b_MLS vector
    do a=1,ntotal            

        if((mod(itype(a),itype_virtual) .le. itype_real_max) .and. (mod(itype(a),itype_virtual) .ge. itype_real_min)) then
            if( gamma_cont(a) .lt. 1.D0) then
                 call inv3d(AA_MLS_inv(:,:,a),AA_MLS(:,:,a))
                 b_MLS(:,a)=MATMUL(AA_MLS_inv(:,:,a),b_MLS(:,a))
            else
                ! we could define b_MLS(1) as gamma_discrt for internal particles 
            endif
            
        endif
    
    enddo    

        
    
    deallocate(AA_MLS, AA_MLS_inv)

end
    