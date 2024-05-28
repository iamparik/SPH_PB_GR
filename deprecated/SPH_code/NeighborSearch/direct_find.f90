!****************************************************************************
!
!  SUBROUTINE: direct_find 
!
!  PURPOSE:  Subroutine to calculate the smoothing funciton for each particle and 
!            the interaction parameters used by the SPH algorithm. Interaction
!            pairs are determined by directly comparing the particle distance
!            with the corresponding smoothing length.
!
!   CREATED:        9/13/2017       by  GR Lab
!   Last Modified:  04/18/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine direct_find
    use config_parameter, only:SPH_dim, itype_virtual, hsml_const
    use particle_data ,   only: niac,ntotal,itype,x,hsml,max_interaction, &     
        & pair_i,pair_j,w,dwdx       
    implicit none

!----------------------------------------------------------------------
!                   i,j: A variable used as a counter in loop
!                   d: A variable used to iterate dimensions
!
!----------------------------------------------------------------------

    integer(4) i,j,d   
    real(8) scale_k,mhsml
    
!  The maximum relative distance is determined, after which the kernel will be zero,
! scale_k is not to be confused with kappa which modifies smoothening length
! scale_k is only the factor used in the kernel to approximate various sections
! as 0<R<1, 1<R<2 will all have different curve approximations.
    call sml_mult_factor(scale_k)

 ! We initailize the variable storing the number of interactions 
    niac=0
    
! Now we check every possible pair of particles in the computation domain 
! and determine particle pairs who lie within each others compact kernel   
! The below is only for cosntant smoothing length
    do i=1,ntotal-1
        ! Only particles or points that are used in SPH integration/summation is to be inlcuded
        ! all itypes 0 are classified as reference points that are not used in simulation or for boundary condition
        if (mod(itype(i),itype_virtual) .ne. 0) then 
            do j = i+1, ntotal
                if (mod(itype(j),itype_virtual) .ne. 0) then
                    ! Since we are considering cosntant smoothing length,
                    ! we can use mean smoothing length as either hsml(i) or hsml(j)
                    mhsml = hsml_const !hsml(i)
            
                    call particlePair(i,j,scale_k, mhsml)
                endif
            
            enddo
        endif
   enddo
      

end