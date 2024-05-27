!****************************************************************************
!
!  SUBROUTINE: 2D_edge_particle_pair 
!
!  PURPOSE:  Subroutine to calculate the smoothing funciton for each particle and 
!            the interaction parameters used by the SPH algorithm. Interaction
!            pairs are determined by directly comparing the particle distance
!            with the corresponding smoothing length.
!
!   CREATED:        04/21/2022       by  PARIKSHIT BOREGOWDA 
!   Last Modified:  04/21/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine edge_particle_pair_2D(khsml,xi,xs1, xs2,ns,interact_edge_particle,error_tol)
    implicit none
    real(8) xi(2),xs1(2), xs2(2), ns(2), khsml , error_tol 
    logical :: interact_edge_particle
    real(8) r_n, dx1(2),dx2(2), r_ab, r_ac, r_bc

    !The normal distance of the particle from the edge is intiaialized to zero
    r_n=0.D0

    ! The vector joining one of the vertice and particle is defined
    dx1(:)= xs1(:)-xi(:)
    dx2(:)= xs2(:)-xi(:)

    !  The dot product of dx.ns gives r_n
    r_n= abs(dot_product(dx1(:),ns(:)))

    r_ab=norm2(xs2(:)-xs1(:))

    r_ac=abs(dot_product(dx1(:),(xs2(:)-xs1(:))/r_ab))  

    r_bc=abs(dot_product(dx2(:),(xs2(:)-xs1(:))/r_ab)) 

    interact_edge_particle= .false.

    ! First we consider only particles within an infinite length rectangular region around the 
    ! edge under consideration, and the wrapping rectangle is twice khsml.
    ! then we restrict the infinite rectangle region, by enclosign with two semi circular regions 
    ! around the end vertices with radius of the circular region equal to khsml
    if(r_n .le. khsml) then
        if((r_ac+r_bc-r_ab) .le. error_tol) then
            interact_edge_particle= .true. 
        elseif((norm2(dx2(:)) .le. khsml) .or. (norm2(dx1(:)) .le. khsml)) then
            interact_edge_particle= .true. 
        endif
    endif
    
end 
    



! Distance between projection



