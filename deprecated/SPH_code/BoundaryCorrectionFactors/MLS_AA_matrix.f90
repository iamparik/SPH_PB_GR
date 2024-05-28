!****************************************************************************
!
!  SUBROUTINE: MLS_AA_matrix_linear
!
!  PURPOSE: This created the Correction matrix Aa for MLS 
!            
!
!   CREATED:        07/20/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  07/20/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine MLS_AA_matrix_linear(Aa, xa, xb, w, massb, rhob, dim)        
    implicit none
    
    integer(4), intent(in) :: dim
    real(8), intent(inout):: Aa(dim+1,dim+1)
    real(8), intent(in)::xa(dim), xb(dim), w, massb, rhob
    integer(4) :: i,j,d
    real(8) :: x_mat, xba(dim+1)

    do d=1,dim+1
        if(d .eq. 1) then            
            xba(d) = 1.D0            
        else
            xba(d) = xb(d-1) - xa(d-1)         
            !xba(d) = xa(d-1) - xb(d-1) 
        endif                
    enddo
    
    do i=1,dim+1
        do j=1,dim+1
            Aa(i,j)= Aa(i,j) + xba(i)*xba(j)*w*massb/rhob
        enddo
    enddo      
    
end
    