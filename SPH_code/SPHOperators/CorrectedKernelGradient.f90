subroutine CorrectedKernelGradient(Cdwdx, scalar_factor, matrix_factor, Scalar0Matrix1, dim)
    implicit None
    integer(4), intent(in) :: dim, Scalar0Matrix1
    real(8), intent(in) :: scalar_factor, matrix_factor(dim,dim)
    real(8), intent(inout) :: Cdwdx(dim)
    integer(4) :: d
    real(8) :: GmaInvtemp(dim,dim)
    
    if(Scalar0Matrix1 .eq. 0) then
        do d=1,dim
            GmaInvtemp(d,d)= 1.D0/scalar_factor
        enddo            
        
        Cdwdx= MATMUL(GmaInvtemp,Cdwdx)
        
    else
        Cdwdx= MATMUL(matrix_factor,Cdwdx)
    endif
    
end 

    
    
            
    