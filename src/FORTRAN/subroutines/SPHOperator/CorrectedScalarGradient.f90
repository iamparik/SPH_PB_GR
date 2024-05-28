subroutine CorrectedScaGradPtoP(delF_a,delF_b,F_a,F_b,dwdx, mass_a, mass_b, rho_a, rho_b, &
                    & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, &
                    & gamma_cont_b, gamma_discrt_b, gamma_mat_b, gamma_mat_inv_b, xi1_mat_inv_b, &
                    & dim, CF_ID, grad_type) 
    !scalar Gradient for two particles is calculated, by first correcting the kernel gradient

    implicit none
    integer(4), intent(in) :: dim, CF_ID, grad_type
    real(8), intent(in) :: F_a, F_b, dwdx(dim), mass_a, mass_b, rho_a, rho_b, &
        & gamma_cont_a, gamma_discrt_a, gamma_mat_a(dim,dim), gamma_mat_inv_a(dim,dim), xi1_mat_inv_a(dim,dim), &
        & gamma_cont_b, gamma_discrt_b, gamma_mat_b(dim,dim), gamma_mat_inv_b(dim,dim), xi1_mat_inv_b(dim,dim)
    real(8), intent(inout) :: delF_a(dim),delF_b(dim)
    integer(4) :: d, Scalar0Matrix1
    real(8) :: Cdwdx_a(dim), Cdwdx_b(dim), matrix_factor(dim,dim), scalar_factor
    
    call CorrectionFactorParsing(scalar_factor,matrix_factor,CF_ID,Scalar0Matrix1, &
        & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, dim)            
    Cdwdx_a(:)=dwdx(:)
    call CorrectedKernelGradient(Cdwdx_a, scalar_factor, matrix_factor, Scalar0Matrix1, dim)
            
    call CorrectionFactorParsing(scalar_factor,matrix_factor,CF_ID,Scalar0Matrix1, &
        & gamma_cont_b, gamma_discrt_b, gamma_mat_b, gamma_mat_inv_b, xi1_mat_inv_b, dim)
    Cdwdx_b(:)=-dwdx(:)
    call CorrectedKernelGradient(Cdwdx_b, scalar_factor, matrix_factor, Scalar0Matrix1, dim)    

    call ScalarGradientPtoP(delF_a,delF_b,F_a,F_b,Cdwdx_a, Cdwdx_b, mass_a, mass_b, rho_a, rho_b, dim, grad_type)
    
endsubroutine
    
    
    
subroutine CorrectedScaGradPtoB(delF_a,F_a,F_b,del_gamma_as,&
                    & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, &
                    & dim, CF_ID, grad_type) 

    !vector divergence for particle itneracting with boundary is calculated, by first correctign the kernel gradient

    implicit none
    integer(4), intent(in) :: dim, CF_ID, grad_type
    real(8), intent(in) :: F_a, F_b, del_gamma_as(dim), &
        & gamma_cont_a, gamma_discrt_a, gamma_mat_a(dim,dim), gamma_mat_inv_a(dim,dim), xi1_mat_inv_a(dim,dim)
    real(8), intent(inout) :: delF_a(dim) 
    integer(4) :: d, Scalar0Matrix1
    real(8) :: Cdgmas(dim), matrix_factor(dim,dim), scalar_factor
    
    call CorrectionFactorParsing(scalar_factor,matrix_factor,CF_ID,Scalar0Matrix1, &
                & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, dim)       
    Cdgmas(:)=del_gamma_as(:)
    call CorrectedKernelGradient(Cdgmas, scalar_factor, matrix_factor, Scalar0Matrix1, dim)  
    
    call ScalarGradientPtoB(delF_a,F_a,F_b,Cdgmas, dim, grad_type)
    
endsubroutine