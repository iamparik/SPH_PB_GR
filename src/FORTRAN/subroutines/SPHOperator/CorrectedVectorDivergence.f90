subroutine CorrectedVecDivPtoP(divF_a,divF_b,F_a,F_b,dwdx, mass_a, mass_b, rho_a, rho_b, &
                    & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, &
                    & gamma_cont_b, gamma_discrt_b, gamma_mat_b, gamma_mat_inv_b, xi1_mat_inv_b, &
                    & FS_a,FS_b, dim, CF_ID, div_type) 
    !vector divergence for two particles is calculated, by first correctign the kernel gradient

    implicit none
    integer(4), intent(in) :: dim, CF_ID, div_type
    real(8), intent(in) :: F_a(dim), F_b(dim), dwdx(dim), mass_a, mass_b, rho_a, rho_b, &
        & gamma_cont_a, gamma_discrt_a, gamma_mat_a(dim,dim), gamma_mat_inv_a(dim,dim), xi1_mat_inv_a(dim,dim), &
        & gamma_cont_b, gamma_discrt_b, gamma_mat_b(dim,dim), gamma_mat_inv_b(dim,dim), xi1_mat_inv_b(dim,dim)
    real(8), intent(inout) :: divF_a,divF_b
    integer(4) :: d, Scalar0Matrix1
    integer(2) :: FS_a,FS_b
    real(8) :: Cdwdx_a(dim), Cdwdx_b(dim), matrix_factor(dim,dim), scalar_factor
    
    call CorrectionFactorParsing(scalar_factor,matrix_factor,CF_ID,Scalar0Matrix1, &
        & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a,FS_a, dim)            
    Cdwdx_a(:)=dwdx(:)
    call CorrectedKernelGradient(Cdwdx_a, scalar_factor, matrix_factor, Scalar0Matrix1, dim)
            
    call CorrectionFactorParsing(scalar_factor,matrix_factor,CF_ID,Scalar0Matrix1, &
        & gamma_cont_b, gamma_discrt_b, gamma_mat_b, gamma_mat_inv_b, xi1_mat_inv_b,FS_b,dim)
    Cdwdx_b(:)=-dwdx(:)
    call CorrectedKernelGradient(Cdwdx_b, scalar_factor, matrix_factor, Scalar0Matrix1, dim)    

    call VectorDivergencePtoP(divF_a,divF_b,F_a,F_b,Cdwdx_a, Cdwdx_b, mass_a, mass_b, rho_a, rho_b, dim, div_type)
    
endsubroutine
    


subroutine CorrectedVecDivPtoB(divF_a,F_a,F_b,del_gamma_as,&
                    & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, &
                    & FS_a, dim, CF_ID, div_type) 

    !vector divergence for particle itneracting with boundary is calculated, by first correctign the kernel gradient

    implicit none
    integer(4), intent(in) :: dim, CF_ID, div_type
    real(8), intent(in) :: F_a(dim), F_b(dim), del_gamma_as(dim), &
        & gamma_cont_a, gamma_discrt_a, gamma_mat_a(dim,dim), gamma_mat_inv_a(dim,dim), xi1_mat_inv_a(dim,dim)
    real(8), intent(inout) :: divF_a 
    integer(4) :: d, Scalar0Matrix1
    integer(2) :: FS_a
    real(8) :: Cdgmas(dim), matrix_factor(dim,dim), scalar_factor
    
    call CorrectionFactorParsing(scalar_factor,matrix_factor,CF_ID,Scalar0Matrix1, &
                & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a,FS_a, dim)       
    Cdgmas(:)=del_gamma_as(:)
    call CorrectedKernelGradient(Cdgmas, scalar_factor, matrix_factor, Scalar0Matrix1, dim)  
    
    call VectorDivergencePtoB(divF_a,F_a,F_b,Cdgmas, dim, div_type)
    
endsubroutine