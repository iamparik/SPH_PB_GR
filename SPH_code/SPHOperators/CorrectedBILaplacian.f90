subroutine CorrectedBILapPtoP(dF_a,dF_b,F_a,F_b,dwdx, mass_a, mass_b, rho_a, rho_b, &
                    & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, &
                    & gamma_cont_b, gamma_discrt_b, gamma_mat_b, gamma_mat_inv_b, xi1_mat_inv_b, &
                    & dim, x_ab, CF_ID) 
    !Morris/Brookshaw's Laplacian for two particles is calculated, by first correctign the kernel gradient
    
    implicit none
    integer(4), intent(in) :: dim, CF_ID 
    real(8), intent(in) :: x_ab(dim),F_a, F_b, dwdx(dim), mass_a, mass_b, rho_a, rho_b, &
        & gamma_cont_a, gamma_discrt_a, gamma_mat_a(dim,dim), gamma_mat_inv_a(dim,dim), xi1_mat_inv_a(dim,dim), &
        & gamma_cont_b, gamma_discrt_b, gamma_mat_b(dim,dim), gamma_mat_inv_b(dim,dim), xi1_mat_inv_b(dim,dim)
    real(8), intent(inout) :: dF_a,dF_b
    integer(4) :: d, Scalar0Matrix1
    real(8) :: Cdwdx_a(dim), Cdwdx_b(dim), matrix_factor(dim,dim), scalar_factor
    
    call CorrectionFactorParsing(CF_ID,Scalar0Matrix1,scalar_factor,matrix_factor, &
        & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, dim)            
    Cdwdx_a(:)=dwdx(:)
    call CorrectedKernelGradient(Cdwdx_a, scalar_factor, matrix_factor, Scalar0Matrix1, dim)
            
    call CorrectionFactorParsing(CF_ID,Scalar0Matrix1,scalar_factor,matrix_factor, &
        & gamma_cont_b, gamma_discrt_b, gamma_mat_b, gamma_mat_inv_b, xi1_mat_inv_b, dim)
    Cdwdx_b(:)=-dwdx(:)
    call CorrectedKernelGradient(Cdwdx_b, scalar_factor, matrix_factor, Scalar0Matrix1, dim)    

    dF_a= dF_a + 2.D0*(F_a-F_b)*(dot_product(x_ab,Cdwdx_a)/norm2(x_ab)**2)*mass_b/rho_b
    dF_b= dF_b + 2.D0*(F_b-F_a)*(dot_product(-x_ab,Cdwdx_b)/norm2(x_ab)**2)*mass_a/rho_a

endsubroutine
    
    
subroutine BILapPtoB_PCG(dF_a, F_s, Dfdn_s, del_gamma_as,&
                    & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, &
                    & dim, CF_ID ,dirich0Neum1) 

    !vector divergence for particle itneracting with boundary is calculated, by first correctign the kernel gradient

    implicit none
    integer(4), intent(in) :: dim, CF_ID, dirich0Neum1
    real(8), intent(in) :: F_s (dim) ,Dfdn_s, del_gamma_as(dim), &
        & gamma_cont_a, gamma_discrt_a, gamma_mat_a(dim,dim), gamma_mat_inv_a(dim,dim), xi1_mat_inv_a(dim,dim)
    real(8), intent(inout) :: dF_a 
    integer(4) :: d, Scalar0Matrix1
    real(8) :: Cdgmas(dim), matrix_factor(dim,dim), scalar_factor
    
    call CorrectionFactorParsing(CF_ID,Scalar0Matrix1,scalar_factor,matrix_factor, &
                & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, dim)       
    Cdgmas(:)=del_gamma_as(:)
    call CorrectedKernelGradient(Cdgmas, scalar_factor, matrix_factor, Scalar0Matrix1, dim)  
    
    if(dirich0Neum1 .eq. 0) then
        dF_a= dF_a - dot_product(F_s,Cdgmas)
    elseif(dirich0Neum1 .eq. 1) then
        dF_a =dF_a - Dfdn_s*norm2(Cdgmas)
    else
        pause
    endif
    
endsubroutine