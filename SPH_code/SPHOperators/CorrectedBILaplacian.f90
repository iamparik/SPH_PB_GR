subroutine CorrectedBILapPtoP(dF_a,dF_b,F_a,F_b,dwdx, mass_a, mass_b, rho_a, rho_b, &
                    & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, &
                    & gamma_cont_b, gamma_discrt_b, gamma_mat_b, gamma_mat_inv_b, xi1_mat_inv_b, &
                    & dim, x_ab,CF_ID) 
    !Morris/Brookshaw's Laplacian for two particles is calculated, by first correctign the kernel gradient
    
    implicit none
    integer(4), intent(in) :: dim, CF_ID
    real(8), intent(in) :: x_ab(dim),F_a, F_b, dwdx(dim), mass_a, mass_b, rho_a, rho_b, &
        & gamma_cont_a, gamma_discrt_a, gamma_mat_a(dim,dim), gamma_mat_inv_a(dim,dim), xi1_mat_inv_a(dim,dim), &
        & gamma_cont_b, gamma_discrt_b, gamma_mat_b(dim,dim), gamma_mat_inv_b(dim,dim), xi1_mat_inv_b(dim,dim)
    real(8), intent(inout) :: dF_a,dF_b
    integer(4) :: d, Scalar0Matrix1 
    real(8) :: Cdwdx_a(dim), Cdwdx_b(dim), matrix_factor(dim,dim), scalar_factor
    
    ! By default only gamma_cont is used for BILaplacian formulation
    
    call CorrectionFactorParsing(scalar_factor,matrix_factor,CF_ID,Scalar0Matrix1, &
        & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, dim)            
    Cdwdx_a(:)=dwdx(:)
    call CorrectedKernelGradient(Cdwdx_a, scalar_factor, matrix_factor, Scalar0Matrix1, dim)
            
    call CorrectionFactorParsing(scalar_factor,matrix_factor,CF_ID,Scalar0Matrix1, &
        & gamma_cont_b, gamma_discrt_b, gamma_mat_b, gamma_mat_inv_b, xi1_mat_inv_b, dim)
    Cdwdx_b(:)=-dwdx(:)
    call CorrectedKernelGradient(Cdwdx_b, scalar_factor, matrix_factor, Scalar0Matrix1, dim)    

    dF_a= dF_a + 2.D0*(F_a-F_b)*(dot_product(x_ab,Cdwdx_a)/norm2(x_ab)**2)*mass_b/rho_b
    dF_b= dF_b + 2.D0*(F_b-F_a)*(dot_product(-x_ab,Cdwdx_b)/norm2(x_ab)**2)*mass_a/rho_a

endsubroutine
    
    
    
    
    
    
subroutine CorrectedBILapPtoB(dF_a, vec_val, scalar_val, del_gamma_as,&
                    & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, &
                    & dim, dirich0Neum1, CF_ID) 

    !vector divergence for particle itneracting with boundary is calculated, by first correcting the kernel gradient

    implicit none
    integer(4), intent(in) :: dim, dirich0Neum1, CF_ID
    real(8), intent(in) :: vec_val(dim), scalar_val , del_gamma_as(dim), &
        & gamma_cont_a, gamma_discrt_a, gamma_mat_a(dim,dim), gamma_mat_inv_a(dim,dim), xi1_mat_inv_a(dim,dim)
    real(8), intent(inout) :: dF_a 
    integer(4) :: d, Scalar0Matrix1
    real(8) :: Cdgmas(dim), matrix_factor(dim,dim), scalar_factor
    
    call CorrectionFactorParsing(scalar_factor,matrix_factor,CF_ID,Scalar0Matrix1, &
                & gamma_cont_a, gamma_discrt_a, gamma_mat_a, gamma_mat_inv_a, xi1_mat_inv_a, dim)       
    Cdgmas(:)=del_gamma_as(:)
    call CorrectedKernelGradient(Cdgmas, scalar_factor, matrix_factor, Scalar0Matrix1, dim)  
    

    if(dirich0Neum1 .eq. 0) then
        dF_a= dF_a - dot_product(vec_val,Cdgmas)
    elseif(dirich0Neum1 .eq. 1) then
        dF_a =dF_a - scalar_val*norm2(Cdgmas)
    endif
        

    
endsubroutine
    

    
    
    
subroutine BILViscousBdry(vec_val,scalar_val,dirich0Neum1, grad_vel_a, grad_vel_s, vx_a, vx_s, delx_ab, dx_r, &
                    & surf_norm_s, d, dim, BIL_type)

    implicit none
    integer(4), intent(in) :: d, dim, BIL_type
    integer(4), intent(out) :: dirich0Neum1
    real(8), intent(out) :: vec_val(dim), scalar_val 
    real(8), intent(in) :: grad_vel_a(dim), grad_vel_s(dim), vx_a(dim), vx_s(dim), delx_ab(dim), surf_norm_s(dim), dx_r

    
    if(BIL_type .eq. 1) then    !BIL-PCG Formulation (using ∇v_a + ∇v_s)
        vec_val(:)=grad_vel_a(:) + grad_vel_s(:) 
        dirich0Neum1=0
        
    elseif(BIL_type .eq. 2) then    !BIL-PCG Formulation (using 2*∇v_a)
        vec_val(:)=2.D0*grad_vel_a(:)
        dirich0Neum1=0
     
    elseif(BIL_type .eq. 3) then    !BIL-PCG Formulation (using 2*∇v_s)
        vec_val(:)=2.D0*grad_vel_s(:)
        dirich0Neum1=0
     
    elseif(BIL_type .eq. 4) then    !BIL-Macia Formulation
        scalar_val=2.D0*(vx_a(d)-vx_s(d))*dot_product(delx_ab, surf_norm_s)/norm2(delx_ab)**2
        dirich0Neum1=1
    
    elseif(BIL_type .eq. 5) then    !BIL-NTG formulation
        scalar_val=2.D0*(vx_a(d)-vx_s(d))/dot_product(delx_ab, surf_norm_s) 
        dirich0Neum1=1
     
    elseif(BIL_type .eq. 6) then    !BIL-BLT (USAW) formulation

        vec_val(:)=(vx_a - dot_product(vx_a,surf_norm_s) * surf_norm_s) &
                    & /(max(dot_product(delx_ab,surf_norm_s), dx_r))
        scalar_val = vec_val(d)
        dirich0Neum1=1
    else
        write(*,*) 'BILtype not defined'
        pause
        
        
    endif   
    
    
endsubroutine