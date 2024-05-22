subroutine CorrectionFactorParsing(scalar_factor,matrix_factor,id,Scalar0Matrix1,g_c, g_d, g_m, g_m_i, xi_m_i, dim)
    !This subroutine helps identify the correction factor for a given id
    implicit None
    integer(4), intent(in) :: dim,id
    integer(4), intent(out) :: Scalar0Matrix1
    real(8), intent(out) :: scalar_factor, matrix_factor(dim,dim)
    real(8), intent(in) :: g_c, g_d, g_m(dim,dim),g_m_i(dim,dim), xi_m_i(dim,dim)
    integer(4) :: d

    
    matrix_factor = 0.D0
    scalar_factor = 0.D0
    
    if( id .eq. 1) then
        !This uses gamma_cont as a factor
        scalar_factor = g_c
        Scalar0Matrix1 = 0
     
    elseif( id .eq. 2 ) then
        !This uses gamma_discrt as a factor
        scalar_factor = g_d
        Scalar0Matrix1 = 0
        
    elseif( id .eq. 3 ) then
        !This uses PCBI matrix inverse as a factor
        matrix_factor=g_m_i
        Scalar0Matrix1 = 1
    
    elseif( id .eq. 4 ) then
        !This uses diagonal terms of PCBI matrix as a factor
        do d=1,dim
            matrix_factor(d,d)= 1.D0/g_m(d,d)
        enddo
        Scalar0Matrix1 = 1
    
    elseif( id .eq. 5 ) then
        !This uses L matrix inverse (essentially PCBI without bdry) as a factor
        matrix_factor=xi_m_i
        Scalar0Matrix1 = 1
        
    else
        scalar_factor = 1.D0
        Scalar0Matrix1 = 0
    endif
    
endsubroutine
    
    
        
        
    