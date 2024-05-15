subroutine surfaceNormalBdry(temp_surf_norm, temp_x, dim)
    implicit none
    integer(4) :: dim
    real(8), intent(out) :: temp_surf_norm(dim)
    real(8), intent(in) :: temp_x(dim,dim)
    real(8) :: mn
    
    if(dim .eq. 2) then
        temp_surf_norm(1)= -(temp_x(2,2)- temp_x(2,1))
        temp_surf_norm(2)= (temp_x(1,2)- temp_x(1,1))
        mn= norm2(temp_surf_norm(:))
        temp_surf_norm(1)= temp_surf_norm(1)/mn
        temp_surf_norm(2)= temp_surf_norm(2)/mn
        
    else
        write(*,*) "dimensn not defined"
        
    endif
    
endsubroutine