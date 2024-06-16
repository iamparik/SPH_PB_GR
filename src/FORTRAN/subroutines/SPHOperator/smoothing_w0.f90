!****************************************************************************
!
!  SUBROUTINE: smoothing_w0
!
!  PURPOSE: This function allows calculation of smoothing weight at radius zero
!           That is W(x_a-x_a, kh)
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/18/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine smoothing_w0
    use config_parameter, only:SPH_dim
    use particle_data ,   only: ntotal,itype,hsml,w_aa
    implicit none
    
    real(8) mhsml,a,dx_dw_0(SPH_dim)
    
    if (.NOT. Allocated(w_aa)) allocate(w_aa(ntotal))
    w_aa=0.D0
    dx_dw_0=0.D0
    
    
    do a=1,ntotal
        mhsml=hsml(a)
        call kernel(0.D0,dx_dw_0,mhsml,w_aa(a),dx_dw_0)
    enddo
    
end
    
    