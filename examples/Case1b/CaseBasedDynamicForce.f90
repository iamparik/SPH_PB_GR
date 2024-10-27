subroutine CaseBasedDynamicForce(F_ext,del_t)
    use particle_data, only: vx, ntotal, itype
    use config_parameter, only: SPH_dim, itype_periodic, itype_real_max, itype_real_min
    use case1b_mod, only: long_v_n_1,long_v_n_2, long_v_req
    
    implicit none
    
    real(8), dimension(SPH_dim), intent(inOut):: F_ext
    real(8), intent(in) :: del_t
    real(8) :: vx_current
    integer(4):: nperiodic,a
    
    
    nperiodic=0
    vx_current=0.D0
    
    do a=1,ntotal
       if( ((itype(a)-itype_periodic) .le. itype_real_max) .and. ((itype(a)-itype_periodic) .ge. itype_real_min)  )  then
           vx_current=vx_current+vx(1,a)
           nperiodic=nperiodic+1
       endif
    enddo
    vx_current=vx_current/dble(nperiodic)
    F_ext(1)= F_ext(1) + (long_v_req - 2.D0*long_v_n_1 + long_v_n_2)/del_t
    
    if (F_ext(1)*del_t .gt. long_v_req ) F_ext(1)=long_v_req/(1.D1*del_t)
    if( F_ext(1)*del_t .lt. 0.D0 ) F_ext(1)=0.D0
    
    long_v_n_2 = long_v_n_1
    long_v_n_1 = vx_current
    
    
    write(*,*) " vx_current = ", vx_current, " and F_ext = ", F_ext
    
    
    
    
end subroutine