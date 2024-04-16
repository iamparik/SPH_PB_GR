subroutine marker_point_value(f)
    use config_parameter, only: itype_real_min
    use particle_data, only: itype, w, ntotal, niac, pair_i, pair_j, gamma_discrt, rho, mass
    implicit none
    
    integer(4) k,d,a,b
    real(8) f(ntotal)
    
    do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)

        if((itype(a) .eq. itype_real_min) .and. (rho(b) .ne. 0))  call fncnApproxOperator(f(a),f(b),mass(b),rho(b),w(k)/gamma_discrt(a))
        if((itype(b) .eq. itype_real_min) .and. (rho(a) .ne. 0))  call fncnApproxOperator(f(b),f(a),mass(a),rho(a),w(k)/gamma_discrt(b))

     enddo

    
end
    