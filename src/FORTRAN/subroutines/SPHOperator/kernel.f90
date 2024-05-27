      subroutine kernel(r,dx,hsml,w,dwdx)
    !----------------------------------------------------------------------
    !   Subroutine to calculate the smoothing kernel wij and its 
    !   derivatives dwdxij.
    !     if skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
    !            = 2, Gauss kernel   (Gingold and Monaghan 1981) 
    !            = 3, Quintic kernel (Morris 1997)
    !     r    : Distance between particles i and j                     [in]
    !     dx   : x-, y- and z-distance between i and j                  [in]  
    !     hsml : Smoothing length                                       [in]
    !     w    : Kernel for all interaction pairs                      [out]
    !     dwdx : Derivative of kernel with respect to x, y and z       [out]
          use config_parameter, only:SPH_dim,skf,pi
          implicit none
    !      
          real(8) r, dx(SPH_dim), hsml, w, dwdx(SPH_dim)
          integer(4) d
          real(8) q, dw, factor
          real(8) a, k
    !
          q = r/hsml 
          w = 0.D0
          do d=1,SPH_dim         
                dwdx(d) = 0.D0
          enddo
    !----
          if (skf.eq.1) then  !Cubic spline kernel
    !      
                if (SPH_dim.eq.1) then
                    factor = 1.D0/hsml
                elseif (SPH_dim.eq.2) then
                    factor = 15.D0/(7.D0*pi*hsml*hsml)
                elseif (SPH_dim.eq.3) then
                    factor = 3.D0/(2.D0*pi*hsml*hsml*hsml)
                else
                     print *,' >>> Error <<< : Wrong dimension: SPH_dim =',SPH_dim
                     pause
                endif     
                                              
                if (q.ge.0.and.q.le.1.e0) then          
                      w = factor * (2.D0/3.D0 - q*q + q**3 / 2.D0)
                      do d = 1, SPH_dim
                            dwdx(d) = factor * (-2.D0+3.D0/2.D0*q)/hsml**2 * dx(d)       
                      enddo   
                else if (q.gt.1.e0.and.q.le.2) then          
                     w = factor * 1.D0/6.D0 * (2.D0-q)**3 
                     do d = 1, SPH_dim
                            dwdx(d) =-factor * 1.D0/6.D0 * 3.D0*(2.D0-q)**2/hsml * (dx(d)/r)
                     enddo              
	         else
	             w=0.
                    do d= 1, SPH_dim
                            dwdx(d) = 0.
                    enddo             
                endif     
    !---                                    
          else if (skf.eq.2) then  !Gauss kernel
    !      
               factor = 1.D0 / (hsml**SPH_dim * pi**(SPH_dim/2.D0))      
	        if(q.ge.0.and.q.le.3) then
	             w = factor * exp(-q*q)
                    do d = 1, SPH_dim
                        dwdx(d) = w * (-2.D0*dx(d)/hsml/hsml)
                    enddo 
	        else
	             w = 0.
                    do d = 1, SPH_dim
                        dwdx(d) = 0.
                    enddo 	   
	        endif	       
    !----	
          else if (skf.eq.3) then  !Quintic kernel
    !      
                if (SPH_dim.eq.1) then
                    factor = 1.0D0 / (120.0D0*hsml)
                elseif (SPH_dim.eq.2) then
                    factor = 7.0D0 / (478.0D0*pi*hsml*hsml)
                elseif (SPH_dim.eq.3) then
                    factor = 1.0D0 / (120.0D0*pi*hsml*hsml*hsml)
                else
                    write(*,'(A,2x,I2)') ' >>> Error <<< : Wrong dimension: SPH_dim =', SPH_dim
                    pause

                endif 
                             
	        if(q.ge.0.and.q.le.1) then
                      w = factor * ( (3.0D0-q)**5 - 6.0D0*(2.D0-q)**5 + 15.0D0*(1.0D0-q)**5 )
                      do d= 1, SPH_dim
                            dwdx(d) = factor * ( -5.0D0*(3.0D0-q)**4 + 30.0D0*(2.0D0-q)**4 - 75.0D0*(1.0D0-q)**4 ) / hsml * (dx(d)/r) 
                            if(r.le.1D-5) then
                              ! write(*,'(A)') 'kernel: two particles are too close!!'
                            endif        
                     enddo 
	        else if(q.gt.1.and.q.le.2) then
                      w = factor * ( (3.0D0-q)**5 - 6.0D0*(2.0D0-q)**5 )
                      do d= 1, SPH_dim
                        dwdx(d) = factor * (-5.0D0*(3.0D0-q)**4 + 30.0D0*(2.0D0-q)**4) / hsml * (dx(d)/r)
                        if(r.le.1D-5) then
                           write(*,'(A)') 'kernel: two particles are too close!!'
                        endif  
                      enddo 
                else if(q.gt.2.and.q.le.3) then
                      w = factor * (3.0D0-q)**5 
                      do d= 1, SPH_dim
                        dwdx(d) = factor * (-5.0D0*(3.0D0-q)**4) / hsml * (dx(d)/r)
                        if(r.le.1D-5) then
                           write(*,'(A)') 'kernel: two particles are too close!!'
                        endif  
                      enddo 
                else   
	              w = 0.
                      do d = 1, SPH_dim
                        dwdx(d) = 0.0D0
                      enddo  
                endif
    !---
          else if (skf.eq.4) then  !Yang kernel
                a = -10.D0/((k+31.D0)*pi*hsml*hsml)
                k = -1.D0
                if(q<1) then
                     w = a*(k*q**3 - 3.D0*(k+1.D0)*q*q + 3.D0*(k+3.D0)*q - k - 7.D0)
                     a = a*(3.D0*k*q*q - 6.D0*(k+1.D0)*q + 3.D0*(k+3.D0))/hsml/r
                     dwdx = a*dx
                else
                     w = a*(q-2.D0)**3
                     a = a*3.D0*(q-2.D0)**2/hsml/r
                     dwdx = a*dx
                end if
                
          elseif (skf.eq.5) then ! Wendland Quintic Kernel https://pysph.readthedocs.io/en/latest/reference/kernels.html
                if (SPH_dim.eq.1) then
                    print *,' Wendlan Quintic not defined for dimension: SPH_dim =',SPH_dim
                    pause
                elseif (SPH_dim.eq.2) then
                    factor = 7.D0/(4.D0*pi*hsml*hsml)
                elseif (SPH_dim.eq.3) then
                    factor = 21.D0/(16.D0*pi*hsml*hsml*hsml)
                else
                     print *,' >>> Error <<< : Wrong dimension: SPH_dim =',SPH_dim
                     pause
                endif 
                
                if ((q.ge.0.D0).and. (q.le. 2.D0)) then          
                    w = factor * (1.D0 - q/2.D0)**4 *(2.D0*q+1.D0)
                    do d = 1, SPH_dim
                        !dwdx(d) = factor *((1.D0-q/2.D0)**3)* (-5.D0*q)/ hsml * (dx(d)/r) 
                        dwdx(d) = factor *((1.D0-q/2.D0)**3)* (-5.D0/ hsml**2) * dx(d) 
                    enddo   
                else
                    w=0.D0
                    do d= 1, SPH_dim
                        dwdx(d) = 0.
                    enddo             
                endif
    !                
          endif
    end