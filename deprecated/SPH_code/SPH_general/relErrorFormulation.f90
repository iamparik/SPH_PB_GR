subroutine relErrorFormulation(F_analytical,F_SPH,RE,REB,REI,errorType)
    use config_parameter, only: SPH_dim,itype_real_max,itype_real_min
    use particle_data, only: gamma_cont,ntotal,itype, mass, rho
    
implicit none

integer(4)::  a,n
logical::errorType
real(8), intent(in):: F_analytical(ntotal), F_SPH(ntotal)
real(8)::RE(2),REB(2),REI(2)
real(8):: relError, relErrorBoundary, relErrorInside, Error, ErrorBoundary, ErrorInside
real(8):: maxrelError, maxrelErrorBoundary, maxrelErrorInside, maxError, maxErrorBoundary, maxErrorInside

    relError=0.D0
    relErrorBoundary= 0.D0   
    relErrorInside=0.D0
    Error=0.D0
    ErrorBoundary=0.D0
    ErrorInside=0.D0
    
    maxrelError=0.D0
    maxrelErrorBoundary=0.D0
    maxrelErrorInside=0.D0
    maxError=0.D0
    maxErrorBoundary=0.D0
    maxErrorInside=0.D0
    
    do a=1, ntotal
        
        if ((itype(a) .gt. itype_real_min) .and. (itype(a) .le. itype_real_max)) then 
            
            
            if(errorType) then
                if(F_analytical(a) .ne. 0.D0) then 
                    relError= relError + (mass(a)/rho(a))*((F_analytical(a)-F_SPH(a))/F_analytical(a))**2.D0
                    maxrelError= max(maxrelError,abs((F_analytical(a)-F_SPH(a))/F_analytical(a)))
                
                    if(gamma_cont(a) .lt. 1.D0) then
                        relErrorBoundary=relErrorBoundary + (mass(a)/rho(a))*((F_analytical(a)-F_SPH(a))/F_analytical(a))**2.D0
                        maxrelErrorBoundary= max(maxrelErrorBoundary,abs((F_analytical(a)-F_SPH(a))/F_analytical(a)))                    
                    endif
                
                    if(gamma_cont(a) .eq. 1.D0 ) then
                        relErrorInside=relErrorInside + (mass(a)/rho(a))*((F_analytical(a)-F_SPH(a))/F_analytical(a))**2.D0
                        maxrelErrorInside= max(maxrelErrorInside,abs((F_analytical(a)-F_SPH(a))/F_analytical(a)))
                    endif
                else
                    write(*,*) "relError will be NAN hence fanalytical is used as 1"
                    relError= relError + (mass(a)/rho(a))*((F_analytical(a)-F_SPH(a))/1.D0)**2.D0
                    maxrelError= max(maxrelError,abs((F_analytical(a)-F_SPH(a))/1.D0))
                
                    if(gamma_cont(a) .lt. 1.D0) then
                        relErrorBoundary=relErrorBoundary + (mass(a)/rho(a))*((F_analytical(a)-F_SPH(a))/1.D0)**2.D0
                        maxrelErrorBoundary= max(maxrelErrorBoundary,abs((F_analytical(a)-F_SPH(a))/1.D0))                    
                    endif
                
                    if(gamma_cont(a) .eq. 1.D0 ) then
                        relErrorInside=relErrorInside + (mass(a)/rho(a))*((F_analytical(a)-F_SPH(a))/1.D0)**2.D0
                        maxrelErrorInside= max(maxrelErrorInside,abs((F_analytical(a)-F_SPH(a))/1.D0))
                    endif
                    
                endif
                
            else            
                
                Error= Error + (mass(a)/rho(a))*(F_analytical(a)-F_SPH(a))**2.D0
                maxError= max(maxError,abs(F_analytical(a)-F_SPH(a)))
                
                if(gamma_cont(a) .lt. 1.D0) then
                    ErrorBoundary= ErrorBoundary + (mass(a)/rho(a))*(F_analytical(a)-F_SPH(a))**2.D0
                    maxErrorBoundary= max(maxErrorBoundary,abs(F_analytical(a)-F_SPH(a)))                   
                endif
                
                if(gamma_cont(a) .eq. 1.D0 ) then
                    ErrorInside= ErrorInside + (mass(a)/rho(a))*(F_analytical(a)-F_SPH(a))**2.D0
                    maxErrorInside= max(maxErrorInside,abs(F_analytical(a)-F_SPH(a)))
                endif
            endif
             
        endif
        
    enddo
    
    relError = DSQRT(relError)
    relErrorBoundary = DSQRT(relErrorBoundary)
    relErrorInside = DSQRT(relErrorInside)
    
    Error = DSQRT(Error)
    ErrorBoundary = DSQRT(ErrorBoundary)
    ErrorInside = DSQRT(ErrorInside)
    
    !write(*,*)" For case",n," RelError=", relError," relErrorBoundary=", relErrorBoundary," relErrorInside=", relErrorInside
    !write(*,*)" *****"
    !write(*,*)" For case",n," MaxRelError=", maxrelError," MaxErrorBoundary=", maxrelErrorBoundary," MaxrelErrorInside=", maxrelErrorInside
    !write(*,*)" *****"
    !write(*,*)" For case",n," Error=", Error," ErrorBoundary=", ErrorBoundary," ErrorInside=", ErrorInside
    
    !write(*,*) relError," ", relErrorBoundary," ", relErrorInside
    !write(*,*) maxrelError," ", maxrelErrorBoundary," ", maxrelErrorInside
    
    
    if(errorType) then
        
        write(*,*) relError," ", relErrorBoundary," ", relErrorInside
        write(*,*) maxrelError," ", maxrelErrorBoundary," ", maxrelErrorInside
    
        RE(1) = DLOG10(relError)
        RE(2) = DLOG10(maxrelError)
        REB(1)= DLOG10(relErrorBoundary)
        REB(2)= DLOG10(maxrelErrorBoundary)
        REI(1)= DLOG10(relErrorInside)
        REI(2)= DLOG10(maxrelErrorInside)

    else
        
        write(*,*) Error," ", ErrorBoundary," ", ErrorInside
        write(*,*) maxError," ", maxErrorBoundary," ", relErrorInside
        
        
        if(Error .eq. 0) then
            RE(1)= -30
        else            
            RE(1) = DLOG10(Error)
        endif
        
        if(maxError .eq. 0) then
            RE(2)= -30
        else            
            RE(2) = DLOG10(maxError)
        endif
        
        if(ErrorBoundary .eq. 0) then
            REB(1)= -30
        else            
            REB(1) = DLOG10(ErrorBoundary)
        endif
        
        if(maxErrorBoundary .eq. 0) then
            REB(2)= -30
        else            
            REB(2) = DLOG10(maxErrorBoundary)
        endif
        
        if(ErrorInside .eq. 0) then
            REI(1)= -30
        else            
            REI(1) = DLOG10(ErrorInside)
        endif
        
        if(maxErrorInside .eq. 0) then
            REI(2)= -30
        else            
            REI(2) = DLOG10(maxErrorInside)
        endif
        
    endif
    
    
    
end
    