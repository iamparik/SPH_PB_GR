subroutine inputCADtoParticleData(k, nreal_mesh, totVol)
    use config_parameter, only : SPH_dim,ExtInputMeshType, DataConfigPath, dx_r
    use particle_data , only : x, itype, vol
    
    implicit none
    
    integer(4), intent(in) :: nreal_mesh
    real(8), intent(in) :: totVol
    integer(4), intent(inout) :: k
    integer d
    real(8) tempType, avgVol
    
    
    
    ! Read input file containing particle/point information    
    if( ExtInputMeshType .eq. 1) then 
        ! Use cartesian cordinate represent particles with area  
        ! used directly from the CAD file
        open(1,file= DataConfigPath // '/input_particles_cartesianMesh.dat',status='old')

        do while (.not.eof(1))
            k=k+1
            read(1,*) tempType, (x(d, k), d=1,SPH_dim), vol(k)
            call ICinput(itype(k), NINT(tempType))       
        enddo
    
    elseif(ExtInputMeshType .eq. 2) then 
        ! Use cartesian cordinate represent particles with area being  
        ! the area evaluated with particle spacing
        ! for no errors ExtInputMeshType 1 and 2 must be the same
        open(1,file= DataConfigPath // '/input_particles_cartesianMesh.dat',status='old')

        avgVol=dx_r**(SPh_dim)
        
        do while (.not.eof(1))
            k=k+1
            read(1,*) tempType, (x(d, k), d=1,SPH_dim), vol(k)
            vol(k)=avgVol
            call ICinput(itype(k), NINT(tempType))      
        enddo
            
    elseif ( ExtInputMeshType .eq. 4) then 
        ! Use tri/quad mesh centroids with constant avg area for all particles
        ! evaluated with the particle spacing
        open(1,file= DataConfigPath // '/input_particles_bulkMesh.dat',status='old') 

        avgVol=totVol/dble(nreal_mesh)
        do while (.not.eof(1))
            k=k+1
            read(1,*) tempType, (x(d, k), d=1,SPH_dim), vol(k)
            vol(k) =avgVol
            call ICinput(itype(k), NINT(tempType))  
        enddo
    elseif ( ExtInputMeshType .eq. 3) then 
        ! Use tri/quad mesh centroids with area from
        ! CAD file directly used
        open(1,file= DataConfigPath // '/input_particles_bulkMesh.dat',status='old') 
            
        do while (.not.eof(1))
            k=k+1
            read(1,*) tempType, (x(d, k), d=1,SPH_dim), vol(k)
            call ICinput(itype(k), NINT(tempType))  
        enddo
            
    elseif(ExtInputMeshType .eq. 5) then
        ! Use packed particle configuration with area
        ! being the area in particle file
        open(1,file= DataConfigPath // '/input_PP.dat',status='old')

        do while (.not.eof(1))
            k=k+1
            read(1,*) tempType, (x(d, k), d=1,SPH_dim), vol(k) 
            call ICinput(itype(k), NINT(tempType))
        enddo
    endif

end subroutine