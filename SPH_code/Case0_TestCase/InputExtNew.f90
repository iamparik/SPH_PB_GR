!****************************************************************************
!
!  SUBROUTINE: input
!
!  PURPOSE:  Subroutine that accepts any geometry input
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  11/02/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

    subroutine inputExt
        
    use config_parameter, only: SPH_dim,etype_periodic, etype_SolidWall1, etype_SolidWall2, &
        & etype_FreeSurface1, etype_thermal_neumann, etype_thermal_dirichlet, etype_FarWall, &
        & DataConfigPath, itype_virtual, itype_real_min, &
        & dx_r, hsml_const,hydrostaticHeight, rho_init, mu_const, &
        &  g_const, c_sound, F_ext, ExtInputMeshType, packagingIterations
    use particle_data ,   only: maxn, max_interaction, max_e_interaction, maxnv,  &
        & x,vx,mass,rho, vol,p,itype,hsml,mu, temp, nreal, nedge, nflow, nghost, ntotal, edge, &
        & surf_norm, nedge_rel_edge, etype, etotal, maxedge, &
        & pBC_edges,tangent_pBC, &
        &  f0max, dynamicProblem,avgVol, &
        & packed_x,packed_itype, packed_vol, simGridSize, &
        & x_ve
    !use ifport, only:random
    
    implicit none
    real(8) tempType
    integer(4) k,d,s,ke,periodicPairs,i , additionalParticles
    real(8) tempX(SPH_dim,SPH_dim), tempBdry(5), mn, scale_k
    integer(4) nreal_mesh, nrealCartesian, nEbulkBdry, nEdomainBdry
    real(8) totVol, x_ve_temp(SPH_dim,SPH_dim)
    character(40) :: input_file_name
   
    
    ! Confirm if the SPH particles are static/Eulerian or Dynamic/Lagrangian
    dynamicProblem=.true.
    
    ! Retreive the scaling factor for the smoothing length used in the smoothing kernel    
    call sml_mult_factor(scale_k)
    
    ! Read size of the domain from config file 
    write(*,*)'  **************************************************'
    write(*,*)'       Reading size of the SPH domain configuration...   '
    open(1,file= DataConfigPath // '/input_param_SimSize.dat',status='old')
    do while (.not.eof(1))
        ! I might have to modify below to strictly read as integers
          read(1,*) nreal_mesh, nrealCartesian, nEbulkBdry, nEdomainBdry, totVol
    enddo
    close(1)
    
    
    additionalParticles= 100 ! This is necessary only if expect to use periodic particles, or unflow, outflow particles
    if (packagingIterations .eq. 0) then
        maxedge=int(nEdomainBdry)*2
        maxnv=10*maxedge 
        maxn= max(nreal_mesh, nrealCartesian)+ maxedge +additionalParticles
    else
        maxedge=int(nEbulkBdry)*2
        maxnv=10*maxedge 
        maxn= max(nreal_mesh, nrealCartesian)+ maxedge + additionalParticles
    endif
    
    
    !variables required to define geometry are particle propoerties are allocated
    ALLOCATE(x(SPH_dim,maxn))
    ALLOCATE(mass(maxn), rho(maxn), vol(maxn), hsml(maxn), itype(maxn))
    
    ! The allocatable variables are now initialized
    x=0.D0
    mass=0.D0
    rho=0.D0
    vol=0.D0
    hsml=0.D0
    itype=0
    

    ! initialize particle number to 0  
    k=0
        
    ! read approriate data file according to input type
    call inputCADtoParticleData(k, nreal_mesh,totVol)


    nreal=k ! I can distinguish nreal and nedge here if needed
    write(*,*)'      Total number of Real particles : ', nreal    	
    write(*,*)'**************************************************'

    ! initialize edges to 0
    s=0
    
    input_file_name = '/input_bulkBdryEdge.dat'
    ! read bdry data file and store bdry information
    call inputCADtoEdgeData(s, maxedge, input_file_name, 10)
    
    etotal=s
    write(*,*)'      Total number of boundary segments : ', etotal    	
    write(*,*)'**************************************************'      

    ! Define the tangent for the periodic boundary, 
    ! One needs to make sure the normals of the boundary should be pointing awat from
    ! real domain for below to work
    !Now, to determine the tangential sense of the periodic edges we calcualte tangent_pBC
    if (allocated(pBC_edges)) then
        Allocate(tangent_pBC(SPH_dim,2,1))
        tangent_pBC(1,1,1)= x_ve(1,edge(2,pBC_edges(1,1)))- x_ve(1,edge(1,pBC_edges(1,1)))
        tangent_pBC(2,1,1)= x_ve(2,edge(2,pBC_edges(1,1)))- x_ve(2,edge(1,pBC_edges(1,1)))
        mn= norm2(tangent_pBC(:,1,1))
        tangent_pBC(:,1,1)=tangent_pBC(:,1,1)/mn
        
        ! We change the tengent direction of one of the periodic edge of a given pair.
        tangent_pBC(1,2,1)= x_ve(1,edge(1,pBC_edges(2,1)))- x_ve(1,edge(2,pBC_edges(2,1)))
        tangent_pBC(2,2,1)= x_ve(2,edge(1,pBC_edges(2,1)))- x_ve(2,edge(2,pBC_edges(2,1)))
        mn= norm2(tangent_pBC(:,2,1))
        tangent_pBC(:,2,1)=tangent_pBC(:,2,1)/mn
    endif

    !Now create particle per edge (this needs to be deleted when nedge_rel_Edge is decommisioned)
    ALLOCATE(nedge_rel_edge(maxedge))
    ! Now we will find the center of edges, which in Method 1 of Parikshit et. al will work as virtual points
    ! In method 2 it will be used as a real boundary particle 
    ! The below needs to be modified slightly to include more than one edge_temp particle per edge_temp (and to include different dx_v
    nedge=0
    do s=1,etotal
        k= k+1
        nedge_rel_edge(s)=k ! Here nedge_rel_edge(i)=i  is a simple case of one edge_temp particle per edge_temp
        
        do d =1,SPH_dim
            x_ve_temp(:,d)=x_ve(:,edge(d,s))
        enddo        
        call centroidBdrySegment(x(:,k), x_ve_temp, SPH_dim)

       
        itype(k)= itype_virtual+itype_real_min

        nedge=nedge+1
    enddo
    
    ntotal=nreal+nedge
    
    ! Maximum interactions for particle-particle and edge_temp-particle is defined
    max_interaction= 10*maxn*(ceiling((hsml_const/dx_r)*scale_k*2))**SPH_dim
    max_e_interaction= 10*maxedge*(ceiling((hsml_const/dx_r)*scale_k*2))**SPH_dim 

    write(*,*) "maximum particle-particle and particle-point interaction defined : ",max_interaction
    write(*,*) "maximum particle-edge_temp interaction defined : ",max_e_interaction
    write(*,*)'**************************************************'   
    
    ! find size of the simulation
    ALLOCATE( simGridSize(SPH_dim,2))
    simGridSize(:,1)= MINVAL(x_ve,2) 
    simGridSize(:,2)= MAXVAL(x_ve,2) 
    write(*,*)'**************************************************'   
    write(*,*) "simulation domain size is ", simGridSize, "with xmin = ", simGridSize(1,1), "xmax=",simGridSize(1,2)
    write(*,*) "with ymin = ", simGridSize(2,1), "ymax=",simGridSize(2,2)
    write(*,*)'**************************************************'   
    
    do k=1,ntotal
        !For changing density to hydrostatic
        rho(k)=rho_init
        mass(k)=rho(k)*vol(k)
        hsml(k)=hsml_const
    enddo
    

    ! Now perform the particle packing algorithm 
    if(packagingIterations) then
        call particlePackingTimeIntegration(.true.)
    !input : quick_converge_step2C  ! use .true. to enable quickconverge    
   
        deallocate(simGridSize)
    endif

    DEALLOCATE(mass, rho, hsml)
    
    end