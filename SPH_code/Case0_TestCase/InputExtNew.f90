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
        & x,vx,mass,rho, vol,p,itype,hsml,mu, temp, nreal, nflow, nedge,nghost, ntotal, edge, &
        & surf_norm, nedge_rel_edge, etype, etotal, maxedge, &
        & pBC_edges,tangent_pBC, bdryVal_vel, bdryVal_prs, bdryVal_rho, bdryVal_temp, &
        &  f0max, dynamicProblem,avgVol, &
        & packed_x,packed_itype, packed_vol, simGridSize
    !use ifport, only:random
    
    implicit none
    real(8) tempType
    integer(4) k,d,s,ke,periodicPairs,i , additionalParticles
    real(8) tempX(SPH_dim,SPH_dim), tempBdry(5), mn, scale_k
    integer(4) nrealMesh, nrealCartesian, nEbulkBdry, nEdomainBdry
    real(8) totVol
    real(8),DIMENSION(:,:),ALLOCATABLE :: xV_temp, edge_temp
   
    
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
          read(1,*) nrealMesh, nrealCartesian, nEbulkBdry, nEdomainBdry, totVol
    enddo
    close(1)
    
    maxedge=nEdomainBdry*2
    maxnv=3 *maxedge *2
    additionalParticles= 100 ! This is necessary only if expect to use periodic particles, or unflow, outflow particles
    if (packagingIterations .eq. 0) then
        maxn= max(nrealMesh, nrealCartesian) +maxnv +additionalParticles
    else
        maxn= nreal+maxnv + additionalParticles
    endif
    
    
    !vaiuables required to define geometry are particle propoerties are allocated
    ALLOCATE(x(SPH_dim,maxn),vx(SPH_dim,maxn))
    ALLOCATE(mass(maxn), rho(maxn), vol(maxn), p(maxn), hsml(maxn), itype(maxn), mu(maxn), temp(maxn))
    
    ! The allocatable variables are now initialized
    x=0.D0
    vx=0.D0
    mass=0.D0
    rho=0.D0
    p=0.D0
    hsml=0.D0
    itype=0
    mu=0
    temp=0.D0
    vol=0.D0
    
    ! Read input file containing particle/point information  if particle hasnt been packed already  
    if( packagingIterations .eq. 0) then ! Use cartesian cordinate represent particles with area being average area
       
        ! Read input file containing particle/point information    
        if( ExtInputMeshType .eq. 1) then ! Use cartesian cordinate represent particles with area being average area
            open(1,file= DataConfigPath // '/input_particles_cartesianMesh.dat',status='old')
            ! initialize particle number to 0  
            k=0
            do while (.not.eof(1))
                k=k+1
                read(1,*) tempType, (x(d, k), d=1,SPH_dim), vol(k)
                itype(k) = NINT(tempType)        
            enddo
    
        elseif(ExtInputMeshType .eq. 2) then
            open(1,file= DataConfigPath // '/input_particles_cartesianMesh.dat',status='old')
            ! initialize particle number to 0  
            k=0
            avgVol=dx_r**(SPh_dim)
        
            do while (.not.eof(1))
                k=k+1
                read(1,*) tempType, (x(d, k), d=1,SPH_dim), vol(k)
                vol(k)=avgVol
                itype(k) = NINT(tempType)        
            enddo
            
        elseif ( ExtInputMeshType .eq. 3) then ! Use tri/quad mesh centroids with constant avg area for all particles
            open(1,file= DataConfigPath // '/input_particles_bulkMesh.dat',status='old') 
            ! initialize particle number to 0    
            k=0
            avgVol=totVol/dble(nrealMesh)
            do while (.not.eof(1))
                k=k+1
                read(1,*) tempType, (x(d, k), d=1,SPH_dim), vol(k)
                vol(k) =avgVol
                itype(k) = NINT(tempType)        
            enddo
        elseif ( ExtInputMeshType .eq. 4) then ! Use tri/quad mesh centroids with their true 2D area
            open(1,file= DataConfigPath // '/input_particles_bulkMesh.dat',status='old') 
            ! initialize particle number to 0    
            k=0
            
            do while (.not.eof(1))
                k=k+1
                read(1,*) tempType, (x(d, k), d=1,SPH_dim), vol(k)
                itype(k) = NINT(tempType)        
            enddo
            
        elseif(ExtInputMeshType .eq. 5) then
            open(1,file= DataConfigPath // '/input_PP.dat',status='old')
            ! initialize particle number to 0  
            k=0
            avgVol=dx_r**(SPh_dim)
        
            do while (.not.eof(1))
                k=k+1
                read(1,*) tempType, (x(d, k), d=1,SPH_dim), vol(k)
                !vol(k)=avgVol
                itype(k) = NINT(tempType)        
            enddo
        endif
    
    else
        !if packing algorithm is already executed, simply retrieve the information of packed particles
        k=0
        do i=1,nreal
            k=k+1
            x(:,k)=packed_x(:,i)
            itype(k)=packed_itype(i)
            vol(k)=packed_vol(i)
        enddo
        
    endif
    
    nreal=k ! I can distinguish nreal and nedge here if needed
    close(1)
    write(*,*)'      Total number of Real particles : ', nreal    	
    write(*,*)'**************************************************'

    
    
    Allocate(edge_temp(SPH_dim,maxedge))   
    Allocate(etype(maxedge))
    ALLOCATE(xV_temp(SPH_dim,maxnv))
    
    edge_temp=0
    etype=0
    xV_temp=0.D0
    
    ! Read input file containing edge_temp information
    open(1,file= DataConfigPath // '/input_domainBdryEdge.dat',status='old')
    ! initialize particle number to 0    
    s=0
    ke=0
    do while (.not.eof(1))
        s=s+1
        read(1,*) tempType, (tempX(d,1), d=1,SPH_dim), (tempX(d,2), d=1,SPH_dim) !, (tempbdry(d), d=1,5)
        !etype(s) = NINT(tempType)      
        call BCinput(etype(s), NINT(tempType) )
        do d=1,SPH_dim
            ke=ke+1
            xV_temp(:,ke)=tempX(:,d)
            edge_temp(d,s)=ke
        enddo
        
        if(etype(s) .eq. etype_SolidWall1) then
            if( .not. allocated(bdryVal_vel)) then 
                allocate( bdryVal_vel(SPH_dim,maxedge),bdryVal_prs(maxedge), bdryVal_rho(maxedge), bdryVal_temp(maxedge))
            endif
            do d=1,SPH_dim
                bdryVal_vel(d,s)= 0.D0!tempBdry(d)
            enddo
            bdryVal_prs(s)  = 0.D0 !tempBdry(3)
            bdryVal_rho(s)  =rho_init !tempBdry(4)
            bdryVal_temp(s) =300 !tempBdry(5)
             
        elseif(etype(s) .eq. etype_periodic) then  
            write(*,*) "periodic BCS available for only one pair of periodic boundary"
            periodicPairs=1
            
            if( .not. allocated(pBC_edges)) then                       
                allocate(pBC_edges(2,periodicPairs))
                pBC_edges(1,1)=s
            else
                pBC_edges(2,1)=s
            endif
            
            
        endif                  
    enddo
    etotal=s
    close(1)
    write(*,*)'      Total number of boundary segments : ', nedge    	
    write(*,*)'**************************************************'      

        
    ! Define surface normals of walls, with normal per edge_temp
    ! ** this can be made a function/seperate subroutine
    ALLOCATE(surf_norm(SPH_dim, maxedge))    
    do s = 1, etotal   
        surf_norm(1,s)= -(xV_temp(2,edge_temp(2,s))- xV_temp(2,edge_temp(1,s)))
        surf_norm(2,s)= (xV_temp(1,edge_temp(2,s))- xV_temp(1,edge_temp(1,s)))
        mn= norm2(surf_norm(:,s))
        surf_norm(1,s)= surf_norm(1,s)/mn
        surf_norm(2,s)= surf_norm(2,s)/mn
    enddo 
    
    ! Define the tangent for the periodic boundary, 
    ! One needs to make sure the normals of the boundary should be pointing awat from
    ! real domain for below to work
    !Now, to determine the tangential sense of the periodic edges we calcualte tangent_pBC
    if (allocated(pBC_edges)) then
        Allocate(tangent_pBC(SPH_dim,2,1))
        tangent_pBC(1,1,1)= xV_temp(1,edge_temp(2,pBC_edges(1,1)))- xV_temp(1,edge_temp(1,pBC_edges(1,1)))
        tangent_pBC(2,1,1)= xV_temp(2,edge_temp(2,pBC_edges(1,1)))- xV_temp(2,edge_temp(1,pBC_edges(1,1)))
        mn= norm2(tangent_pBC(:,1,1))
        tangent_pBC(:,1,1)=tangent_pBC(:,1,1)/mn
        
        ! We change the tengent direction of one of the periodic edge of a given pair.
        tangent_pBC(1,2,1)= xV_temp(1,edge_temp(1,pBC_edges(2,1)))- xV_temp(1,edge_temp(2,pBC_edges(2,1)))
        tangent_pBC(2,2,1)= xV_temp(2,edge_temp(1,pBC_edges(2,1)))- xV_temp(2,edge_temp(2,pBC_edges(2,1)))
        mn= norm2(tangent_pBC(:,2,1))
        tangent_pBC(:,2,1)=tangent_pBC(:,2,1)/mn
    endif

    
    ALLOCATE(nedge_rel_edge(maxedge))
    ! Now we will find the center of edges, which in Method 1 of Parikshit et. al will work as virtual points
    ! In method 2 it will be used as a real boundary particle 
    ! The below needs to be modified slightly to include more than one edge_temp particle per edge_temp (and to include different dx_v
    do s=1,etotal
        k= k+1
        nedge_rel_edge(s)=k ! Here nedge_rel_edge(i)=i  is a simple case of one edge_temp particle per edge_temp
        x(:,k)=0.D0        
        do d=1,SPH_dim
            x(:,k)= x(:,k)+xV_temp(:,edge_temp(d,s))/dble(SPH_dim) !this defines the location of edge_temp points/particles we need per edge_temp
        enddo
        
         if(etype(s) .eq. etype_SolidWall1) then             
            itype(k)= itype_virtual+itype_real_min
        endif 
 
    enddo
    nedge= k - nreal     
    write(*,*)'      Total number of edge_temp particles : ', nedge    	
    write(*,*)'**************************************************'      

    !Let us define the edge_temp vertices by ghost points to visually demonstrate edges on tecplot
    !This also lets us uniquely label ghost points in a global sense
    do i=1,ke
        k=k+1
        x(:,k)=xV_temp(:,i)  
        itype(k)=0    
    enddo
    
    allocate(edge(SPH_dim,maxedge))
    do s=1,etotal
        do i=1,SPH_dim
            edge(i,s)= edge_temp(i,s) +nreal+nedge
        enddo
    enddo
    
    deallocate(xV_temp,edge_temp)
    
    ! The total ghost points is stored and then printed on screen
    nghost=k-(nreal+nedge)      
    write(*,*)'      Total number of Vertex Points for edges / ghost points ', nghost    	
    write(*,*)'**************************************************'      

    ntotal=nreal+nedge+nghost
    
    ! Maximum interactions for particle-particle and edge_temp-particle is defined
    max_interaction= 10*maxn*(ceiling((hsml_const/dx_r)*scale_k*2))**SPH_dim
    max_e_interaction= 10*maxedge*(ceiling((hsml_const/dx_r)*scale_k*2))**SPH_dim 

    write(*,*) "maximum particle-particle and particle-point interaction defined : ",max_interaction
    write(*,*) "maximum particle-edge_temp interaction defined : ",max_e_interaction
    write(*,*)'**************************************************'   
    
    ! find size of the simulation
    ALLOCATE( simGridSize(SPH_dim,2))
    simGridSize(:,1)= MINVAL(x,2) 
    simGridSize(:,2)= MAXVAL(x,2) 
    write(*,*)'**************************************************'   
    write(*,*) "simulation domain size is ", simGridSize, "with xmin = ", simGridSize(1,1), "xmax=",simGridSize(1,2)
    write(*,*) "with ymin = ", simGridSize(2,1), "ymax=",simGridSize(2,2)
    write(*,*)'**************************************************'   
    
    do k=1,ntotal
        !For changing density to hydrostatic
        rho(k)=rho_init !*((7.D0*9.81D0*(hydroStaticHeight-x(2,k))/(c_sound**2)+1.D0)**(1.D0/7.D0))
        p(k)=(rho_init*c_sound**2/7.D0)*((rho(k)/rho_init)**7.D0-1.D0) 
        mass(k)= rho(k)*vol(k)
        hsml(k)=hsml_const
        mu(k)=mu_const
        p(k)=0.D0
        temp(k)=300
    enddo
    
    
    end