!****************************************************************************
!
!  SUBROUTINE: input
!
!  PURPOSE:  Subroutine that accepts any geometry input
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  11/08/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

    subroutine inputParticlePacking
        
    use config_parameter, only: SPH_dim, &
        & DataConfigPath, itype_virtual, itype_real_min, &
        & dx_r, hsml_const, &
        & rho_init, packagingIterations, ExtInputMeshType
    use particle_data ,   only: maxn, max_interaction, max_e_interaction, maxnv,  &
        & x,mass,rho, vol,itype,hsml, nreal, nflow, nedge,nghost, ntotal, edge, &
        & surf_norm, nedge_rel_edge, etype, etotal, maxedge, &
        & avgVol, packed_x,packed_itype, packed_vol, simGridSize, xstart
    !use ifport, only:random
    
    implicit none
    real(8) tempType
    integer(4) k,d,s,ke,i
    real(8) tempX(SPH_dim,SPH_dim), tempBdry(5), mn, scale_k
    integer(4) nreal_mesh, nrealCartesian, nEbulkBdry, nEdomainBdry
    real(8) totVol
    real(8),DIMENSION(:,:),ALLOCATABLE :: xV_temp, edge_temp
    character (40) :: x2name
           
    
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
    
    maxedge= nEbulkBdry
    maxnv=3 *maxedge
    if (ExtInputMeshType .le. 2) maxn= nrealCartesian +maxnv
    if (ExtInputMeshType .ge. 3) maxn= nreal_mesh +maxnv
    
    !vaiuables required to define geometry are particle propoerties are allocated
    ALLOCATE(x(SPH_dim,maxn))
    ALLOCATE(mass(maxn), rho(maxn), vol(maxn), itype(maxn), hsml(maxn))
    
    ! The allocatable variables are now initialized
    x=0.D0
    mass=0.D0
    rho=0.D0
    itype=0
    vol=0.D0
    hsml=0.D0
    simGridSize=0.D0
    
    ! initialize particle number to 0  
    k=0
        
    ! read approriate data file according to inut type
    call inputCADtoParticleData(k, nreal_mesh,totVol)
    
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
    open(1,file= DataConfigPath // '/input_bulkBdryEdge.dat',status='old')
    ! initialize particle number to 0    
    s=0
    ke=0
    do while (.not.eof(1))
        s=s+1
        read(1,*) tempType, (tempX(d,1), d=1,SPH_dim), (tempX(d,2), d=1,SPH_dim) !, (tempbdry(d), d=1,5)
         
       etype(s)=NINT(tempType)
        
        do d=1,SPH_dim
            ke=ke+1
            xV_temp(:,ke)=tempX(:,d)
            edge_temp(d,s)=ke
        enddo
         
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
                
        itype(k)= itype_virtual+itype_real_min
 
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

    write(*,*) "maximum particle-particle and particle-point interaction defined for particle initialization algorithm: ",max_interaction
    write(*,*) "maximum particle-edge_temp interaction defined for particle initialization algorithm: ",max_e_interaction
    write(*,*)'**************************************************'   
    
    do k=1,ntotal
        !For changing density to hydrostatic
        rho(k)=rho_init
        mass(k)=rho(k)*vol(k)
        hsml(k)=hsml_const
    enddo
    
    
    ! find size of the simulation
    ALLOCATE( simGridSize(SPH_dim,2))
    
    simGridSize(:,1)= MINVAL(x,2) 
    simGridSize(:,2)= MAXVAL(x,2) 

    write(*,*)'**************************************************'   
    write(*,*) "particle packing grid size is ", simGridSize, "with xmin = ", simGridSize(1,1), "xmax=",simGridSize(1,2)
    write(*,*) "with ymin = ", simGridSize(2,1), "ymax=",simGridSize(2,2)
    write(*,*)'**************************************************'   
    
    Allocate(xStart(SPH_dim,ntotal))
    xStart=x        

    ! Now perform the particle packing algorithm 
    call particlePackingTimeIntegration(.true.)
    !input : quick_converge_step2C
    
    deallocate(xstart)
    
    write(x2name, '(A,A)') DataConfigPath,'/input_PP.dat'
    open (1, file = x2name)
    

    ! write particle data to the file input_PP.dat
    do i=1,nreal
        write(1,'(I10,4(e22.10,2x))') itype(i), (x(d, i), d=1,SPH_dim), vol(i)
    enddo
    
    close(1)
        
    
    deallocate(x, itype, vol,mass, rho, etype, surf_norm, nedge_rel_edge, edge, hsml, simGridSize)

    end