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
        
    use config_parameter, only: SPH_dim, DataConfigPath, &
        & dx_r, hsml_const, edge_to_dx_ratio
    use particle_data ,   only: maxn, max_interaction, max_e_interaction, maxnv,  &
        & x,vol,itype,hsml, nreal, nedge, ntotal, edge, &
        & surf_norm, etype, etotal, maxedge, simGridSize, &
        & x_ve, integ_pts_for_edge,val_integ_pts, mid_pt_for_edge
    !use ifport, only:random
    
    implicit none
    real(8) tempType
    integer(4) k,d,s,ke,k_int_pt,i
    real(8)  mn, scale_k
    integer(4) nreal_mesh, nrealCartesian
    real(8) totVol, x_ve_temp(SPH_dim,SPH_dim), len_bulkBdry, len_domainBdry
    character(40) :: input_file_name
    
    ! Retreive the scaling factor for the smoothing length used in the smoothing kernel    
    call sml_mult_factor(scale_k)
    
    ! Read size of the domain from config file 
    write(*,*)'  **************************************************'
    write(*,*)'       Reading size of the SPH domain configuration...   '
    open(1,file= DataConfigPath // '/input_param_SimSize.dat',status='old')
    do while (.not.eof(1))
        ! I might have to modify below to strictly read as integers
          read(1,*) nreal_mesh, nrealCartesian, len_bulkBdry, len_domainBdry, totVol
    enddo
    close(1)
    
    
    ! maximum number of edges and particles
    maxedge=int((len_bulkBdry/dx_r)*(dble(edge_to_dx_ratio)*5.D0+1.D0)) ! *2.D0 is used to account for the 2X vertices of issue #48
    maxnv=10*maxedge 
    maxn= max(nreal_mesh, nrealCartesian)+ maxedge
    
    !variables required to define geometry are particle propoerties are allocated
    ALLOCATE(x(SPH_dim,maxn))
    ALLOCATE( vol(maxn), hsml(maxn), itype(maxn))
    
    ! The allocatable variables are now initialized
    x=0.D0
    vol=0.D0
    hsml=0.D0
    itype=0
    
    
    ! initialize particle number to 0  
    k=0
        
    ! read approriate data file according to input type
    open(1,file= DataConfigPath // '/input_particles_cartesianMesh.dat',status='old')
    
    do while (.not.eof(1))
        k=k+1
        read(1,*) tempType, (x(d, k), d=1,SPH_dim), vol(k)  
        itype(k) = int(tempType)
    enddo

    nreal=k ! I can distinguish nreal and nedge here if needed
    write(*,*)'      Total number of Real particles : ', nreal    	
    write(*,*)'**************************************************'

    ! initialize edges to 0
    s=0
    
    ! read approriate data file according to input type
    input_file_name = '/input_bulkBdryEdge.dat'

    ! read bdry data file and store bdry information
    call inputCADtoEdgeData(s, maxedge, input_file_name, edge_to_dx_ratio)
    
    etotal=s
    write(*,*)'      Total number of boundary segments : ', etotal    	
    write(*,*)'**************************************************'      


    !val_integ_pts initialized using maximum pts per edge times the number of edges
    ! here we use one itnegration pt per edge, hence 1*maxedge
    ! integ_pts_for_edge initialized using the maximum number of edges
    ALLOCATE(integ_pts_for_edge(2,maxedge), val_integ_pts(SPH_dim+1,1*maxedge), mid_pt_for_edge(SPH_dim,maxedge))
    
    ! Now we will find the center of edges, which in Method 1 of Parikshit et. al will work as virtual points
    nedge=0
    
    !initialize kth integration pt to zero
    k_int_pt=0
    do s=1,etotal
               
        do d =1,SPH_dim
            x_ve_temp(:,d)=x_ve(:,edge(d,s))
        enddo      
        
    ! Now define number of integration pts per edge for the below lines. 
    ! This would need to be in some loop for more than one integration pt per edge    
        
        k_int_pt= k_int_pt+1 !update the kth integration point globally numbered
        
        integ_pts_for_edge(1,s)= k_int_pt ! Starting integration pt for a given edge
        
        integ_pts_for_edge(2,s) = 1 ! Total number of integration pts for a given edge
        
        ! Determine the location of integration point
        call centroidBdrySegment(val_integ_pts(1:SPH_dim, k_int_pt), x_ve_temp, SPH_dim)
        
        ! Determine the integration weight assosciated with the integration point
        val_integ_pts(SPH_dim+1, k_int_pt)= norm2(x_ve(:,edge(1,s))-x_ve(:,edge(2,s))) ! in this 2D case, it simply length of 1D bdry edge
       
        
    ! Determine the location of edge mid point
        call centroidBdrySegment(mid_pt_for_edge(:, s), x_ve_temp, SPH_dim)

    enddo
    
    nedge = k_int_pt
    
    ! Add any other particles/ points that would be used in particle particle interatcionts
    ntotal=nreal ! + nvirtual_points
    
    ! Maximum interactions for particle-particle and edge_temp-particle is defined
    max_interaction= edge_to_dx_ratio*maxn*(ceiling((hsml_const/dx_r)*scale_k*2))**SPH_dim
    max_e_interaction= edge_to_dx_ratio*maxedge*(ceiling((hsml_const/dx_r)*scale_k*2))**SPH_dim 

    write(*,*) "maximum particle-particle and particle-point interaction defined : ",max_interaction
    write(*,*) "maximum particle-edge_temp interaction defined : ",max_e_interaction
    write(*,*)'**************************************************'   
    
    ! find size of the simulation
    ALLOCATE( simGridSize(SPH_dim,2))
    simGridSize(:,1)= min(MINVAL(x_ve,2),MINVAL(x,2) ) -dx_r
    simGridSize(:,2)= max(MAXVAL(x_ve,2), MAXVAL(x,2)) +dx_r
    write(*,*)'**************************************************'   
    write(*,*) "simulation domain size is ", simGridSize, "with xmin = ", simGridSize(1,1), "xmax=",simGridSize(1,2)
    write(*,*) "with ymin = ", simGridSize(2,1), "ymax=",simGridSize(2,2)
    write(*,*)'**************************************************'   
    
    do k=1,ntotal
        !For changing density to hydrostatic
        hsml(k)=hsml_const
    enddo
    
    !call OutputSetup
    !needed i=only to debug

    ! Now perform the particle packing algorithm 
     call particlePackingTimeIntegration  
   
    deallocate(simGridSize, x, vol, itype)
    deallocate(x_ve,etype, edge, surf_norm)
    deallocate(integ_pts_for_edge, val_integ_pts, mid_pt_for_edge)

    DEALLOCATE(hsml)
    
    end