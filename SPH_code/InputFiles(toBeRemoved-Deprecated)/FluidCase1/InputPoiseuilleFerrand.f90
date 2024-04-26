!****************************************************************************
!
!  SUBROUTINE: input
!
!  PURPOSE:  Subroutine that creates real particles and edges for 
!           a square geometry [0,1,0,1]  here bottom boundary have dirichlet boundary
!           Right and Left have neuman boundary (as adiabatic), Top is neuman with heat flux.
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/11/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

    subroutine input
        
    use config_parameter, only: SPH_dim,etype_periodic, etype_SolidWall1, etype_SolidWall2, &
        & etype_FreeSurface1, etype_thermal_neumann, etype_thermal_dirichlet
    use config_geometry, only: npw, npl, dx_r, bt_dx_r, length, width, hsml_const, &
        & rho_init, mu_const, g_const, c_sound, F_ext, &
        & BILtype  
    use particle_data ,   only: maxn, max_interaction, max_e_interaction, maxnv,  &
        & x,vx,mass,rho,p,itype,hsml,mu, temp, nreal, nflow, nedge,nghost, ntotal, edge, &
        & surf_norm, nedge_rel_edge, etype, etotal, maxedge, &
        & pBC_edges,tangent_pBC, bdryVal_vel, bdryVal_prs, bdryVal_rho, bdryVal_temp, &
        & WallBoundaryLayer, f0max , dynamicProblem
    use ifport, only:random
    
    implicit none
    
    integer(4) i,j,s, k,ke ,e_divide
    real(8) scale_k, mn, length1, width1, ecf 
    real(8) ke_out, keTop,keBottom, keRight, keLeft!case specific width
    integer(4) npa
    real(8),DIMENSION(:,:),ALLOCATABLE :: xn
    
    write(*,*)'  ***************************************************'
    write(*,*) 'BIL case number'
    write(*,*)'  ***************************************************'
    read (*,*) BILtype
    
    dynamicProblem=.true.
    
    call sml_mult_factor(scale_k)
    
    e_divide=10
    
    ! the below is max number of particles for a rectangle with width*length real domain
    maxn= (npw*npl) + (npw+npl+4)*2*(3)*e_divide ! 3 for 3 points/particles per edge
    maxnv= (npw+npl+4)*(3)*e_divide ! 3 for 3 points/particles per edge 
    
    
    !vaiuables required to define geometry are particle propoerties are allocated
    ALLOCATE(x(SPH_dim,maxn),vx(SPH_dim,maxn))
    ALLOCATE(mass(maxn), rho(maxn), p(maxn), hsml(maxn), itype(maxn), mu(maxn), temp(maxn))
    
    
    !the allocatable variables are now initialized
    x=0.D0
    vx=0.D0
    mass=0.D0
    rho=0.D0
    p=0.D0
    hsml=0.D0
    itype=0
    mu=0
    temp=0.D0
    
    f0max = 1.D0
    
    !Define External Force
    allocate(F_ext(SPH_dim))
    F_ext =0.D0
    F_ext(1)= 8.D-1
    
    !Initialize the local variable
    k=0
    
    !Equation of state factor for 
    ecf=7.D0
    
    !WallBoundaryLayer =.true.
    
    write(*,*) "npw =" , npw, "npl =", npl, "dx_r=", dx_r 
    
    !Create Real Particles
    do i=1,npl
        do j=1,npw
            k=k+1
            x(1,k)= dble(i-1)*dx_r + dx_r/2.D0
            x(2,k)= dx_r*dble(j-1) + dx_r/2.D0                                      
            vx(1,k)=0.D0 
            vx(2,k)=0.D0 
            rho(k)=rho_init!*((ecf*9.81D0*(width-x(2,k))/(c_sound**2)+1.D0)**(1.D0/ecf))
            mass(k)= rho(k)*dx_r*dx_r
            p(k)=0.
            itype(k)=3
            hsml(k)=hsml_const 
            mu(k)=mu_const
            temp(k)= 750+273
    
        enddo 
    enddo
    
    ! This will define length and width of outside edges
    length1= x(1,k) + dx_r/2.D0
    width1 = x(2,k) + dx_r/2.D0 
    
    
! Create particles with no mass to capture value
    ! do j=1,npw
    !    k=k+1
    !    x(1,k)= length1/2.D0 
    !    x(2,k)= dx_r*dble(j-1) + dx_r/2.D0 
    !    vx(1,k)=0.D0 
    !    vx(2,k)=0.D0 
    !    rho(k)=rho_init
    !    mass(k)= 0.D0
    !    p(k)=0.D0
    !    itype(k)=1
    !    hsml(k)=hsml_const 
    !    mu(k)=mu_const
    !enddo
   
    !   The total real particles is stored and then printed on screen    
    nreal=k
    write(*,*)'      Total number of Real particles   ', nreal    	
    write(*,*)'**************************************************'
    
    
    !   The flow particles if applicable are created here
    
    
    !   The total flow particles is stored and then printed on screen     
    nflow=0
    write(*,*)'      Total number of Flow particles   ', nflow    	
    write(*,*)'**************************************************'
    
    
    !   Reference nodes to create edge particles
    ALLOCATE(xn(SPH_dim,maxnv))
    ke=0
    
    
     !   Outside Bottom edge Layer created below
    do i=1,e_divide*npl
        ke=ke+1
        xn(1,ke)=dble(i-1)*dx_r/dble(e_divide) 
        xn(2,ke)=0.D0 
    enddo
    
    keBottom=ke
    
    ! Outside  Right edge corner
    !do i=1,e_divide*npw
        ke=ke+1
        xn(1,ke)=length1 
        xn(2,ke)=0.D0!dble(i-1)*dx_r/dble(e_divide) 
    !enddo
    
    keRight=ke
    
    !  Outside Top edge Layer created below
    do i=1,e_divide*npl
        ke=ke+1
        xn(1,ke)=length1 - dble(i-1)*dx_r/dble(e_divide) 
        xn(2,ke)=width1 
    enddo
    
    keTop=ke
    
    ! Outside Left edge corner
    !do i=1,e_divide*npw
        ke=ke+1
        xn(1,ke)=0.D0 
        xn(2,ke)=width1 !width1 - dble(i-1)*dx_r/dble(e_divide)
    !enddo

    keLeft=ke    
    
    !Store the outside total edge number     
    ke_out=ke
    
   

    ! We create edges using the vertices of an edge. This will be useful to calculate boundary integral and defining surface normal.
    ! we find the total edges present in the intial configuration
    etotal=ke
    
    maxedge= int(2*etotal)
    
    Allocate(edge(SPH_dim,maxedge))   
    Allocate(etype(maxedge))
    etype=0
    do i=1,ke
        !for 2D we only have 2 vertices per edge, but an edge in 3d is a plane,
        ! which is atleast a triangle
        edge(1,i) = i  
        edge(2,i) = i+1

        if(i .eq. ke_out) then
            edge(2,i)=edge(1,1)
        endif
        
        if(i .le. keBottom) etype(i)=etype_SolidWall1!etype_thermal_dirichlet
            
        if((i .gt. keBottom) .and. (i .le. keRight)) etype(i)=etype_periodic!etype_thermal_neumann    
        
        if((i .gt. keRight) .and. (i .le. keTop)) etype(i)=etype_SolidWall1!etype_thermal_neumann
        
        if((i .gt. keTop) .and. (i .le. keLeft)) etype(i)=etype_periodic!etype_thermal_neumann    
            

        
    enddo 
    
    ! Define surface normals of walls, with normal per edge
    ! ** this can be made a function/seperate subroutine
    ALLOCATE(surf_norm(SPH_dim, maxedge))
    do s = 1, etotal   
        surf_norm(1,s)= -(xn(2,edge(2,s))- xn(2,edge(1,s)))
        surf_norm(2,s)= (xn(1,edge(2,s))- xn(1,edge(1,s)))
        mn= norm2(surf_norm(:,s))
        surf_norm(1,s)= surf_norm(1,s)/mn
        surf_norm(2,s)= surf_norm(2,s)/mn
    enddo 
    
    ! Initialize the variable nedge_rel_edge, that will store edges related to a boundary particle
    ALLOCATE(nedge_rel_edge(maxedge))
    
    ! Now we will find the center of edges, which in Method 1 of Parikshit et. al will work as virtual points
    ! In method 2 it will be used as a real boundary particle 
    ! The below needs to be modified slightly to include more than one edge particle per edge (and to include different dx_v
    do s=1,etotal
        k=k+1
        nedge_rel_edge(s)=k ! Here nedge_rel_edge(i)=i  is a simple case of one edge particle per edge
        x(:,k)= (xn(:,edge(1,s))+ xn(:,edge(2,s)))/2.D0 !this defines the location of edge points/particles we need per edge
        vx(:,k)=0.D0

        itype(k)=101!This will be 100 if these are to be treated as virtual points in simulation  
        hsml(k)=hsml_const
        
        ! For Method 1 boundary poitns, they have no mass, and the below are really interpolated values,
        ! or imposed boundary values. They are not particles, but just points
        ! For particles, mass below is updated
        rho(k)= rho_init
        mass(k)=0.D0                                                    
        p(k)=0.D0
        
        ! mu needs to be carefully incorporated in 
        ! multiphase flows while using certian algorithm 
        !like ferrands as mu of wall changes depending 
        !on particle in it viccinity        
        mu(k)=mu_const 
    enddo  
        
    ! The total edge particles is stored and then printed on screen
    nedge=k-(nreal+nflow)      
    write(*,*)'      Total number of Edge particles   ', nedge    	
    write(*,*)'**************************************************'      
    
    !Let us define the edge vertices by ghost points to visually demonstrate edges on tecplot
    !This also lets us uniquely label ghost points in a global sense
    do i=1,ke
        k=k+1
        x(:,k)=xn(:,i)  
        itype(k)=0
    enddo
    
    ! The total ghost points is stored and then printed on screen
    nghost=k-(nreal+nflow+nedge)      
    write(*,*)'      Total number of Vertex Points for edges / ghost points ', nghost    	
    write(*,*)'**************************************************'      
    
    ! The edge vertices need to also be globally labelled
    do i=1,ke
        !for 2D we only have 2 vertices per edge, but an edge in 3d is a plane,
        ! which is atleast a triangle
        edge(1,i) = i + (nreal+nflow+nedge)  
        edge(2,i) = i+1 + (nreal+nflow+nedge) 
        
        if(i .eq. ke_out) then
            edge(2,i)=edge(1,1)
        endif
        
    enddo 
    

    ! Select periodic edge pairs
    Allocate(pBC_edges(2,1))    !Here we have just one pair of periodic edges
    pBC_edges(1,1)=int(ke/2)
    pBC_edges(2,1)=ke
    
    ! Determine the tangential sense of the periodic edges
    Allocate(tangent_pBC(SPH_dim,2,1))
    tangent_pBC(1,1,1)= 0.D0 !x(1,edge(2,pBC_edges(1,1)))- x(1,edge(1,pBC_edges(1,1)))
    tangent_pBC(2,1,1)= 1.D0 !x(2,edge(2,pBC_edges(1,1)))- x(2,edge(1,pBC_edges(1,1)))
    !mn= norm2(tangent_pBC(:,1,1))
    !tangent_pBC(:,1,1)=tangent_pBC(:,1,1)/mn
    
    tangent_pBC(1,2,1)= 0.D0 !x(1,edge(1,pBC_edges(2,1)))- x(1,edge(2,pBC_edges(2,1)))
    tangent_pBC(2,2,1)= 1.D0 !x(2,edge(1,pBC_edges(2,1)))- x(2,edge(2,pBC_edges(2,1)))
    !mn= norm2(tangent_pBC(:,2,1))
    !tangent_pBC(:,2,1)=tangent_pBC(:,1,1)/mn
    
    
    
    ! The total number of all particle types and point representations  are printed on screen
    ntotal=nreal+nflow+nedge+nghost
    write(*,*)'      Total number of Edge particles   ', ntotal    	
    write(*,*)'**************************************************'  
    
  
    
! Apply boundary values
    allocate( bdryVal_vel(SPH_dim,etotal),bdryVal_prs(etotal), bdryVal_rho(etotal), bdryVal_temp(etotal))
    
    do s=1,etotal
        
         ! Below gives Dirichlet boundary condition value
        if(etype(s) .eq. etype_thermal_dirichlet) then
            
            !if(s .le. keBottom)     bdryVal_temp(s)=1000+273 
            
            !if((s .gt. keBottom) .and. (s .le. keRight))     
        
            !if((s .gt. keRight) .and. (s .le. keTop))   bdryVal_temp(s)=500+273 
        
            !if((s .gt. keTop) .and. (s .le. keLeft))  
       
        endif
        
        ! Below gives Neumann boundary condition value along normal
        if(etype(s) .eq. etype_thermal_neumann) then
            
            !if(s .le. keBottom) 
            
            !if((s .gt. keBottom) .and. (s .le. keRight))    bdryVal_temp(s)= 0.D0    
        
            !if((s .gt. keRight) .and. (s .le. keTop))       bdryVal_temp(s)= 500+273
        
            !if((s .gt. keTop) .and. (s .le. keLeft))    bdryVal_temp(s)= 0.D0 

        endif
        
        ! Below gives Dirichlet boundary for Velocity(and temp?) and Neumann Boundary for Pressure and Density
        if(etype(s) .eq. etype_SolidWall1) then
            ! here v, dp/dn, drho/dn are calculated
            !bdryVal_vel(s)=0.D0
            !bdryVal_prs(s)=0.D0
            !bdryVal_rho(s)=0.D0
            if(s .le. keBottom) then
                bdryVal_vel(:,s)=0.D0
                bdryVal_prs(s)=0.D0
                bdryVal_rho(s)=0.D0
                bdryVal_temp(s)=273
            
            elseif((s .gt. keBottom) .and. (s .le. keRight)) then  
                bdryVal_vel(:,s)=0.D0
                bdryVal_prs(s)=0.D0
                bdryVal_rho(s)=0.D0
                bdryVal_temp(s)= 273  
        
            elseif((s .gt. keRight) .and. (s .le. keTop)) then
                bdryVal_vel(:,s)=0.D0 !1.D-3
                bdryVal_prs(s)=0.D0
                bdryVal_rho(s)=0.D0
                bdryVal_temp(s)= 273
        
            elseif((s .gt. keTop) .and. (s .le. keLeft)) then
                bdryVal_vel(:,s)=0.D0
                bdryVal_prs(s)=0.D0
                bdryVal_rho(s)=0.D0
                bdryVal_temp(s)= 273
            endif            
            
        endif
        
        
        if(etype(s) .eq. etype_FreeSurface1) then
            ! here dv/dn, p, rho are calulated
            !bdryVal_vel(s)=0.D0
            !bdryVal_prs(s)=0.D0
            !bdryVal_rho(s)=rho_init
        endif
        
    enddo
    
    
    
    ! Maximum interactions for particle-particle and edge-particle is defined
    max_interaction= maxn*(ceiling((hsml_const/dx_r)*scale_k*2))**SPH_dim*e_divide 
    max_e_interaction= ke*(ceiling((hsml_const/dx_r)*scale_k*2))**SPH_dim*e_divide 
    
    write(*,*) "max_interaction",max_interaction, "max_e_interaction",max_e_interaction
    
    ! randomise the particles a bit
    !do  k=1,nreal
    !    if(itype(k) .ne. 1) then
    !        x(1,k)= x(1,k) + (random(0)-0.5D0)*dx_r/2.D0  
    !        x(2,k)= x(2,k) + (random(0)-0.5D0)*dx_r/2.D0 
    !    endif
    !enddo
    !
        
    deallocate(xn)
    
    end