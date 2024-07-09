subroutine inputCADtoEdgeData(s, max_edge, input_file_name, avgEdgeSize_ratio )
    use config_parameter, only : SPH_dim, ExtInputMeshType, DataConfigPath, dx_r, etype_periodic, periodic_pairs, packagingIterations
    use particle_data, only : x_ve, etype, edge, surf_norm, pBC_edges, ve_total, maxnv
    
    implicit none
    integer(4), intent(inout) :: s, max_edge
    character (40), intent(in) :: input_file_name
    integer(4), intent(in) :: avgEdgeSize_ratio
    real(8) :: tempType, ve_point_ID(SPH_dim), tempX(SPH_dim,SPH_dim), dx_edge_reqd, max_edge_size, temp_surf_norm(SPH_dim)
    integer(4) :: d, num_divisions, i, num_vertex_pts, temp_edgetype_int
    real(8),allocatable :: divided_tempX(:,:,:)
    
    
    ! find minimum edge size, by using the preferred rtio of dx_r/dx_edge
    dx_edge_reqd = dx_r/dble(avgEdgeSize_ratio)
    
    Allocate(x_ve(SPH_dim,maxnv),etype(max_edge), edge(SPH_dim, max_edge), surf_norm(SPH_dim, max_edge))
    num_vertex_pts=0
    x_ve=0
    etype=0
    edge =0
    
    ! Read input file containing vertices of edge
    open(1,file= DataConfigPath // input_file_name ,status='old')
    
    do while (.not. eof(1))
        
        
        
        read(1,*) tempType, (tempX(d,1), d=1,SPH_dim), (tempX(d,2), d=1,SPH_dim) !, (temp_surf_norm(d), d=1,SPH_dim) ! this can be generalized with x_temp(d1,d2)
        
        !obtain the surface normal for the bdry type from the given order of vertices
        ! the below wouldnt be required if surf_norm is read directly from bdry file
        call surfaceNormalBdry(temp_surf_norm, tempX, SPH_dim)
        
        !find edgetype 
        call BCinput(temp_edgetype_int,NINT(tempType), packagingIterations) 
        
        
        if(temp_edgetype_int .eq. etype_periodic) then ! this is to treat general open type boundaries like periodic, inlet/outlet flow etc.
                
            ! each division above is a boundary element
            s=s+1
            
            ! make sure the array size is big enough
            if(max_edge .eq. num_vertex_pts) then
                write(*,*) "Increase value of max edge and rerun simulation"
                write(*,*) "current maxedge =", max_edge 
                write(*,*) "current num_vertex_pts =", num_vertex_pts 
                write(*,*) "Paused in the conditional statement if(temp_edgetype_int .eq. etype_periodic) in inputCADtoEdgeData"
                pause
            endif
                
            do d= 1,SPH_dim
                
                num_vertex_pts=num_vertex_pts+1
                
                !the below can be optimized to avoid duplicating the vertices
                x_ve(:,num_vertex_pts)=tempX(:,d)
                
                ! assign the vertex numbers to the edge counterparts
                edge(d,s)= num_vertex_pts      
            enddo
            
            ! assign the same surfacenormal to the public variable
            surf_norm(:,s) = temp_surf_norm(:)
            
            ! assign the edge type to the public variable
            etype(s) = temp_edgetype_int
                

            !the below needs to be generalized for more than one periodic pairs
            if( .not. allocated(pBC_edges)) then                       
                allocate(pBC_edges(2,periodic_pairs))
                pBC_edges(1,1)=s
            else
                pBC_edges(2,1)=s
            endif
            
            
        else ! for non open type bc we need to see if the bdry element can be divided further
            
            ! check if the distance between any two vertices representing the edge too big.
            call MaxEdgeSize(max_edge_size,tempX,SPH_dim)
        
            ! check if the current maximum edge size connecting any two vertices is bigger than
            ! required edge size.
            ! we use a tolerance of 0.01%, ie. (max_edge_size - dx_edge_reqd)/max_edge_size = 10^-2
            if(max_edge_size .gt. dx_edge_reqd/(1.D0-1.D-2)) then            
                
                ! Find how many divisions of hte element are required
                if(SPH_dim .eq. 2) then
                    num_divisions=ceiling(max_edge_size/dx_edge_reqd)
                    allocate(divided_tempX(SPH_dim,SPH_dim,num_divisions))

                elseif(SPH_dim .eq. 3) then
                    write(*,*) "need to add for 3D"
                    pause
                endif
                
                ! divide a boundary element into multiple elements to ensure they are the size of the bdry elements
                call divideBdryElement(divided_tempX,tempX, num_divisions, SPH_dim)        
            
            else
                ! no divisions needed if size of elements are small enough
                num_divisions = 1                
                allocate(divided_tempX(SPH_dim,SPH_dim, num_divisions))                
                divided_tempX(:,:,num_divisions)=tempX(:,:)
            endif
            
            ! go through all boundary elements (if divided) assosciated with the read element
            do i = 1, num_divisions                 
                ! each division above is a boudanry element
                s=s+1

                ! make sure the array size is big enough
                if(max_edge .eq. num_vertex_pts) then
                    write(*,*) "Increase value of max edge and rerun simulation"
                    write(*,*) "current maxedge =", max_edge 
                    write(*,*) "current num_vertex_pts =", num_vertex_pts 
                    write(*,*) "Paused in do i = 1, num_divisions loop inside inputCADtoEdgeData"
                    pause
                endif
                
                do d= 1,SPH_dim
                    !the below can be optimized to avoid duplicating the vertices
                    num_vertex_pts=num_vertex_pts+1
                    x_ve(:,num_vertex_pts)=divided_tempX(:,d,i)    
                
                    ! assign the vertex numbers to the edge counterparts
                    edge(d,s)= num_vertex_pts      
                enddo
            
                ! assign the same surfacenorm to the divisions
                surf_norm(:,s) = temp_surf_norm(:)
            
                ! assign the same edgetype to the divisions
                etype(s) = temp_edgetype_int
                
            enddo
            
            deallocate(divided_tempX)
        endif
    enddo

    ve_total = num_vertex_pts
        
        !(ve_point_ID(d), d=1,SPH_dim) ! not used yet, but needs be implemented in next update
        
        
endsubroutine
        
    
subroutine MaxEdgeSize(edge_size, tempX,dim)
!This checks for the maximum length of any two vertices
    implicit none
    
    integer(4), intent(in) :: dim
    real(8), intent(in) :: tempx(dim,dim)
    real(8), intent(out) :: edge_size
    real(8):: edge_len(dim)
    integer(4) :: d
    
    edge_size=0
        
    if(dim .eq. 2) then
        edge_size=norm2(tempX(:,1)-tempX(:,2))
    
    elseif(dim .eq. 3) then
        do d= 1,dim
            edge_len(d)=norm2(tempX(:,d)-tempX(:,mod(d, dim) + 1))
            edge_size = max(edge_size,edge_len(d))
        enddo
    endif
    
endsubroutine



subroutine divideBdryElement(divided_tempX,tempX, num_divisions, dim)
! This divides a bdry element into additional elements
    implicit none
    
    integer(4), intent(in) :: num_divisions
    integer(4), intent(in) :: dim
    real(8), intent(out) :: divided_tempX(dim,dim, num_divisions)
    real(8), intent(in) :: tempx(dim,dim)    
    real(8) :: vector_dir_2D(2)
    integer(4):: i
    
    
    if(dim .eq. 2) then
        vector_dir_2D(:)= (tempX(:,2)-tempX(:,1))/num_divisions

        divided_tempX(:,1,1) = tempX(:,1) 
        divided_tempX(:,2,1) = tempX(:,1)+ vector_dir_2D(:) 
        do i = 2,num_divisions-1
            divided_tempX(:,1,i) = divided_tempX(:,2,i-1)
            divided_tempX(:,2,i) = divided_tempX(:,2,i-1) + vector_dir_2D(:)
        enddo
        divided_tempX(:,1,num_divisions) = divided_tempX(:,2,num_divisions-1)
        divided_tempX(:,2,num_divisions) = tempX(:,2)       
        
        
        
    elseif(dim .eq. 3) then
        
    endif
    
endsubroutine
