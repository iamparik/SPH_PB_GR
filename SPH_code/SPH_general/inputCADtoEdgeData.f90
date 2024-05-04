subroutine inputCADtoEdgeData(s, max_edge, input_file_name, avgEdgeSize_ratio )
    use config_parameter, only : SPH_dim, ExtInputMeshType, DataConfigPath, dx_r
    use particle_data, only : x_ve, etype
    
    implicit none
    integer(4), intent(inout) :: s, max_edge
    character (40), intent(in) :: input_file_name
    integer(4), intent(in) :: avgEdgeSize_ratio
    real(8) :: tempType, ve_point_ID(SPH_dim), tempX(SPH_dim,SPH_dim), dx_edge_reqd, max_edge_size, temp_surf_norm(SPH_dim)
    integer(4) :: num_divisions, i, num_vertex_pts, temp_edgetype_int
    real(8), dimension(:,:,:), allocatable :: divided_tempX
    
    ! start with zero edges
    s=0
    
    ! find minimum edge size
    dx_edge_reqd = dx_r/dble(input_file_name)
    
    Allocate(x_ve(SPH_dim,max_edge),etype(max_edge), edge(SPH_dim, max_edge))
    num_vertex_pts=0
    x_ve=0
    etype=0
    edge =0
    
    ! Read input file containing vertices of edge
    open(1,file= DataConfigPath // input_file_name ,status='old')
    
    do while (.not. eof(1))
        read(1,*) tempType, (tempX(d,1), d=1,SPH_dim), (tempX(d,2), d=1,SPH_dim), (temp_surf_norm(d), d=1,SPH_dim) ! this can be generalized with x_temp(d1,d2)
        
        !find edgetype 
        call BCinput(temp_edgetype_int,NINT(tempType)) 
        
        
        if(temp_edgetype_int .eq. etype_periodic) then ! this is to treat general open type boundaries like periodic, inlet/outlet flow etc.
            
            
            
        else ! for non open type bc we need to see if the bdry element can be divided further
            
            ! check if the distance between any two vertices representing the edge too big.
            call MaxEdgeSize(max_edge_size,tempX,SPH_dm)
        
            if(max_edge_size .gt. dx_edge_reqd) then            
                ! divide a boundary element into multiple elements to ensure they are the size of the bdry elements
                call divideBdryElement(num_divisions,divided_tempX,tempX, max_edge_size,dx_edge_reqd, SPH_dim)        
            
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
        
        
        
        
        !(ve_point_ID(d), d=1,SPH_dim) ! not used yet, but needs be implemented in next update
        
        
endsubroutine
        
    
subroutine MaxEdgeSize(edge_size, tempX,dim)
!This checks for the maximum length of any two vertices
    implicit none
    
    integer(4), intent(in) :: dim
    real(8), intent(in) :: tempx(dim,dim)
    real(8), itnent(out) :: edge_size
    real(8):: edge_len(dim)
    integer(4) :: d
    
    edge_size=0
        
    if(dim .eq. 2) then
        edge_size=norm2(tempX(:,1),tempX(:,2))
    
    elseif(dim .eq. 3) then
        do d= 1,dim
            edge_len(d)=norm2(tempX(:,d),tempX(:,mod(d, dim) + 1))
            edge_size = max(edge_size,edge_len(d))
        enddo
    endif
    
endsubroutine



subroutine divideBdryElement(num_divisions,divided_tempX,tempX, max_edge_size,dx_edge_reqd, dim)
! This divides a bdry element into additional elements
    implicit none
    
    integer(4), intent(inout) :: num_divisions
    real(8), intent(inout), dimension (:,:,:), allocatable :: divided_tempX
    real(8), intent(in) :: tempx(dim,dim), max_edge_size, dx_edge_reqd
    integer(4), intent(in) :: dim
    real(8) :: vector_dir_2D(2)
    integer(4):: i
    
    
    if(dim .eq. 2) then
        !Find how many divisions are needed
        num_divisions=ceiling(max_edge_size/dx_edge_reqd)

        vector_dir_2D(:)= (tempX(:,2)-tempX(:,1))/num_divisions
        
        allocate(divided_tempX(2,2,num_divisions))
        
        divided_tempX(:,1,1) = tempX(:,1) 
        divided_tempX(:,2,1) = tempX(:,1)+ vector_dir_2D(:) 
        do i = 2,num_divisions-1
            divided_tempX(:,1,i) = divided_tempX(:,2,i-1)
            divided_tempX(:,2,i) = divided_tempX(:,2,i-1) + vector_dir_2D(:)
        enddo
        divided_tempX(:,1,num_divisions) = divided_tempX(:,2,num_divisions)
        divided_tempX(:,2,num_divisions) = tempX(:,2)       
        
        
        
    elseif(dim .eq. 3) then
        
    endif
    

endsubroutine
        
    
    