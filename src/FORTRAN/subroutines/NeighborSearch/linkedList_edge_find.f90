!****************************************************************************
!
!  SUBROUTINE: link_list 
!
!  PURPOSE: Linked-list algorithm is a temporary mesh of cells
!           overlaid on the problem domain. Particles are assigned
!           to cells and identified through grid linked lists
!
!   CREATED:        4/1/2020         by  PARIKSHIT BOREGOWDA
!   Last Modified:  11/27/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************  
subroutine linkedList_edge_find
    use config_parameter, only:SPH_dim, itype_virtual, hsml_const
    use particle_data ,   only: niac,ntotal,itype,x,hsml,max_interaction, &     
        & simGridSize,pair_i,pair_j,w,dwdx, x_ve ,etotal, eniac, pbc_eniac, mid_pt_for_edge, edge
    implicit none
    
    integer(4) gridMax(3), adjCells, d, dd, numCells , nc(3), cell, &
        & a,b,c, i,j,k, ii,jj,kk, cell2, cl, p, q, r, coord_min(3),coord_max(3)
    real(8) cellwidth, scale_k, khsml, x_mid
    integer(4),DIMENSION(:),ALLOCATABLE:: firstPinCell, lastPinCell, linkedPlist, startCellofEdge, lastCellofEdge
    
    integer(4) :: s
    
    ! we need to find k of kh
    call sml_mult_factor(scale_k)
    
    !do d=1,13
    !    write(*,*) directions(:,d)
    if(SPH_dim .eq. 3) write(*,*) "Note 3D linkedlist not yet tested!!"
    !enddo
    
    !Here we define the width of each background cell. This can be further
    !multiplied by a factor to decrease or increase the width
    cellwidth = scale_k*hsml_const !max(2.D0*scale_k*hsml_const,scale_k*hsml_const)  
    
    !initialize number of cells in the domain as 0
    numCells=0
     
    ! In a 3D space, there can be maximum one cell
    gridMax(:)=1
    
    !First we find number of cells in each direction, using simulation bounds of grid size
    do d=1,SPH_dim
        gridMax(d) = ceiling((simGridSize(d,2) - simGridSize(d,1))/cellwidth)
        if( gridMax(d) .eq. 0) gridMax(d) = 1 ! this to ensure if a dimension input with no thickness is accounted accurately
    enddo
    
    !total cells, is simple i*j*k
    numCells=product(gridMax)
    
    
    allocate(firstPinCell(numCells),lastPinCell(numCells), linkedPlist(ntotal), startCellofEdge(etotal), lastCellofEdge(etotal))
    firstPinCell=0    
    lastPinCell=0
    linkedPlist=0
    startCellofEdge=1
    lastCellofEdge=1
    
    ! Create a linked list of particles assosicated to a given cell. If a cell is empty then its firstPinCell will remain 0.    
    do a=1,ntotal
        nc(:)=1
        if(mod(itype(a),itype_virtual) .ne. 0) then
            do d=1,SPH_dim
                !Find the i,j,k position of cell the particle belongs to            
                nc(d) = ceiling((x(d,a) - simGridSize(d,1))/cellwidth)
                if (nc(d) .eq. 0) nc(d) =1
                if (nc(d) .lt. 0) then
                    write(*,*) "value of cell number negative in linkedlist NNPS algorithm for x at ", x(:,a), "for d =", d   
                    pause
                    nc(d)=1
                endif            
            enddo   
        
            cell= (nc(3)-1)*gridMax(2)*gridMax(1)+(nc(2)-1)*gridMax(1)+ nc(1)
        
            if(cell .gt. numcells) then
                write(*,*) "A particle has left the boundary"    
                pause
            endif
            
            if (lastPinCell(cell) .eq. 0) then
                firstPinCell(cell) = a
                lastPinCell(cell)=a
            else
                linkedPlist(lastPinCell(cell))= a
                lastPinCell(cell)=a
            endif
        endif
            
    enddo
    
    
    
    do s = 1, etotal
        nc(:)=1
        coord_min(:)=gridMax(:)
        coord_max(:) =1
        do dd=1,SPH_dim
            do d=1,SPH_dim
                !Find the i,j,k position of cell the particle belongs to            
                nc(d) = ceiling((x_ve(d, edge(dd,s)) - simGridSize(d,1))/cellwidth)
                if (nc(d) .eq. 0) nc(d) =1
                if (nc(d) .lt. 0) then
                    write(*,*) "value of cell number negative in linkedlist NNPS algorithm for x at ", x_ve(:, edge(dd,s)), "for d =", d   
                    pause
                    nc(d)=1
                endif   
                coord_min(d)= max(min(nc(d),coord_min(d))-1,1)
                coord_max(d)= min(max(nc(d),coord_max(d))+1,gridMax(d))
            enddo                     
        enddo
        cell= (coord_min(3)-1)*gridMax(2)*gridMax(1)+(coord_min(2)-1)*gridMax(1) + coord_min(1)
        startCellofEdge(s) = cell
        cell= (coord_max(3)-1)*gridMax(2)*gridMax(1)+(coord_max(2)-1)*gridMax(1) + coord_max(1)
        lastCellofEdge(s) = cell
    enddo
    
    ! We initailize the variable storing the number of interactions 
    eniac=0
    
    pbc_eniac=0
    
    
    ! We now go through all particles of a cell with its adjacent cell.
    ! start from origin cell, and then update aong x direction, then y and then finally along z
    do s=1,etotal
        c=startCellofEdge(s)
        i=mod(c-1,gridMax(1))+1
        j=mod((c-1 -(i-1))/gridMax(1),gridMax(2))+1
        k= ((c-1)-(i-1) - (j-1)*gridMax(1))/(gridMax(2)*gridMax(1)) +1  
        
        cl=lastCellofEdge(s)
        p=mod(cl-1,gridMax(1))+1
        q=mod((cl-1 -(i-1))/gridMax(1),gridMax(2))+1
        r= ((cl-1)-(i-1) - (j-1)*gridMax(1))/(gridMax(2)*gridMax(1)) +1  
 
        do ii=i-1,p+1
            do jj=j-1,q+1
                do kk=k-1,r+1                       
                    if((kk .le. gridMax(3)) .and.(jj .le. gridMax(2)) .and. (ii .le. gridMax(1)) .and. (kk .ge. 1) .and.(jj .ge. 1) .and. (ii .ge. 1)) then
                    ! Determine the cell number of adjacent cell
                        cell2=(kk-1)*gridMax(2)*gridMax(1)+(jj-1)*gridMax(1)+(ii-1)+1
                        if(lastPinCell(cell2) .gt. 0) then ! Ensures neighboring cell has atleast one particle in it 
                            a=0 ! Initialize particle number of neibhoring cell to zero                                    
                            do while (a .ne. lastPinCell(cell2))
                                if (a .eq. 0) then
                                    a=firstPinCell(cell2) ! start with first particle in cell
                                else
                                    a=linkedPlist(a) ! after first particle use linkedlist to update
                                endif
                                
                                khsml = scale_k*hsml(a)
            
                                call particleEdgePair(a,s,scale_k,khsml)
                            enddo                         
                        endif
                    endif                      
                enddo                 
            enddo
        enddo     
    enddo
    
    deallocate(firstPinCell,lastPinCell, linkedPlist,startCellofEdge, lastCellofEdge)
    
end subroutine
