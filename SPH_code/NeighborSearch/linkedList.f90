!****************************************************************************
!
!  SUBROUTINE: link_list 
!
!  PURPOSE: Linked-list algorithm is a temporary mesh of cells
!           overlaid on the problem domain. Particles are assigned
!           to cells and identified through grid linked lists
!
!   CREATED:        4/1/2020       by  PARIKSHIT BOREGOWDA
!   Last Modified:  11/27/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************  
subroutine link_list
    use config_parameter, only:SPH_dim, itype_virtual
    use config_geometry, only: hsml_const
    use particle_data ,   only: niac,ntotal,itype,x,hsml,max_interaction, &     
        & simGridSize,pair_i,pair_j,w,dwdx       
    implicit none
    
    integer(4) gridMax(3), adjCells, d, dd, numCells , nc(3), cell, &
        & a,b,c, i,j,k, ii,jj,kk, cell2
    real(8) cellwidth, scale_k, mhsml
    integer(4),DIMENSION(:),ALLOCATABLE:: firstPinCell, lastPinCell, linkedPlist
    integer :: directions(3, 13) = reshape([ &
        & 0, 1, 0, &
        & 0, 0, 1, &
        & 0, 1, 1, &
        & 1, 0, 0, &
        & 1, 1, 0, &
        & 1, 0, 1, &
        & 1, 1, 1, &
        & -1, 1, 0, &
        & -1, 0, 1, &
        & -1, 1, 1, &
        & 0, -1, 1, &
        & 1, -1, 1, &
        & 1, 1, -1 ], [3, 13])
    
    ! we need to find k of kh
    call sml_mult_factor(scale_k)
    
    !do d=1,13
    !    write(*,*) directions(:,d)
    if(SPH_dim .eq. 3) write(*,*) "Note 3D linkedlist not eyt tested!!"
    !enddo
    
    !Here we define the width of each background cell. This can be further
    !multiplied by a factor to decrease or increase the width
    cellwidth = max(1.D0*scale_k*hsml_const,scale_k*hsml_const)
    
    adjCells=1   !max(1, ceiling((scale_k*hsml_const)/cellwidth)) !-1
    
    

    !initialize number of cells in the domain as 0
    numCells=0
     
    ! In a 3D space, there can be maximum one cell
    gridMax(:)=1
    
    !First we find number of cells in each direction, using simulation bounds of grid size
    do d=1,SPH_dim
        gridMax(d) = ceiling((simGridSize(d,2) - simGridSize(d,1))/cellwidth)
        if( gridMax(d) .eq. 0) gridMax(d) = 1 ! this to ensure if a dimension input with no thickness is accoutned accurately
    enddo
    
    !total cells, is simple i*j*k
    numCells=product(gridMax)
    
    
    allocate(firstPinCell(numCells),lastPinCell(numCells), linkedPlist(ntotal))
    firstPinCell=0    
    lastPinCell=0
    linkedPlist=0
    
    ! Create a linked list of particles assosicated to a given cell. If a cell is empty then its firstPinCell will remain 0.    
    do a=1,ntotal
        nc(:)=1
        if(mod(itype(a),itype_virtual) .ne. 0) then
            do d=1,SPH_dim
                !Find the i,j,k position of cell the particle belongs to            
                nc(d) = ceiling((x(d,a) - simGridSize(d,1))/cellwidth)
                if (nc(d) .eq. 0) nc(d) =1
                if (nc(d) .lt. 0) then
                    write(*,*) "value of cell number negative in linkedlist NNPS algorithm"    
                    pause
                    nc(d)=1
                endif            
            enddo   
        
            cell= (nc(3)-1)*gridMax(2)*gridMax(1)+(nc(2)-1)*gridMax(1)+ (nc(1)-1)+1
        
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
    
    ! We initailize the variable storing the number of interactions 
    niac=0
    
    ! We go through all the particles in a given cell
    do c=1, numcells
        
        if(lastPinCell(c) .gt. 0) then
            a=0 ! Initialize particle number of cell to zero
            b=0
            do while (a .ne. lastPinCell(c))
                if (a .eq. 0) then
                    a=firstPinCell(c) ! start with first particle in cell
                else
                    a=linkedPlist(a) ! after first particle use linkedlist to update
                endif
                
                b=linkedPlist(a)! This ensures we start from the particle next to a in linked list
                
                ! Now we need to go over all particles after a, not including a. Also, we need to ensure we dont
                ! continue if there is only 1 particle in cell or on the last particle in cell already at a,
                ! in which case b=0
                do while (b .ne. 0) 
                    
                    !The average smoothinglength is calculated, but we can simply use hsml_const
                    mhsml = hsml_const
                    call particlePair(a,b,scale_k, mhsml) ! now check if the two particles form a particle pair. If yes then assign the particle pair
                    
                    b=linkedPlist(b) ! find next particle in linked list after a, and repeat till last in list
                enddo
            enddo
        endif
        
    enddo

    
    ! We now go through all particles of a cell with its adjacent cell.
   ! start from origin cell, and then update aong x direction, then y and then finally along z
    do c=1,numcells
        
        i=mod(c-1,gridMax(1))+1
        j=mod((c-1 -(i-1))/gridMax(1),gridMax(2))+1
        k= ((c-1)-(i-1) - (j-1)*gridMax(1))/(gridMax(2)*gridMax(1)) +1
        
        if(lastPinCell(c) .gt. 0) then ! this ensure we consider a cell that has atleast one particle in it
            
            do d=1,13
                
                ii=directions(1,d)
                jj=directions(2,d)
                kk=directions(3,d)
                
                if((k+kk .le. gridMax(3)) .and.(j+jj .le. gridMax(2)) .and. (i+ii .le. gridMax(1)) .and. (k+kk .ge. 1) .and.(j+jj .ge. 1) .and. (i+ii .ge. 1)) then
                    ! Determine the cell number of adjacent cell
                    cell2=(k+kk-1)*gridMax(2)*gridMax(1)+(j+jj-1)*gridMax(1)+(i+ii-1)+1
                    if(lastPinCell(cell2) .gt. 0) then ! Ensures neighboring cell has atleast one particle in it
                        a=0 ! Initialize particle number of cell to zero
                        do while (a .ne. lastPinCell(c)) !go over all particles in linkedlist till the last particle in cell
                                    
                            if (a .eq. 0) then
                                a=firstPinCell(c) ! start with first particle in cell
                            else
                                a=linkedPlist(a) ! after first particle use linkedlist to update
                            endif
                                    
                            b=0 ! Initialize particle number of neibhoring cell to zero
                                    
                            do while (b .ne. lastPinCell(cell2))
                                    
                                if (b .eq. 0) then
                                    b=firstPinCell(cell2) ! start with first particle in cell
                                else
                                    b=linkedPlist(b) ! after first particle use linkedlist to update
                                endif
                                        
                                !The average smoothinglength is calculated, but we can simply use hsml_const
                                mhsml = hsml_const
                                call particlePair(a,b,scale_k, mhsml) ! now check if the two particles form a particle pair. If yes then assign the particle pair
                                    
                            enddo
                        enddo                            
                                
                    endif
                endif
            enddo    
             
             
             
        endif
        
    enddo
end subroutine
    
  