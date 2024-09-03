!****************************************************************************
!
!  SUBROUTINE: nnps_algorithm
!
!  PURPOSE:  nnps algorithm finds the nearest particles to a given particle. We 
!            can use either direct find, link list or tree search
!
!   CREATED:        12/23/2020       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/18/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine nnps_algorithm(adnl_srch_rds)
    use config_parameter, only: SPH_dim, nnps, print_step 
    use particle_data ,   only: itimestep, max_interaction, &
        & pair_i,pair_j,w,dwdx, niac
    implicit none
!-------------------------------------------------------------------
![output]           itimestep:The timestep of the SPH loop
![output]           ntotal: Total number of all particles in the 
!                           simulation
![Input]            hsml:   smoothing length of particle
![input/output]     x:      coordinate of particles
![output]           niac:   number of active particle interataction
!                           pairs
![output]           pair_i: ith particle of a given interaction pair
!                           for all pairs
![output]           pair_j: jth particle of a given interaction pair
!                           for all pairs
![output]           w:      smoothening kernel value for a given a
!                           pair
![output]           dwdx:   Smoothening kernel derivative for given
!                           pair
![output]           ns:     number of neighbouring particles for a 
!                           give particle
![input/output]     itype:  determines the type of particle like 
!                           water, soil, etc and also particle 
!                           type real/virtual/flow

!---------------------------------------------------------------------------------
    integer(4) i
    character (40) :: xname
    real(8), intent(in) :: adnl_srch_rds
    
if (.NOT. Allocated(pair_i)) allocate(pair_i(max_interaction)) 
if (.NOT. Allocated(pair_j)) allocate(pair_j(max_interaction)) 
if (.NOT. Allocated(w)) allocate(w(max_interaction)) 
if (.NOT. Allocated(dwdx)) allocate(dwdx(SPH_dim,max_interaction)) 

! We ensure the following variables are zero for every tiem the algorithm is called
pair_i=0
pair_j=0
w=0
dwdx=0


! Depending on value of nnps, different algorithms for 
! neighbor search can be used. 
! When one of the search algorithms is called the first time
! this code also displays the algorithm used
if (nnps.eq.1) then 
    call direct_find(adnl_srch_rds)
    if (mod(itimestep,print_step).eq.0)   then
        write(*,'(A)') 'direct_find for particle particle pair has been called!'
        write(*,*) ' Total number of particle particle itneractions =', niac
    endif
elseif (nnps.eq.2) then
    call linkedList(adnl_srch_rds)
    if (mod(itimestep,print_step).eq.0)   then
        write(*,'(A)') 'linkedList for particle particle pair has been called!'
        write(*,*) ' Total number of particle particle itneractions =', niac
    endif
endif


!write(*,*) ' Total number of particle particle itneractions =', niac
    !write(xname, '(A)') 'PPInteract_directfind.dat'
    !write(xname, '(A)') 'PPInteract_linklist.dat'
    !open (1, file = xname)
    !write(1,*) 'variables = nnps, a, b'
    !do i=1,niac
    !    write(1,1002) i, pair_i(i), pair_j(i)    
    !enddo
    !write(*,*) ' Total number of particle particle itneractions =', niac
    !pause
    !
    !1002  FORMAT(3(I10))
end