!****************************************************************************
!
!  SUBROUTINE: nnes_algorithm
!
!  PURPOSE:  nnps algorithm finds the nearest particles to a given particle. We 
!            can use either direct find, link list or tree search
!
!   CREATED:        12/23/2020       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/18/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine nnes_algorithm
    use config_parameter, only: SPH_dim, nnes
    use particle_data ,   only: itimestep,epair_a, epair_s, max_e_interaction, &
            &   etotal, pBC_epair_a, pBC_epair_s
    
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
integer(4) max_pBC_e_interaction

max_pBC_e_interaction=int(40*max_e_interaction/etotal)
    
if (.NOT. Allocated(epair_a)) allocate(epair_a(max_e_interaction)) 
if (.NOT. Allocated(epair_s)) allocate(epair_s(max_e_interaction)) 

if (.NOT. Allocated(pBC_epair_a)) allocate(pbc_epair_a(max_pBC_e_interaction)) 
if (.NOT. Allocated(pBC_epair_s)) allocate(pbc_epair_s(max_pBC_e_interaction)) 

epair_a=0
epair_s=0

pBC_epair_a=0
pBC_epair_s=0

! Depending on value of nnps, different algorithms for 
! neighbor search can be used. 
! When one of the search algorithms is called the first time
! this code also displays the algorithm used
if (nnes.eq.1) then 
    call direct_edge_find
    if (itimestep.eq.1)   write(*,'(A)') 'direct_find for edge + particle pair has been called!'
endif

end