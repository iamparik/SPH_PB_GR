module particle_data
implicit none

!-------------------------------------------------------
!   Module:   particle_data 
!            
!   PURPOSE:  Creates and store global Variables that are required by
!               all kinds of simulations
!    
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/17/2022       by  PARIKSHIT BOREGOWDA 
!-------------------------------------------------------

![Output]           itimestep:      Tiemstep iterations for which the simulation has been run
public itimestep
integer(4):: itimestep
    
![Output]           maxn:      Maximum number of all particles possible
![Output]           max_interaction:      Maximum number of particle-particle interactions
![Output]           maxnv:      maximum number of boundaries, vertices and ghost points (all non real particles)
![Output]           max_e_interaction:      Maximum number of edge-particle interactions
![Output]           maxedge:      Maximum number of edge elements possible
    public maxn,max_interaction,maxnv,max_e_interaction, maxedge, simGridSize
    integer(4):: maxn
    integer(4):: max_interaction
    integer(4):: maxnv
    integer(4):: max_e_interaction
    integer(4):: maxedge
    real(8),dimension(:,:),ALLOCATABLE::simGridSize

    
![input/output]     x:      coordinate of particles  
![input/output]     vx:     velocity of particles
![input/output]     mass:   mass of particles
![input/output]     rho:    density of particles
![Output]           p:      Pressure value of particle
![Input]            hsml:   smoothing length of particle    
![Input]            mu:     viscosity of particle
![Input]            itype:     type of particle (real, edge, ghost etc)
![Input]            etype:     type of edge (real, periodic, etc)

    public x,vx,mass,rho,vol,p,hsml,mu,temp, itype,etype, avgVol, KE, PE, TE, max_vel, KE_prev,free_surf_particle, free_surf_val
    real(8),DIMENSION(:,:),ALLOCATABLE :: x
    real(8),DIMENSION(:,:),ALLOCATABLE :: x_prev
    real(8),DIMENSION(:,:),ALLOCATABLE :: vx
    real(8),DIMENSION(:),ALLOCATABLE:: mass
    real(8),DIMENSION(:),ALLOCATABLE:: rho
    real(8),DIMENSION(:),ALLOCATABLE:: vol
    real(8),DIMENSION(:),ALLOCATABLE:: p
    real(8),DIMENSION(:),ALLOCATABLE:: hsml
    real(8),DIMENSION(:),ALLOCATABLE:: mu
    real(8), dimension(:),allocatable:: temp 
    integer(2),DIMENSION(:),ALLOCATABLE::itype
    integer(4),DIMENSION(:),ALLOCATABLE::etype
    real(8) :: avgVol
    real(8),DIMENSION(:),ALLOCATABLE :: KE 
    real(8),DIMENSION(:),ALLOCATABLE :: PE
    real(8),DIMENSION(:),ALLOCATABLE :: TE
    real(8),DIMENSION(:),ALLOCATABLE :: max_vel
    real(8):: KE_prev
    integer(2),DIMENSION(:),ALLOCATABLE :: free_surf_particle 
    real(8), dimension(:), ALLOCATABLE :: free_surf_val
    

    
    
    
    
    public packed_x,packed_itype, packed_vol, delC,delCMax,delCAvg, &
        & delCL2,ga_Max,ga_Avg,ga_L2, xStart, g_a_min, g_a_max
    real(8),DIMENSION(:,:),ALLOCATABLE :: packed_x
    integer(2),DIMENSION(:),ALLOCATABLE::packed_itype
    real(8),DIMENSION(:),ALLOCATABLE:: packed_vol
    real(8),DIMENSION(:,:),ALLOCATABLE :: delC
    real(8),DIMENSION(:),ALLOCATABLE :: delCAvg
    real(8),DIMENSION(:),ALLOCATABLE :: delCMax
    real(8),DIMENSION(:),ALLOCATABLE :: delCL2
    real(8),DIMENSION(:),ALLOCATABLE :: ga_Avg
    real(8),DIMENSION(:),ALLOCATABLE :: ga_Max
    real(8),DIMENSION(:),ALLOCATABLE :: ga_L2
    real(8),DIMENSION(:),ALLOCATABLE :: g_a_min
    real(8),DIMENSION(:),ALLOCATABLE :: g_a_max
    real(8),DIMENSION(:,:),ALLOCATABLE :: xStart
    
![Input]            edge:               stores vertices of solid wall boudnary (edge in 2d and plane in 3d)    
![Output]           integ_pts_for_edge: first index returns the starting index of integration points for a given edge 
!                                       and second index returns the number of itnegration pts for the edge
![Output]           val_integ_pts:      first set of index determines the coordinates of the integration pt and 
!                                       the last index returns the weight of the numerical integration assosciated with the pt
![Output]           mid_pt_for_edge:    coordinate values of edge mid point.
![Output]           surf_norm:          surface normals for wall edges
![input/output]     x_ve:               coordinate of vertices  
![input/output]     vx_ve:              velocity veoctor of vertices 
    public edge,integ_pts_for_edge, val_integ_pts, surf_norm, x_ve
    integer(4), DIMENSION(:,:), ALLOCATABLE:: edge
    integer(4), DIMENSION(:,:), ALLOCATABLE:: integ_pts_for_edge
    real(8), DIMENSION(:,:), ALLOCATABLE:: val_integ_pts
    real(8), DIMENSION(:,:), ALLOCATABLE:: mid_pt_for_edge
    real(8),DIMENSION(:,:),ALLOCATABLE:: surf_norm  
    real(8),DIMENSION(:,:),ALLOCATABLE :: x_ve
    real(8),DIMENSION(:,:),ALLOCATABLE :: vx_ve


![input/output]     nreal:  total number of real particles
![input/output]     nflow:  total number of flow boundary particles
![input/output]     nedge:  total number of wall boundary particles
![input/output]     nghost: total number of reference points for edges/ ghost points
![input/output]     ntotal: total number of all particles (and referecne points to be printed) 
![input/output]     etotal: total number of all edges  
![input/output]     ereal: total number of real edges  
![input/output]     ve_total: total number of vertices
    public nreal,nflow,nedge,nghost, ntotal, etotal, ve_total
    integer(4):: nreal
    integer(4):: nflow
    integer(4):: nedge
    integer(4):: nghost
    integer(4):: ntotal
    integer(4):: etotal
    integer(4)::ve_total
    
![Output]           niac:   Total number of particle interactions
![Output]           pair_i: i'th particle of a given particle pair 
![Output]           pair_j: j'th particle of a given particle pair 
![Output]           w:      Value of the interpolation weight (value of smoothing kernel for a given particle distance), w_ab
![Output]           w_aa:     Value of the interpolation weight for r=0, for a give particle w_aa
![Output]           dwdx:   Derivative of smoothing Kernel for a given particle distance, with respect to one particle    
    public niac, pair_i, pair_j, w, dwdx
    integer(4):: niac
    integer(4), dimension(:),allocatable::  pair_i
    integer(4), dimension(:),allocatable::  pair_j 
    real(8), dimension(:),allocatable:: w
    real(8), dimension(:),allocatable:: w_aa
    real(8), dimension(:,:),allocatable::  dwdx 

    
![Output]           eniac:   Total number of particle edge interactions
![Output]           epair_a: ath particle of a given particle-edge pair 
![Output]           epair_s: sth edge of a given particle-edge pair 
    public epair_a, epair_s, eniac 
    integer(4), dimension(:),allocatable::  epair_a
    integer(4), dimension(:),allocatable::  epair_s
    integer(4):: eniac
    
    
!defining all normalization and correction factors
 ![Output]           gamma_discrt: Σ(W_ab*m_b/rho_b), shephards renormalization factor
   public gamma_discrt, gamma_cont, gamma_cont_prev, del_gamma_as, del_gamma_as_prev, &
       & del_gamma, xi_cont_mat, xi_cont_mat_inv, xi1_mat, xi1_mat_inv, beta_mat, gamma_mat, &
       & gamma_mat_inv, dynamicProblem, b_MLS, gamma_density_cont
    real(8), dimension(:), Allocatable :: gamma_discrt
    real(8), dimension(:), Allocatable :: gamma_cont    
    real(8), dimension(:), Allocatable :: gamma_cont_prev
    real(8), dimension(:,:), Allocatable :: del_gamma_as
    real(8), dimension(:), Allocatable :: del_gamma_as_prev
    real(8), dimension(:,:),allocatable:: del_gamma
    real(8), dimension(:,:,:),allocatable:: xi_cont_mat
    real(8), dimension(:,:,:),allocatable:: xi_cont_mat_inv
    real(8), dimension(:,:,:),allocatable:: xi1_mat
    real(8), dimension(:,:,:),allocatable:: xi1_mat_inv
    real(8), dimension(:,:,:),allocatable:: beta_mat
    real(8), dimension(:,:,:),allocatable:: gamma_mat
    real(8), dimension(:,:,:),allocatable:: gamma_mat_inv
    real(8),DIMENSION(:,:),ALLOCATABLE :: b_MLS
    logical :: dynamicProblem
    real(8), dimension(:), Allocatable :: gamma_density_cont
    
!the below parameters are used in the context of periodic BC
    public ntotal_prev,etotal_prev, ve_total_prev, pBC_edges, pBC_epair_a,pBC_epair_s, &
        & pBC_eniac,tangent_pBC, edge_pBC_pairs, pBC_duplicate_pair, nperiodic
    integer(4):: ntotal_prev
    integer(4):: etotal_prev
    integer(4):: ve_total_prev
    integer(4), dimension(:,:),allocatable:: pBC_edges
    integer(4), dimension(:),allocatable::  pBC_epair_a
    integer(4), dimension(:),allocatable::  pBC_epair_s
    integer(4):: pBC_eniac
    real(8), dimension(:,:,:), Allocatable ::tangent_pBC
    integer(4), dimension(:,:), allocatable:: edge_pBC_pairs
    integer(4), dimension(:,:), allocatable:: pBC_duplicate_pair
    integer(4):: nperiodic
    
! physics based paramters for flow problems
    public drho,rho_prev, dstress, dtemp, temp_prev, temp_analytical, dgrho_prev, rho_s,prsr_bdry_val
    real(8), dimension(:),allocatable:: drho 
    real(8), dimension(:),allocatable:: rho_prev
    real(8), dimension(:,:),allocatable:: dstress
    real(8), dimension(:),allocatable:: dtemp
    real(8), dimension(:),allocatable:: temp_prev 
    real(8), dimension(:),allocatable:: temp_analytical  
    real(8), dimension(:),allocatable:: dgrho_prev
    real(8), dimension(:),allocatable:: rho_s
    real(8), dimension(:),allocatable:: prsr_bdry_val
    
! paramter/s for boundary values
    public bdryVal_seg, p_counter
    real(8), dimension(:,:), allocatable:: bdryVal_seg
    real(8), dimension(:), allocatable:: p_counter
    
    
! paramter/s used to find error/relative error
    public ErrorL2norm, hsml_factor_name , hsml_factor_code, f0max
    real(8), dimension(:,:), allocatable :: ErrorL2norm
    character(9) :: hsml_factor_name = "_default_"
    integer(4) :: hsml_factor_code = 0
    real(8) :: f0max=1.D0
        
    
    end module
    