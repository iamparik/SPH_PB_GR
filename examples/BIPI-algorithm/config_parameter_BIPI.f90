    
!****************************************************************************
!
!  MODULE: config_parameter
!
!  PURPOSE:  Set configuration of SPH simulation
!           
!
!   CREATED:        08/24/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/23/2024       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
module config_parameter
    implicit none 
        
!     save_step set to 1 for instant record of virtual particle changes      
    public pi
    real(8),parameter :: pi=3.14159265358979323846D0    

        
!   itype_virtual   = used in the context of modulus with 100, so 0,100,200,300... are virtual points,
!                    itype_virtual also indicates the edge point        
!   itype_real_max  = all real particles are currently defined to be between 1-99, so max is 99
!   itype_real_min  = all real particles are currently defined to be between 1-99, so min is 1        
    public itype_virtual, itype_real_max, itype_real_min, itype_periodic
    integer(2), parameter :: itype_virtual=100  
    integer(4), parameter :: itype_real_max=99
    integer(4), parameter :: itype_real_min=1
    integer(4), parameter :: itype_periodic=200
        
!   etype_virtual   = used in the context of modulus with 100, so 0,100,200,300... are virtual edges with different boundary types
!   etype_periodic  =  periodic boundary
!   etype_real_max  = all real particles are currently defined to be between 1-99, so max is 99
!   etype_real_min  = all real particles are currently defined to be between 1-99, so min is 1  
!   etype_SolidWall1 = Solid wall with velocity Dirichlet BC, and neumann Pressure and density
!   etype_SolidWall2 = Solid wall with Velocity Neumann BC, and Neumann Pressure and density
!   etype_FreeSurface1= Free surface
!   periodic_pairs   = number of periodic pairs - only 1 is tested
    public etype_virtual, etype_periodic, etype_real_max, etype_real_min, &
        & etype_thermal_dirichlet, etype_thermal_neumann, &
        & etype_SolidWall1,etype_SolidWall2,etype_FarWall, etype_FreeSurface1, &
        & periodic_pairs
        
    integer(4), parameter :: etype_virtual=100
    integer(4), parameter :: etype_real_max=99
    integer(4), parameter :: etype_real_min=1
    integer(4), parameter :: etype_SolidWall1=2
    integer(4), parameter :: etype_SolidWall2=4
    integer(4), parameter :: etype_FarWall=5
    integer(4), parameter :: etype_FreeSurface1=3
    integer(4), parameter :: etype_thermal_dirichlet=51
    integer(4), parameter :: etype_thermal_neumann=52
    integer(4), parameter :: etype_periodic=200
    integer(4), parameter :: periodic_pairs=1


!dataOutputPath: Directory to store the output file
!dataOutputPath3: Directory to store an alternative output file
    public DataConfigPath, dataOutputPath, dataOutputErrors ,dataPackingPath
    character(LEN=*) , parameter :: DataConfigPath='data_geo_config' 
    character(LEN=*) , parameter :: dataOutputPath='data_output' 
    character(LEN=*) , parameter :: dataOutputErrors='data_output_errors'        
    character(LEN=*) , parameter :: dataPackingPath='data_packing'
  
        
! Below are variables that can be defiend at run time
    public  SPH_dim, skf, eos, solver_type, &
        & time_step,print_step,save_step,save_packing, backup_step, dx_r, ExtInputMeshType, packagingIterations,  &
        & hsml_const, rho_init, mu_const, ref_vel_max, c_sound, F_ext,  &
        & NumericalSimCase,  timeIntegrationScheme, SumDenstype, summationDensity, artViscType, MLS_density_bound,  &
        & MLS_step,BILtype, PrsrGradtype,PSTCoeff, ConDivtype, densDiffType, delta_SPH, prsrBdryType, HG_density_correction,PSTtype, &
        & nnps, nnes,WallBoundaryLayer, FScutoff, edge_to_dx_ratio, pack_step2a, shorten_step2a, pack_step2b, pack_step2c, shorten_step2c, PST_step, &
        & time_ev_par_op, external_input_InitialCondition, P_EOS_extra
    
!SPH_dim : Dimension of the problem (1, 2 or 3)    
    integer(4), parameter :: SPH_dim = 2
    
!Smoothing kernel function 
!skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
!      2, Gauss kernel   (Gingold and Monaghan 1981) 
!      3, Quintic kernel (Morris 1997)
!	   4, yang
!      5. Wendland Quintic Kernel
    integer(4), parameter :: skf=5
    
!Equation of State 
!eos = 1, Gamma law for ideal gas
!      2, Artificial EOS by Monaghan (1994) 
!      3, Artificial EOS by Morris (1997)
!      4, Mie-Gruneisen Equation of State 
!      5, Tillotson Equation of State
!      6, JWL & Mie-Gruneisen
!      7, Soil (From H.H.Bui')
    integer(4), parameter :: eos =2
    
! Solver for linear system of equations (of pressure poisson)
! solver_type = 0 for Direct solution of Poisson's equation
!               1 for BICGSTAB2 method
    integer(4), parameter :: solver_type =1
    
    real(8) :: time_step = 2.5D-4
    
    integer(4) :: print_step
    
    integer(4) :: save_step
    
    integer(4) :: save_packing
    
    integer(4) :: backup_step
    
    real(8) :: dx_r
    
    integer(4), parameter :: ExtInputMeshType =2
    
    logical :: packagingIterations
    
    real(8) :: hsml_const
    
    real(8) , parameter :: rho_init = 1000.D0
    
    real(8), parameter  :: mu_const = 0.D0
    
    real(8), parameter :: ref_vel_max = 2.5D0
    
    real(8) :: c_sound = 20.D0
    
    real(8), dimension(2), parameter  :: F_ext = 0.D0 
    
    integer(4) , parameter :: NumericalSimCase = 4
    
    integer(4), parameter  :: timeIntegrationScheme =1
    
    integer(4) , parameter :: BILtype= 23
    
    integer(4), parameter  :: PrsrGradtype = 31
    
    integer(4), parameter  :: ConDivtype = 21
    
    integer(4) , parameter :: SumDenstype = 6
    
    logical, parameter  :: summationDensity = .false.
    
    integer(4), parameter  :: artViscType = 0
    
    logical, parameter  :: MLS_density_bound = .false.
    
    integer(4), parameter :: MLS_step = 1
    
    integer(4), parameter :: densDiffType = 0
    
    real(8) :: delta_SPH =0.1D0
    
    integer(4), parameter  :: prsrBdryType = 12
    
    logical, parameter  :: HG_density_correction = .false.
    
    logical, parameter :: FS_density_correction = .false.
    
    integer(4), parameter  :: PSTtype = 0
        
    integer(4)  :: PST_step = 1
    
    real(8) :: PSTCoeff
    
    real(8) :: delC_Cap
    
    logical, parameter :: WallBoundaryLayer = 0
    
!Nearest neighbor particle searching (nnps) method
!nnps = 1 : Simplest and direct searching
!       2 : Sorting grid linked list
!       3 : Tree algorithm
    integer(4) :: nnps
    
    integer(4) :: nnes
    
    real(8), parameter :: FScutoff = 1.5D0
    
    integer(4) :: edge_to_dx_ratio
    
    integer(4) :: pack_step2a
    
    logical :: shorten_step2a =.true.
    
    real(8) :: pack_step2b
    
    integer(4) :: pack_step2c
    
    logical :: shorten_step2c =.true.
   
    logical :: time_ev_par_op =.false.
    
    logical :: external_input_InitialCondition = .false.
    real(8) :: P_EOS_extra=0.D0
    
    !integer(4) :: calc_prsr_bdry_IDs(4) = (/6, 7, 8, 9/)
    
    contains

    
    subroutine read_config_file
        character(100) :: filename
        integer :: io_status
        character(100) :: line
        character(50) :: param_name
        real(8) :: param_value

        filename = 'data_geo_config/config_BIPI.txt'
        open(unit=10, file=filename, status='old', action='read', iostat=io_status)

        do
            read(10, '(A)', iostat=io_status) line
            if (io_status /= 0) exit

            ! Parse the line and extract parameter name and value
            read(line, *) param_name,param_value

            ! Set the corresponding parameter in the module
            call set_parameter(param_name, param_value)
        end do

        close(10)        
        
    end subroutine read_config_file
    
    
    
    subroutine set_parameter(param_name, param_value)
        character(len=*), intent(in) :: param_name
        real(8), intent(in) :: param_value

        ! Set the corresponding parameter based on its name
        select case (trim(param_name))
            case ('save_packing')
                save_packing = param_value
                save_step= param_value
                print_step= param_value
            case ('dx_r')
                dx_r = param_value
            case ('packagingIterations')
                if(param_value .eq. 0) then
                    packagingIterations = .false.
                else
                    packagingIterations = .true.    
                endif  
            case ('pack_step2a')
                pack_step2a =param_value
            case ('pack_step2b')
                pack_step2b =param_value
            case ('pack_step2c')
                pack_step2c =param_value
            case ('shorten_step2a')
                if(param_value .eq. 0) then
                    shorten_step2a = .false.
                else
                    shorten_step2a = .true.
                endif 
            case ('shorten_step2c')
                if(param_value .eq. 0) then
                    shorten_step2c = .false.
                else
                    shorten_step2c = .true.
                endif 
            case ('hsml_const')
                hsml_const = param_value
            case ('delta_SPH')
                delta_SPH = param_value
            case('PST_step')
                PST_step = param_value
            case('PSTCoeff')
                PSTCoeff = param_value
            case('delC_Cap')
                delC_Cap = param_value
            case ('nnps') 
                nnps = param_value
            case ('nnes') 
                nnes = param_value
            case ('edge_to_dx_ratio')
                edge_to_dx_ratio =param_value
            case('time_ev_par_op')            
                if(param_value .eq. 1) then
                    time_ev_par_op = .true.
                else
                    time_ev_par_op = .false.    
                endif
            case default
                ! Handle unknown parameter name
                write(*, *) 'Unknown parameter: ', trim(param_name)
            end select   
            
            
           
    end subroutine set_parameter

end module config_parameter