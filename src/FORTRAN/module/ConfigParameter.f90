    
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
        & time_step,print_step,save_step,backup_step, dx_r, ExtInputMeshType, packagingIterations,  &
        & hsml_const, rho_init, mu_const, hydroStaticHeight, g_const, ref_vel_max, c_sound, F_ext,  &
        & NumericalSimCase,  timeIntegrationScheme, SumDenstype, summationDensity, artViscType, MLS_density_bound,  &
        & MLS_step,BILtype, PrsrGradtype,PSTCoeff, ConDivtype, densDiffType, delta_SPH, prsrBdryType, HG_density_correction,PSTtype, &
        & nnps, nnes,WallBoundaryLayer, FScutoff
    
!SPH_dim : Dimension of the problem (1, 2 or 3)    
    integer(4) :: SPH_dim
    
!Smoothing kernel function 
!skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
!      2, Gauss kernel   (Gingold and Monaghan 1981) 
!      3, Quintic kernel (Morris 1997)
!	   4, yang
!      5. Wendland Quintic Kernel
    integer(4) :: skf
    
!Equation of State 
!eos = 1, Gamma law for ideal gas
!      2, Artificial EOS by Monaghan (1994) 
!      3, Artificial EOS by Morris (1997)
!      4, Mie-Gruneisen Equation of State 
!      5, Tillotson Equation of State
!      6, JWL & Mie-Gruneisen
!      7, Soil (From H.H.Bui')
    integer(4) :: eos
    
! Solver for linear system of equations (of pressure poisson)
! solver_type = 0 for Direct solution of Poisson's equation
!               1 for BICGSTAB2 method
    integer(4) :: solver_type
    
    real(8) :: time_step
    
    integer(4) :: print_step
    
    integer(4) :: save_step
    
    integer(4) :: backup_step
    
    real(8) :: dx_r
    
    integer(4) :: ExtInputMeshType
    
    logical :: packagingIterations
    
    real(8) :: hsml_const
    
    real(8) :: rho_init
    
    real(8) :: mu_const
    
    real(8) :: hydroStaticHeight
    
    real(8) :: g_const
    
    real(8) :: ref_vel_max
    
    real(8) :: c_sound
    
    real(8), dimension(2) :: F_ext    
    
    integer(4) :: NumericalSimCase
    
    integer(4) :: timeIntegrationScheme
    
    integer(4) :: BILtype
    
    integer(4) :: PrsrGradtype
    
    integer(4) :: ConDivtype
    
    integer(4) :: SumDenstype
    
    logical :: summationDensity
    
    integer(4) :: artViscType
    
    logical :: MLS_density_bound
    
    integer(4) :: MLS_step
    
    integer(4) :: densDiffType
    
    real(8) :: delta_SPH
    
    integer(4) :: prsrBdryType 
    
    logical :: HG_density_correction
    
    integer(4) :: PSTtype
    
    real(8) :: PSTCoeff
    
    real(8) :: delC_Cap
    
    logical :: WallBoundaryLayer
    
!Nearest neighbor particle searching (nnps) method
!nnps = 1 : Simplest and direct searching
!       2 : Sorting grid linked list
!       3 : Tree algorithm
    integer(4) :: nnps
    
    integer(4) :: nnes
    
    real(8) :: FScutoff
    
    integer(4) :: edge_to_dx_ratio
    
    integer(4) :: pack_step2a
    
    real(8) :: pack_step2b
    
    integer(4) :: pack_step2c
    
    integer(4) :: PST_step
    
    contains

    
    subroutine read_config_file
        character(100) :: filename
        integer :: io_status
        character(100) :: line
        character(50) :: param_name
        real(8) :: param_value

        filename = 'data_geo_config/config_data.txt'
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
        
        ! Now you can use the parameters in the config_parameters module
        write(*, *) 'SPH_dim = ', SPH_dim
        write(*, *) 'skf = ', skf
        write(*, *) 'eos = ', eos
        write(*, *) 'poisson eqn solver_type = ', solver_type
        write(*, *) 'time_step =  ', time_step
        write(*, *) 'dx_r = ', dx_r
        write(*,*)  'ExtInputMeshType = ', ExtInputMeshType
        write(*, *) 'packagingIterations = ', packagingIterations
        write(*, *) 'pack_step2a = ',  pack_step2a
        write(*, *) 'pack_step2b = ',  pack_step2b
        write(*, *) 'pack_step2c = ',  pack_step2c
        write(*, *) 'hsml_const = ', hsml_const
        write(*, *) 'rho_init = ', rho_init
        write(*, *) 'mu_const = ', mu_const
        write(*, *) 'hydroStaticHeight = ', hydroStaticHeight
        write(*, *) 'g_const = ', g_const
        write(*, *) 'ref_vel_max = ', ref_vel_max
        write(*, *) 'c_sound = ', c_sound
        write(*, *) 'F_ext = ', F_ext
        write(*, *) 'NumericalSimCase = ', NumericalSimCase
        write(*, *) 'timeIntegrationScheme = ', timeIntegrationScheme
        write(*, *) 'BILtype = ', BILtype
        write(*, *) 'PrsrGradtype = ', PrsrGradtype
        write(*, *) 'ConDivtype = ', ConDivtype
        write(*, *) 'SumDenstype = ', SumDenstype
        write(*, *) 'summationDensity = ', summationDensity
        write(*, *) 'artViscType = ', artViscType
        write(*, *) 'MLS_density_bound = ', MLS_density_bound
        write(*, *) 'MLS_step = ', MLS_step
        write(*, *) 'densDiffType = ', densDiffType
        write(*, *) 'delta_SPH = ', delta_SPH
        write(*, *) 'prsrBdryType = ', prsrBdryType
        write(*, *) 'HG_density_correction = ', HG_density_correction
        write(*, *) 'PSTtype = ',PSTtype
        write(*, *) 'PST_step = ',PST_step
        write(*, *) 'PSTCoeff = ', PSTCoeff
        write(*, *) 'WallBoundaryLayer = ', WallBoundaryLayer
        write(*, *) 'nnps = ', nnps
        write(*, *) 'nnes = ', nnes 
        write(*,*) 'FScutoff = ' , FScutoff
        write(*,*) 'edge_to_dx_ratio = ' , edge_to_dx_ratio
        
    end subroutine read_config_file
    
    
    
    subroutine set_parameter(param_name, param_value)
        character(len=*), intent(in) :: param_name
        real(8), intent(in) :: param_value

        ! Set the corresponding parameter based on its name
        select case (trim(param_name))
            case ('SPH_dim')
                SPH_dim = param_value
            case ('skf')
                skf = param_value
            case ('eos')
                eos = param_value
            case ('solver_type')
                solver_type = param_value
            case ('time_step')
                time_step = param_value
            case ('print_step')
                print_step= param_value
            case ('save_step')
                save_step= param_value
            case ('backup_step')
                backup_step= param_value
            case ('dx_r')
                dx_r = param_value
            case ('ExtInputMeshType')   
                ExtInputMeshType = param_value
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
            case ('hsml_const')
                hsml_const = param_value
            case ('rho_init')
                rho_init = param_value
            case ('mu_const')
                mu_const = param_value
            case ('hydroStaticHeight')
                hydroStaticHeight = param_value
            case ('g_const')
                g_const = param_value
            case ('ref_vel_max')
                ref_vel_max= param_value
            case ('c_sound')
                c_sound = param_value
            case ('F_ext_x')
                F_ext(1) = param_value
            case ('F_ext_y')
                F_ext(2) = param_value
            case ('NumericalSimCase')
                NumericalSimCase = param_value
            case ('timeIntegrationScheme')
                timeIntegrationScheme = param_value
            case ('BILtype')
                BILtype = param_value
            case ('PrsrGradtype')
                PrsrGradtype = param_value
            case ('ConDivtype')
                ConDivtype = param_value
            case ('SumDenstype')
                SumDenstype = param_value
            case ('summationDensity')
                if(param_value .eq. 0) then
                    summationDensity = .false.
                else
                    summationDensity = .true.    
                endif              
            case ('artViscType')
                artViscType = param_value
            case ('MLS_density_bound')
                if(param_value .eq. 0) then
                    MLS_density_bound = .false.
                else
                    MLS_density_bound = .true.    
                endif  
            case ('MLS_step')
                MLS_step = param_value
            case ('densDiffType')
                densDiffType = param_value
            case ('delta_SPH')
                delta_SPH = param_value
            case ('prsrBdryType')
                prsrBdryType = param_value
            case ('HG_density_correction')
                if(param_value .eq. 0) then
                    HG_density_correction = .false.
                else
                    HG_density_correction = .true.    
                endif  
            case ('PSTtype') 
                PSTtype = param_value
            case('PST_step')
                PST_step = param_value
            case('PSTCoeff')
                PSTCoeff = param_value
            case('delC_Cap')
                delC_Cap = param_value
                
            case ('WallBoundaryLayer')
                if(param_value .eq. 0) then
                    WallBoundaryLayer = .false.
                else
                    WallBoundaryLayer = .true.    
                endif   
            case ('nnps') 
                nnps = param_value
            case ('nnes') 
                nnes = param_value
            case ('FScutoff')
                FScutoff = param_value
            case ('edge_to_dx_ratio')
                edge_to_dx_ratio =param_value
            case default
                ! Handle unknown parameter name
                write(*, *) 'Unknown parameter: ', trim(param_name)
            end select           
           
    end subroutine set_parameter

end module config_parameter