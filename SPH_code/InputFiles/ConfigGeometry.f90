    
!****************************************************************************
!
!  MODULE: config_geometry
!
!  PURPOSE:  Geemetric Configuration of input file is set here
!           
!
!   CREATED:        08/24/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  11/16/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
module config_geometry
    implicit none

    public  time_step,print_step,save_step,backup_step, dx_r, ExtInputMeshType, packagingIterations,  &
            & hsml_const, rho_init, mu_const, hydroStaticHeight, g_const, ref_vel_max, c_sound, F_ext,  &
            & NumericalSimCase,  timeIntegrationScheme, SumDenstype, summationDensity, artViscType, MLS_density_bound,  &
            & MLS_step,BILtype, PrsrGradtype,PSTCoeff, ConDivtype, densDiffType, delta_SPH, prsrBdryType, HG_density_correction,PSTtype, &
            & nnps, nnes,WallBoundaryLayer, FScutoff
    
    
    real(8) :: time_step
    integer(4) :: print_step
    integer(4) :: save_step
    integer(4) :: backup_step
    real(8) :: dx_r
    integer(4) :: ExtInputMeshType
    integer(4) :: packagingIterations
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
    integer(4) :: nnps
    integer(4) :: nnes
    real(8) :: FScutoff
    
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
        
        
          ! Now you can use the parameters in the config_geometry module
            write(*, *) 'time_step =  ', time_step
            write(*, *) 'dx_r = ', dx_r
            write(*,*) 'ExtInputMeshType = ', ExtInputMeshType
            write(*, *) 'packagingIterations = ', packagingIterations
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
            write(*, *) 'PSTCoeff = ', PSTCoeff
            write(*, *) 'WallBoundaryLayer = ', WallBoundaryLayer
            write(*, *) 'nnps = ', nnps
            write(*, *) 'nnes = ', nnes 
            write(*,*) 'FScutoff = ' , FScutoff
        
    end subroutine read_config_file
    
    
    
    subroutine set_parameter(param_name, param_value)
        character(len=*), intent(in) :: param_name
        real(8), intent(in) :: param_value

        ! Set the corresponding parameter based on its name
        select case (trim(param_name))
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
                packagingIterations = param_value
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
            case default
                ! Handle unknown parameter name
                write(*, *) 'Unknown parameter: ', trim(param_name)
            end select
            
            
            
           
    end subroutine set_parameter

end module config_geometry