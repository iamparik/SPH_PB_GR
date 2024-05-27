!****************************************************************************
!
!  MODULE: config_parameter
!
!  PURPOSE:  General Configuration of file is set here
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/21/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
    module config_parameter
        implicit none
        !SPH_dim : Dimension of the problem (1, 2 or 3)
        public SPH_dim
        integer(4),parameter :: SPH_dim=2

        !dataOutputPath: Directory to store the output file
        !dataOutputPath3: Directory to store an alternative output file
        public DataConfigOutputPath,dataOutputPath, DataFactorsOutputPath, dataOutputErrors
        character(LEN=*) , parameter :: DataConfigOutputPath='data_geo_config' 
        character(LEN=*) , parameter :: dataOutputPath='data_output' 
        character(LEN=*) , parameter :: DataFactorsOutputPath='data_output_factors'
        character(LEN=*) , parameter :: dataOutputErrors='data_output_errors'


        !Nearest neighbor particle searching (nnps) method
        !nnps = 1 : Simplest and direct searching
        !       2 : Sorting grid linked list
        !       3 : Tree algorithm
        public nnps
        integer(4),parameter :: nnps=1

        
        !Smoothing kernel function 
        !skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
        !      2, Gauss kernel   (Gingold and Monaghan 1981) 
        !      3, Quintic kernel (Morris 1997)
        !	   4, yang
        !      5. Wendland Quintic Kernel
        public skf
        integer(4),parameter :: skf=5   

        !Equation of State 
        !eos = 1, Gamma law for ideal gas
        !      2, Artificial EOS by Monaghan (1994) 
        !      3, Artificial EOS by Morris (1997)
        !      4, Mie-Gruneisen Equation of State 
        !      5, Tillotson Equation of State
        !      6, JWL & Mie-Gruneisen
        !      7, Soil (From H.H.Bui')
        public eos
        integer(4),parameter :: eos=2
        


        !     Control parameters for output 
        !     int_stat = .true. : Print statistics about SPH particle interactions.
        !                         including virtual particle information.
        !     print_step: Print Timestep (On Screen)
        !     save_step : Save Timestep    (To Disk File)
        !	  backup_step: used for restart or continue of analysis or as backup data
        
        public int_stat,print_step,save_step,backup_step

        logical,parameter :: int_stat=.true.
        integer(4),parameter :: print_step = 100
        integer(4),parameter :: save_step = 100 
        integer(4),parameter :: backup_step = 10000000



        !     save_step set to 1 for instant record of virtual particle changes      
        public pi
        real(8),parameter :: pi=3.14159265358979323846D0    


        public solver_type
        ! solver_type = 0 for Direct solution of Poisson's equation
        !               1 for BICGSTAB2 method
        integer(4),parameter :: solver_type=1
        
        
        !   itype_virtual   = used in the context of modulus with 100, so 0,100,200,300... are virtual points,
        !                    itype_virtual also indicates the edge point        
        !   itype_real_max  = all real particles are currently defined to be between 1-99, so max is 99
        !   itype_real_min  = all real particles are currently defined to be between 1-99, so min is 1        
        public itype_virtual, itype_real_max, itype_real_min, itype_periodic
        integer(4), parameter :: itype_virtual=100  
        integer(4), parameter :: itype_real_max=99
        integer(4), parameter :: itype_real_min=1
        integer(4), parameter :: itype_periodic=200
        
        !   etype_virtual   = used in the context of modulus with 100, so 0,100,200,300... are virtual edges with different boundary types
        !   etype_periodic  =  periodic boundary
        !   etype_real_max  = all real particles are currently defined to be between 1-99, so max is 99
        !   etype_real_min  = all real particles are currently defined to be between 1-99, so min is 1  
        !   etype_SolidWall1 = Solid wall
        !   etype_FreeSurface1= Free surface
        public etype_virtual, etype_periodic, etype_real_max, etype_real_min, &
            & etype_thermal_dirichlet, etype_thermal_neumann, &
            & etype_SolidWall1,etype_SolidWall2, etype_FreeSurface1
        
        integer(4), parameter :: etype_virtual=100
        integer(4), parameter :: etype_real_max=99
        integer(4), parameter :: etype_real_min=1
        integer(4), parameter :: etype_SolidWall1=2
        integer(4), parameter :: etype_SolidWall2=4
        integer(4), parameter :: etype_FreeSurface1=3
        integer(4), parameter :: etype_thermal_dirichlet=51
        integer(4), parameter :: etype_thermal_neumann=52
        integer(4), parameter :: etype_periodic=200
        

    end module

    
    
!****************************************************************************
!
!  MODULE: config_geometry
!
!  PURPOSE:  Geemetric Configuration of input file is set here
!           
!
!   CREATED:        08/24/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/21/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
    module config_geometry
        implicit none

        public time_step
        real(8),parameter ::  time_step =1.D-4 !2.5D-4

        
        public npw,npl,width,length,dx_r, bt_dx_r
        real(8), parameter :: width = 1.D0
        real(8), parameter :: length= 0.5D0
        integer(4), parameter ::  npw= 50!20
        real(8), parameter :: bt_dx_r=0.D0 ! ratio:  (boundary thickness/dx_r)
        !dx_r= (width - b_thickness*2)/npw !we can now define thickness in terms of dx_r
        real(8), parameter ::  dx_r= width/(dble(npw) + bt_dx_r*2.D0)
        integer(4), parameter ::  npl=nint(length/dx_r)        


        !create constant hsml and smoothening length multiplier
        public hsml_const,rho_init, mu_const, c_sound, g_const, F_ext, C_shift, &
            & NumericalSimCase, BILtype, PrsrGradtype, ConDivtype
        real(8), parameter ::   hsml_const=dx_r*2.D0
        real(8), parameter ::   rho_init=1000.D0
        real(8), parameter ::   mu_const = 1.D-3!1.D-3   
        real(8), parameter ::   c_sound = 1.D-3 !20.D0
        real(8), parameter ::   g_const = 0.D0!9.81D0 
        real(8), dimension(:), allocatable ::   F_ext !5.D-5   
        real(8), parameter ::   C_shift = 0.7D0      
        integer(4), parameter ::   NumericalSimCase=1 !if unknown or relError Calcualtion is to be avoided type 0
        integer(4) ::    BILtype=0
        integer(4), parameter ::    PrsrGradtype=0
        integer(4), parameter ::    ConDivtype=0
      

    end module