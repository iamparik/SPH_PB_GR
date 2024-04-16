!****************************************************************************
!
!  MODULE: config_parameter
!
!  PURPOSE:  General Configuration of file is set here
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/07/2024       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
    module config_parameter
        implicit none

        !dataOutputPath: Directory to store the output file
        !dataOutputPath3: Directory to store an alternative output file
        public DataConfigPath, dataOutputPath, dataOutputErrors ,dataPackingPath
        character(LEN=*) , parameter :: DataConfigPath='data_geo_config' 
        character(LEN=*) , parameter :: dataOutputPath='data_output' 
        character(LEN=*) , parameter :: dataOutputErrors='data_output_errors'        
        character(LEN=*) , parameter :: dataPackingPath='data_packing'

        !     save_step set to 1 for instant record of virtual particle changes      
        public pi
        real(8),parameter :: pi=3.14159265358979323846D0    
        
        
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
        !   etype_SolidWall1 = Solid wall with velocity Dirichlet BC, and neumann Pressure and density
        !   etype_SolidWall2 = Solid wall with Velocity Neumann BC, and Neumann Pressure and density
        !   etype_FreeSurface1= Free surface
        public etype_virtual, etype_periodic, etype_real_max, etype_real_min, &
            & etype_thermal_dirichlet, etype_thermal_neumann, &
            & etype_SolidWall1,etype_SolidWall2,etype_FarWall, etype_FreeSurface1
        
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
        
        
    end module
