        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:35 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BDRYLAYERVEL__genmod
          INTERFACE 
            SUBROUTINE BDRYLAYERVEL(SPH_DIM,BDRYVALWALLVEL)
              USE PARTICLE_DATA, ONLY :                                 &
     &          EPAIR_A,                                                &
     &          EPAIR_S,                                                &
     &          ENIAC,                                                  &
     &          SURF_NORM,                                              &
     &          NEDGE_REL_EDGE,                                         &
     &          MU,                                                     &
     &          X,                                                      &
     &          VX
              INTEGER(KIND=4), INTENT(IN) :: SPH_DIM
              REAL(KIND=8) :: BDRYVALWALLVEL(SPH_DIM,ENIAC)
            END SUBROUTINE BDRYLAYERVEL
          END INTERFACE 
        END MODULE BDRYLAYERVEL__genmod
