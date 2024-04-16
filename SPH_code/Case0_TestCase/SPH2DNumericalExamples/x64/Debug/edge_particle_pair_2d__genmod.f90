        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  4 12:45:26 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EDGE_PARTICLE_PAIR_2D__genmod
          INTERFACE 
            SUBROUTINE EDGE_PARTICLE_PAIR_2D(KHSML,XI,XS1,XS2,NS,       &
     &INTERACT_EDGE_PARTICLE,ERROR_TOL)
              REAL(KIND=8) :: KHSML
              REAL(KIND=8) :: XI(2)
              REAL(KIND=8) :: XS1(2)
              REAL(KIND=8) :: XS2(2)
              REAL(KIND=8) :: NS(2)
              LOGICAL(KIND=4) :: INTERACT_EDGE_PARTICLE
              REAL(KIND=8) :: ERROR_TOL
            END SUBROUTINE EDGE_PARTICLE_PAIR_2D
          END INTERFACE 
        END MODULE EDGE_PARTICLE_PAIR_2D__genmod
