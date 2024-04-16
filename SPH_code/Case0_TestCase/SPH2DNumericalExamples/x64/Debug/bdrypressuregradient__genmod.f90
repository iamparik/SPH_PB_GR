        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:28 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BDRYPRESSUREGRADIENT__genmod
          INTERFACE 
            SUBROUTINE BDRYPRESSUREGRADIENT(D1F,F0,GMA)
              USE PARTICLE_DATA, ONLY :                                 &
     &          EPAIR_A,                                                &
     &          EPAIR_S,                                                &
     &          ENIAC,                                                  &
     &          ETOTAL,                                                 &
     &          MASS,                                                   &
     &          RHO,                                                    &
     &          ETYPE,                                                  &
     &          DEL_GAMMA_AS,                                           &
     &          NEDGE_REL_EDGE,                                         &
     &          NTOTAL,                                                 &
     &          X,                                                      &
     &          VX,                                                     &
     &          SURF_NORM,                                              &
     &          GAMMA_DENSITY_CONT
              REAL(KIND=8), INTENT(INOUT) :: D1F(2,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: F0(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: GMA(NTOTAL)
            END SUBROUTINE BDRYPRESSUREGRADIENT
          END INTERFACE 
        END MODULE BDRYPRESSUREGRADIENT__genmod
