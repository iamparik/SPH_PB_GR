        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr  1 17:04:34 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PERIODICPARAMETER__genmod
          INTERFACE 
            SUBROUTINE PERIODICPARAMETER(FX)
              USE PARTICLE_DATA, ONLY :                                 &
     &          NPERIODIC,                                              &
     &          PBC_DUPLICATE_PAIR,                                     &
     &          NTOTAL
              REAL(KIND=8) :: FX(NTOTAL)
            END SUBROUTINE PERIODICPARAMETER
          END INTERFACE 
        END MODULE PERIODICPARAMETER__genmod
