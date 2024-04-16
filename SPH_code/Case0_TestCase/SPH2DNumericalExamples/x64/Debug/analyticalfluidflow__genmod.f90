        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:32 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ANALYTICALFLUIDFLOW__genmod
          INTERFACE 
            SUBROUTINE ANALYTICALFLUIDFLOW(ITIMESTEP,DT,VEL,BDRYVALS)
              USE PARTICLE_DATA, ONLY :                                 &
     &          NTOTAL,                                                 &
     &          ETOTAL,                                                 &
     &          ITYPE,                                                  &
     &          X,                                                      &
     &          RHO,                                                    &
     &          MU
              INTEGER(KIND=4), INTENT(IN) :: ITIMESTEP
              REAL(KIND=8), INTENT(IN) :: DT
              REAL(KIND=8) :: VEL(NTOTAL)
              LOGICAL(KIND=4), INTENT(IN) :: BDRYVALS
            END SUBROUTINE ANALYTICALFLUIDFLOW
          END INTERFACE 
        END MODULE ANALYTICALFLUIDFLOW__genmod
