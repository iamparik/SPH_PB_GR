        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:31 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ANALYTICALTHERMALDIFFUSION__genmod
          INTERFACE 
            SUBROUTINE ANALYTICALTHERMALDIFFUSION(ITIMESTEP,DT,TEMP,    &
     &BDRYVALS)
              USE PARTICLE_DATA, ONLY :                                 &
     &          NTOTAL,                                                 &
     &          ETOTAL,                                                 &
     &          ITYPE,                                                  &
     &          X
              INTEGER(KIND=4), INTENT(IN) :: ITIMESTEP
              REAL(KIND=8), INTENT(IN) :: DT
              REAL(KIND=8) :: TEMP(NTOTAL)
              LOGICAL(KIND=4), INTENT(IN) :: BDRYVALS
            END SUBROUTINE ANALYTICALTHERMALDIFFUSION
          END INTERFACE 
        END MODULE ANALYTICALTHERMALDIFFUSION__genmod
