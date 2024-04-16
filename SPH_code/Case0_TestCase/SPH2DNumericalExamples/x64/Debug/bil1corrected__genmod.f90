        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  4 12:45:27 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BIL1CORRECTED__genmod
          INTERFACE 
            SUBROUTINE BIL1CORRECTED(LAPF,DF0,GMA_INV,SPH_DIM,ENIAC,    &
     &EPAIR_A,EPAIR_S,DGMAS,ETYPE,ETOTAL,NTOTAL,COEFF)
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              INTEGER(KIND=4), INTENT(IN) :: ETOTAL
              INTEGER(KIND=4), INTENT(IN) :: ENIAC
              INTEGER(KIND=4), INTENT(IN) :: SPH_DIM
              REAL(KIND=8), INTENT(INOUT) :: LAPF(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: DF0(SPH_DIM,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: GMA_INV(SPH_DIM,SPH_DIM,NTOTAL&
     &)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_A(ENIAC)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_S(ENIAC)
              REAL(KIND=8), INTENT(IN) :: DGMAS(SPH_DIM,ENIAC)
              INTEGER(KIND=4), INTENT(IN) :: ETYPE(ETOTAL)
              REAL(KIND=8), INTENT(IN) :: COEFF
            END SUBROUTINE BIL1CORRECTED
          END INTERFACE 
        END MODULE BIL1CORRECTED__genmod
