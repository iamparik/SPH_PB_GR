        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  4 12:45:21 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIVERGENCEKCBI_CSPM__genmod
          INTERFACE 
            SUBROUTINE DIVERGENCEKCBI_CSPM(D1F,F0,XI_INV,SPH_DIM,NIAC,  &
     &PAIR_I,PAIR_J,DWDX,ENIAC,EPAIR_A,DGMAS,MASS,RHO,ITYPE,NTOTAL,     &
     &NUM_EDGES)
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              INTEGER(KIND=4), INTENT(IN) :: ENIAC
              INTEGER(KIND=4), INTENT(IN) :: NIAC
              INTEGER(KIND=4), INTENT(IN) :: SPH_DIM
              REAL(KIND=8), INTENT(INOUT) :: D1F(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: F0(SPH_DIM,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: XI_INV(SPH_DIM,SPH_DIM,NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_I(NIAC)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_J(NIAC)
              REAL(KIND=8), INTENT(IN) :: DWDX(SPH_DIM,NIAC)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_A(ENIAC)
              REAL(KIND=8), INTENT(IN) :: DGMAS(SPH_DIM,ENIAC)
              REAL(KIND=8), INTENT(IN) :: MASS(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: RHO(NTOTAL)
              INTEGER(KIND=2), INTENT(IN) :: ITYPE(NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: NUM_EDGES
            END SUBROUTINE DIVERGENCEKCBI_CSPM
          END INTERFACE 
        END MODULE DIVERGENCEKCBI_CSPM__genmod
