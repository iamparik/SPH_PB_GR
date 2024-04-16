        !COMPILER-GENERATED INTERFACE MODULE: Wed May  3 13:46:59 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GRADIENTFERRANDPCBI__genmod
          INTERFACE 
            SUBROUTINE GRADIENTFERRANDPCBI(D1F,F0,FS0,GMA,GMA_INV,      &
     &SPH_DIM,NIAC,PAIR_I,PAIR_J,DWDX,ENIAC,EPAIR_A,EPAIR_S,DGMAS,MASS, &
     &RHO,RHOS,ITYPE,NTOTAL,NUM_EDGES)
              INTEGER(KIND=4), INTENT(IN) :: NUM_EDGES
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              INTEGER(KIND=4), INTENT(IN) :: ENIAC
              INTEGER(KIND=4), INTENT(IN) :: NIAC
              INTEGER(KIND=4), INTENT(IN) :: SPH_DIM
              REAL(KIND=8), INTENT(INOUT) :: D1F(SPH_DIM,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: F0(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: FS0(NUM_EDGES)
              REAL(KIND=8), INTENT(IN) :: GMA(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: GMA_INV(SPH_DIM,SPH_DIM,NTOTAL&
     &)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_I(NIAC)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_J(NIAC)
              REAL(KIND=8), INTENT(IN) :: DWDX(SPH_DIM,NIAC)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_A(ENIAC)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_S(ENIAC)
              REAL(KIND=8), INTENT(IN) :: DGMAS(SPH_DIM,ENIAC)
              REAL(KIND=8), INTENT(IN) :: MASS(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: RHO(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: RHOS(NUM_EDGES)
              INTEGER(KIND=2), INTENT(IN) :: ITYPE(NTOTAL)
            END SUBROUTINE GRADIENTFERRANDPCBI
          END INTERFACE 
        END MODULE GRADIENTFERRANDPCBI__genmod
