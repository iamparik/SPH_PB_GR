        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:38 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GRADIENTCSPM__genmod
          INTERFACE 
            SUBROUTINE GRADIENTCSPM(D1F,F0,XI_INV,SPH_DIM,NIAC,PAIR_I,  &
     &PAIR_J,DWDX,MASS,RHO,ITYPE,NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              INTEGER(KIND=4), INTENT(IN) :: NIAC
              INTEGER(KIND=4), INTENT(IN) :: SPH_DIM
              REAL(KIND=8), INTENT(INOUT) :: D1F(SPH_DIM,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: F0(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: XI_INV(SPH_DIM,SPH_DIM,NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_I(NIAC)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_J(NIAC)
              REAL(KIND=8), INTENT(IN) :: DWDX(SPH_DIM,NIAC)
              REAL(KIND=8), INTENT(IN) :: MASS(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: RHO(NTOTAL)
              INTEGER(KIND=2), INTENT(IN) :: ITYPE(NTOTAL)
            END SUBROUTINE GRADIENTCSPM
          END INTERFACE 
        END MODULE GRADIENTCSPM__genmod
