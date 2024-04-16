        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:31 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZEROTHORDERFUNCTION__genmod
          INTERFACE 
            SUBROUTINE ZEROTHORDERFUNCTION(F0,F,GMA,NIAC,PAIR_I,PAIR_J,W&
     &,W_AA,MASS,RHO,ITYPE,NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              INTEGER(KIND=4), INTENT(IN) :: NIAC
              REAL(KIND=8), INTENT(IN) :: F0(NTOTAL)
              REAL(KIND=8), INTENT(INOUT) :: F(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: GMA(NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_I(NIAC)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_J(NIAC)
              REAL(KIND=8), INTENT(IN) :: W(NIAC)
              REAL(KIND=8), INTENT(IN) :: W_AA(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: MASS(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: RHO(NTOTAL)
              INTEGER(KIND=2), INTENT(IN) :: ITYPE(NTOTAL)
            END SUBROUTINE ZEROTHORDERFUNCTION
          END INTERFACE 
        END MODULE ZEROTHORDERFUNCTION__genmod
