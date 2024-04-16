        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:33 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MLSFUNCTIONAPPROX__genmod
          INTERFACE 
            SUBROUTINE MLSFUNCTIONAPPROX(F,F0,B_MLS,X,NIAC,PAIR_I,PAIR_J&
     &,W,W_AA,MASS,RHO,ITYPE,GAMMA_CONT,NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              INTEGER(KIND=4), INTENT(IN) :: NIAC
              REAL(KIND=8), INTENT(INOUT) :: F(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: F0(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: B_MLS(3,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: X(2,NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_I(NIAC)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_J(NIAC)
              REAL(KIND=8), INTENT(IN) :: W(NIAC)
              REAL(KIND=8), INTENT(IN) :: W_AA(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: MASS(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: RHO(NTOTAL)
              INTEGER(KIND=2), INTENT(IN) :: ITYPE(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: GAMMA_CONT(NTOTAL)
            END SUBROUTINE MLSFUNCTIONAPPROX
          END INTERFACE 
        END MODULE MLSFUNCTIONAPPROX__genmod
