        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  4 12:45:24 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MLS_AA_MATRIX_LINEAR__genmod
          INTERFACE 
            SUBROUTINE MLS_AA_MATRIX_LINEAR(AA,XA,XB,W,MASSB,RHOB,DIM)
              INTEGER(KIND=4), INTENT(IN) :: DIM
              REAL(KIND=8), INTENT(INOUT) :: AA(DIM+1,DIM+1)
              REAL(KIND=8), INTENT(IN) :: XA(DIM)
              REAL(KIND=8), INTENT(IN) :: XB(DIM)
              REAL(KIND=8), INTENT(IN) :: W
              REAL(KIND=8), INTENT(IN) :: MASSB
              REAL(KIND=8), INTENT(IN) :: RHOB
            END SUBROUTINE MLS_AA_MATRIX_LINEAR
          END INTERFACE 
        END MODULE MLS_AA_MATRIX_LINEAR__genmod
