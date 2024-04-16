        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  4 12:45:28 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FNCNAPPROXOPERATOR__genmod
          INTERFACE 
            SUBROUTINE FNCNAPPROXOPERATOR(FNCN_A,FNCN_B,MASS_B,RHO_B,   &
     &WEIGHT)
              REAL(KIND=8) :: FNCN_A
              REAL(KIND=8), INTENT(IN) :: FNCN_B
              REAL(KIND=8), INTENT(IN) :: MASS_B
              REAL(KIND=8), INTENT(IN) :: RHO_B
              REAL(KIND=8), INTENT(IN) :: WEIGHT
            END SUBROUTINE FNCNAPPROXOPERATOR
          END INTERFACE 
        END MODULE FNCNAPPROXOPERATOR__genmod
