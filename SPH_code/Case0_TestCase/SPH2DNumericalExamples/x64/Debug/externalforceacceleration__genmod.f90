        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:28 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EXTERNALFORCEACCELERATION__genmod
          INTERFACE 
            SUBROUTINE EXTERNALFORCEACCELERATION(DSTRESS,RHO,ITYPE,     &
     &NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              REAL(KIND=8) :: DSTRESS(2,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: RHO(NTOTAL)
              INTEGER(KIND=2), INTENT(IN) :: ITYPE(NTOTAL)
            END SUBROUTINE EXTERNALFORCEACCELERATION
          END INTERFACE 
        END MODULE EXTERNALFORCEACCELERATION__genmod
