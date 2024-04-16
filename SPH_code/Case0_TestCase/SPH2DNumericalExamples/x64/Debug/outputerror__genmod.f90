        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  4 12:45:20 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE OUTPUTERROR__genmod
          INTERFACE 
            SUBROUTINE OUTPUTERROR(N,RE,REB,REI,S,DATAOUTPUTPATH,NREAL, &
     &HSML_FACTOR_NAME,HFC)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: RE(2,N)
              REAL(KIND=8), INTENT(IN) :: REB(2,N)
              REAL(KIND=8), INTENT(IN) :: REI(2,N)
              REAL(KIND=8), INTENT(IN) :: S
              CHARACTER(*) :: DATAOUTPUTPATH
              INTEGER(KIND=4) :: NREAL
              CHARACTER(LEN=9), INTENT(IN) :: HSML_FACTOR_NAME
              INTEGER(KIND=4), INTENT(IN) :: HFC
            END SUBROUTINE OUTPUTERROR
          END INTERFACE 
        END MODULE OUTPUTERROR__genmod
