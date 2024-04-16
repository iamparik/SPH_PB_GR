        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:36 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BILMACIA__genmod
          INTERFACE 
            SUBROUTINE BILMACIA(LAPF,F0,FS0,X,XS,GMA,SPH_DIM,ENIAC,     &
     &EPAIR_A,EPAIR_S,DGMAS,ETYPE,ETOTAL,NTOTAL,KK_S)
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              INTEGER(KIND=4), INTENT(IN) :: ETOTAL
              INTEGER(KIND=4), INTENT(IN) :: ENIAC
              INTEGER(KIND=4), INTENT(IN) :: SPH_DIM
              REAL(KIND=8), INTENT(INOUT) :: LAPF(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: F0(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: FS0(ETOTAL)
              REAL(KIND=8), INTENT(IN) :: X(SPH_DIM,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: XS(SPH_DIM,ETOTAL)
              REAL(KIND=8), INTENT(IN) :: GMA(NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_A(ENIAC)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_S(ENIAC)
              REAL(KIND=8), INTENT(IN) :: DGMAS(SPH_DIM,ENIAC)
              INTEGER(KIND=4), INTENT(IN) :: ETYPE(ETOTAL)
              REAL(KIND=8), INTENT(IN) :: KK_S(ETOTAL)
            END SUBROUTINE BILMACIA
          END INTERFACE 
        END MODULE BILMACIA__genmod
