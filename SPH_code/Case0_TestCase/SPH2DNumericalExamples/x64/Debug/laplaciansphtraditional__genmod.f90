        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  4 12:45:23 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LAPLACIANSPHTRADITIONAL__genmod
          INTERFACE 
            SUBROUTINE LAPLACIANSPHTRADITIONAL(LAPF,F0,G,X,GMA,SPH_DIM, &
     &NIAC,PAIR_I,PAIR_J,DWDX,MASS,RHO,ITYPE,NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              INTEGER(KIND=4), INTENT(IN) :: NIAC
              INTEGER(KIND=4), INTENT(IN) :: SPH_DIM
              REAL(KIND=8), INTENT(INOUT) :: LAPF(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: F0(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: G(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: X(SPH_DIM,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: GMA(NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_I(NIAC)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_J(NIAC)
              REAL(KIND=8), INTENT(IN) :: DWDX(SPH_DIM,NIAC)
              REAL(KIND=8), INTENT(IN) :: MASS(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: RHO(NTOTAL)
              INTEGER(KIND=2), INTENT(IN) :: ITYPE(NTOTAL)
            END SUBROUTINE LAPLACIANSPHTRADITIONAL
          END INTERFACE 
        END MODULE LAPLACIANSPHTRADITIONAL__genmod
