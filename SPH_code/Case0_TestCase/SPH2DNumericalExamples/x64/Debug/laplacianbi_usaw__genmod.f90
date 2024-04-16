        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  4 12:45:19 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LAPLACIANBI_USAW__genmod
          INTERFACE 
            SUBROUTINE LAPLACIANBI_USAW(LAPF,F0,DF0,DFS0,X,KK,KK_S,GMA, &
     &SPH_DIM,NIAC,PAIR_I,PAIR_J,DWDX,ENIAC,EPAIR_A,EPAIR_S,DGMAS,MASS, &
     &RHO,ITYPE,NTOTAL,ETOTAL)
              INTEGER(KIND=4), INTENT(IN) :: ETOTAL
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              INTEGER(KIND=4), INTENT(IN) :: ENIAC
              INTEGER(KIND=4), INTENT(IN) :: NIAC
              INTEGER(KIND=4), INTENT(IN) :: SPH_DIM
              REAL(KIND=8), INTENT(INOUT) :: LAPF(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: F0(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: DF0(SPH_DIM,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: DFS0(SPH_DIM,ETOTAL)
              REAL(KIND=8), INTENT(IN) :: X(SPH_DIM,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: KK(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: KK_S(ETOTAL)
              REAL(KIND=8), INTENT(IN) :: GMA(NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_I(NIAC)
              INTEGER(KIND=4), INTENT(IN) :: PAIR_J(NIAC)
              REAL(KIND=8), INTENT(IN) :: DWDX(SPH_DIM,NIAC)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_A(ENIAC)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_S(ENIAC)
              REAL(KIND=8), INTENT(IN) :: DGMAS(SPH_DIM,ENIAC)
              REAL(KIND=8), INTENT(IN) :: MASS(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: RHO(NTOTAL)
              INTEGER(KIND=2), INTENT(IN) :: ITYPE(NTOTAL)
            END SUBROUTINE LAPLACIANBI_USAW
          END INTERFACE 
        END MODULE LAPLACIANBI_USAW__genmod
