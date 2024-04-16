        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:34 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BILUSAW__genmod
          INTERFACE 
            SUBROUTINE BILUSAW(LAPF,DF0,DFS0,GMA,SPH_DIM,ENIAC,EPAIR_A, &
     &EPAIR_S,DGMAS,ETYPE,ETOTAL,NTOTAL,ETYPE_PERIODIC,NEUMANN,DIRICHLET&
     &,BDRYVAL,BDRYVALDIM)
              INTEGER(KIND=4), INTENT(IN) :: BDRYVALDIM
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              INTEGER(KIND=4), INTENT(IN) :: ETOTAL
              INTEGER(KIND=4), INTENT(IN) :: ENIAC
              INTEGER(KIND=4), INTENT(IN) :: SPH_DIM
              REAL(KIND=8), INTENT(INOUT) :: LAPF(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: DF0(SPH_DIM,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: DFS0(SPH_DIM,ETOTAL)
              REAL(KIND=8), INTENT(IN) :: GMA(NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_A(ENIAC)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_S(ENIAC)
              REAL(KIND=8), INTENT(IN) :: DGMAS(SPH_DIM,ENIAC)
              INTEGER(KIND=4), INTENT(IN) :: ETYPE(ETOTAL)
              INTEGER(KIND=4), INTENT(IN) :: ETYPE_PERIODIC
              INTEGER(KIND=4), INTENT(IN) :: NEUMANN
              INTEGER(KIND=4), INTENT(IN) :: DIRICHLET
              REAL(KIND=8), INTENT(IN) :: BDRYVAL(BDRYVALDIM)
            END SUBROUTINE BILUSAW
          END INTERFACE 
        END MODULE BILUSAW__genmod
