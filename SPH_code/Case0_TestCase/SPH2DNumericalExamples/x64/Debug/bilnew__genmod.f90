        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 16 13:59:08 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BILNEW__genmod
          INTERFACE 
            SUBROUTINE BILNEW(LAPF,DFS0,GMA,SPH_DIM,ENIAC,EPAIR_A,      &
     &EPAIR_S,DGMAS,ETYPE,ETOTAL,NTOTAL,ETYPE_PERIODIC,NEUMANN,DIRICHLET&
     &,BDRYVAL)
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              INTEGER(KIND=4), INTENT(IN) :: ETOTAL
              INTEGER(KIND=4), INTENT(IN) :: ENIAC
              INTEGER(KIND=4), INTENT(IN) :: SPH_DIM
              REAL(KIND=8), INTENT(INOUT) :: LAPF(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: DFS0(SPH_DIM,ETOTAL)
              REAL(KIND=8), INTENT(IN) :: GMA(NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_A(ENIAC)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_S(ENIAC)
              REAL(KIND=8), INTENT(IN) :: DGMAS(SPH_DIM,ENIAC)
              INTEGER(KIND=4), INTENT(IN) :: ETYPE(ETOTAL)
              INTEGER(KIND=4), INTENT(IN) :: ETYPE_PERIODIC
              INTEGER(KIND=4), INTENT(IN) :: NEUMANN
              INTEGER(KIND=4), INTENT(IN) :: DIRICHLET
              REAL(KIND=8), INTENT(IN) :: BDRYVAL(ETOTAL)
            END SUBROUTINE BILNEW
          END INTERFACE 
        END MODULE BILNEW__genmod
