        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  4 12:45:20 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BIL1__genmod
          INTERFACE 
            SUBROUTINE BIL1(LAPF,DF0,GMA,SPH_DIM,ENIAC,EPAIR_A,EPAIR_S, &
     &DGMAS,ETYPE,ETOTAL,NTOTAL,ETYPE_PERIODIC,NEUMANN,DIRICHLET,BDRYVAL&
     &,COEFF)
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              INTEGER(KIND=4), INTENT(IN) :: ETOTAL
              INTEGER(KIND=4), INTENT(IN) :: ENIAC
              INTEGER(KIND=4), INTENT(IN) :: SPH_DIM
              REAL(KIND=8), INTENT(INOUT) :: LAPF(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: DF0(SPH_DIM,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: GMA(NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_A(ENIAC)
              INTEGER(KIND=4), INTENT(IN) :: EPAIR_S(ENIAC)
              REAL(KIND=8), INTENT(IN) :: DGMAS(SPH_DIM,ENIAC)
              INTEGER(KIND=4), INTENT(IN) :: ETYPE(ETOTAL)
              INTEGER(KIND=4), INTENT(IN) :: ETYPE_PERIODIC
              INTEGER(KIND=4), INTENT(IN) :: NEUMANN
              INTEGER(KIND=4), INTENT(IN) :: DIRICHLET
              REAL(KIND=8), INTENT(IN) :: BDRYVAL(ETOTAL)
              REAL(KIND=8), INTENT(IN) :: COEFF
            END SUBROUTINE BIL1
          END INTERFACE 
        END MODULE BIL1__genmod
