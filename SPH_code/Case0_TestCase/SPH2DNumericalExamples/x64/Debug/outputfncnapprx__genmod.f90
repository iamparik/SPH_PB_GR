        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:31 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE OUTPUTFNCNAPPRX__genmod
          INTERFACE 
            SUBROUTINE OUTPUTFNCNAPPRX(N,X,FNCN0,FNCN,DFNCN0,DFNCN,     &
     &LAP_FNCN0,LAP_FNCN,NEDGE_REL_EDGE,EDGE,ETOTAL,SPH_DIM,NTOTAL,ITYPE&
     &,ETYPE,DATAOUTPUTPATH)
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              INTEGER(KIND=4), INTENT(IN) :: SPH_DIM
              INTEGER(KIND=4), INTENT(IN) :: ETOTAL
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: X(SPH_DIM,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: FNCN0(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: FNCN(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: DFNCN0(SPH_DIM,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: DFNCN(SPH_DIM,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: LAP_FNCN0(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: LAP_FNCN(NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: NEDGE_REL_EDGE(ETOTAL)
              INTEGER(KIND=4), INTENT(IN) :: EDGE(SPH_DIM,ETOTAL)
              INTEGER(KIND=2), INTENT(IN) :: ITYPE(NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: ETYPE(ETOTAL)
              CHARACTER(*) :: DATAOUTPUTPATH
            END SUBROUTINE OUTPUTFNCNAPPRX
          END INTERFACE 
        END MODULE OUTPUTFNCNAPPRX__genmod
