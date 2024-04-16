        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:39 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FLUIDMOMENTUM__genmod
          INTERFACE 
            SUBROUTINE FLUIDMOMENTUM(DSTRESS,X,VX,RHO,P,ITYPE,NTOTAL)
              USE PARTICLE_DATA, ONLY :                                 &
     &          MAXN,                                                   &
     &          PBC_EDGES
              INTEGER(KIND=4), INTENT(IN) :: NTOTAL
              REAL(KIND=8), INTENT(INOUT) :: DSTRESS(2,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: X(2,MAXN)
              REAL(KIND=8), INTENT(IN) :: VX(2,MAXN)
              REAL(KIND=8), INTENT(IN) :: RHO(MAXN)
              REAL(KIND=8), INTENT(IN) :: P(MAXN)
              INTEGER(KIND=2), INTENT(IN) :: ITYPE(MAXN)
            END SUBROUTINE FLUIDMOMENTUM
          END INTERFACE 
        END MODULE FLUIDMOMENTUM__genmod
