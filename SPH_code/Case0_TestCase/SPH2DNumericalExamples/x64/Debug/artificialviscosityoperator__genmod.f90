        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:34 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ARTIFICIALVISCOSITYOPERATOR__genmod
          INTERFACE 
            SUBROUTINE ARTIFICIALVISCOSITYOPERATOR(DSTRESS,VX,X,RHO,    &
     &OPRTRTYPE)
              USE PARTICLE_DATA, ONLY :                                 &
     &          PAIR_I,                                                 &
     &          PAIR_J,                                                 &
     &          NIAC,                                                   &
     &          NTOTAL,                                                 &
     &          ITYPE,                                                  &
     &          DWDX,                                                   &
     &          MASS
              REAL(KIND=8) :: DSTRESS(2,NTOTAL)
              REAL(KIND=8) :: VX(2,NTOTAL)
              REAL(KIND=8) :: X(2,NTOTAL)
              REAL(KIND=8) :: RHO(NTOTAL)
              INTEGER(KIND=4) :: OPRTRTYPE
            END SUBROUTINE ARTIFICIALVISCOSITYOPERATOR
          END INTERFACE 
        END MODULE ARTIFICIALVISCOSITYOPERATOR__genmod
