        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:30 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SUMMATIONDENSITYOPERATOR__genmod
          INTERFACE 
            SUBROUTINE SUMMATIONDENSITYOPERATOR(RHO,OPRTTYPE)
              USE PARTICLE_DATA, ONLY :                                 &
     &          PAIR_I,                                                 &
     &          PAIR_J,                                                 &
     &          NIAC,                                                   &
     &          NTOTAL,                                                 &
     &          W,                                                      &
     &          W_AA,                                                   &
     &          MASS,                                                   &
     &          GAMMA_CONT,                                             &
     &          GAMMA_DISCRT,                                           &
     &          ITYPE,                                                  &
     &          MAXN,                                                   &
     &          DGRHO_PREV,                                             &
     &          GAMMA_DENSITY_CONT,                                     &
     &          FREESURFACEVAR
              REAL(KIND=8) :: RHO(MAXN)
              INTEGER(KIND=4), INTENT(IN) :: OPRTTYPE
            END SUBROUTINE SUMMATIONDENSITYOPERATOR
          END INTERFACE 
        END MODULE SUMMATIONDENSITYOPERATOR__genmod
