        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:32 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIVERGENCESUMMATIONBDRY__genmod
          INTERFACE 
            SUBROUTINE DIVERGENCESUMMATIONBDRY(D1F,F0,RHO_IN)
              USE PARTICLE_DATA, ONLY :                                 &
     &          PAIR_I,                                                 &
     &          PAIR_J,                                                 &
     &          NIAC,                                                   &
     &          NTOTAL,                                                 &
     &          W,                                                      &
     &          W_AA,                                                   &
     &          DWDX,                                                   &
     &          MASS,                                                   &
     &          GAMMA_CONT,                                             &
     &          GAMMA_DISCRT,                                           &
     &          ITYPE,                                                  &
     &          MAXN,                                                   &
     &          DGRHO_PREV
              REAL(KIND=8), INTENT(INOUT) :: D1F(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: F0(2,NTOTAL)
              REAL(KIND=8), INTENT(INOUT) :: RHO_IN(NTOTAL)
            END SUBROUTINE DIVERGENCESUMMATIONBDRY
          END INTERFACE 
        END MODULE DIVERGENCESUMMATIONBDRY__genmod
