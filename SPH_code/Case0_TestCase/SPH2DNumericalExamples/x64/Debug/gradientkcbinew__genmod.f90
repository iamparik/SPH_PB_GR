        !COMPILER-GENERATED INTERFACE MODULE: Sat Nov 25 15:01:31 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GRADIENTKCBINEW__genmod
          INTERFACE 
            SUBROUTINE GRADIENTKCBINEW(D1F,F0,FS0)
              USE PARTICLE_DATA, ONLY :                                 &
     &          GAMMA_CONT,                                             &
     &          NIAC,                                                   &
     &          PAIR_I,                                                 &
     &          PAIR_J,                                                 &
     &          DWDX,                                                   &
     &          ENIAC,                                                  &
     &          EPAIR_A,                                                &
     &          EPAIR_S,                                                &
     &          DEL_GAMMA_AS,                                           &
     &          MASS,                                                   &
     &          RHO,                                                    &
     &          ITYPE,                                                  &
     &          NTOTAL,                                                 &
     &          ETOTAL
              REAL(KIND=8), INTENT(INOUT) :: D1F(2,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: F0(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: FS0(ETOTAL)
            END SUBROUTINE GRADIENTKCBINEW
          END INTERFACE 
        END MODULE GRADIENTKCBINEW__genmod
