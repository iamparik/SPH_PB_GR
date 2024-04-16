        !COMPILER-GENERATED INTERFACE MODULE: Wed Nov 29 16:47:18 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GRADIENTBIF__genmod
          INTERFACE 
            SUBROUTINE GRADIENTBIF(D1F,F0,FS0)
              USE PARTICLE_DATA, ONLY :                                 &
     &          NIAC,                                                   &
     &          PAIR_I,                                                 &
     &          PAIR_J,                                                 &
     &          DWDX,                                                   &
     &          X,                                                      &
     &          ENIAC,                                                  &
     &          EPAIR_A,                                                &
     &          EPAIR_S,                                                &
     &          DEL_GAMMA_AS,                                           &
     &          MASS,                                                   &
     &          RHO,                                                    &
     &          ITYPE,                                                  &
     &          NTOTAL,                                                 &
     &          ETOTAL,                                                 &
     &          NEDGE_REL_EDGE,                                         &
     &          SURF_NORM,                                              &
     &          GAMMA_CONT
              REAL(KIND=8), INTENT(INOUT) :: D1F(2,NTOTAL)
              REAL(KIND=8), INTENT(IN) :: F0(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: FS0(ETOTAL)
            END SUBROUTINE GRADIENTBIF
          END INTERFACE 
        END MODULE GRADIENTBIF__genmod
