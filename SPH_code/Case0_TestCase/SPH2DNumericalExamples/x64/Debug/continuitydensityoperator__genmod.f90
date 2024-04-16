        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:31 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CONTINUITYDENSITYOPERATOR__genmod
          INTERFACE 
            SUBROUTINE CONTINUITYDENSITYOPERATOR(DIVFNCN,FNCN,RHO_IN,   &
     &OPRTRTYYPE)
              USE PARTICLE_DATA, ONLY :                                 &
     &          PAIR_I,                                                 &
     &          PAIR_J,                                                 &
     &          NIAC,                                                   &
     &          NTOTAL,                                                 &
     &          ITYPE,                                                  &
     &          DWDX,                                                   &
     &          EPAIR_A,                                                &
     &          EPAIR_S,                                                &
     &          ENIAC,                                                  &
     &          ETOTAL,                                                 &
     &          ETYPE,                                                  &
     &          DEL_GAMMA_AS,                                           &
     &          NEDGE_REL_EDGE,                                         &
     &          GAMMA_CONT,                                             &
     &          GAMMA_DISCRT,                                           &
     &          GAMMA_MAT,                                              &
     &          GAMMA_MAT_INV,                                          &
     &          XI1_MAT_INV,                                            &
     &          XI_CONT_MAT_INV,                                        &
     &          MASS,                                                   &
     &          FREESURFACEVAR,                                         &
     &          GAMMA_DENSITY_CONT
              REAL(KIND=8) :: DIVFNCN(NTOTAL)
              REAL(KIND=8) :: FNCN(2,NTOTAL)
              REAL(KIND=8) :: RHO_IN(NTOTAL)
              INTEGER(KIND=4) :: OPRTRTYYPE
            END SUBROUTINE CONTINUITYDENSITYOPERATOR
          END INTERFACE 
        END MODULE CONTINUITYDENSITYOPERATOR__genmod
