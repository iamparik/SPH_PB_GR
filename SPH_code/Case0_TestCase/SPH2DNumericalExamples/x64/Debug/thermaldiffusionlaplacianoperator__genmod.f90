        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:36 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE THERMALDIFFUSIONLAPLACIANOPERATOR__genmod
          INTERFACE 
            SUBROUTINE THERMALDIFFUSIONLAPLACIANOPERATOR(FNCN,LAP_FNCN, &
     &OPRTRTYYPE)
              USE PARTICLE_DATA, ONLY :                                 &
     &          RHO,                                                    &
     &          MASS,                                                   &
     &          ITYPE,                                                  &
     &          X,                                                      &
     &          NREAL,                                                  &
     &          NEDGE,                                                  &
     &          NGHOST,                                                 &
     &          NTOTAL,                                                 &
     &          NIAC,                                                   &
     &          PAIR_I,                                                 &
     &          PAIR_J,                                                 &
     &          W,                                                      &
     &          W_AA,                                                   &
     &          DWDX,                                                   &
     &          ENIAC,                                                  &
     &          EPAIR_A,                                                &
     &          EPAIR_S,                                                &
     &          NEDGE_REL_EDGE,                                         &
     &          EDGE,                                                   &
     &          SURF_NORM,                                              &
     &          GAMMA_CONT,                                             &
     &          GAMMA_MAT,                                              &
     &          GAMMA_DISCRT,                                           &
     &          DEL_GAMMA_AS,                                           &
     &          GAMMA_MAT_INV,                                          &
     &          XI1_MAT_INV,                                            &
     &          XI_CONT_MAT_INV,                                        &
     &          ETOTAL,                                                 &
     &          ETYPE,                                                  &
     &          NPERIODIC,                                              &
     &          PBC_DUPLICATE_PAIR,                                     &
     &          PBC_EDGES,                                              &
     &          BDRYVAL_TEMP
              REAL(KIND=8) :: FNCN(NTOTAL)
              REAL(KIND=8) :: LAP_FNCN(NTOTAL)
              INTEGER(KIND=4) :: OPRTRTYYPE
            END SUBROUTINE THERMALDIFFUSIONLAPLACIANOPERATOR
          END INTERFACE 
        END MODULE THERMALDIFFUSIONLAPLACIANOPERATOR__genmod
