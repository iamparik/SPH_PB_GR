        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:38 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DENSITYDIFFUSIONOPERATOR__genmod
          INTERFACE 
            SUBROUTINE DENSITYDIFFUSIONOPERATOR(DRHO,RHO,DENSDIFFTYPE,  &
     &DELTA_SPH)
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
     &          X,                                                      &
     &          PBC_EDGES
              REAL(KIND=8), INTENT(INOUT) :: DRHO(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: RHO(NTOTAL)
              INTEGER(KIND=4), INTENT(IN) :: DENSDIFFTYPE
              REAL(KIND=8), INTENT(IN) :: DELTA_SPH
            END SUBROUTINE DENSITYDIFFUSIONOPERATOR
          END INTERFACE 
        END MODULE DENSITYDIFFUSIONOPERATOR__genmod
