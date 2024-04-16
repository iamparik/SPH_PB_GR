        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:34 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE VISCOUSSTRESSOPERATOR__genmod
          INTERFACE 
            SUBROUTINE VISCOUSSTRESSOPERATOR(DSTRESS,FNCN,RHO,OPRTRTYYPE&
     &)
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
     &          SURF_NORM,                                              &
     &          GAMMA_CONT,                                             &
     &          GAMMA_DISCRT,                                           &
     &          GAMMA_MAT,                                              &
     &          GAMMA_MAT_INV,                                          &
     &          XI1_MAT_INV,                                            &
     &          XI_CONT_MAT_INV,                                        &
     &          BDRYVAL_VEL,                                            &
     &          MASS,                                                   &
     &          MU,                                                     &
     &          X,                                                      &
     &          PBC_EDGES
              REAL(KIND=8) :: DSTRESS(2,NTOTAL)
              REAL(KIND=8) :: FNCN(2,NTOTAL)
              REAL(KIND=8) :: RHO(NTOTAL)
              INTEGER(KIND=4) :: OPRTRTYYPE
            END SUBROUTINE VISCOUSSTRESSOPERATOR
          END INTERFACE 
        END MODULE VISCOUSSTRESSOPERATOR__genmod
