        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:27 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BOUNDARYCONDITIONNTOD__genmod
          INTERFACE 
            SUBROUTINE BOUNDARYCONDITIONNTOD(FNCN_S,BDRYVAL,FNCN)
              USE PARTICLE_DATA, ONLY :                                 &
     &          NTOTAL,                                                 &
     &          ETOTAL,                                                 &
     &          ETYPE,                                                  &
     &          EPAIR_A,                                                &
     &          EPAIR_S,                                                &
     &          ENIAC,                                                  &
     &          HSML,                                                   &
     &          MASS,                                                   &
     &          RHO,                                                    &
     &          NEDGE_REL_EDGE,                                         &
     &          X,                                                      &
     &          SURF_NORM
              REAL(KIND=8) :: FNCN_S(ETOTAL)
              REAL(KIND=8) :: BDRYVAL(ETOTAL)
              REAL(KIND=8) :: FNCN(NTOTAL)
            END SUBROUTINE BOUNDARYCONDITIONNTOD
          END INTERFACE 
        END MODULE BOUNDARYCONDITIONNTOD__genmod
