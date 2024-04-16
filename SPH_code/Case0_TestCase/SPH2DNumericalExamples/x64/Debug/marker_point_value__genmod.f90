        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:38 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MARKER_POINT_VALUE__genmod
          INTERFACE 
            SUBROUTINE MARKER_POINT_VALUE(F)
              USE PARTICLE_DATA, ONLY :                                 &
     &          ITYPE,                                                  &
     &          W,                                                      &
     &          NTOTAL,                                                 &
     &          NIAC,                                                   &
     &          PAIR_I,                                                 &
     &          PAIR_J,                                                 &
     &          GAMMA_DISCRT,                                           &
     &          RHO,                                                    &
     &          MASS
              REAL(KIND=8) :: F(NTOTAL)
            END SUBROUTINE MARKER_POINT_VALUE
          END INTERFACE 
        END MODULE MARKER_POINT_VALUE__genmod
