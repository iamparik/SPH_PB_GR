        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:33 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RELERRORFORMULATION__genmod
          INTERFACE 
            SUBROUTINE RELERRORFORMULATION(F_ANALYTICAL,F_SPH,RE,REB,REI&
     &,ERRORTYPE)
              USE PARTICLE_DATA, ONLY :                                 &
     &          GAMMA_CONT,                                             &
     &          NTOTAL,                                                 &
     &          ITYPE,                                                  &
     &          MASS,                                                   &
     &          RHO
              REAL(KIND=8), INTENT(IN) :: F_ANALYTICAL(NTOTAL)
              REAL(KIND=8), INTENT(IN) :: F_SPH(NTOTAL)
              REAL(KIND=8) :: RE(2)
              REAL(KIND=8) :: REB(2)
              REAL(KIND=8) :: REI(2)
              LOGICAL(KIND=4) :: ERRORTYPE
            END SUBROUTINE RELERRORFORMULATION
          END INTERFACE 
        END MODULE RELERRORFORMULATION__genmod
