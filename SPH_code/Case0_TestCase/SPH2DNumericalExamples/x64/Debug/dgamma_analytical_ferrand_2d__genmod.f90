        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  4 12:45:26 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DGAMMA_ANALYTICAL_FERRAND_2D__genmod
          INTERFACE 
            SUBROUTINE DGAMMA_ANALYTICAL_FERRAND_2D(XA,XEI,XE_NORM,     &
     &DGAMMA_EDGE,H,PI,SCALE_K)
              REAL(KIND=8) :: XA(2)
              REAL(KIND=8) :: XEI(2,2)
              REAL(KIND=8), INTENT(IN) :: XE_NORM(2)
              REAL(KIND=8), INTENT(OUT) :: DGAMMA_EDGE(2)
              REAL(KIND=8), INTENT(IN) :: H
              REAL(KIND=8), INTENT(IN) :: PI
              REAL(KIND=8), INTENT(IN) :: SCALE_K
            END SUBROUTINE DGAMMA_ANALYTICAL_FERRAND_2D
          END INTERFACE 
        END MODULE DGAMMA_ANALYTICAL_FERRAND_2D__genmod
