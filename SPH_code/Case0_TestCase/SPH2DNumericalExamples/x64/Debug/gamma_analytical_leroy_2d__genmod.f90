        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  4 12:45:26 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GAMMA_ANALYTICAL_LEROY_2D__genmod
          INTERFACE 
            SUBROUTINE GAMMA_ANALYTICAL_LEROY_2D(XA,XV,S_NORM,GAMMA_AS,H&
     &,PI)
              REAL(KIND=8), INTENT(IN) :: XA(2)
              REAL(KIND=8), INTENT(IN) :: XV(2,2)
              REAL(KIND=8), INTENT(IN) :: S_NORM(2)
              REAL(KIND=8), INTENT(OUT) :: GAMMA_AS
              REAL(KIND=8), INTENT(IN) :: H
              REAL(KIND=8), INTENT(IN) :: PI
            END SUBROUTINE GAMMA_ANALYTICAL_LEROY_2D
          END INTERFACE 
        END MODULE GAMMA_ANALYTICAL_LEROY_2D__genmod
