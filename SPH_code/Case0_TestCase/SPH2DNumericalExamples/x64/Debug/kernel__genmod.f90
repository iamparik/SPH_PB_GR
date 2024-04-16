        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  9 18:25:37 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE KERNEL__genmod
          INTERFACE 
            SUBROUTINE KERNEL(R,DX,HSML,W,DWDX)
              REAL(KIND=8) :: R
              REAL(KIND=8) :: DX(2)
              REAL(KIND=8) :: HSML
              REAL(KIND=8) :: W
              REAL(KIND=8) :: DWDX(2)
            END SUBROUTINE KERNEL
          END INTERFACE 
        END MODULE KERNEL__genmod
