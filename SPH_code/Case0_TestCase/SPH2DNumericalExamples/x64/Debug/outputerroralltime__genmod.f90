        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  4 12:45:22 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE OUTPUTERRORALLTIME__genmod
          INTERFACE 
            SUBROUTINE OUTPUTERRORALLTIME(MAXTIMESTEP,CURRENTSTEP,FX,   &
     &ISIZE,DATAOUTPUTPATH,RUNTYPE,DT)
              INTEGER(KIND=4), INTENT(IN) :: ISIZE
              INTEGER(KIND=4), INTENT(IN) :: MAXTIMESTEP
              INTEGER(KIND=4), INTENT(IN) :: CURRENTSTEP
              REAL(KIND=8), INTENT(IN) :: FX(MAXTIMESTEP,ISIZE)
              CHARACTER(*) :: DATAOUTPUTPATH
              INTEGER(KIND=4), INTENT(IN) :: RUNTYPE
              REAL(KIND=8), INTENT(IN) :: DT
            END SUBROUTINE OUTPUTERRORALLTIME
          END INTERFACE 
        END MODULE OUTPUTERRORALLTIME__genmod
