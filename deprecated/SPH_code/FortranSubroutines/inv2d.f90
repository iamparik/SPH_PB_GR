!****************************************************************************
!
!  SUBROUTINE: inv2d(A,A-1)
!			   A simple 2D matrix inverter
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        05/05/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  05/05/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine inv2d(B,A)
    implicit none

    real(8), intent(in), dimension(2,2):: A
    real(8), dimension(2,2):: B, A_cofactor
    real(8):: A_det
    integer i,j
    
    A_cofactor(1,1)=A(2,2)
    A_cofactor(2,2)=A(1,1)
    A_cofactor(1,2)=-A(1,2)
    A_cofactor(2,1)=-A(2,1)
    
    A_det=A(2,2)*A(1,1)-A(1,2)*A(2,1)
    
    if (isNAN(1/A_det)) then
        write(*,*) 'Warning : matrix is non invertible, caution while using the inverse directly'
    endif
    
    do i=1,2
        do j=1,2
            B(i,j)=A_cofactor(i,j)/A_det
		enddo
	enddo
 end   
            
    
    
    
    
    