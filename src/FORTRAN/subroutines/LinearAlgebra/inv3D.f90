!****************************************************************************
!
!  SUBROUTINE: inv3d(A-1,A)
!			   A simple 3D matrix inverter
!
!  PURPOSE:  calculation of the inverse of a 3×3 matrix.
!
!   CREATED:        07/24/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  07/24/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine inv3d(B,A)
    implicit none

    real(8), intent(in), dimension(3,3):: A
    real(8), dimension(3,3):: B
    real(8):: A_det
    integer i,j
    
    
    ! Calculate the inverse determinant of the matrix
    A_det = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
       
    if (isNAN(1/A_det))     write(*,*) 'Warning : matrix is non invertible, caution while using the inverse directly'
    
    ! Calculate the inverse of the matrix
    B(1,1) = +A_det * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -A_det * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +A_det * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -A_det * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +A_det * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -A_det * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +A_det * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -A_det * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +A_det * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    
end   