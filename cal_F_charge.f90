program main
implicit none
    interface 
        subroutine inverse(matrix, matrixInverse)
        implicit none
            complex(kind = 8), intent(in) :: matrix
            complex(kind = 8), intent(out) :: matrixInverse
        end subroutine inverse
    end interface
    integer :: i, j, k, l
    integer :: row, col
    integer :: matrixDimension
    complex(kind = 8), dimension(:, :), allocatable :: gamma_charge
end program main
