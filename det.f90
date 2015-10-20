subroutine det(matrix, determinant)
    implicit none
    complex(kind = 8), intent(in) :: matrix(:, :)
    complex(kind = 8), intent(out) :: determinant
    complex(kind = 8), dimension(:, :), allocatable :: tempMatrix
    integer :: matrixDimension
    integer :: i, j
    complex(kind = 8) :: prod
    integer :: info
    integer, dimension(:), allocatable :: IPIV
    
    matrixDimension = size(matrix, 1)
    if (matrixDimension .ne. size(matrix, 2)) then
        write(*, *) "Only square matrix is accepted. "
        stop
    endif
    allocate(tempMatrix(matrixDimension, matrixDimension))
    allocate(IPIV(matrixDimension))

    do i = 1, matrixDimension
    do j = 1, matrixDimension
        tempMatrix(i, j) = matrix(i, j)
    enddo
    enddo

    call zgetrf(matrixDimension, matrixDimension, tempMatrix, matrixDimension, IPIV, info)
    
    prod = dcmplx(1, 0)
    do i = 1, matrixDimension
        prod = prod*tempMatrix(i, i)
    enddo
    determinant = prod

    deallocate(IPIV)
    deallocate(tempMatrix)
end subroutine det
