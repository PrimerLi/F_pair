subroutine eigenvalues(matrix, eigenvalue)
    implicit none
    complex(kind = 8), intent(in) :: matrix(:, :)
    complex(kind = 8), intent(out) :: eigenvalue(:)
    integer :: matrixDimension
    complex(kind = 8), dimension(:, :), allocatable :: tempMatrix
    complex(kind = 8), dimension(:), allocatable :: tempVector, VL
    complex(kind = 8), dimension(:, :), allocatable :: rightEigenvectors
    integer :: info, lwork
    complex(kind = 8), dimension(:), allocatable :: work
    real(kind = 8), dimension(:), allocatable :: rwork
    integer :: i, j, k, l
    
    matrixDimension = size(matrix, 1)
    if (size(matrix, 2) .ne. matrixDimension) then
        write(*, *) "Only square matrix is accepted. "
        stop
    endif
    if (size(matrix, 1) .ne. size(eigenvalue, 1)) then
        write(*, *) "Dimension incompatible. "
        stop
    endif

    lwork = 2*matrixDimension
    allocate(tempMatrix(matrixDimension, matrixDimension))
    allocate(tempVector(matrixDimension), VL(matrixDimension))
    allocate(work(lwork))
    allocate(rwork(2*matrixDimension))
    allocate(rightEigenvectors(matrixDimension, matrixDimension))

    do i = 1, matrixDimension
    do j = 1, matrixDimension
        tempMatrix(i, j) = matrix(i, j)
    enddo
    enddo
    call zgeev("N", "V", matrixDimension, tempMatrix, matrixDimension, &
            &tempVector, VL, matrixDimension, rightEigenvectors, &
            &matrixDimension, work, lwork, rwork, info)

    do i = 1, matrixDimension
        eigenvalue(i) = tempVector(i)
    enddo

    deallocate(tempMatrix)
    deallocate(tempVector, VL)
    deallocate(work, rwork)
    deallocate(rightEigenvectors)
end subroutine eigenvalues
