program main
    implicit none
    interface
        subroutine inverse(A, AInverse)
        implicit none
            complex(kind = 8), intent(in) :: A(:, :)
            complex(kind = 8), intent(out) :: AInverse(:, :)
        end subroutine inverse
        subroutine eigenvalues(matrix, eigenvalue)
        implicit none
            complex(kind = 8), intent(in) :: matrix(:, :)
            complex(kind = 8), intent(out) :: eigenvalue(:)
        end subroutine eigenvalues
        subroutine lev(vector, leading, leadingIndex)
        implicit none
            complex(kind = 8), intent(in) :: vector(:)
            complex(kind = 8), intent(out) :: leading
            integer, intent(out) :: leadingIndex
        end subroutine lev
    end interface
    integer :: Niom
    real(kind = 8) :: beta, Ueff, omega0
    integer :: Niom_factor
    complex(kind = 8), dimension(:, :, :), allocatable :: giom_average
    complex(kind = 8), dimension(:, :, :), allocatable :: g_temp
    real(kind = 8), dimension(:), allocatable :: xiom
    integer :: nw, nw1, i, j, k, l, row, col
    integer :: Niom2
    integer :: Niom2_factor
    complex(kind = 8), dimension(:, :, :), allocatable :: Chi_p_0_i
    complex(kind = 8), dimension(:, :, :, :), allocatable :: Chi_p_i
    integer :: Ndim
    real(kind = 8) :: x, x1
    complex(kind = 8) :: z
    complex(kind = 8), dimension(:, :), allocatable :: Z_hlp1, Z_hlp2
    complex(kind = 8), dimension(:, :), allocatable :: Z_hlp1_inv, Z_hlp2_inv
    complex(kind = 8), dimension(:, :), allocatable :: gamma_pair, F_pair
    integer :: delta
    integer :: L_dim, norb
    integer :: matrixDimension
    integer :: nu, mu, ku, lu
    integer :: narg
    character *100 :: buffer
    complex(kind = 8), dimension(:, :, :, :), allocatable :: Chi_0_l_pair
    integer :: ixq, nxq
    real(kind = 8), dimension(:), allocatable :: xq
    complex(kind = 8), dimension(:, :), allocatable :: Chi_0_l_p_matrix
    complex(kind = 8), dimension(:, :), allocatable :: resultMatrix
    complex(kind = 8), dimension(:), allocatable :: F_pair_eigenvalues
    complex(kind = 8) :: F_pair_lev
    integer :: F_pair_lev_index
    complex(kind = 8), dimension(:, :), allocatable :: intermediate
    complex(kind = 8), dimension(:, :), allocatable :: intermediate_two
    complex(kind = 8) :: s

    narg = command_argument_count()
    if (narg .ne. 1) then
        write(*, *) "beta is the first argument. "
        stop
    endif

    call getarg(1, buffer)
    read(buffer, *) beta

    Niom_factor = 5
    Niom2_factor = 2
    Niom = max(100, int(Niom_factor*beta))
    Niom2 = max(20, int(Niom2_factor*beta))
    Ndim = 2*Niom2
    norb = 2
    L_dim = 4
    matrixDimension = Ndim*L_dim
    nxq = 2

    allocate(xiom(Niom))
    allocate(giom_average(norb, norb, Niom))
    allocate(g_temp(norb, norb, -Niom+1:Niom))
    allocate(Chi_p_i(L_dim, L_dim, -Niom2+1:Niom2, -Niom2+1:Niom2))
    allocate(Chi_p_0_i(L_dim, L_dim, -Niom2+1:Niom2))
    allocate(Z_hlp1(Ndim*L_dim, Ndim*L_dim), Z_hlp2(matrixDimension, matrixDimension))
    allocate(Z_hlp1_inv(matrixDimension, matrixDimension), Z_hlp2_inv(matrixDimension, matrixDimension))
    allocate(gamma_pair(matrixDimension, matrixDimension))
    allocate(F_pair(matrixDimension, matrixDimension))
    allocate(Chi_0_l_pair(nxq, -Niom2+1:Niom2, L_dim, L_dim))
    allocate(xq(nxq))
    allocate(Chi_0_l_p_matrix(matrixDimension, matrixDimension))
    allocate(resultMatrix(matrixDimension, matrixDimension))
    allocate(F_pair_eigenvalues(matrixDimension))
    allocate(intermediate(matrixDimension, matrixDimension))
    allocate(intermediate_two(matrixDimension, matrixDimension))

    open(unit = 23, file = "Chi_0_l_pair", action = "read")
    do nw = -Niom2+1, Niom2
        do ixq = 1, nxq
            read(23, *) xq(ixq)
            do i = 1, L_dim
            do j = 1, L_dim
                read(23, *) row, col
                read(23, *) Chi_0_l_pair(ixq, nw, row, col)
            enddo
            enddo
        enddo
    enddo
    close(23)
    
    open(unit = 13, file = "giom_average", action = "read")
    do i = 1, norb
    do j = 1, norb
        read(13, *) row, col
        do nw = 1, Niom
            read(13, *) xiom(nw), giom_average(row, col, nw)
        enddo
    enddo
    enddo
    close(13)

    do i = 1, norb
    do j = 1, norb
        do nw = 1, Niom
            g_temp(i, j, nw) = giom_average(i, j, nw)
            g_temp(i, j, -nw + 1) = conjg(giom_average(i, j, nw))
        enddo
    enddo
    enddo

    do nw = -Niom2 + 1, Niom2
        i = 0
        do nu = 1, norb
        do mu = 1, norb
            i = i + 1
            j = 0
            do ku = 1, norb
            do lu = 1, norb
                j = j + 1
                Chi_p_0_i(i, j, nw) = g_temp(nu, ku, -nw+1)*g_temp(mu, lu, nw)
            enddo
            enddo
        enddo
        enddo
    enddo

    open(unit = 13, file = "Chi_p_average", action = "read")
    do i = 1, L_dim
    do j = 1, L_dim
        read(13, *) row, col
        do nw = -Niom2+1, Niom2
        do nw1 = 1, Niom2
            read(13, *) x, x1
            read(13, *) z
            Chi_p_i(i, j, nw, nw1) = z
            Chi_p_i(i, j, -nw + 1, -nw1+1) = conjg(z)
        enddo
        enddo
    enddo
    enddo
    close(13)

    Z_hlp1 = dcmplx(0, 0)
    Z_hlp2 = dcmplx(0, 0)
    Z_hlp1_inv = dcmplx(0, 0)
    Z_hlp2_inv = dcmplx(0, 0)

    do i = 1, L_dim
    do j = 1, L_dim
        do nw = -Niom2+1, Niom2
            k = (i - 1)*Ndim + nw + Niom2
            l = (j - 1)*Ndim + nw + Niom2
            Z_hlp1(k, l) = conjg(Chi_p_0_i(i, j, nw))
        enddo
        do nw = -Niom2+1, Niom2
        do nw1 = -Niom2+1, Niom2
            k = (i - 1)*Ndim + nw + Niom2
            l = (j - 1)*Ndim + nw1 + Niom2
            Z_hlp2(k, l) = Chi_p_i(i, j, nw, nw1)
        enddo
        enddo
    enddo
    enddo

    call inverse(Z_hlp1, Z_hlp1_inv)
    call inverse(Z_hlp2, Z_hlp2_inv)

    gamma_pair = Z_hlp1_inv - Z_hlp2_inv

    do i = 1, matrixDimension
    do j = 1, matrixDimension
        s = dcmplx(0, 0)
    	do k = 1, matrixDimension
	        s = s + Z_hlp1_inv(i, k)*Z_hlp2(k, j)
	    enddo
	    intermediate(i, j) = s
    enddo
    enddo

    do i = 1, matrixDimension
    do j = 1, matrixDimension
    	s = dcmplx(0, 0)
    	do k = 1, matrixDimension
	        s = s + intermediate(i, k) * Z_hlp1_inv(k, j)
	    enddo
	    intermediate_two(i, j) = s
    enddo
    enddo

    do i = 1, matrixDimension
    do j = 1, matrixDimension
    	F_pair(i, j) = intermediate_two(i, j) - Z_hlp1_inv(i, j)
    enddo
    enddo

    if (.false.) then
        open(unit = 23, file = "F_pair", action = "write")
        do i = 1, matrixDimension
        do j = 1, matrixDimension
            write(23, *) i, j, F_pair(i, j)
        enddo
        enddo
        close(23)

        Z_hlp2 = dcmplx(0, 0)
        intermediate = dcmplx(0, 0)
        do i = 1, matrixDimension
        do j = 1, matrixDimension
            do k = 1, matrixDimension
                Z_hlp2(i, j) = Z_hlp2(i, j) + gamma_pair(i, k)*Z_hlp1(k, j)
            enddo
            if (i == j) then
                intermediate(i, j) = dcmplx(1, 0) - Z_hlp2(i, j)
            else 
                intermediate(i, j) = -Z_hlp2(i, j)
            endif
        enddo
        enddo
        intermediate_two = dcmplx(0, 0)
        call inverse(intermediate, intermediate_two)
        intermediate = dcmplx(0, 0)
        do i = 1, matrixDimension
        do j = 1, matrixDimension
            do k = 1, matrixDimension
                intermediate(i, j) = intermediate(i, j) + &
                &intermediate_two(i, k)*gamma_pair(k, j)
            enddo
        enddo
        enddo
        open(unit = 12, file = "F_pair_gamma", action = "write")
        do i = 1, matrixDimension
        do j = 1, matrixDimension
            write(12, *) i, j, intermediate(i, j)
        enddo
        enddo
        close(12)
    endif

    Chi_0_l_p_matrix = dcmplx(0, 0)
    do i = 1, L_dim
    do j = 1, L_dim
        do nw = -Niom2+1, Niom2
            k = (i-1)*Ndim + nw + Niom2
            l = (j-1)*Ndim + nw + Niom2
            Chi_0_l_p_matrix(k, l) = conjg(Chi_0_l_pair(1, nw, i, j))
        enddo
    enddo
    enddo

    Z_hlp2 = Chi_0_l_p_matrix - Z_hlp1 !Z_hlp1 is Chi_p_0_i

    do i = 1, matrixDimension
    do j = 1, matrixDimension
        resultMatrix(i, j) = dcmplx(0, 0)
        do k = 1, matrixDimension
            resultMatrix(i, j) = resultMatrix(i, j) + F_pair(i, k) * Z_hlp2(k, j)
        enddo
    enddo
    enddo

    call eigenvalues(resultMatrix, F_pair_eigenvalues)
    call lev(F_pair_eigenvalues, F_pair_lev, F_pair_lev_index)

    open(unit = 12, file = "F_pair_eigenvalues", action = "write")
    do i = 1, matrixDimension
        write(12, *) i, F_pair_eigenvalues(i)
    enddo
    close(12)

    write(*, *) beta, F_pair_lev

    deallocate(xiom)
    deallocate(giom_average)
    deallocate(g_temp)
    deallocate(Z_hlp1, Z_hlp2)
    deallocate(Z_hlp1_inv, Z_hlp2_inv)
    deallocate(Chi_p_i)
    deallocate(Chi_p_0_i)
    deallocate(gamma_pair, F_pair)
    deallocate(Chi_0_l_pair)
    deallocate(xq)
    deallocate(Chi_0_l_p_matrix, resultMatrix)
    deallocate(F_pair_eigenvalues)
    deallocate(intermediate)
    deallocate(intermediate_two)
end program main

integer function delta(i, j)
    implicit none
    integer, intent(in) :: i, j
    if (i == j) then
        delta = 1
    else
        delta = 0
    endif
end function delta
