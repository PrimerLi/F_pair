subroutine lev(vector, leading, leadingIndex)
    implicit none
    complex(kind = 8), intent(in) :: vector(:)
    complex(kind = 8), intent(out) :: leading
    integer, intent(out) :: leadingIndex
    integer :: i
    integer :: length
    real(kind = 8) :: minValue
    complex(kind = 8) :: temp
    real(kind = 8) :: eps
    integer, parameter :: dp = 8
    integer :: countLength 
    complex(kind = 8) :: one

    one = dcmplx(1, 0)

    length = size(vector)
    
    minValue = abs(vector(1) - one)
    leading = vector(1)
    leadingIndex = 1
    do i = 1, length
        if (abs(vector(i) - one) < minValue) then
            leading = vector(i)
            leadingIndex = i
	    minValue = abs(vector(i) - one)
        endif
    enddo

end subroutine lev
