subroutine Dx_mult(m, length, x, y)
    implicit none
    integer, intent(in) :: m, length
    real, dimension(m**2), intent(in) :: x
    real, dimension(m**2), intent(out) :: y

    integer :: i, j, k, start, ending

    call sparsegather(x, length)

    do k = 1, m
        start = (k-1)*m + 1
        ending = k*m
        y(start) = -x(k+2) + x(ending)

        do i = start+1, ending-1
            j = i - (k-1)*m
            y(i) = -x(j) + x(i)
        end do
        y(ending) = -x(start) + x(ending-1)
    end do
end subroutine Dx_mult
