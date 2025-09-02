subroutine Dy_mult(m, length, x, y)
    implicit none
    integer, intent(in) :: m, length
    real, dimension(m**2), intent(in) :: x
    real, dimension(m**2), intent(out) :: y

    integer :: i, j, k

    call sparsegather(x, length)

    do i = 1, m
        k = i + m
        j = m**2 - m + i
        y(i) = -x(k) + x(j)
    end do

    do i = m+1, m**2-m
        k = i + m
        j = i - m
        y(j) = -x(k) + x(j)
    end do

    do i = m**2-m+1, m**2
        k = i - m**2 + m
        j = i - m
        y(j) = -x(k) + x(j)
    end do

    ! call sparsegather(y, m**2)

    print*, "Dy_mult", y
end subroutine Dy_mult
