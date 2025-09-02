subroutine Mat_Mult(A,u,b,length)
    use header
    implicit none
    type(Matrix), intent(in) :: A    ! Matrix to multiply with
    type(Vector), intent(in) :: u    ! Input vector u
    type(Vector), intent(inout) :: b ! Output vector b = A*u
    integer, intent(in) :: length
    integer :: i, i_j, j

    ! Calculate each component of b by taking the scalar product of the i-th row of A and the vector u
    call sparsegather(u,length)

    do i = b%ibeg, b%iend
        b%xx(i) = 0.0_rk
        do i_j = int(A%ii(i)), int(A%ii(i+1)-1)
            j = A%jj(i_j)
            b%xx(i) = b%xx(i) + A%aa(i_j) * u%xx(j)
        end do
    end do
end subroutine Mat_Mult
