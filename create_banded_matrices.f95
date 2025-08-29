subroutine create_banded_matrices(n,alpha,Dx_hat,H_hat,w_vec)
implicit none
integer, intent(in) :: n
real(kind=8), dimension(2,n), intent(inout) :: Dx_hat
real(kind=8), dimension(4,n), intent(inout) :: H_hat
real(kind=8), intent(in) :: alpha
real(kind=8), dimension(n), intent(inout) :: w_vec
real(kind=8), dimension(n) :: h, h_vec
real(kind=8) :: rho
integer :: info
integer, dimension(n) :: ipiv
real(kind=8), external :: ddot
H_hat(1,:) = 0.0_8
H_hat(2,:) = alpha**2.0_8
H_hat(2,1) = 0.0_8
H_hat(3,:) = 1.0_8 + 2.0_8 * alpha**2.0_8
H_hat(3,1) = 1.0_8 + 3.0_8 * alpha**2.0_8
H_hat(3,n) = 1.0_8 + 3.0_8 * alpha**2.0_8
H_hat(4,:) = alpha**2.0_8
H_hat(4,n) = 0.0_8
Dx_hat(1,:) = 1.0_8
Dx_hat(1,1) = 0.0_8
Dx_hat(2,:) = -1.0_8
h(:) = 0.0_8
h(1) = alpha
h(n) = alpha
h_vec = h
call dgbsv(n, 1, 1, 1, H_hat, 4, ipiv, h_vec, n, info)
if (info /= 0) then
print *, 'dgbsv in create_banded_matrices is wrong'
print *, info
end if
rho = ddot(n, h_vec, 1, h, 1)
w_vec = h_vec / (1.0_8 - rho)
end subroutine create_banded_matrices
