subroutine timestepping_improved(n,theta,nu,mmax,Dx_hat,H_hat,w_vec,u,eta)
implicit none
integer, intent(in) :: n, mmax
real(kind=8), intent(in) :: theta, nu
real(kind=8), dimension(2,n), intent(in) :: Dx_hat
real(kind=8), dimension(4,n), intent(in) :: H_hat
real(kind=8), dimension(n), intent(in) :: w_vec
real(kind=8), dimension(n), intent(inout) :: u, eta
integer :: i
real(kind=8), dimension(n) :: r, eta_slash, eta_hat
do i = 1, mmax
call dgbmv('N', n, n, 0, 1, -nu*theta, Dx_hat, 2, eta, 1, 1.0_8, u, 1)
u(n) = u(n) - nu*theta*eta(1)
call dgbmv('T', n, n, 0, 1, nu, Dx_hat, 2, u, 1, 0.0_8, r, 1)
r(1) = r(1) + nu*u(n)
eta_hat = eta
call banded_helmholtz_solve(n, theta*nu, H_hat, w_vec, eta_hat, r)
eta = eta + eta_hat
eta_slash = (2*theta-1)*eta_hat + (1-theta)*eta
call dgbmv('N', n, n, 0, 1, -nu, Dx_hat, 2, eta_slash, 1, 1.0_8, u, 1)
u(n) = u(n) - nu*eta_slash(1)
end do
end subroutine timestepping_improved

