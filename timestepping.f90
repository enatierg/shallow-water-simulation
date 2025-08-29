subroutine timestepping(n,theta,nu,mmax,Dx,H,u,eta)
implicit none
integer, intent(in) :: n, mmax
real(kind=8), intent(in) :: theta, nu
real(kind=8), dimension(n,n), intent(in) :: Dx, H
real(kind=8), dimension(n), intent(inout) :: u
real(kind=8), dimension(n), intent(inout) :: eta
real(kind=8), dimension(n) :: eta_hat, eta_slash, r
integer :: i
do i = 1, mmax
call dgemv('N', n, n, -nu*theta, Dx, n, eta, 1, 1.0_8, u, 1)
call dgemv('T', n, n, nu, Dx, n, u, 1, 0.0_8, r, 1)
call helmholtz_solve(n, H, eta_hat, r)
eta = eta + eta_hat
eta_slash = (2*theta-1)*eta_hat + (1-theta)*eta
call dgemv('N', n, n, -nu, Dx, n, eta_slash, 1, 1.0_8, u, 1)
end do
end subroutine timestepping
