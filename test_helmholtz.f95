subroutine test_helmholtz(n,H)
implicit none
integer, intent(in) :: n
real(kind=8), dimension(n,n), intent(in) :: H
real(kind=8), dimension(n) :: eta_in, eta, r
real(kind=8) :: norm
real(kind=8), external :: dnrm2
call random_number(eta_in(:))
eta = eta_in
call dgemv('N', n, n, 1.0_8, H, n, eta_in, 1, 0.0_8, r, 1)
call helmholtz_solve(n, H, eta, r)
norm = dnrm2(n, eta - eta_in, 1)
print *, '||eta-eta_in||2=', norm
end subroutine test_helmholtz
