function field_delta(n,u1,eta1,u2,eta2) result(delta)
implicit none
integer, intent(in) :: n
real(kind=8), dimension(n), intent(in) :: u1, eta1, u2, eta2
real(kind=8) :: delta
real(kind=8), external :: dnrm2
delta = 1.0_8/n * (dnrm2(n, u1-u2, 1) + dnrm2(n, eta1-eta2, 1))
end function field_delta
