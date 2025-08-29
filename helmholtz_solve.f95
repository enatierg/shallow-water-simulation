subroutine helmholtz_solve(n,H,eta,r)
implicit none
integer, intent(in) :: n
real(kind=8), dimension(n,n), intent(in) :: H
real(kind=8), dimension(n), intent(inout) :: eta
real(kind=8), dimension(n), intent(in) :: r
integer :: info
real(kind=8), dimension(n,n) :: H_kopija
H_kopija = H
eta = r
call dposv('L',n,1,H_kopija,n,eta,n,info)
if (info /= 0) then
print, 'something is wrong with dposv in helmholtz_solve'
print, info
return
end if
end subroutine helmholtz_solve
