subroutine create_matrices(n,alpha,Dx,H)
implicit none
integer, intent(in) :: n
real(kind=8), intent(in) :: alpha
real(kind=8), dimension(n,n) :: Dx
real(kind=8), dimension(n,n) :: H
integer :: i
Dx(:,:) = 0.0_8
do i = 1, n
Dx(i,i) = -1.0_8
if (i == n) then
Dx(i,1) = +1.0_8
else
Dx(i,i+1) = +1.0_8
end if
end do
call dgemm('T','N',n,n,n,alpha**2.0_8,Dx,n,Dx,n,0.0_8,H,n)
do i=1,n
H(i,i)=H(i,i)+1.0_8
end do
end subroutine create_matrices
