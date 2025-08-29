program shallow_water
implicit none
real(kind=8), dimension(:), allocatable :: u, eta
real(kind=8), dimension(:,:), allocatable :: Dx, H
integer :: n
real(kind=8) :: cg
real(kind=8) :: deltat
real(kind=8) :: theta
real(kind=8) :: tfinal
real(kind=8) :: alpha, nu
integer :: mmax
real(kind=8) :: t_start, t_finish
open(unit=69, file='input.dat')
read(69,) n, deltat, cg, theta, tfinal
mmax = nint(tfinal/deltat)
nu = cg * deltat * n
alpha = nu * theta
write(,) "Shallow water equations"
write(,*) "=============="
write(*,)
write(,'("n = ", I12)') n
write(*,'("dx = ", F12.4)') 1.0_8/n
write(*,'("deltat = ", F12.4)') deltat
write(*,'("cg = ", F12.4)') cg
write(*,'("theta = ", F12.4)') theta
write(*,'("tfinal = ", F12.4)') tfinal
write(*,'("mmax = ", I12)') mmax
write(*,'("nu = ", F12.4)') nu
write(*,'("alpha = ", F12.4)') alpha
write(,)
allocate(u(n))
allocate(eta(n))
allocate(Dx(n,n))
allocate(H(n,n))
call initial_condition(n, u, eta)
call create_matrices(n, alpha, Dx, H)
call cpu_time(t_start)
call timestepping(n, theta, nu, mmax, Dx, H, u, eta)
call cpu_time(t_finish)
write(,) "Elapsed time = ", t_finish - t_start, " s"
write(,) "Time per timestep = ", (t_finish - t_start)/mmax, " s"
call save_fields(n, u, eta)
call test_helmholtz(n, H)
deallocate(u)
deallocate(eta)
deallocate(Dx)
deallocate(H)
end program shallow_water

