program shallow_water_improved
! Main program for solving the shallow water equations
implicit none

! Vectors for height perturbation and velocity
real(kind=8), dimension(:), allocatable :: u, eta
! banded matrices Dx_hat and H_hat
real(kind=8), dimension(:,:), allocatable :: Dx_hat, H_hat
! vector w_vec
real(kind=8), dimension(:), allocatable :: w_vec
! Numerical constants
integer :: n ! problem size n
real(kind=8) :: cg ! gravity wave speed
real(kind=8) :: deltat ! time step size Delta t
real(kind=8) :: theta ! off-centering parameter
real(kind=8) :: tfinal ! final time
real(kind=8) :: alpha, nu ! alpha and Courant number
integer :: mmax ! Number of time steps
! Start and finish times
real(kind=8) :: t_start, t_finish

! Read parameters from file
open(unit=69, file='input.dat')
read(69,*) n, deltat, cg, theta, tfinal

mmax = nint(tfinal/deltat)
nu = cg * deltat * n
alpha = nu * theta

write(*,*) "Shallow water equations"
write(*,*) "---"
write(*,*)
write(*,'("n      = ",I6)') n
write(*,'("dx     = ",F12.4)') 1.0_8/n
write(*,'("dt     = ",F12.4)') deltat
write(*,'("cg     = ",F12.4)') cg
write(*,'("theta  = ",F12.4)') theta
write(*,'("tfinal = ",F12.4)') tfinal
write(*,'("mmax   = ",I12)') mmax
write(*,'("nu     = ",F12.4)') nu
write(*,'("alpha  = ",F12.4)') alpha

! allocate memory
allocate(u(n))
allocate(eta(n))
allocate(Dx_hat(2,n))
allocate(H_hat(4,n))
allocate(w_vec(n))

! Initialize velocity and height perturbation
call initial_condition(n,u,eta)
! create banded matrices Dx_hat, H_hat, w_vec
call create_banded_matrices(n,alpha,Dx_hat,H_hat,w_vec)

! Integrate shallow water equations
call cpu_time(t_start)
call timestepping_improved(n,theta,nu,mmax,Dx_hat,H_hat,w_vec,u,eta)
call cpu_time(t_finish)
write(*,*) "Elapsed time = ",t_finish-t_start," s"
write(*,*) "Time per timestep = ",(t_finish-t_start)/mmax," s"
! save fields to disk
call save_fields(n,u,eta)

! adding a line that runs test_helmholtz to test banded_helmholtz_solve
call test_helmholtz(n,alpha,H_hat,w_vec)

! Deallocate memory
deallocate(u)
deallocate(eta)
deallocate(Dx_hat)
deallocate(H_hat)
deallocate(w_vec)

end program shallow_water_improved
