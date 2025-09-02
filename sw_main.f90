program sw_main

use header
implicit none
include "mpif.h"

! Model parameters
! ---
integer :: m = 256    ! Horizontal grid size
real(kind=rk) :: dt = 0.0005_rk    ! Size of leapfrog time step
integer :: kmax = 1000    ! Number of time steps
real(kind=rk) :: cg = 1.0_rk    ! Phase velocity
real(kind=rk) :: sigma = 0.2_rk    ! ) Width and position of peak
real(kind=rk) :: x0 = 0.25_rk    ! ) in initial condition
real(kind=rk) :: y0 = 0.25_rk    !

! Other variables
integer :: n    ! Problem size n = m^2
real(kind=rk) :: nu    ! Courant number nu = cg*dt/dx
real(kind=rk) :: T    ! Final time, T = kmax*dt
real(kind=rk) :: s    ! average over eta

! Model fields
type(Vector) :: u    ! ) 2d velocity field
type(Vector) :: v    ! )
type(Vector) :: eta    ! Discrete solution
type(Vector) :: eike

! Discretisation matrices
type(Matrix) :: Dx, Dy    ! Matrices for derivatives

integer :: ibeg, iend, myid, nproc, ierr, size

real(kind=rk) :: t0, tkmax, ttotal    ! Timing variables

call MPI_Init(ierr)
call MPI_Comm_Rank(MPI_COMM_WORLD, myid, ierr)
call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr)

open(69, file='inputs.dat')
read(69, *) m, dt, kmax

! Calculate matrix size and final time
n = m**2
T = kmax*dt
nu = cg*dt*m
size = n/nproc
ibeg = (myid*size)+1
iend = (myid+1)*size
Dx%ibeg = ibeg
Dy%ibeg = ibeg
u%ibeg = ibeg
v%ibeg = ibeg
Dx%iend = iend
Dy%iend = iend
u%iend = iend
v%iend = iend

! Print out model parameters
if (myid == 0) then
    write(*,'("Shallow water equations")')
    write(*,'("---")')
    write(*,'("m    = ",I10)') m
    write(*,'("matrix size  = ",I10)') n
    write(*,'("dt    = ",F10.4)') dt
    write(*,'("dx    = ",F10.4)') 1.0_rk/m
    write(*,'("cg    = ",F10.4)') cg
    write(*,'("nu    = ",F10.4)') nu
    write(*,'("T    = ",F10.4)') T
    write(*,'("timesteps   = ",I10)') kmax
    write(*,'("sigma    = ",F10.4)') sigma
end if

! Allocate memory for fields
u%n = n
v%n = n
eta%n = n
allocate(u%xx(n))
allocate(v%xx(n))
allocate(eta%xx(n))

! Setup matrices D_x and D_y
call create_matrices(m, n, Dx, Dy)

! ==== Test case: peak in initial conditions ====
! Initialise fields
call eta_peak(m, sigma, x0, y0, eta)
! u = 0
! v = 0
u%xx(:) = 0.0_rk
v%xx(:) = 0.0_rk

! Calculate average over eta at t=0
call average(n, eta, s, myid, nproc, ierr, ibeg, iend)
if (myid == 0) write(*,'("Initial average <eta> = ",E24.18)') s

t0 = MPI_Wtime()
! Forward timestepping with leapfrog algorithm
call leapfrog(m, kmax, nu, Dx, Dy, u, v, eta, myid, nproc)
tkmax = MPI_Wtime()

ttotal = tkmax - t0
if (myid == 0) write(*,'("Total time it takes to run leapfrog() = ",E24.18)') 1E6*(ttotal)/real(kmax)

eike = eta
call MPI_Gather(eta%xx(ibeg:iend), size, MPI_DOUBLE_PRECISION, eike%xx, size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

! Calculate average over eta at t=T
call average(n, eta, s, myid, nproc, ierr, ibeg, iend)

do m = ibeg, iend
    eta%xx(m) = (cg**2)*(u%xx(m)**2 + v%xx(m)**2) / 2.0_rk
end do
call MPI_AllGather(eta%xx(ibeg:iend), size, MPI_DOUBLE_PRECISION, eta%xx, size, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

if (myid == 0) write(*,'("Average at t=T <eta> = ",E24.18)') s

call average(n, eta, s, myid, nproc, ierr, ibeg, iend)
if (myid == 0) write(*,'("Average kinetic energy = ",E24.18)') s

! Save solution
if (myid == 0) call save_field(m, eike, "eta.dat")

! Deallocate all memory
! Vectors
deallocate(u%xx)
deallocate(v%xx)
deallocate(eta%xx)

! Sparse matrices
deallocate(Dx%ii)
deallocate(Dx%jj)
deallocate(Dx%aa)
deallocate(Dy%ii)
deallocate(Dy%jj)
deallocate(Dy%aa)

call MPI_Finalize(ierr)

end program sw_main
