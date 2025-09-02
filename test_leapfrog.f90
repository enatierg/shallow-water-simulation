program test_leapfrog

use header
implicit none
include "mpif.h"
integer :: m, kmax, myid, nproc, ierr
real(kind=rk) :: dt, cg, sigma, x0, y0
integer :: n, size, ibeg, iend
real(kind=rk) :: nu, T, I_k, I_prev
type(Vector) :: u, v, eta
type(Matrix) :: Dx, Dy
integer :: k, numerr

! setting parameters
m = 64
dt = 0.004_rk
kmax = 1000
n = m**2
cg = 1.0_rk
sigma = 0.2_rk
x0 = 0.25_rk
y0 = 0.25_rk

call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

! calculating parameters
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

! allocating memory for fields
u%n = n
v%n = n
eta%n = n
allocate(u%xx(n))
allocate(v%xx(n))
allocate(eta%xx(n))

! testing leapfrog()
! setting up Dx and Dy
call create_matrices(m, n, Dx, Dy)

! initializing fields and numerr
call eta_peak(m, sigma, x0, y0, eta)
u%xx(:) = 0.0_rk
v%xx(:) = 0.0_rk
numerr = 0

! calculating initial average I(0)
call average(n, eta, I_prev, myid, nproc, ierr, ibeg, iend)
if (myid == 0) print*, 'initial average I(0)=', I_prev

do k = 1, kmax
    call leapfrog(m, 1, nu, Dx, Dy, u, v, eta, myid, nproc)
    call average(n, eta, I_k, myid, nproc, ierr, ibeg, iend)
    if (myid == 0) then
        ! checking whether I(k)-I(k-1)=0
        if (abs(I_k - I_prev) > 1e-10_rk) then
            print*, 'I(k) is not constant'
            print*, '|I_k-I_prev|=', abs(I_k - I_prev)
            numerr = numerr + 1
        end if
        I_prev = I_k
    end if
end do

if (myid == 0) then
    print*, 'exited loop with ', numerr, 'discrepancies in the average I(k)'
    if (numerr == 0) print*, 'leapfrog() runs successfully. good job!'
end if

! testing average()
u%xx(:) = 1.0_rk
call average(n, u, I_k, myid, nproc, ierr, ibeg, iend)

if (myid == 0) then
    if (abs(I_k - 1.0_rk) < 1e-10_rk) print*, 'average() runs successfully. good job!'
end if

deallocate(u%xx)
deallocate(v%xx)
deallocate(eta%xx)
deallocate(Dx%ii)
deallocate(Dx%aa)
deallocate(Dx%jj)
deallocate(Dy%aa)
deallocate(Dy%ii)
deallocate(Dy%jj)

call MPI_Finalize(ierr)

end program test_leapfrog
