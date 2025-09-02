subroutine leapfrog(m,kmax,nu,Dx,Dy,u,v,eta,rank,nproc)

use header
implicit none
include "mpif.h"

integer, intent(in) :: m    ! Number of grid cells in one direction
integer, intent(in) :: kmax    ! Number of time steps
real(kind=rk), intent(in) :: nu    ! Courant number
type(Matrix), intent(in) :: Dx    !) Discretisation matrices for derivatives
type(Matrix), intent(in) :: Dy    !) in x- and y-direction
type(Vector), intent(inout) :: u    !) Components of velocity field
type(Vector), intent(inout) :: v    !)
type(Vector), intent(inout) :: eta ! Height perturbation field
integer, intent(in) :: rank, nproc

! *******************************
integer :: k    ! for loop variable
integer :: ibeg, iend    ! mpi variables
real(kind=rk) :: kappa
type(Vector) :: u_, v_, eta_, me    ! temporary vectors

! initializing the fields at t=0:
ibeg = u%ibeg
iend = u%iend

! Allocate and initialize previous time step vectors
u_%n = u%n
u_%ibeg = ibeg
u_%iend = iend
allocate(u_%xx(u_%n))
u_%xx = u%xx

v_%n = v%n
v_%ibeg = ibeg
v_%iend = iend
allocate(v_%xx(v_%n))
v_%xx = v%xx

eta_%n = eta%n
eta_%ibeg = ibeg
eta_%iend = iend
allocate(eta_%xx(eta_%n))
eta_%xx = eta%xx

! Temporary vector
me%n = m*m
me%ibeg = ibeg
me%iend = iend
allocate(me%xx(me%n))

! running the for loop:
do k = 1, kmax
    ! setting the kappa values
    if (mod(k,2) == 0) then
        kappa = -0.5_rk * nu
    else
        kappa = -nu
    end if

    ! calculating the fields at the next timestep using eqns. (15)
    call Mat_Mult(Dx, eta, me, m-1)
    u%xx = u%xx + kappa * me%xx
    call Mat_Mult(Dy, eta, me, m)
    v%xx = v%xx + kappa * me%xx
    call Mat_Mult(Dx, u, me, m-1)
    eta%xx = eta%xx + kappa * me%xx
    call Mat_Mult(Dy, v, me, m)
    eta%xx = eta%xx + kappa * me%xx

    ! swapping fields by rewriting the temporary vector
    me%xx = eta%xx
    eta%xx = eta_%xx
    eta_%xx = me%xx
    me%xx = u%xx
    u%xx = u_%xx
    u_%xx = me%xx
    me%xx = v%xx
    v%xx = v_%xx
    v_%xx = me%xx
end do

deallocate(me%xx)
deallocate(u_%xx)
deallocate(v_%xx)
deallocate(eta_%xx)
! ***********************************************************************
end subroutine leapfrog
