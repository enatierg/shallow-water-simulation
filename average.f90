subroutine average(n,phi,s,myid,nproc,ierr,ibeg,iend)

use header
implicit none
include "mpif.h"

integer, intent(in) :: n, nproc    ! Dimension of vector
type(Vector), intent(in) :: phi    ! Vector to average
real(kind=rk), intent(out) :: s    ! Resulting average
real(kind=rk) :: suml, global_sum
integer, intent(in) :: myid, ierr, ibeg, iend
integer :: iter

suml = 0.0_rk
global_sum = 0.0_rk

! calculating local sum
do iter = ibeg, iend
    suml = suml + phi%xx(iter)
end do

! calculating global sum
call MPI_AllReduce(suml, global_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

! calculating the average
s = global_sum / real(n, rk)

end subroutine average
