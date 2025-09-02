subroutine sparsegather(x,length)
    use header
    implicit none
    include "mpif.h"
    type(vector), intent(inout) :: x
    integer :: myid,numprocs,stat(MPI_STATUS_SIZE),ierr
    integer :: id_left, id_right,length
    call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
    id_right = myid+1
    id_left = myid-1
    if (numprocs > 1) then
        if (mod(myid,2) .eq. 0) then
            if (id_right < numprocs) then
                ! Send to Right neighbour
                call MPI_Send(x%xx(x%iend-length+1:x%iend),length,&
                    MPI_DOUBLE_PRECISION,id_right,id_right, MPI_COMM_WORLD,ierr)
            end if
            if (id_left >= 0) then
                ! Receive from Left neighbour
                call MPI_Recv(x%xx(modulo(x%ibeg-2,x%n)+1-length+1:modulo(x%ibeg-2,x%n)+1),length,&
                    MPI_DOUBLE_PRECISION,id_left, &
                    MPI_ANY_TAG,MPI_COMM_WORLD, stat,ierr)
                ! Send to Left neighbour
                call MPI_Send(x%xx(x%ibeg:x%ibeg+length-1),length,&
                    MPI_DOUBLE_PRECISION,id_left,id_left, MPI_COMM_WORLD,ierr)
            end if
            if (id_right < numprocs) then
                ! Receive from Right neighbour
                call MPI_Recv(x%xx(modulo(x%iend,x%n)+1:modulo(x%iend,x%n)+length),length,&
                    MPI_DOUBLE_PRECISION,id_right, &
                    MPI_ANY_TAG,MPI_COMM_WORLD, stat,ierr)
            end if
        else
            if (id_left >= 0) then
                ! Receive from left neighbour
                call MPI_Recv(x%xx(modulo(x%ibeg-2,x%n)+1-length+1:modulo(x%ibeg-2,x%n)+1),length,&
                    MPI_DOUBLE_PRECISION,id_left, &
                    MPI_ANY_TAG,MPI_COMM_WORLD, stat,ierr)
            end if
            if (id_right < numprocs) then
                ! Send to Right neighbour
                call MPI_Send(x%xx(x%iend-length+1:x%iend),length,&
                    MPI_DOUBLE_PRECISION,id_right,id_right, MPI_COMM_WORLD,ierr)
                ! Receive from Right neighbour
                call MPI_Recv(x%xx(modulo(x%iend,x%n)+1:modulo(x%iend,x%n)+length),length,&
                    MPI_DOUBLE_PRECISION,id_right, &
                    MPI_ANY_TAG,MPI_COMM_WORLD, stat,ierr)
            end if
            if (id_left >= 0) then
                ! Send to Left neighbour
                call MPI_Send(x%xx(x%ibeg:x%ibeg+length-1),length,&
                    MPI_DOUBLE_PRECISION,id_left,id_left, MPI_COMM_WORLD,ierr)
            end if
            if (myid==0) then
                call MPI_Send(x%xx(1:length),length,MPI_DOUBLE_PRECISION,&
                    numprocs-1,numprocs-1,MPI_COMM_WORLD,ierr)
                call MPI_Recv(x%xx(x%n-length+1:x%n),length,MPI_DOUBLE_PRECISION,&
                    numprocs-1,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
            else if (myid==numprocs-1) then
                call MPI_Send(x%xx(x%n-length+1:x%n),length,MPI_DOUBLE_PRECISION,&
                    0,0,MPI_COMM_WORLD,ierr)
                call MPI_Recv(x%xx(1:length),length,MPI_DOUBLE_PRECISION,&
                    0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
            end if
        end if
    end if
end subroutine sparsegather
