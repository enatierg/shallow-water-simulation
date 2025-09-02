subroutine create_matrices(m,n,Dx,Dy)

use header

implicit none

integer, intent(in) :: m    ! Number of gridpoints in one direction  
integer, intent(in) :: n    ! Size of vectors  
type(Matrix), intent(inout) :: Dx    ! ) Discretisation matrices for  
type(Matrix), intent(inout) :: Dy    ! ) derivatives in x- and y-direction  

integer :: nrow, ip1, im1, jp1, jm1, i, j ! Loop variables  
integer :: inzDx    ! Counter for ii (D_x)  
integer :: inzDy    ! Counter for ii (D_y)  

! Set up discretisation matrices  
! Set size of matrix and total numbers of non-zeros  

Dx%n = n  
Dx%nnz = 2*n  
Dy%n = n  
Dy%nnz = 2*n  

! ... and allocate memory (need to deallocate this before exiting)  
allocate(Dx%aa(Dx%nnz))  
allocate(Dx%ii(Dx%n+1))  
allocate(Dx%jj(Dx%nnz))  

allocate(Dy%aa(Dy%nnz))  
allocate(Dy%ii(Dy%n+1))  
allocate(Dy%jj(Dy%nnz))  

! Set up matrix Dx  

inzDx = 1  

! Loop over matrix rows  
do nrow = 1, n  
! Work out 2d coordinates i and j  
j = (nrow - 1)/m + 1
i = nrow - (j - 1) * m
! === Dx ===
Dx%ii(nrow) = inzDx
! Work out coordinates of neighbours in x-direction
ipl = i+1
iml = i-1
! Use periodic boundary conditions
if (ipl > m) ipl = ipl - m
if (iml < 1) iml = iml + m
! Set entries in matrix Dx
Dx%jj(inzDx) = iml + (j-1)*m
Dx%aa(inzDx) = -1.0_rk
inzDx = inzDx + 1
Dx%jj(inzDx) = ipl + (j-1)*m
Dx%aa(inzDx) = +1.0_rk
inzDx = inzDx + 1
end do
Dx%ii(n+1) = inzDx

*********************************************** setting up Dy *********************
inzDy=1
do nrow=1,n
! working out 2d coordinates i and j
j=(nrow-1)/m+1
i=nrow-(j-1)*m
! === Dy ===
Dy%ii(nrow)=inzDy
! working out coordinates of neighbours in x-direction
jpl=j+1
jml=j-1
! using periodic boundary conditions
if (jpl>m) jpl=jpl-m
if (jml<1) jml=jml+m
! allocating entries in matrix Dy
Dy%jj(inzDy)=(jpl-1)*m+i
Dy%aa(inzDy)=+1.0_rk
inzDy=inzDy+1
Dy%jj(inzDy)=(jml-1)*m+i
Dy%aa(inzDy)=-1.0_rk
inzDy=inzDy+1
end do
Dy%ii(n+1)=inzDy
! ***********************
end subroutine create_matrices
