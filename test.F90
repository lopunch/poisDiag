program driver

    use poisDiag, poisson_solver_3d => poisson_solver_3d_sp
!   use poisDiag, poisson_solver_3d => poisson_solver_3d_dp

    implicit none

    type(poisson_solver_3d)               :: solver
    integer                               :: nx, ny, nz
    real,   dimension(:,:,:), allocatable :: tmp
!   real*8, dimension(:,:,:), allocatable :: f, p
    real  , dimension(:,:,:), allocatable :: f, p

    integer             :: i

    open(11,file='poisson_test.dat',form='unformatted',access='sequential')
    read(11) nx
    read(11) ny
    read(11) nz

    allocate(tmp(nx,ny,nz), &
               f(nx,ny,nz), &
               p(nx,ny,nz))

    read(11) tmp
    close(11)
    
!   f = dble(tmp)
    f = tmp
    deallocate(tmp)

!   call solver % init( [nx, ny, nz], [3000.d0, 3000.d0, 500.d0], [(neumann, i = 1,6)] )
    call solver % init( [nx, ny, nz], [3000.e0, 3000.e0, 500.e0], [(neumann, i = 1,6)] )
    call solver % solve(f, p)

    open(11,file='p.dat',form='unformatted',access='sequential')
    write(11) nx
    write(11) ny
    write(11) nz
    write(11) f
    write(11) p
    close(11)

    call solver % destroy

    deallocate(f, p)
end program driver
