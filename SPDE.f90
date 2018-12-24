module SPDE
implicit none
private
    public :: dw, exit_traj
    real(8), parameter :: pi = 3.14159265358979323846_8
contains

function dw(ndim, dt) result(rvs)
implicit none
    integer, intent(in) :: ndim
    real(8), intent(in) :: dt
    real(8), dimension(ndim) :: rvs
    real(8) :: rand_u, tmp1
    integer :: i
    !-- Box-Muller algorithm (cos side)
    do i=1,ndim
        call random_number(rand_u)
        tmp1 = dsqrt(-2*dlog(rand_u))
        call random_number(rand_u)
        rvs(i) = tmp1*dcos(2*pi*rand_u)
    enddo
    rvs = rvs * dsqrt(dt)
end function dw

subroutine exit_traj(x, dt, esti)
    real(8) :: x(2), dt
    real(8), intent(out) :: esti
    real(8) :: xc(2)
    integer :: i
    
    if(x(1)**2+x(2)**2 > 1) stop 'initial point error'
    xc = x
    esti = 0.0
    do while (sum(xc**2) <= 1)
        !-- EM format
        xc = xc + xc * dt + dw(2,dt)
        !-- nonhomogeneous term
        esti = esti - ( sum(x**2) + 1.0 ) * dt
    enddo
    !-- boundary term
    esti = esti + 0.5
end subroutine exit_traj

end module SPDE

program elliptic_PDE
use SPDE
implicit none
    real(8), allocatable :: u_mesh(:,:)
    real(8) :: dt, x(2), esti_val
    integer :: i,j,k,N,NMC
    N = 20
    NMC = 1000 
    allocate(u_mesh(N,N))
    dt = 0.0005

    u_mesh = 0
    do i=1,N
        do j=1,N
            x(1) = real(i)/real(N) * sqrt(0.5)
            x(2) = real(j)/real(N) * sqrt(0.5)
            if(sum(x**2) > 1.0) exit
            do k=1,NMC
                call exit_traj(x, dt, esti_val)
                u_mesh(i,j) = u_mesh(i,j) + esti_val
            enddo
            u_mesh(i,j) = u_mesh(i,j) / real(NMC)
        enddo
    enddo
    
    open(unit=222, file='u_mesh.dat')
    do i=1,N
        write(222,*) u_mesh(i,:)
    enddo
    close(unit=222)
end program elliptic_PDE
