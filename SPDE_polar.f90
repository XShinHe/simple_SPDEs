module SPDE
implicit none
private
    public :: init_seed, dw, exit_traj, pi
    real(8), parameter :: pi = 3.14159265358979323846_8
contains

subroutine init_seed()
 	integer :: n, ival(8), v(3), i
	integer, allocatable :: seed(:)
    call date_and_time(values=ival)
    v(1) = ival(8) + 2048*ival(7)
    v(2) = ival(6) + 64*ival(5)     ! value(4) isn't really 'random'
    v(3) = ival(3) + 32*ival(2) + 32*8*ival(1)
    call random_seed(size=n)
    allocate(seed(n))
    call random_seed()   ! Give the seed an implementation-dependent kick
    call random_seed(get=seed)
    do i=1, n
        seed(i) = seed(i) + v(mod(i-1, 3) + 1)
    enddo
    call random_seed(put=seed)
    deallocate(seed)
end subroutine init_seed

function dw(ndim, dt) result(rvs)
implicit none
    integer, intent(in) :: ndim
    real(8), intent(in) :: dt
    real(8), dimension(ndim) :: rvs
    real(8) :: rand0, tmp1
    integer :: i
    !-- Box-Muller algorithm (cos side)
    do i=1,ndim
        call random_number(rand0)
        tmp1 = dsqrt(-2*dlog(rand0))
        call random_number(rand0)
        rvs(i) = tmp1*dcos(2*pi*rand0)
    enddo
    !-- scale by sqrt(dt)
    rvs = rvs * dsqrt(dt)
end function dw

subroutine exit_traj(x, dt, esti)
    real(8) :: x(2), dt
    real(8), intent(out) :: esti
    real(8) :: xc(2)
    integer :: i
    
    if(x(1)**2+x(2)**2 > 1.0) stop 'initial point error'
    xc = x
    esti = 0.0
    do while (sum(xc**2) <= 1)
        !-- EM format
        xc = xc + xc * dt + dw(2,dt)
        !-- nonhomogeneous term, g(xc) = 1 + xc(1)**2 + xc(2)**2
        esti = esti - ( sum(xc**2) + 1.0 ) * dt
    enddo
    !-- boundary term, f(x_exit) = 0.5
    esti = esti + 0.5_8
end subroutine exit_traj

end module SPDE

program elliptic_PDE
use SPDE
implicit none
    real(8), allocatable :: u_mesh(:), u_mesh_exact(:)
    real(8) :: dt, x(2), du2, esti_val, dr
    character(len=100) :: buffer
	integer :: i,j,k,N,NMC
    
	open(unit=111, file='sde.parm')
	read(111,*) buffer
	read(111,*) N, dt, NMC
    allocate(u_mesh(N), u_mesh_exact(N))

	call init_seed()

    u_mesh = 0
	u_mesh_exact = 0
	du2 = 0
	dr = 1.0/real(N)
	
	!-- choose a polar axis
    do i=1,N
        x(1) = real(i-1)/real(N)
        x(2) = 0.0_8
        do k=1,NMC
            call exit_traj(x, dt, esti_val)
            u_mesh(i) = u_mesh(i) + esti_val
        enddo
        u_mesh(i) = u_mesh(i) / real(NMC)
        u_mesh_exact(i) = 0.5 * sum(x**2)
        
        !-- L2 norm of u-u_exact in D
		du2 = du2 + 2*pi*x(1)*dr*(u_mesh(i)-u_mesh_exact(i))**2
		
		print *, 'complete ', real(100*i)/real(N), '%'
    enddo
    
	print *, 'dtime=', dt, ', with error2 = ', du2
	print *, sum(2*pi*0.5*0.1*(u_mesh-u_mesh_exact)**2)

    open(unit=222, file='u_mesh.dat')
    open(unit=333, file='u_mesh_exact.dat')
    do i=1,N
    write(222,*) u_mesh(i)
    write(333,*) u_mesh_exact(i)
    enddo
    close(unit=222)
	close(unit=333)
	
end program elliptic_PDE
