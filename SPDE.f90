module SPDE
implicit none
private
    public :: init_seed, dw, exit_traj
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
    
    if(x(1)**2+x(2)**2 > 1.0) stop 'initial point error'
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
    real(8), allocatable :: u_mesh(:,:), u_mesh_exact(:,:)
    real(8) :: dt, x(2), du2, esti_val, test_val
    character(len=100) :: buffer
	integer :: i,j,k,N,NMC
    
	open(unit=111, file='sde.parm')
	read(111,*) buffer
	read(111,*) N, dt, NMC
    allocate(u_mesh(N,N), u_mesh_exact(N,N))

	call init_seed()

    u_mesh = 0
	u_mesh_exact = 0
	du2 = 0

	test_val = 0
	do i=1, NMC
		x = 0.5_8
		call exit_traj(x,dt,esti_val)
		test_val = test_val + esti_val
	enddo
	test_val = test_val / real(NMC)
	print *, '(0.5,0.5) exact 0.25  exp ', test_val
	stop 'debug'
	
    do i=1,N
        do j=1,N
            x(1) = real(i-1)/real(N)
            x(2) = real(j-1)/real(N)
            if(sum(x**2) > 1.0) exit
            do k=1,NMC
                call exit_traj(x, dt, esti_val)
                u_mesh(i,j) = u_mesh(i,j) + esti_val
            enddo
            u_mesh(i,j) = u_mesh(i,j) / real(NMC)
        	u_mesh_exact(i,j) = 0.5 * sum(x**2)
			du2 = du2 + (u_mesh(i,j)-u_mesh_exact(i,j))**2
		enddo
    enddo
    
	print *, 'dtime=', dt, ', with error2 = ', du2/real(N*N)

    open(unit=222, file='u_mesh.dat')
    open(unit=333, file='u_mesh_exact.dat')
	do i=1,N
        write(222,*) u_mesh(i,:)
		write(333,*) u_mesh_exact(i,:)
    enddo
    close(unit=222)
	close(unit=333)
end program elliptic_PDE
