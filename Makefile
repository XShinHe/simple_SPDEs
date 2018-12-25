FC=gfortran

# default only calculate along the radius
default:
	${FC} SPDE_polar.f90 -o SDE
origin:
	${FC} SPDE.f90 -o SDE
clean:
	rm *.mod SDE
