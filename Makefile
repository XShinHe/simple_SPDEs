FC=gfortran

default:
	${FC} SPDE.f90 -o test
