FC = gfortran
FLAGS = --fixed-line-length-none --free-line-length-none -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace

susy_qpi : susy_qpi.f90 dcuhre.o invert4x4.o inc/susy_qpi.f90 inc/write_data.f90 inc/sG0.f90
	$(FC) $(FLAGS) -o susy_qpi susy_qpi.f90 dcuhre.o invert4x4.o

dcuhre.o : dcuhre.f
	$(FC) $(FLAGS) -O -c dcuhre.f

invert4x4.o : invert4x4.f
	$(FC) $(FLAGS) -O -c invert4x4.f